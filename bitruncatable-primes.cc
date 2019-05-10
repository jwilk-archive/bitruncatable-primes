// Copyright © 2019 Jakub Wilk <jwilk@jwilk.net>
// SPDX-License-Identifier: MIT

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gmpxx.h>

static mpz_class record;

static constexpr intmax_t est_max_c = 331780864;
static intmax_t cnt = 4;

static constexpr int max_half_width = 63;

static intmax_t stats[max_half_width] = {4};

double timer()
{
    static auto t_start = std::chrono::steady_clock::now();
    auto t_now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(t_now - t_start);
    return duration.count();
}

std::div_t prev_eta = {};

static bool operator ==(std::div_t x, std::div_t y)
{
    return (x.quot == y.quot) && (x.rem == y.rem);
}
static bool operator !=(std::div_t x, std::div_t y)
{
    return !(x == y);
}

class ETA
{
private:
    std::div_t prev;
public:
    void update(double time)
    {
        std::div_t eta = std::div(time / 60.0, 60);
        if (eta != this->prev) {
            std::cerr << "ETA: " << eta.quot << "h " << eta.rem << "m" << std::endl;
            this->prev = eta;
        }
    }
}

static eta = ETA();

static void explore(char *s, int w)
{
    intmax_t c;
    #pragma omp atomic
    stats[w]++;
    #pragma omp atomic capture
    c = ++cnt;
    if ((c & 0xFFF) == 0 && c < est_max_c)
    {
        double dt = timer();
        double eta_sec = dt * (est_max_c - c) / c;
        #pragma omp critical
        eta.update(eta_sec);
    }
    bool dead_end = true;
    w += 1;
    assert(w <= max_half_width);
    for (int d1 = 1; d1 <= 9; d1++) {
        s[-w] = d1 + '0';
        s[+w] = '1';
        mpz_class n(s - w, 10);
        for (int d2 = 1; d2 <= 9; d2 += 2, n += 2) {
            if (d2 == 5)
                continue;
            if (mpz_probab_prime_p(n.get_mpz_t(), 16)) {
                s[+w] = '0' + d2;
                explore(s, w);
                dead_end = false;
            }
        }
    }
    s[-w] = '\0';
    s[+w] = '\0';
    w--;
    if (dead_end)
    #pragma omp critical
    {
        mpz_class n(s - w, 10);
        if (n > record) {
             std::cout << n << " (" << (2 * w + 1) << " digits)" << std::endl;
             record = n;
        }
    }
}

int main(int argc, char **argv)
{
    std::vector<mpz_class> initial_primes;
    for (int i = 100; i <= 999; i++) {
        switch ((i / 10) % 10) {
        case 2:
        case 3:
        case 5:
        case 7:
            mpz_class j(i);
            if (mpz_probab_prime_p(j.get_mpz_t(), 0))
                initial_primes.push_back(j);
        }
    }
    assert(initial_primes.size() == 59);
    #pragma omp parallel for
    for (size_t i = 0; i < initial_primes.size(); i++) {
        char s[2 * max_half_width + 2] = {};
        mpz_get_str(s + max_half_width - 1, 10, initial_primes[i].get_mpz_t());
        explore(s + max_half_width, 1);
    }
    // double-check with higher number of Miller-Rabin iterations
    mpz_class m = record;
    while (true) {
        assert(mpz_probab_prime_p(m.get_mpz_t(), 40));
        std::string s = m.get_str();
        if (s.length() == 1)
            break;
        m = mpz_class(s.substr(1, s.length() - 2), 10);
    }
    std::cout << "# Statistics (length, count):" << std::endl;
    for (int w = 0; stats[w]; w++)
        std::cout << "#\t" << 2 * w + 1 << "\t" << stats[w] << std::endl;
}

// vim:ts=4 sts=4 sw=4 et
