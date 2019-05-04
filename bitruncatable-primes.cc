// Copyright © 2019 Jakub Wilk <jwilk@jwilk.net>
// SPDX-License-Identifier: MIT

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gmpxx.h>

static mpz_class record;

static const int est_max_c = 332579483;
static int c = 0;

auto t_start = std::chrono::steady_clock::now();

void explore(std::string s)
{
    #pragma omp atomic
    c++;
    if ((c & 0xFFF) == 0 && c < est_max_c)
    {
        auto t_now = std::chrono::steady_clock::now();
        double dt = std::chrono::duration_cast<std::chrono::duration<double>>(t_now - t_start).count();
        double eta_sec = dt * (est_max_c - c) / c;
        auto eta = std::div(eta_sec / 60.0, 60);
        #pragma omp critical
        std::cerr << "ETA: " << eta.quot << "h " << eta.rem << "m" << std::endl;
    }
    size_t len = s.length();
    std::vector<char> z(len + 3);
    std::copy(s.begin(), s.end(), z.begin() + 1);
    bool dead_end = true;
    for (int d1 = 1; d1 <= 9; d1++) {
        z[0] = d1 + '0';
        z[len + 1] = '0';
        mpz_class n0(z.data());
        for (int d2 = 1; d2 <= 9; d2 += 2) {
            if (d2 == 5)
                continue;
            mpz_class n = n0 + d2;
            if (mpz_probab_prime_p(n.get_mpz_t(), 16)) {
                z[len + 1] = '0' + d2;
                explore(z.data());
                dead_end = false;
            }
        }
    }
    if (dead_end)
    #pragma omp critical
    {
        mpz_class n(s);
        if (n > record) {
             std::cout << n << " (" << s.length() << " digits)" << std::endl;
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
    for (size_t i = 0; i < initial_primes.size(); i++)
        explore(initial_primes[i].get_str());
    // double-check with higher number of Miller-Rabin iterations
    mpz_class m = record;
    while (true) {
        assert(mpz_probab_prime_p(m.get_mpz_t(), 40));
        std::string s = m.get_str();
        if (s.length() == 1)
            break;
        m = mpz_class(s.substr(1, s.length() - 2));
    }
}

// vim:ts=4 sts=4 sw=4 et
