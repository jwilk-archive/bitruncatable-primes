// Copyright © 2019 Jakub Wilk <jwilk@jwilk.net>
// SPDX-License-Identifier: MIT

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gmpxx.h>

int main(int argc, char **argv)
{
    std::vector<mpz_class> primes = {2, 3, 5, 7};
    int len = 1;
    mpz_class record;
    const int est_max_c = 332579483;
    int c = 0;
    auto t_start = std::chrono::steady_clock::now();
    while (primes.size()) {
        std::vector<mpz_class> new_primes;
        #pragma omp parallel for
        for (size_t i = 0; i < primes.size(); i++) {
            std::vector<char> s(len + 3);
            mpz_get_str(&s[1], 10, primes[i].get_mpz_t());
            for (int d1 = 1; d1 <= 9; d1++) {
                s[0] = d1 + '0';
                s[len + 1] = '0';
                mpz_class n0(s.data());
                for (int d2 = 1; d2 <= 9; d2 += 2) {
                    if (d2 == 5)
                        continue;
                    mpz_class n = n0 + d2;
                    if (mpz_probab_prime_p(n.get_mpz_t(), 16))
                        #pragma omp critical
                        new_primes.push_back(n);
                }
            }
            #pragma omp atomic
            c++;
            if ((c & 0xFFFF) == 0 && c < est_max_c)
            {
                auto t_now = std::chrono::steady_clock::now();
                double dt = std::chrono::duration_cast<std::chrono::duration<double>>(t_now - t_start).count();
                double eta_sec = dt * (est_max_c - c) / c;
                auto eta = std::div(eta_sec / 60.0, 60);
                #pragma omp critical
                std::cerr << "ETA: " << eta.quot << "h " << eta.rem << "m" << std::endl;
            }
        }
        if (new_primes.size() == 0)
            record = *std::max_element(primes.begin(), primes.end());
        primes = std::move(new_primes);
        len += 2;
        std::cerr << "len=" << len << ", #=" << primes.size() << std::endl;
    }
    // double-check with higher number of Miller-Rabin iterations
    mpz_class m = record;
    while (true) {
        assert(mpz_probab_prime_p(m.get_mpz_t(), 40));
        std::string s = m.get_str();
        std::cout << s << std::endl;
        if (s.length() == 1)
            break;
        m = mpz_class(s.substr(1, s.length() - 2));
    }
}

// vim:ts=4 sts=4 sw=4 et
