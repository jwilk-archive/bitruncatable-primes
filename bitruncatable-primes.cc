// Copyright © 2019 Jakub Wilk <jwilk@jwilk.net>
// SPDX-License-Identifier: MIT

#include <cassert>
#include <iostream>
#include <vector>

#include <gmpxx.h>

int main(int argc, char **argv)
{
    std::vector<mpz_class> primes = {2, 3, 5, 7};
    int len = 1;
    mpz_class record;
    while (primes.size()) {
        std::vector<mpz_class> new_primes;
        #pragma omp parallel for
        for (int d1 = 1; d1 <= 9; d1++) {
            std::vector<char> s(len + 3);
            s[0] = d1 + '0';
            for (mpz_class &p : primes) {
                mpz_get_str(&s[1], 10, p.get_mpz_t());
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
