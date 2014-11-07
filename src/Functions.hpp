#include <gmpxx.h>
#include <iostream>

mpz_class BinaryGCD(mpz_class& u, mpz_class v) {
    if (u == 0) return v;
    if (v == 0) return u;
    mpz_class g = 1;
    while ((u % 2) == 0 && (v % 2) == 0) {
        u >>= 1; // u / 2
        v >>= 1;
        g <<= 1; // 2 * g
    }
    while (u > 0) {
        if (u % 2 == 0) {
            u >>= 1;
        } else if (v % 2 == 0) {
            v >>= 1;
        } else {
            mpz_class t = abs(u - v) / 2;
            if (u < v) {
                v = t;
            } else {
                u = t;
            }
        }
    }
    return v * g;
}

mpz_class ModularExponentiation(mpz_class& base, mpz_class& exponent, const mpz_class& modulus) {
    mpz_class result = 1;
    base = base % modulus;

    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modulus;
        }
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}

bool MillerRabinPrime(const mpz_class& n, size_t k, gmp_randclass& rand) {
    if (n < 4) { return true; }

    mpz_class d  = n - 1;
    mpz_class s  = 0;

    // S is the maximal power of dividing n - 1
    // n - 1 = 2^s d
    while (d % 2 == 0) {
        d >>= 1;
        ++s;
    }

    // k times
    for (size_t i = 0; i < k; ++i) {
        mpz_class a = 2 + rand.get_z_range(n-3); // Random number [2, n-2]
        mpz_class x = ModularExponentiation(a, d, n); // Debugging: mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

        if (x == 1 || x == (n - 1)) {
            continue;
        }

        // s - 1 times
        for (mpz_class j = 0; j < s; ++j) {
            x = (x * x) % n;
            if (x == 1) {
                return false;
            } else if (x == (n - 1)) {
                break;
            }
        }
        if (x != (n - 1)) {
            return false;
        }
    }
    return true;
}

mpz_class PollardsRhoAlgorithm(const mpz_class& n, gmp_randclass& rand) {
    if (n % 2 == 0) return 2; // TODO
    if (n % 3 == 0) return 3; // TODO
    if (n % 5 == 0) return 5; // TODO
    if (n % 7 == 0) return 7; // TODO

    mpz_class x = rand.get_z_range(n+1); // [0, n]
    mpz_class y = x;
    mpz_class c = rand.get_z_range(n+1);
    mpz_class d = 1;
    mpz_class diff_abs = 0;

    auto g = [] (mpz_class& x, mpz_class& c, const mpz_class& n) -> mpz_class { return ((x * x ) % n + c) % n; };

    while (d == 1) {
        x = g(x, c, n);
        y = g(y, c, n);
        y = g(y, c, n);
        diff_abs = abs(x-y);
        d = BinaryGCD(diff_abs, n); // Debug: mpz_gcd (d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
        if (d == n) { // TODO
            break;
        }
    }
    return d;
}

