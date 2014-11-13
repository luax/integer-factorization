#include "Functions.hpp"
#include <vector>
#include <cmath>
#include <cassert>

#define pr(a) std::cout << a << std::endl
// // Checks if the number n is B-smooth
// bool IsSmooth() {

// }
mpz_class FactorOut(mpz_class N, mpz_class p) {
    if (N % p == 0) {
        std::cout << "HAH" << std::endl;
    }
    while (N % p == 0) {
        N = N / p;
    }
    return N;
}

double Log(const mpz_class& a){
    int bits = mpz_sizeinbase(a.get_mpz_t(), 2);
    mpf_class factor = 1;
    factor <<= bits; // 2^bits
    mpf_class t = (a / factor);
    return bits + std::log2(t.get_d());
}

mpz_class GetS(const mpz_class& N, mpz_class p) {
    // TODO Actual tonelli shanlk
    mpz_class s = 0;
    while (true) {
        ++s;
        if ((s * s % p) == (N % p)) {
            return s;
        }
    }
}


class QuadraticSieve {

    private:
    const mpz_class& N;
    const long B;  // Let B be a bound on primes

    std::vector<long> primes;

    mpz_class Legendre(long p) {
        return ModularExponentiation(N, (p-1)/2, p);
    }

    mpz_class TonelliShank(mpz_class p) {
        // TODO Actual tonelli shanlk
        mpz_class s = 0;
        while (true) {
            ++s;
            if (((s * s) - N) % p == 0) { // (N % p)
                return s;
            }
        }
    }

    long GetSmoothnessBound() {
        double l = Log(N) * Log(Log(N));
        return std::exp(0.5 * std::sqrt(l));
    }

    mpz_class SqrtCeil(const mpz_class& n) {
        mpf_class tmp(n);
        tmp = sqrt(tmp);
        tmp = ceil(tmp);
        mpz_class tmp2(tmp);
        return tmp2;
    }

    mpz_class Mod(const mpz_class& x, long p) {
        return ((x % p) + p ) % p;
    }

    public:

    void GetFactors() {
        // http://www.cs.virginia.edu/crab/QFS_Simple.pdf
        // Retrieve factor base
        std::vector<long> base;
        //base.push_back(-1);
        for (long prime : SieveOfEratosthenes(B)) {
            //if (prime == 2) { 
            //    // Ignore the "oddest prime"
            //    continue;
            //}
            if (Legendre(prime) == 1) {
                base.push_back(prime);
            }
        }

        size_t M = 100;
        std::vector<mpz_class> Q(M);

        mpz_class a = SqrtCeil(N); // √N

        // q(x) = (X - √N)² - N
        auto q = [&] (const mpz_class& i) -> mpz_class { return ((i + a) * (i + a)) - N; };

        for (size_t i = 0; i < Q.size(); ++i) {
           Q[i] = q(i);
        }


        for (long p : base) {
            double logP = std::log2(p);
            mpz_class x0 = TonelliShank(p);
            mpz_class x1 = p - x0;

            mpz_class x[2] = {Mod(x0 - a, p), Mod(x1 - a, p)};
            for (int j = 0; j < 2; ++j) {
                for (size_t index = x[j].get_si(); index < M; index += p) {
                     while (Q[index] % p == 0) {
                        Q[index] /= p;
                    }
                }
                if (p == 2) break;
            }
        }

        std::vector<mpz_class> smoothNumbers;
        double threshold = std::log2(base.back());
        for (size_t i = 0; i < Q.size(); ++i) {
            if (abs(Q[i]) < threshold) {
                pr(Q[i]);
                mpz_class b = q(i);
                pr(b);
                for (long p : base) {
                    while (b % p == 0) {
                        b /= p;
                    }
                }
                if (b == 1) {
                    smoothNumbers.push_back(q(i));
                }
            }
        }
        pr("");
        for(auto a : smoothNumbers)
            pr(a);

    }

    QuadraticSieve(const mpz_class& N) : N(N), B(GetSmoothnessBound()) {}
};