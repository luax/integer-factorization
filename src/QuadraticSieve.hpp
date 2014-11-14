#include "Functions.hpp"
#include <vector>
#include <cmath>
#include <cassert>

#define pr(a) std::cout << a << std::endl

class QuadraticSieve {

    public:
    static const size_t M = 10000000; // Interval size
    const mpz_class& N;
    const long B;  // Let B be a bound on primes

    std::vector<long> primes;

    mpz_class Legendre(const mpz_class& n, long p) {
        // ändrad
        if (n % p == 0) return 0;
        mpz_class i = ModularExponentiation(n, (p-1)/2, p);
        if(i == 1 ) return 1;
        return -1;
    }

    mpz_class TonelliShank(long p) {
        mpz_class s = 0;
        mpz_class q = p - 1;
        while ((q & 1) == 0) {
            q >>= 1;
            ++s; 
        }

        if (s == 1)
            return ModularExponentiation(N, (p+1)/4, p);

        mpz_class z = 2;
        while (ModularExponentiation(z, (p-1)/2, p) != p - 1)
            ++z;

        mpz_class c = ModularExponentiation(z, q, p);
        mpz_class r = ModularExponentiation(N, (q + 1) / 2, p);
        mpz_class t = ModularExponentiation(N, q, p);
        mpz_class m = s;
        while (t != 1) {
            mpz_class tmp = t;
            mpz_class i = 0;
            while (tmp != 1) {
              tmp = (tmp * tmp) % p;
              ++i;
            }
            mpz_class b = ModularExponentiation(c, ModularExponentiation(2, (m - i - 1), p - 1), p);
            mpz_class b2 = (b * b) % p;
            r = (r * b) % p;
            t = (t * b2) % p;
            c = b2;
            m = i;
        }
        return r;
    }

    double Log(const mpz_class& a){
        int bits = mpz_sizeinbase(a.get_mpz_t(), 2);
        mpf_class factor = 1;
        factor <<= bits; // 2^bits
        mpf_class t = (a / factor);
        return bits + std::log2(t.get_d());
    }

    long GetSmoothnessBound() {
        double l = Log(N) * Log(Log(N));
        double b = std::exp(0.5 * std::sqrt(l));
        //return b < 1000 ? 1000 : b;
        return 10000000;
    }

    mpz_class SqrtCeil(const mpz_class& n) {
        mpf_class tmp(n);
        tmp = sqrt(tmp);
        tmp = ceil(tmp);
        mpz_class tmp2(tmp);
        return tmp2;
    }

    long Mod(const mpz_class& x, long p) {
        mpz_class tmp = x % p;
        return (tmp.get_si() + p ) % p;
    }

    public:

    void GetFactors() {
        std::cout << "B: " << B << std::endl;
        // http://www.cs.virginia.edu/crab/QFS_Simple.pdf
        // Retrieve factor base
        std::vector<long> base;
        //base.push_back(-1);
        for (long prime : SieveOfEratosthenes(B)) {
            if (Legendre(N, prime) == 1) {
                base.push_back(prime);
            }
        }

        std::cout << "primes: " << base.back() << std::endl;
        
        mpz_class sqrtN = SqrtCeil(N); // √N

        // q(x) = (X - √N)² - N
        auto q = [&] (const mpz_class& i) -> mpz_class { return ((i + sqrtN) * (i + sqrtN)) - N; };


        struct PrimeData{
            long p;
            double log;
            long x[2];

            PrimeData(long p, double log, long x1, long x2) : p(p), log(log) {
                x[0] = x1;
                x[1] = x2;
            }
        };

        // precalculate prime data
        std::vector<PrimeData> data;
        //for(size_t t = 1; t < 4; ++t){
            for(long p : base){
                long pt = std::pow(p, 1);
                double logP = std::log2(p);
                mpz_class x0 = TonelliShank(pt);
                assert((x0*x0) % p == N % p);
                mpz_class x1 = pt - x0;
                data.push_back(PrimeData(pt, logP, Mod(x0 - sqrtN, pt), Mod(x1 - sqrtN, pt)));
            }
        //}

        std::cout << "prime data" << std::endl;

        int a = 0;
        int b = M;
        double threshold = std::log2(base.back());
        std::vector<double> Q(M);
        std::vector<mpz_class> smoothNumbers;
        while(smoothNumbers.size() < base.size() + 1){
            // init Q array
            for (size_t i = 0; i < Q.size(); ++i) {
               Q[i] = Log(abs(q(i + a)));
            }
            for (PrimeData& pd : data) {
                for (int j = 0; j < 2; ++j) {
                    while(pd.x[j] < b) {
                        Q[pd.x[j] - a] -= pd.log;
                        pd.x[j] += pd.p;
                    }
                    if (pd.p == 2) break;
                }
            }
            std::cout << "refine " << std::endl;

            for (size_t i = 0; i < Q.size(); ++i) {
                if (abs(Q[i]) < threshold) {
                    mpz_class val = q(i + a);
                    for (long p : base) {
                        while (val % p == 0) {
                            val /= p;
                        }
                    }
                    if (val == 1) {
                        smoothNumbers.push_back(q(i + a));
                    };
                }
            }
            pr("");
            pr("a: " << a);
            pr("b: " << b);
            pr("size: " << smoothNumbers.size() << " " << base.size() + 1);
            a = b;
            b += M;
        } 
    }

    QuadraticSieve(const mpz_class& N) : N(N), B(GetSmoothnessBound()) {}
};