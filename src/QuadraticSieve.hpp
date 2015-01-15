#include "Functions.hpp"
#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>
#include <queue>
#include <boost/dynamic_bitset.hpp>
#include <chrono>

#define pr(a) std::clog << a << std::endl
#define prn(a) std::clog << #a << ": " << a << std::endl;

#ifndef QUADRATICSIEVE_H
#define QUADRATICSIEVE_H

class QuadraticSieve {

    public:
    mpz_class N;
    long B;    // Let B be a bound on primes
    mpz_class start_x;
    const unsigned long interval; // Interval size
    const double T; // Interval size

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

    float Log2(const mpz_class& a){
        return mpz_sizeinbase(a.get_mpz_t(), 2);
    }

    long GetSmoothnessBound() {
        double l = Log(N) * Log(Log(N));
        double b = std::exp(0.5 * std::sqrt(l));
        return (b + 1000) > 1000000 ? 1000000 : b;
        // return 1000000;
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
    
    std::vector<long> GetPrimeBase() {
        std::vector<long> base;
        // base.push_back(-1);
        for (long prime : SieveOfEratosthenes(B)) {
            if (Legendre(N, prime) == 1) {
                base.push_back(prime);
            }
        }
        return base;
    }

    void FastGauss(boost::dynamic_bitset<>* M, unsigned long rows, unsigned long columns) {
        // Gauss
        // https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
        for (unsigned long i = 0; i < rows; ++i) {
           int p = -1;
           for (unsigned long j = 0; j < columns; ++j) {
               if (M[i][j] == 1) {
                   p = j;
                   break;
               }
           }
           if (p == -1) {
               continue;
           }
 
           for (unsigned long k = 0; k < rows; ++k) {
               if (k == i) {
                   continue;
               }
               if (M[k][p] == 1) {
                   M[k] ^= M[i];
               }
           }
        }
    }

    mpz_class GetFactor(const std::atomic<bool>& interrupt) {
        B = GetSmoothnessBound();
        pr("Getting prime base");
        std::vector<long> base = GetPrimeBase();
        unsigned long needed_smooth_numbers = base.size() + 20;
        double threshold = Log2(base.back()); // T * std::log2(base.back());
        unsigned long rows = needed_smooth_numbers;
        unsigned long columns = base.size() + needed_smooth_numbers;
        pr("Threshold: " << threshold);
        pr("Needed smooth numbers: " << needed_smooth_numbers);

        pr("Creating matrix");
        boost::dynamic_bitset<>* M = new boost::dynamic_bitset<>[rows];

        mpz_class sqrtN = SqrtCeil(N); 
        start_x = 0;

        // q(x) = (X - √N)² - N
        auto q = [&] (const mpz_class& i) -> mpz_class { return ((i + sqrtN) * (i + sqrtN)) - N; };

        pr("Calculating prime data");
        struct PrimeData{
            long p;
            float log;
            mpz_class x[2];

            PrimeData(long p, double log, mpz_class x1, mpz_class x2) : p(p), log(log) {
                x[0] = x1;
                x[1] = x2;
            }
        };

        mpz_class a = start_x;
        mpz_class b = a + interval;

        // precalculate prime data
        std::vector<PrimeData> data;
        for(long p : base){
            float logP = Log2(p);
            mpz_class x0 = TonelliShank(p);
            mpz_class x1 = p - x0;
            data.push_back(PrimeData(p, logP, Mod(x0 - sqrtN, p) + a, Mod(x1 - sqrtN, p) + a));
        }

        pr("Finding smooth numbers");

        std::vector<float> Q(interval);
        std::vector<std::tuple<mpz_class, mpz_class>> smooth_numbers;
        smooth_numbers.reserve(needed_smooth_numbers);

        // Clock stuff
        typedef std::chrono::system_clock Clock;
        Clock::system_clock::time_point start = Clock::now();
        Clock::system_clock::time_point end;
        std::chrono::duration<double> elapsed_seconds;
        int counter = 0;
        int zero = 0;
        while(smooth_numbers.size() < needed_smooth_numbers && !interrupt){
            // init Q array
            for (unsigned long i = 0; i < Q.size(); ++i) {
               Q[i] = Log2(abs(q(i + a)));
            }
            for (PrimeData& pd : data) {
                for (int j = 0; j < 2; ++j) {
                    while(pd.x[j] < b) {
                        mpz_class index = pd.x[j] - a;
                        Q[index.get_si()] -= pd.log;
                        pd.x[j] += pd.p;
                    }
                    if (pd.p == 2) break; // Only one solution for 2
                }
            }

            for (unsigned long i = 0; i < Q.size(); ++i) {
                if (Q[i] < threshold) {
                    mpz_class y = q(i + a);
                    mpz_class y_copy = y;
                    boost::dynamic_bitset<> matrix_row(columns);
                    for(unsigned long j = 0; j < base.size(); ++j){
                        while(y_copy % base[j] == 0){
                            y_copy /= base[j];
                            matrix_row.flip(j);
                        }   
                    }
                    if (y_copy == 1) {
                        std::tuple<mpz_class, mpz_class> t(y, (i + a + sqrtN));
                        smooth_numbers.push_back(t);
                        unsigned long n = smooth_numbers.size() - 1;
                        M[n] = matrix_row;
                        M[n][n + base.size()] = 1;
                        if(smooth_numbers.size() >= needed_smooth_numbers)
                            break;
                    };
                }
            }

            end = Clock::now();
            elapsed_seconds  = end - start;
            float sns = ((float)smooth_numbers.size() / elapsed_seconds.count());
            if (++counter % 5 == 0) {
                pr("LOG: " << N << " " << sns << " sn/s, size: " << smooth_numbers.size() << " / " << needed_smooth_numbers << " interval: " << a << " - " << b);
            }
            if (sns < 0.5) {
                ++zero;
                pr("LOG: too slow " << zero);
            } else {
                if (zero > 0) {
                    --zero;
                }
            }
            if (zero > 100) {
                break;
            }
            a = b;
            b += interval;
        }

        end = Clock::now();
        elapsed_seconds = end-start;

        if (interrupt) {
            pr("Sieving was interrupted after " << elapsed_seconds.count() << " seconds");
        } else if (zero > 100) {
            pr("Sieving was too bad after " << elapsed_seconds.count() << " seconds");
        } else {
            pr("Sieving took " << elapsed_seconds.count() << " seconds");
        }

        pr("");
        start = Clock::now();
        pr("Gaussing");
        FastGauss(M, smooth_numbers.size(), base.size());
        end = Clock::now();
        elapsed_seconds = end-start;
        pr("Gaussing took " << elapsed_seconds.count() << " seconds");

        pr("Looking for linear dependency");
        mpz_class factor = -1;
        for(unsigned long i = 0; i < smooth_numbers.size(); ++i){
            bool is_zero = true;
            for(unsigned long j = 0; j < base.size(); ++j){
                if(M[i][j] == 1){
                    is_zero = false;
                    break;
                }
            }
            if(is_zero){
                std::vector<unsigned long> linear_dependency;
                for(unsigned long j = base.size(); j < columns; ++j){
                    if(M[i][j] == 1){
                        linear_dependency.push_back(j - base.size());
                    }
                }

                pr("Found linear dependency: " << linear_dependency.size());
                if (linear_dependency.size() < 2)
                    continue;
                mpz_class sum_y = 1;
                mpz_class sum_x = 1;
                for(unsigned long j : linear_dependency) {
                    sum_y = (sum_y * std::get<0>(smooth_numbers[j]));
                    sum_x = (sum_x * std::get<1>(smooth_numbers[j]));
                }
                sum_y = sqrt(sum_y);
                sum_y %= N;
                sum_x %= N;
                mpz_class f = BinaryGCD(sum_x - sum_y, N);
                if(f != 1 && f != N){
                    factor = f;
                    break;
                }
            }
        }
        delete[] M;
        return factor;
    }

    QuadraticSieve(const mpz_class& N, unsigned long interval, double T) : N(N), interval(interval), T(T) {}
};

#endif /* QUADRATICSIEVE_H */
