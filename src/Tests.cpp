#include <cxxtest/TestSuite.h>
#include <gmpxx.h>
#include "Functions.hpp"
#include "QuadraticSieve.hpp"

class TestSuite : public CxxTest::TestSuite {
  public:

    void TestRabinMiller() {
        gmp_randclass rand(gmp_randinit_default);

        int k = 25;
        mpz_class prime = 1;
        mpz_class next_prime = 1;

        for(int i = 0; i < 1000; ++i) {
            mpz_nextprime(next_prime.get_mpz_t(), prime.get_mpz_t());
            TSM_ASSERT_EQUALS(next_prime.get_str() + " was false", MillerRabinPrime(next_prime, k, rand), true);
            while (++prime < next_prime) {
                TSM_ASSERT_EQUALS(prime.get_str()  + " was true", MillerRabinPrime(prime, k, rand), false);
            }
            prime = next_prime; 
        }
    }

    void TestPrimeGeneration() {
        gmp_randclass rand(gmp_randinit_default);
        std::vector<long> primes = SieveOfEratosthenes(1000);
        for(size_t i = 0; i < primes.size(); ++i) {
            TS_ASSERT_EQUALS(MillerRabinPrime(primes[i], 25, rand), true);
        }
    }

    void TestQuadraticSieve() {
        mpz_class test;
        // test = "1212121234000000000000000000000000000000000000000000000000000000000001";
        test = "24961";
        test = "15347";
        QuadraticSieve q(test);
        q.GetFactors();
    }
};
