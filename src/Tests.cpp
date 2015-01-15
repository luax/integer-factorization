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
        test = "9855429793132908353190283108657761069678331880366954624985078217911109";
        test = "522622270112106451781453377997875144414706787346078596434190435827619639239";
        test = "32407920105165278953559997155995882803747990815437685798227121";
        //test = "3101584057286468";
        double T = 0.5;
        long M = 1000000;
        QuadraticSieve q(test, M, T);
        std::cout << q.GetFactor(false) << std::endl;     
    }

};
