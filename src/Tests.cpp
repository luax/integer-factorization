#include <cxxtest/TestSuite.h>
#include <gmpxx.h>
#include "Functions.hpp"

class TestSuite : public CxxTest::TestSuite {
  public:

    void TestRabinMiller() {
        gmp_randclass rand(gmp_randinit_default);

        int k = 25;
        mpz_class prime = 1;
        mpz_class next_prime = 1;

        for(int i = 0; i < 10000; ++i) {
            mpz_nextprime(next_prime.get_mpz_t(), prime.get_mpz_t());
            TSM_ASSERT_EQUALS(next_prime.get_str() + " was false", MillerRabinPrime(next_prime, k, rand), true);
            while (++prime < next_prime) {
                TSM_ASSERT_EQUALS(prime.get_str()  + " was true", MillerRabinPrime(prime, k, rand), false);
            }
            prime = next_prime; 
        }
    }
};
