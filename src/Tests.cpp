#include <cxxtest/TestSuite.h>
#include <gmpxx.h>
#include "Functions.hpp"

class TestSuite : public CxxTest::TestSuite {
  public:

    void TestRabinMiller() {
      gmp_randclass rand(gmp_randinit_default);
      int k = 25;
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(5), k, rand), true);
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(7), k, rand), true);
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(11), k, rand), true);
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(169743212279), k, rand), true);
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(169743212280), k, rand), false);
      TS_ASSERT_EQUALS(MillerRabinPrime(mpz_class(169743212280007), k, rand), true);
    }

};
