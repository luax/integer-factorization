#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "Functions.hpp"

std::vector<mpz_class> Factorize(const mpz_class& number) {
    // Create random generator
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));

    std::vector<mpz_class> queue;
    std::vector<mpz_class> factors;

    queue.push_back(number);

    while (!queue.empty()) {
        mpz_class factor = queue.back();
        queue.pop_back();

        if (MillerRabinPrime(factor, 20, rand)) {
            std::cout << factor << " ";
            factors.push_back(factor);
            continue;
        }

        mpz_class divisor = PollardsRhoAlgorithm(factor, rand);
        queue.push_back(divisor);
        queue.push_back(factor / divisor);
    }
    std::cout << std::endl;

    return factors;
}

void PrintFactors(std::vector<mpz_class>& factors) {
    std::sort(factors.begin(), factors.end());
    for (auto f : factors) {
       std::cout << f << " ";
    }
    std::cout << std::endl;
}

std::vector<mpz_class> ReadInput(int argc, char** argv) {
    std::vector<mpz_class> numbers;
    if (argc >= 2) {
        mpz_class num;
        num = argv[1];
        numbers.push_back(num);
    } else {
        for (std::string line; std::getline(std::cin, line);) {
            numbers.push_back(mpz_class(line));
        }
    }
    return numbers;
}

int main(int argc, char **argv) {
    std::vector<mpz_class> numbers = ReadInput(argc, argv);
    for (size_t i = 0; i < numbers.size(); ++i) {
        std::cout << "--- " << (i+1) << ") Factoring number: " << numbers[i] << " ---" << std::endl;
        std::vector<mpz_class> factors = Factorize(numbers[i]);
        std::cout << "--- Done ---" << std::endl;
    }
    return 0;
}
