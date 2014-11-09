#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <queue>
#include <gmpxx.h>
#include <ctime>

#include <map>
#include <sstream>

#include <thread>
#include <chrono>
#include <atomic>

#include "Functions.hpp"

std::vector<mpz_class> Factorize(const mpz_class& number, const std::atomic<bool>& interrupt) {
    // Create random generator
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));

    std::queue<mpz_class> queue;
    std::vector<mpz_class> factors;

    queue.push(number);

    while (!queue.empty() && !interrupt) {
        mpz_class factor = queue.front();
        queue.pop();
        if (factor == 1) {
            continue;
        }
        if (MillerRabinPrime(factor, 25, rand)) { // mpz_probab_prime_p(factor.get_mpz_t(), 25)
            std::cout << "  - " << factor << std::endl;
            factors.push_back(factor);
            continue;   
        }
        mpz_class divisor = BrentsRhoAlgorithm(factor, rand, interrupt);
        queue.push(divisor);
        queue.push(factor / divisor);
    }
    return factors;
}

std::string Output(std::vector<mpz_class>& factors) {
    std::sort(factors.begin(), factors.end());
    std::map <mpz_class, int> m;
    auto inc = [&m] (mpz_class value) { ++m[value]; }; 
    std::for_each(factors.begin(), factors.end(), inc);
    std::stringstream s;
    for (auto f : m) {
        // Prime and occurrence
        s << "  - - " << f.first << std::endl << "    - " << f.second << std::endl;
    }
    return s.str();
}

std::vector<mpz_class> ReadNumbers(int argc, char** argv) {
    std::vector<mpz_class> numbers;
    if (argc >= 3) {
        mpz_class num;
        num = argv[2];
        numbers.push_back(num);
    } else {
        for (std::string line; std::getline(std::cin, line);) {
            numbers.push_back(mpz_class(line));
        }
    }
    return numbers;
}

int main(int argc, char **argv) {
    const int sleep_for = 10;
    int timeout = std::atoi(argv[1]);
    std::vector<mpz_class> numbers = ReadNumbers(argc, argv);

    std::cout << "---" << std::endl;
    std::cout << "- :timeout: " << timeout        << std::endl;

    //std::time_t end_time = Clock::to_time_t(end);
    //std::cout << std::ctime(&end_time);

    typedef std::chrono::system_clock Clock;

    for (size_t i = 0; i < numbers.size(); ++i) {
        
        std::atomic<bool> done(false);
        std::atomic<bool> running(false);
        std::atomic<bool> interrupt(false);

        Clock::system_clock::time_point thread_start;

        std::thread t = std::thread([&] {
            Clock::system_clock::time_point end;
            thread_start = Clock::now();

            running = true;

            std::cout << "- :index: "   << (i+1)      << std::endl;
            std::cout << "  :number: "  << numbers[i] << std::endl;
            std::cout << "  :factors: " << std::endl;

            std::vector<mpz_class> factors = Factorize(numbers[i], interrupt);
            if (interrupt) { return; }

            end = Clock::now();
            std::chrono::duration<double> elapsed_seconds  = end - thread_start;
            
            std::cout << "  :output: " << std::endl;
            std::cout << Output(factors);       
            std::cout << "  :time_elapsed: " << elapsed_seconds.count() << std::endl;
            done = true;
        });


        Clock::system_clock::time_point now;
        while (!done) {
            if(!running) { continue; }

            now = Clock::now();
            std::chrono::duration<double> elapsed_seconds = now - thread_start;

            if (elapsed_seconds.count() > timeout) {
                std::cout << "  :timed_out: true" << std::endl;
                interrupt = true;
                break;
            } else {
                std::this_thread::sleep_for(std::chrono::seconds(sleep_for));
            }
        }
        t.join();
    }
    return 0;
}
