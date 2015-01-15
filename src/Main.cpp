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
#include "QuadraticSieve.hpp"

void FactorFound(std::vector<mpz_class>& factors, const mpz_class& factor, mpz_class& current_number) {
    factors.push_back(factor);
    current_number /= factor;
    std::clog << "LOG: found factor -> "  << factor << std::endl;
    std::clog << "LOG: current number: " << current_number << std::endl;
    std::cout << "  - " << factor << std::endl;
} 

void TrialDivision(std::vector<mpz_class>& factors, mpz_class& n, const std::atomic<bool>& interrupt) {
    const static std::vector<long> primes = SieveOfEratosthenes(1000000000);
    for (size_t i = 0; i < primes.size(); ++i) {
        if (primes[i] * primes[i] > n) break;
        while (n % primes[i] == 0) {
            FactorFound(factors, primes[i], n); // TODO
        }
    }
}

std::vector<mpz_class> Factorize(const mpz_class& number, const std::atomic<bool>& interrupt) {
    // Create random generator
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));

    std::queue<mpz_class> queue;
    std::vector<mpz_class> factors;

    mpz_class current_number = number;

    std::clog << "LOG: Doing trial divison " << std::endl;

    TrialDivision(factors, current_number, interrupt);

    queue.push(current_number);

    std::clog << "LOG: Doing Pollard's rho " << std::endl;

    // Clock stuff
    typedef std::chrono::system_clock Clock;
    Clock::system_clock::time_point thread_start;
    Clock::system_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    std::atomic<bool> done(false);
    std::atomic<bool> running(false);
    std::atomic<bool> interrupt2(false);

    std::thread t = std::thread([&] {
        thread_start = Clock::now();
        running = true;
        while (!queue.empty() && !interrupt2) {
            mpz_class factor = queue.front();
            queue.pop();

            if (factor == 1) {
                continue;
            }

            if (MillerRabinPrime(factor, 25, rand)) {
                FactorFound(factors, factor, current_number);
                continue;   
            }
            mpz_class divisor = BrentsRhoAlgorithm(factor, rand, interrupt2);
            queue.push(divisor);
            queue.push(factor / divisor);
        }
        done = true;
    });

    Clock::system_clock::time_point now;
    while (!done) {
        if(!running) { continue; }
        now = Clock::now();
        std::chrono::duration<double> elapsed_seconds = now - thread_start;
        if (elapsed_seconds.count() > 3000) {
            interrupt2 = true;
            std::clog << "LOG: Pollard's timed out" << std::endl;
            break;
        } else {
            std::this_thread::sleep_for(std::chrono::seconds(100));
        }
    }
    
    t.join();

    std::clog << "LOG: Doing Quadratic sieve " << std::endl;

    while (!queue.empty()) {
        mpz_class factor = queue.front();
        queue.pop();

        if (factor == 1) {
            continue;
        }

        if (factor != -1 && MillerRabinPrime(factor, 25, rand)) {
            FactorFound(factors, factor, current_number);
            continue;   
        }

        if (interrupt) {
            break;
        }
        QuadraticSieve q(factor, 1000000, 0.5);
        mpz_class divisor = q.GetFactor(interrupt);
 
        if(divisor == -1){
            std::clog << "DID NOT FIND ANY FACTORS!" << std::endl;
            std::clog << "Moving on..." << std::endl;
        }else{
            queue.push(divisor);
            queue.push(factor / divisor);
        }
    }
    return factors;
}

std::string OutputFormat(std::vector<mpz_class>& factors) {
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
    for (std::string line; std::getline(std::cin, line);) {
        numbers.push_back(mpz_class(line));
    }
    return numbers;
}

int main(int argc, char **argv) {
    int timeout = std::atoi(argv[1]);
    int sleep_for = timeout > 10 ? 10 : timeout;

    std::vector<mpz_class> numbers = ReadNumbers(argc, argv);

    std::cout << "---" << std::endl;
    std::cout << "- :timeout: " << timeout << std::endl;

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

            std::clog << "LOG: factoring "  << numbers[i] << std::endl;
            std::cout << "- :number: "  << numbers[i] << std::endl;
            std::cout << "  :factors: " << std::endl;

            std::vector<mpz_class> factors = Factorize(numbers[i], interrupt);
            if (interrupt) { return; }

            end = Clock::now();
            std::chrono::duration<double> elapsed_seconds  = end - thread_start;
            
            std::cout << "  :output: " << std::endl;
            std::cout << OutputFormat(factors);       
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
                std::clog << "LOG: timed out" << std::endl;
                interrupt = true;
                break;
            } else {
                std::this_thread::sleep_for(std::chrono::seconds(sleep_for));
            }
        }
        std::clog << "LOG: waiting..." << std::endl;
        t.join();
        std::clog << "LOG: joined" << std::endl;
    }

    std::clog << "LOG: Done!" << std::endl;
    return 0;
}
