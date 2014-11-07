#!/usr/bin/ruby
# encoding: UTF-8
large_integer = "121212123400000000000000000000000000000000000000000000000000000000000000000000000000000001"
large_integer = ARGV[0] unless ARGV[0].nil?
fail unless system("make")
puts "compiled"
system("./bin/Main.out #{large_integer}")
