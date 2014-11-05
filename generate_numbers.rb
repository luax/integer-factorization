#!/usr/bin/ruby
# encoding: UTF-8
j = ARGV[1].to_i
puts (1..100).map {|i| ARGV[0].to_i * 10 ** (60 + j) + i }
