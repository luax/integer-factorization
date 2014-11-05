#include <iostream>
#include <gmp.h>

int main(int argc, char **argv) 
{
    if(argc != 2) { return 1; }
    mpz_t large_integer;
    mpz_init_set_str(large_integer, argv[1], 10);
    gmp_printf ("%Zd\n", large_integer);
    return 0;
}
