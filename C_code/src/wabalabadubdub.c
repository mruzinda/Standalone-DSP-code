#include <stdio.h>
#include <math.h>

int main()
{
    printf("Wabalabadubduuuuub! I turned myself into a pickle, Morty!\n");

    /*
    // Mod test
    float a, b;
    int c;
    a = 3;
    b = 26;
    c = 24;
    printf("Remainder of %f / %d is %lf\n", a, c, fmod(a,c));
    printf("Remainder of %f / %d is %lf\n", b, c, fmod(b,c));
    */
    // Bitwise shift
    int a, b; 
    float c;
    a = 1;
    b = 17;
    c = (a<<b);
    printf("Bitwise shift %d right by %d or %d<<%d = %lf", a,b,a,b,c);
    
    return 0;
}