#include <stdio.h>

/* This method returns a boolean value indicating wether the machine where
 * the code is running is Big Endian or not */
int am_big_endian (void) {
  long one = 1;
  return !(*((char *)(&one)));
}

/* This method swaps the bytes in a 32-bit word, changing it's effective
 * endianness between Big Endian and Little Endian */
unsigned int swap_endianness (unsigned int a) {
	unsigned int a1, a2, a3, a4, b;

	a1 = a & 0xFF000000;
	a2 = a & 0x00FF0000;
	a3 = a & 0x0000FF00;
	a4 = a & 0x000000FF;

	a1 = a1 >> (8+8+8);
	a2 = a2 >> 8;
	a3 = a3 << 8;
	a4 = a4 << (8+8+8);

	b = a1 | a2 | a3 | a4;
	
	return b;
}

/* Test program to check for the correctness of the previous routines */
int main (void) {
  unsigned int f2;
  int a;
  float f;
  
  if (am_big_endian()) {
    printf ("This is Big Endian\n");
  } else {
    printf ("This is Little Endian\n");
  }

  a = 0x12345678;
  printf ("Original: 0x%x, Intercanviat: 0x%x\n", a, swap_endianness (a));
  
  f = 0.0;
  f2 = 0;
  f = 0.1;
  f2 = *(unsigned int *)(void *)&f;
  printf ("Float: %f, 0x%x, Swapped: 0x%x\n", (double)f, f, f2);
  //memcpy (&f2, &f, sizeof (float));
  printf ("Float: %f, 0x%x, Swapped: 0x%x\n", (double)f, f, f2);
  f2 = swap_endianness (f2);
  printf ("Float: %f, ", f);
  printf ("%x, ", f);
  printf ("Swapped: %x\n", f2);

  return 0;
}
