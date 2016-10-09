#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "SOURCE/alloc.h"

int chomp(char *str);

int nucval(char cc)
{
  switch(cc)
    {
    case 'A':
      return(0);
    case 'C':
      return(1);
    case 'G':
      return(2);
    case 'T':
      return(3);
    }
  return(-1);
}

int main(int argc, char *argv[])
{

  int i;
  unsigned long power4[30];
  unsigned long h1, d1, r1;
  int kmer_len;

  power4[0] = 1;
  for(i = 1; i < 30; i++) {
    power4[i] = power4[i-1] * 4;
  }

  kmer_len = strlen(argv[1]);
  h1 = 0;
  for(i = 0; i < kmer_len; i++) {
    h1 += power4[i] * nucval(argv[1][i]);
  }
  d1 = h1 / atoi(argv[2]);
  r1 = h1 % atoi(argv[2]);

  fprintf(stdout,"hh=%lu r1=%lu d1=%lu\n",h1,r1,d1);
  return 0;

}

/*******************************************************************************/
/*                                                                             */
/*                               function chomp                                */
/*                                                                             */
/*******************************************************************************/
int chomp(char *str)
{

/*
 * Get rid of the last '\n' character and return the new string length
 */

  int ll;

  ll = strlen(str);

  if(str[ll-1] != '\n') {
    fprintf(stdout,"Function chomp warning: character string (%s) does not end with a newline character\n\n",str);
    return(ll);
  } else {
    str[ll-1] = '\0';
    return(ll-1);
  }

}
