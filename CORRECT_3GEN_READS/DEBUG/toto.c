#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

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
  char str[100];
  int hts;
  int kmer;
  int hh;
  int power4[100];
  int base;

  strcpy(str,argv[1]);
  hts = atoi(argv[2]);
  kmer = strlen(str);

  power4[0] = 1;
  for(i = 1; i < kmer; i++) {
    power4[i] = power4[i-1] * 4;
  }
  hh = 0;
  for(i = 0; i < kmer; i++) {
    base = nucval(str[i]);
    hh += base * power4[i];
  }
  fprintf(stdout,"hh=%d dividend=%d remainder=%d\n",hh,hh/hts,hh%hts);
  
  return 0;
}

