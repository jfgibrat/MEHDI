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

  int nerr = 0;
  int ii = 0;
  int i, j;
  char buffer[1000];
  char read[1000];
  char kmer[1000];
  unsigned long power4[30];
  int idum, len, ll;
  char dum[100];
  int kmer_len;
  unsigned long hh, dd, rr;
  unsigned long h1, d1, r1;
  int HashMaxSize;
  FILE *fp;

  power4[0] = 1;
  for(i = 1; i < 30; i++) {
    power4[i] = power4[i-1] * 4;
  }

  FOPEN(fp,"/Users/gibrat/MEHDI/SOURCE/tutu","r");
  fgets(buffer,1000,fp);
  sscanf(buffer,"%d",&HashMaxSize);
  while(fgets(buffer,1000,fp) != NULL) {
    if(buffer[0] == '\n') {
      continue;
    }
    if(buffer[0] == '>') {
      sscanf(buffer+1,"%s %d %d",dum,&idum,&len);
      fgets(buffer,1000,fp);
      ll = chomp(buffer);
      if(ll != len) {
	fprintf(stderr,"Error: ll=%d and len= %d\n\n",ll,len);
	exit(1);
      }
      strcpy(read,buffer);
      ii = 0;
    } else {
      sscanf(buffer,"%s %ld %s %ld %s %ld %s",dum,&hh,dum,&dd,dum,&rr,kmer);
      kmer_len = strlen(kmer);
      for(i = 0, j = ii; i < kmer_len; i++, j++) {
	if(read[j] != kmer[i]) {
	  nerr++;
	  fprintf(stderr,"for kmer %d the nucleotide %d is different '%c' vs' %c'\n",ii+1,i,read[j],kmer[i]);
	}
      }
      h1 = 0;
      for(i = 0; i < kmer_len; i++) {
	h1 += power4[i] * nucval(kmer[i]);
      }
      d1 = h1 / HashMaxSize;
      r1 = h1 % HashMaxSize;
      if(hh != h1) {
	nerr++;
	fprintf(stderr,"for kmer %d hh=%lu not equal to h1=%lu\n",ii+1,hh,h1);
      }
      if(dd != d1) {
	nerr++;
	fprintf(stderr,"for kmer %d dd=%lu not equal to d1=%lu\n",ii+1,dd,d1);
      }
      if(rr != r1) {
	nerr++;
	fprintf(stderr,"for kmer %d rr=%lu not equal to r1=%lu\n",ii+1,rr,r1);
      }
      ii++;
    }
  }
  if(nerr == 0) {
    fprintf(stdout,"No error\n");
  } else {
    fprintf(stdout,"\n\nThere were %d errors\n",nerr);
  }
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
