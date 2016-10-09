#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "alloc.h"

// Prototypes of functions called by main
void read_genome(char **genome, uint64_t *glen);

char *get_arg_value(char *val)
{
  return("/Users/gibrat/MEHDI/GEN_3RD_READS/human_chromosomeIV.fa");
}


int main(int argc, char *argv[])
{

  FILE *fp;
  char *genome;
  uint64_t glen;
  char buffer[BUFSIZE];
  char dum[100];
  int lim1=0, lim2=0;
  int nr;
  int nerr = 0;

  /*
   * Read the genome from which the reads will be extracted
   */
  read_genome(&genome,&glen);

#ifdef DEBUG
  fprintf(stdout,"glen=%llu\n",glen);
  fprintf(stdout,"%s\n",genome);
#endif
  FOPEN(fp,argv[1],"r");

  nr = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '>') {
      sscanf(buffer,"%s %s %s %d %s %d",dum,dum,dum,&lim1,dum,&lim2);
      nr++;
    } else {
      if(strncmp(buffer,genome+lim1,lim2-lim1) != 0) {
	fprintf(stdout,"Read %d does not match the genome\n",nr);
	fprintf(stdout,"%s\n",genome+lim1);
	fprintf(stdout,"%s",buffer);
	nerr++;
      }
    }
  }
  if(nerr == 0) {
    fprintf(stdout,"No error\n");
  }

  FCLOSE(fp);
  FREE(genome);

  return(0);
}
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
int return_indx(double *cumproba, int Nbins)
{
  int i;
  double rdm;

  rdm = drand48();
  for(i = 1; i <= Nbins; i++) {
    if(rdm < cumproba[i]) {
      return(i);
    }
  }
  if(i > Nbins) {
    fprintf(stderr,"You should never see this message. Something went seriously wrong\n");
    exit(1);
  }

  return(-1);
}
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
int return_read_length(int lim1, int lim2)
{
  double rdm;
  int len;

  rdm = drand48();
  len = lim1 + floor((lim2 - lim1 + 1) * rdm);

  return(len);
}
 
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
uint64_t return_start_pos(uint64_t glen)
{
  double rdm;
  uint64_t pos;

  rdm = drand48();
  pos = floor(glen * rdm);

  return(pos);
}
 
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
void read_genome(char **genome, uint64_t *glen)
{
  /*
   * Read a fasta file in mode multi lines or single line
   */
  FILE *fp;
  uint32_t ll, lg;
  char buffer[BUFSIZE];

  FOPEN(fp,get_arg_value("genome_seq"),"r");
  
  *glen = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] != '>') {
      ll = strlen(buffer);
      if(ll == BUFSIZE-1 && buffer[ll-1] != '\n') { // When the buffer size is smaller than the length of the line to be read
	*glen += BUFSIZE - 1;                       // fgets reads BUFSIZE-1 characters and adds '\0' at the end of the buffer 
      } else {                
	*glen += ll-1;
      }
    } 
  }

  ALLOC(*genome, *glen + 1);

  rewind(fp);

  lg = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] != '>') {
      ll = strlen(buffer);
      if(ll != BUFSIZE-1 || buffer[ll-1] == '\n') {
	ll--;
	buffer[ll] = '\0';
      }
      strcpy(*genome+lg,buffer);
      lg += ll;
    } 
  }
  (*genome)[lg] = '\0';

  FCLOSE(fp);

}
