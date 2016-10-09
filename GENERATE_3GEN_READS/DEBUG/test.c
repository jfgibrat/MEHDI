#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "alloc_test.h"

char *get_arg_value(char *val)
{
  return("/Users/gibrat/MEHDI/GEN_3RD_READS/human_chromosomeIV.fa");
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
  (*genome)[lg+1] = '\n';

  FCLOSE(fp);

}

int main(int argc, char *argv[])
{
  char *genome;
  uint64_t glen;
  read_genome(&genome,&glen);
  fprintf(stdout,"glen=%llu\n",glen);
  fprintf(stdout,"%s\n",genome);
}
