/*
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * You can contact the main authors via email at:
 * Jean-François Gibrat (jean-francois ! gibrat () inra ! fr)
 * Mathématiques et Informatique Appliquées du Génome à l'Environnement (MaIAGE)
 * INRA - Domaine de Vilvert
 * 78350 Jouy-en-Josas cedex
 * France
 * Copyright (C) 2016
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "alloc.h"

// Prototypes of functions called by main
void process_line_args(int argc, char *argv[]);
char *get_arg_value(char *arg_name);
int return_indx(double *cumproba, int Nbins);
int return_read_length(int lim1, int lim2);
void read_genome(char **genome, uint64_t *glen);
uint64_t return_start_pos(uint64_t glen);
void read_error_rates(float *err_rates);
void generate_errors(char *buffer, float *cum_err_rates);

int main(int argc, char *argv[])
{

  FILE *fp;
  int ii, jj;
  int Nbins;
  char buffer[BUFSIZE];
  int *Nreads;
  int **bins = NULL;
  int sum;
  float *proba;
  double *cumproba;
  int Ngener;
  int *gener_bins;
  long int seedval;
  int len;
  char *genome;
  uint64_t glen;
  uint64_t pos;
  float err_rates[3];   // First position is mismatch error rate then insertion then deletion
  float cum_err_rates[4];

  process_line_args(argc,argv);

  /*
   * Initialize the random number generator
   */
  seedval = atol(get_arg_value("SeedVal"));
  srand48(seedval);

  /*
   *  Reads 3rd generation sequencing technology error rates.
   */
  read_error_rates(err_rates);

  cum_err_rates[0] = err_rates[0];  
  for(ii = 1; ii < 3; ii++) {
    cum_err_rates[ii] = cum_err_rates[ii-1] + err_rates[ii];
  } 
  cum_err_rates[ii] = 1.0;

  /*
   * Read the empirical PacBio distribution of reads
   */  
  FOPEN(fp,get_arg_value("PacBioDistriFile"),"r");

  Nbins = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    Nbins++;
  }

  rewind(fp);

  ALLOC(Nreads,Nbins);
  CALLOC2(bins,Nbins,2);

  ii = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    sscanf(buffer,"%d %d %d",bins[ii],bins[ii]+1,Nreads+ii);
    ii++;
  }

  FCLOSE(fp);

  /*
   * Compute the cumulative distribution of probability
   */
  sum = 0;
  for(ii = 0; ii < Nbins; ii++) {
    sum += Nreads[ii];
  }

  ALLOC(proba,Nbins);
  for(ii = 0; ii < Nbins; ii++) {
    proba[ii] = (float) Nreads[ii] / (float) sum;
  }

  ALLOC(cumproba, Nbins+1);
  cumproba[0] = 0.0;
  for(ii = 1; ii <= Nbins; ii++) {
    cumproba[ii] = cumproba[ii-1] + proba[ii-1];
  }
  cumproba[Nbins] = 1.0;

#ifdef DEBUG
  for(ii = 1; ii <= Nbins; ii++) {
    fprintf(stdout,"ii=%d %.7f\n",ii,cumproba[ii]);
  }
#endif

  /*
   * Generate the read lengths according to the empirical PacBio distribution
   */
  ALLOC(gener_bins, Nbins);
  for(ii = 0; ii < Nbins; ii++) {
    gener_bins[ii] = 0;
  }

  Ngener = atoi(get_arg_value("NreadsToGenerate"));
  for(ii = 0; ii < Ngener; ii++) {
    gener_bins[return_indx(cumproba,Nbins)-1] += 1;
  }

  for(ii = 0; ii < Nbins; ii++) {
    fprintf(stdout,"In interval [%5d - %5d] number of reads= %4d\t diff= %d\n",bins[ii][0],bins[ii][1],gener_bins[ii],gener_bins[ii]-Nreads[ii]);
  }

  /*
   * Read the genome from which the reads will be extracted
   */
  read_genome(&genome,&glen);

  fprintf(stdout,"glen=%llu\n",glen);
#ifdef DEBUG
  fprintf(stdout,"%s\n",genome);
#endif

  /*
   * Write the extracted reads to a file
   */
  if(BUFSIZE < bins[Nbins-1][1]) { // Make sure array buffer will not overflow below
    fprintf(stderr,"You should increase the value of BUFSIZE in alloc.h and recompile the program. The largest sequence is %d and BUFSIZE is %d\n",bins[ii-1][1], BUFSIZE);
    exit(1);
  }

  FOPEN(fp,get_arg_value("ouputFile"),"w");

  int kk = 0;
  for(ii = 0; ii < Nbins; ii++) {
#ifdef DEBUG
    fprintf(stdout,"In interval [%5d - %5d] %d reads are generated\n",bins[ii][0],bins[ii][1],gener_bins[ii]);
#endif
    jj = 1;
    while(jj <= gener_bins[ii]) {
      len = return_read_length(bins[ii][0],bins[ii][1]);
      pos = return_start_pos(glen);
#ifdef DEBUG
      fprintf(stdout,"\ttrying to generate a read of length %d between [%llu %llu] ... ",len,pos,pos+len-1);
#endif
      if(pos+len-1 < glen) {  // Make sure we do not go beyond the genome array limit
	kk++;
	strncpy(buffer,genome+pos,len);
	buffer[len] = '\0';
	generate_errors(buffer,cum_err_rates);
#ifdef KLUDGE
	if(strlen(buffer) > 65535) {
	  fprintf(stdout,"Read %d has length %lu that is greater than 65535. Skipping it...\n",kk,strlen(buffer));
	  jj++;
	  continue;    // Kludge to avoid generating reads length greater than 65535 that would not fit in uint16_t variable
	}
#endif
	fprintf(fp,">Read %d [ %llu - %llu ]\n",kk,pos,pos+len);
	fprintf(fp,"%s\n",buffer);
	jj++;
#ifdef DEBUG
	fprintf(stdout,"OK\n");
      } else {
	fprintf(stdout,"Failed\n");
#endif
      }
    }
  }

  FCLOSE(fp);
  FREE(genome);
  FREE(Nreads);
  FREE(proba);
  FREE(cumproba);
  FREE(gener_bins);
  FREE2(bins);

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
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
void read_error_rates(float *err_rates)
{
  /*
   * Reads 3rd generation sequencing technology error rates.
   * There are 3 types of errors : mismatch, insertion and deletion
   * Error rates must appear in THIS ORDER in the input file! 
   */
  FILE *fp;
  char buffer[BUFSIZE];
  char err_types[3][10] = {"mismatch ","insertion","deletion "};
  int ii = 0;
  
  FOPEN(fp,get_arg_value("read_err_rates"),"r");
  
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '#') continue;
    sscanf(buffer,"%f",err_rates+ii);
    if(err_rates[ii] > 1.0 || err_rates[ii] < 0.0) {
      fprintf(stderr,"Error in read_error_rates: the values of the error rates must lie between 0.0 and 1.0 (%f)\n",err_rates[ii]);
      exit(1);
    }
    ii++;
  }
  
  for(ii = 0; ii < 3; ii++) {
    fprintf(stdout,"%s error rate= %.4f\n",err_types[ii],err_rates[ii]);
  }

  FCLOSE(fp);

}
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
char new_nucleotide(char nuc);

void generate_errors(char *buffer, float *cum_err_rates)
{
  double rdm;
  int ii, jj;
  char tmp[BUFSIZE];

  for(ii = 0, jj = 0; buffer[ii] != '\0'; ii++) {
    rdm = drand48();
    if(rdm < cum_err_rates[0]) {        // mismatch
      tmp[jj] = new_nucleotide(buffer[ii]);
      jj++;
    } else if(rdm < cum_err_rates[1]) { // insertion
      tmp[jj] = new_nucleotide('I');
      jj++;
      tmp[jj] = buffer[ii];
      jj++;
    } else if(rdm < cum_err_rates[2]) { // deletion
      // skip the current nuclotide, i.e., do nothing with tmp
    } else {                            // no modification of the nucleotide
      tmp[jj] = buffer[ii];
      jj++;
    }
  }
  tmp[jj] = '\0';
  strcpy(buffer,tmp);

}
/****************************************************************
 *                                                              *
 *                                                              *
 ****************************************************************/
char new_nucleotide(char nuc)
{
  float cum_proba[] = {0.25,0.50,0.75,1.0};
  double rdm;
  
  if(nuc == 'I') {
    rdm = drand48();
    if(rdm < cum_proba[0]) {
      return('A');
    } else if(rdm < cum_proba[1]) {
      return('C');
    } else if(rdm < cum_proba[2]) {
      return('G');
    } else {
      return('T');
    }
  } else {
    while(1) { // For a mismatch the letter returned must be different from the original one
      rdm = drand48();
      if(rdm < cum_proba[0]) {
	if(nuc != 'A') return('A');
      } else if(rdm < cum_proba[1]) {
	if(nuc != 'C') return('C');
      } else if(rdm < cum_proba[2]) {
	if(nuc != 'G') return('G');
      } else {
	if(nuc != 'T') return('T');
      }
    }
  }
  
  return('Z'); // to avoid compiler's complain
}
