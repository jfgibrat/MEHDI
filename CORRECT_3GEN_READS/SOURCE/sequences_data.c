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
#include <stdint.h>
#include <stdlib.h>
#include "alloc.h"

static uint8_t **seq_bits;   // Array containing the compressed sequence of the reads
static int *seq_lengths;     // Vector containing the length of the reads
static char **seq_titles;    // Array containing the titles associated with the reads
static int Nreads = -1;      // Number of reads
static int *seq_len_indx;    // Index to sort the sequence lengths
static int Nbases = -1;      // Total number of nucleotides in the reads
static int first_LR = -1;    // Index of the shortest long read *when the read lengths are sorted*

/*==================================================================================================+
 +==================================================================================================*/
uint8_t **get_seq_bits()
{
  return(seq_bits);
}
/*==================================================================================================+
 +==================================================================================================*/
void set_seq_bits(uint8_t **value)
{
  seq_bits = value;
}
/*==================================================================================================+
 +==================================================================================================*/
int *get_seq_lengths()
{
  return(seq_lengths);
}
/*==================================================================================================+
 +==================================================================================================*/
void set_seq_lengths(int *value)
{
  seq_lengths = value;
}
/*==================================================================================================+
 +==================================================================================================*/
char **get_seq_titles()
{
  return(seq_titles);
}
/*==================================================================================================+
 +==================================================================================================*/
void set_seq_titles(char **value)
{
  seq_titles = value;
}
/*==================================================================================================+
 +==================================================================================================*/
void set_Nreads(int value)
{
  Nreads = value;
}
/*==================================================================================================+
 +==================================================================================================*/
int get_Nreads()
{
  if(Nreads == -1) {
    fprintf(stderr,"Error: Nreads has not been initialized yet\n\n");
    exit(1);
  }
  return(Nreads);
}

/*==================================================================================================+
 +==================================================================================================*/
int get_Nbases()
{
  if(Nbases == -1) {
    fprintf(stderr,"Error: Nbases has not been initialized yet\n\n");
    exit(1);
  }
  return(Nbases);
}

/*==================================================================================================+
 +==================================================================================================*/
int get_first_LR()
{
  if(first_LR == -1) {
    fprintf(stderr,"Error: first_LR has not been initialized yet\n\n");
    exit(1);
  }
  return(first_LR);
}

/*==================================================================================================+
 +==================================================================================================*/
int get_sorted_len_indx(int i)
{

  /*
   * This function returns the new index of the reads when they are sorted according to their length (from shortest to longest)
   */

  if(i < 0 || i >= Nreads){
    fprintf(stderr,"Error in get_sorted_seq_indx: invalid value of index i (%d). Its value should be betwween 0 and %d\n\n",i, Nreads-1);
    exit(1);
  }

  return(seq_len_indx[i]);

}
/*==================================================================================================+
 +==================================================================================================*/
int get_back_original_indx(int sorted_i)
{

  /*
   * This function does the opposite of the previous function. It returns the original index of a read from its sorted index
   */
  if(sorted_i < 0 || sorted_i >= Nreads){
    fprintf(stderr,"Error in get_back_original_indx: invalid value of sorted index i (%d). Its value should be betwween 0 and %d\n\n",sorted_i, Nreads-1);
    exit(1);
  }

  for(int i = 0; i < Nreads; i++) {
    if(seq_len_indx[sorted_i] == i) {
      return(i);
    }
  }

  return(-1);

}

/*==================================================================================================+
 +==================================================================================================*/
int get_sorted_len(int i)
{

  if(i < 0 || i >= Nreads){
    fprintf(stderr,"Error in get_sorted_seq_indx: invalid value of index i (%d). Its value should be betwween 0 and %d\n\n",i, Nreads-1);
    exit(1);
  }

  return(seq_lengths[seq_len_indx[i]]);

}

/*******************************************************************************/
/*                                                                             */
/*                      routine sort_seq_length                                */
/*                                                                             */
/*******************************************************************************/
void iQsort(int v[], int indx[], int left, int right);
char *get_arg_value(char *str);

void split_seqs_in_LRs_SRs()
{

  /*
   * Determine the size that separates long reads (LRs) from short reads (SRs).
   * One must have: total number of nucleotides in SRs = coverage * total number of nucleotides in LRs
   */

  int i;
  int coverage;
  int quot;
  int tmp;

  /*
   * Sort lengths by increasing value
   */

  ALLOC(seq_len_indx,Nreads);
  for(i = 0 ; i < Nreads; i++) {
    seq_len_indx[i] = i;
  }
  iQsort(seq_lengths, seq_len_indx, 0, Nreads-1);

#ifdef DEBUG
  fprintf(stdout,"\n");
  for(i = 0 ; i < Nreads; i++) {
    if(i == 0) {
      fprintf(stdout,"%2dst shortest reads is number %d\n",i+1, seq_len_indx[i]+1);
    } else if(i == 1) {
      fprintf(stdout,"%2dnd shortest reads is number %d\n",i+1, seq_len_indx[i]+1);
    } else if(i == 2) {
      fprintf(stdout,"%2drd shortest reads is number %d\n",i+1, seq_len_indx[i]+1);
    } else {
      fprintf(stdout,"%2dth shortest reads is number %d\n",i+1, seq_len_indx[i]+1);
    }
  }
#endif
  /*
   * Determine the limit between short reads and long reads
   */

  Nbases = 0;
  for(i = 0; i < Nreads; i++) {
    Nbases += seq_lengths[i]; 
  }
  fprintf(stdout,"\nThere is a grand total of %d sequenced nucleotides in the %d reads\n",Nbases, Nreads);

  coverage = atoi(get_arg_value("coverage"));
  quot = Nbases / (coverage + 1);

  tmp = 0;
  for(i = Nreads-1; i >= 0; i--) {
    tmp+= seq_lengths[seq_len_indx[i]];
    if(tmp > quot) {
      break;
    }
  }
  first_LR = i;

  fprintf(stdout,"There are %d 'long reads' totaling %d nucleotides\n",Nreads-first_LR,tmp);
  fprintf(stdout,"\t- The shortest 'long read' contains %d nucleotides. Its index is %d (the index starts at zero)\n",seq_lengths[seq_len_indx[first_LR]], first_LR);
  fprintf(stdout,"\t- The longest  'long read' contains %d nucleotides.\n",seq_lengths[seq_len_indx[Nreads-1]]);
  fprintf(stdout,"There are %d 'short reads' totaling %d nucleotides\n",first_LR,Nbases-tmp);
  fprintf(stdout,"\t- The shortest 'short read' contains %d nucleotides.\n",seq_lengths[seq_len_indx[0]]);
  fprintf(stdout,"\t- The longest  'short read' contains %d nucleotides.\n",seq_lengths[seq_len_indx[first_LR-1]]);

}

