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
uint32_t determine_hash_table_size();
void read_sequence_file();
void initialise_hash_table(uint32_t HashMaxSize);
int *get_seq_lengths();
char **get_seq_titles();
uint8_t **get_seq_bits();
int get_Nreads();
void split_seqs_in_LRs_SRs();
void hash_long_reads(uint32_t HashTableSize);
int get_sorted_len_indx(int i);
int get_sorted_len(int i);
int get_first_LR();
void sanity_check(int NberLongReads, int LongestReadLength, int Kmer_len, uint32_t HashMaxSize);
int explore_lkd_list(uint32_t hash_val);
void dna_bit_encode(uint8_t **m_data, char *dna_str, int dna_len);
void dna_bit_decode(char **dna_str, uint8_t *m_data, int dna_len, int skip_alloc);
void align_short_reads(uint32_t HashTableSize);

int main(int argc, char *argv[])
{

  uint32_t HashMaxSize;
  int Kmer_size;
  int nn;

  process_line_args(argc,argv);

  Kmer_size = atoi(get_arg_value("Kmer_len"));
  fprintf(stdout,"Size of k-mer= %d\n",Kmer_size);
  
  /*
   * Read the sequences of the '3rd generation reads' 
   */  
  read_sequence_file();

#ifdef DEBUG
  char **p1;
  int *p2;
  uint8_t **p3;
  char *p4;
  uint8_t *p5;

  p1 = get_seq_titles();
  p2 = get_seq_lengths();
  p3 = get_seq_bits();
  nn = get_Nreads();
  fprintf(stdout,"Number of reads=%d\n",nn);
  for(int i = 0; i < nn; i++) {
    fprintf(stdout,"%s --- length= %d\n",p1[i],p2[i]);
    dna_bit_decode(&p4,p3[i],p2[i],0);
    fprintf(stdout,"%s\n",p4);
    dna_bit_encode(&p5,p4,p2[i]);
    FREE(p4);
    FREE(p5);
  }
#endif

  /*
   * Divide the reads into long reads (LRs) and short reads (SRs) according to the specified coverage
   */
  split_seqs_in_LRs_SRs();

  /*
   * Initialize the hash table. Its size is HashMaxSize (a prime number)
   */
  if(atoi(get_arg_value("HashTabSiz_tmp")) < 0) {
    HashMaxSize = -atoi(get_arg_value("HashTabSiz_tmp"));
  } else {
    HashMaxSize = determine_hash_table_size();
  }

  nn = get_Nreads();
  sanity_check(nn-get_first_LR(),get_sorted_len(nn-1),atoi(get_arg_value("Kmer_len")),HashMaxSize);
  
  if(atoi(get_arg_value("try"))) {
    return(0);
  }

  initialise_hash_table(HashMaxSize);

  /*
   * Hash the long reads
   */
  hash_long_reads(HashMaxSize);

#ifdef DEBUG
  int *tally;
  int *nbers;
  int mxl;
  uint32_t ii;
  int i;

  ALLOC(tally,HashMaxSize);
  for(ii = 0; ii < HashMaxSize; ii++) {
    tally[ii] = explore_lkd_list(ii);    
  }
  mxl = tally[0];
  for(i = 1; i < HashMaxSize; i++) {
    if(tally[i] > mxl) mxl = tally[i];
  }
  ALLOC(nbers,mxl+1);
  for(i = 0; i <= mxl; i++) {
    nbers[i] = 0;
  }
  for(i = 0; i < HashMaxSize; i++) {
    nbers[tally[i]] += 1;;
  }
  for(i = 0; i <= mxl; i++) {
    fprintf(stdout,"There are %d positions in the hash table with %d nodes\n",nbers[i],i);
  }
  FREE(tally);
  FREE(nbers);
#endif

  /*
   * Align the short reads on the long reads
   */
  align_short_reads(HashMaxSize);

  return(0);

}
