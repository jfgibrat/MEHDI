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
#include "struct.h"

// Prototypes of functions called by align_short_reads
struct node *get_HashTable();
uint8_t **get_seq_bits();
int *get_seq_lengths();
int get_Nreads();
int get_first_LR();
char *get_arg_value(char *arg_name);
void band_dyn_prog(uint32_t sorted_LRnum, int beg, int end, uint32_t sorted_SRnum);
void initialize_band_dyn_prog();
int get_sorted_len_indx(int i);
int get_sorted_len(int i);
struct hash_return_value *hash(uint8_t *v, int pos, int kmer_len, uint32_t HashMaxSize) ;
void free_hash_table(uint32_t HashMaxSize);
void LR_excision(int min, int32_t pos, int SR_len, int LR_len, int Kmer_len, int margin, int *beg, int *end);

static struct node1 **heads;    // Contains pointers to the head of linked lists

void align_short_reads(uint32_t HashMaxSize)
{
  
  uint8_t **seq_bits;           // Array containing the compressed sequence of the reads
  int *seq_lengths;             // Vector containing the length of the reads
  int Nreads;                   // Number of reads
  int first_LR;                 // Index of the shortest long read *when the read lengths are sorted by increasing values*
  uint8_t *work;                // Vector containing a read in compressed form
  int Kmer_len;                 // Size of K-mer
  struct node *HashTable;       // Long reads hash table
  int lr;
  int minNkmer;
  uint32_t pos_minNkmer;
  int beg, end;
  int margin;
  struct node1 *zz, *hd;
  struct node *pn;
  struct hash_return_value *pp;
  int jj, kn, jj_start;
  uint32_t pos, SRnum; 
  uint32_t Nkmer;
  uint32_t Nmatches;            // Number of short reads k-mers that map on the long read
  uint32_t NTkmers;             // Number of non overlapping k-mers in the short read
  uint32_t pos_diff;            // Position difference between K-mer x and K-mer y when mapped on long reads
  uint32_t pos_kmer;            // Position difference between K-mer x and K-mer y in short reads
  uint32_t significant;         // Minimum number of k-mer matches to decide that the short read indeed maps to the long read
  uint32_t tolerance;           // Tolerance on the difference between pos_diff and pos_kmer due to the fact
                                // that there are insertions and deletions both in short and long reads

  /*
   * Get the required data
   */
  HashTable = get_HashTable();
  seq_bits = get_seq_bits();
  Nreads = get_Nreads();
  first_LR = get_first_LR();
  seq_lengths = get_seq_lengths();
  Kmer_len = atoi(get_arg_value("Kmer_len"));
  tolerance = atoi(get_arg_value("tolerance"));
  significant = atoi(get_arg_value("significant"));
  margin = atoi(get_arg_value("margin"));

  /*
   * Allocate and initialize the heads of the linked lists for each long read
   */
  ALLOC(heads,Nreads-first_LR);
  for(int i = 0; i < Nreads-first_LR; i++) {
    heads[i] = NULL;
  }

  /*
   * Map the short read k-mers to the long reads using the hash table
   */

  for(int ii = 0, nSR = 1; ii < first_LR; ii++, nSR++) {   // Note: short read numbering (nSR) starts at 1 
    lr = get_sorted_len_indx(ii);
    work = seq_bits[lr];
    NTkmers = seq_lengths[lr] / Kmer_len;                  // Number of non overlapping k-mers in the short read
    jj_start = (NTkmers - 1) * Kmer_len;                   // Starting position of the last k-mer in the short read 

    // K-mer are generated in reverse order to have them stored in increasing order in the linked lists below
    for(jj = jj_start, kn = NTkmers; jj >= 0 ; jj -= Kmer_len, kn--) { // Note: Kmer numbering (kn) starts at 1
      pp = hash(work,jj,Kmer_len,HashMaxSize);

      for(pn = &HashTable[pp->hash_val]; pn != NULL; pn = pn->next) {
	if(pp->kmer_tag != pn->kmer_tag) continue;  // Select the right k-mer (i.e., take care of collisions in the linked list)
	ALLOC(zz,1);
	zz->ShortReadNumber = nSR;
	zz->position = pn->position;
	zz->Nkmer = kn;
	zz->next = heads[pn->LongReadNumber-1];     // The new node is inserted right after the head, thus
	heads[pn->LongReadNumber-1] = zz;           // the short read K-mers are stored in increasing order
      }
    }
  }

  /*
   * The hash table is no longer needed
   */
  free_hash_table(HashMaxSize);

  /*
   * Align the short reads with the long reads using banded dynamic programming
   */

  initialize_band_dyn_prog();

  for(int LRnum = first_LR, jj = 0; LRnum < Nreads; LRnum++, jj++) {
    hd = heads[jj];   // *hd is the first node of current long read linked list
    SRnum = hd->ShortReadNumber;   
    Nkmer = hd->Nkmer;
    pos = hd->position;
    Nmatches = 0;
    minNkmer = Nkmer;
    pos_minNkmer = pos;

    for(zz = hd->next; zz != NULL; zz = zz->next) { 
      if(SRnum == zz->ShortReadNumber) {
	if(zz->Nkmer < minNkmer) {
	  minNkmer = zz->Nkmer;
	  pos_minNkmer = zz->position;
	} 
	pos_diff = zz->position - pos;
	pos_kmer = (zz->Nkmer - Nkmer - 1) * Kmer_len;
	if(abs((int) (pos_diff-pos_kmer)) <= tolerance) {
	  Nmatches++;
	}
      } else {                           // New short read found in linked list
	if(Nmatches >= significant) {    // Process the previous one
	  LR_excision(minNkmer,pos_minNkmer,get_sorted_len(SRnum-1),get_sorted_len(LRnum),Kmer_len,margin,&beg,&end);
	  band_dyn_prog(LRnum,beg,end,SRnum-1);
	}
	SRnum = zz->ShortReadNumber;
	Nkmer = zz->Nkmer;
	pos = zz->position;
	Nmatches = 0;
	minNkmer = Nkmer;
	pos_minNkmer = pos;
      }
    }
    if(Nmatches >= significant) {        // Process the last short read found
      LR_excision(minNkmer,pos_minNkmer,get_sorted_len(SRnum-1),get_sorted_len(LRnum),Kmer_len,margin,&beg,&end);
      band_dyn_prog(LRnum,beg,end,SRnum-1);
    }
    
  }
  
}
/*******************************************************************************/
/*                                                                             */
/*                                                                             */
/*******************************************************************************/
void LR_excision(int min, int32_t pos, int SR_len, int LR_len, int Kmer_len, int margin, int *beg, int *end)
{
  /*
   * This function returns the part of the long read that matches the short read, adding a
   * safety margin on 3' and 5' ends to cope with indels caused by sequencing errors.
   *
   * min is the first short read k-mer that has been found matching the long read, for instance, the 5th k-mer
   * pos is the starting position of this k-mer on the long read
   * SR_len is the short read length
   * LR_len is the long read length
   * Kmer_len is the K-mer length
   * *beg and *end are the starting and ending positions of the long read excised part
   */

  *beg = pos - (min - 1) * Kmer_len;
  *end = *beg + SR_len - 1;
  // Taking indels into account
  *beg = *beg - margin;
  if(*beg < 0) {
    *beg = 0;
  }
  *end = *end + margin;
  if(*end >= LR_len) {
    *end = LR_len - 1;
  }
}
