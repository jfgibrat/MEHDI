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

// Prototypes of functions called by main
struct node *get_HashTable();
uint8_t **get_seq_bits();
int *get_seq_lengths();
int get_Nreads();
int get_first_LR();
char *get_arg_value(char *arg_name);
void band_dyn_prog(uint32_t LRnum, uint32_t SRnum);
void initialize_band_dyn_prog();

static struct node1 *heads;

void align_short_reads(uint32_t HashTableSize)
{

  uint8_t **seq_bits;           // Array containing the compressed sequence of the reads
  int *seq_lengths;             // Vector containing the length of the reads
  int Nreads;                   // Number of reads
  int first_LR;                 // Index of the shortest long read *when the read lengths are sorted by increasing values*
  uint8_t *work;                // Vector containing a read in compressed form
  int Kmer_len;                 // Size of K-mer
  struct node *HashTable;       // Long reads hash table
  int lr;
  struct node1 *zz;
  int jj, kn, jj_start;
  uint32_t pos, SRnum; 
  uint32_t Nkmer, Nmatches;
  uint32_t pos_diff;            // Position difference between K-mer x and K-mer y when mapped on long reads
  uint32_t pos_kmer;            // Position difference between K-mer x and K-mer y in short reads
  uint32_t significant;         // Minimum number of matches to decide that the short read indeed maps to the long read
  uint32_t tolerance;           // Tolerance on the difference between pos_diff and pos_kmer due to the fact
                                // that there are insertions and deletions both in short and long reads

  /*
   * Get the required data
   */
  seq_bits = get_seq_bits();
  Nreads = get_Nreads();
  first_LR = get_first_LR();
  seq_lengths = get_seq_lengths();
  Kmer_len = atoi(get_arg_value("Kmer_len"));
  tolerance = atoi(get_arg_value("tolerance"));
  significant = atoi(get_arg_value("significant"));

  /*
   * Allocate and initialize the heads of the linked lists for each long read
   */
  ALLOC(heads,Nreads-first_LR);
  for(int i = 0; i < Nreads-first_LR; i++) {
    heads[i].next = NULL;
    heads[i].position = 0;
    heads[i].ShortReadNumber = 0;
    heads[i].Nkmer = 0;
  }

  /*
   * Map the short read k-mers to the long reads using the hash table
   */

  for(int ii = 0, nSR = 1; ii < first_LR; ii++, nSR++) {   // Note: short read numbering (nSR) starts at 1 
    lr = get_sorted_len_indx(ii);
    work = seq_bits[lr];
    jj_start = (seq_lengths[lr] / Kmer_len - 1) * Kmer_len;

    // K-mer are generated in reverse order to have them stored in increasing order in the linked lists below
    for(jj = jj_start, kn = 1; jj >= 0 ; jj -= Kmer_len, kn++) { // Note: Kmer numbering (kn) starts at 1
      pp = hash(work,jj,Kmer_len,HashTableSize);

      for(pn = &HashTable[pp->hash_val]; pn != NULL; pn = pn->next) {
	if(pp->kmer_tag == pn->kmer_tag) {
	  if(heads[pn->LongReadNumber-1].ShortReadNumber == 0) {
	    heads[pn->LongReadNumber-1].position = pn->position;
	    heads[pn->LongReadNumber-1].ShortReadNumber = nSR;	
	    heads[pn->LongReadNumber-1].Nkmer = kn;
	  } else {
	    ALLOC(zz,1);
	    zz.ShortReadNumber = nSR;
	    zz.position = pn->position;
	    zz.Nkmer = kn;
	    zz.next = heads[pn->LongReadNumber-1].next; // The new node is inserted right after the head, thus
	    heads[pn->LongReadNumber-1].next = &zz;     // the short read K-mers are stored in increasing order
	  }
	}
      }
    }
  }

  /*
   * Align the short reads with the long reads using banded dynamic programming
   */

  initialize_band_dyn_prog();

  for(int LRnum = first_LR, jj = 0; LRnum < Nreads; LRnum++, jj++) {
    SRnum = heads[jj].ShortReadNumber;   // In heads[], long reads are stored starting at 0
    Nkmer = heads[jj].Nkmer;
    pos = heads[jj].position;
    Nmatches = 0;
    for(zz = heads[jj]->next; zz != NULL; zz = zz->next) { 
      if(SRnum == zz->ShortReadNumber) { 
	pos_diff = zz->position - pos;
	pos_kmer = (zz->Nkmer - Nkmer - 1) * Kmer_len;
	if(abs(pos_diff-pos_kmer) <= tolerance) {
	  Nmatches++;
	}
      } else {                           // New short read found in linked list
	if(Nmatches >= significant) {    // Process the previous one
	  band_dyn_prog(LRnum,SRnum-1);
	}
	Nmatches = 0;
      }
      SRnum = zz->ShortReadNumber;
      Nkmer = zz->Nkmer;
      pos = zz->position;
    }
    if(Nmatches >= significant) {        // Process the last short read found
      band_dyn_prog(LRnum,SRnum-1);
    }

  }

}
