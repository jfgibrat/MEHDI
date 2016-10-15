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

static struct node *HashTable;

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
void initialise_hash_table(uint32_t HashMaxSize)
{
  int i;
  ALLOC(HashTable, HashMaxSize);

  for(i = 0; i < HashMaxSize; i++) {
    HashTable[i].LongReadNumber = 0;
    HashTable[i].position = 0;
    HashTable[i].kmer_tag = 0;
    HashTable[i].next = NULL;
  }
}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
void free_hash_table(uint32_t HashMaxSize)
{
  struct node *keep;

  for(int i = 0; i < HashMaxSize; i++) {
    while((keep = HashTable[i].next) != NULL) {
      HashTable[i].next = keep->next;
      FREE(keep);
    }
  }

  FREE(HashTable);

}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
struct node *get_HashTable()
{
  return(HashTable);
}

/************************************************************************************************************
 *                                                                                                          *
 *                                                                                                          *
 ************************************************************************************************************/
char *get_arg_value(char *arg_name);
FILE *Get_File_Ptr(char *fname, char *access_mode);

uint32_t determine_hash_table_size()
{
  /* 
   *  This function returns the closest prime number to HashTabSiz_tmp among those contained in the file
   *  defined by the argument value 'PN_infile' 
   *  HashTabSiz_tmp is the order of magnitude chosen by the user for the size of the hash table.
   */

  FILE *fp;
  int i;
  uint32_t HashMaxSize;
  uint32_t tmp;
  char *fname;
  char buffer[BUFSIZE];
  uint32_t *PrimeNumbers;
  int Nprimes;

  fname = get_arg_value("PN_infile");
  tmp = atoi(get_arg_value("HashTabSiz_tmp"));

  fp = Get_File_Ptr(fname,"r");
  Nprimes = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    Nprimes++;
  }
  rewind(fp);
  ALLOC(PrimeNumbers, Nprimes);
  for(i = 0; i < Nprimes; i++) {
    fscanf(fp,"%u",PrimeNumbers+i);
  }
  FCLOSE(fp);

  for(i = 0; i < Nprimes; i++) {
    if(tmp < PrimeNumbers[i]) {
      break;
    }
  }

  if(i == 0) {
    printf("\nWarning the value of the hash table size (%u) is less than the first prime number available (%u)\n",tmp,PrimeNumbers[0]);
    printf("Setting hash table size to the latter value\n");
    HashMaxSize = PrimeNumbers[0];
  } else if(i == Nprimes) {
    printf("\nWarning the value of the hash table size (%u) is larger than the last prime number available (%u)\n",tmp,PrimeNumbers[Nprimes-1]);
    printf("Setting hash table size to the latter value\n");
    HashMaxSize = PrimeNumbers[Nprimes-1];
  } else {
    printf("\nSetting hash table size to %u\n",PrimeNumbers[i]);
    HashMaxSize = PrimeNumbers[i];
  }

  return(HashMaxSize);
}
/*----------------------------------------------------------*/
/*|                     DEBUG FUNCTION                     |*/
/*----------------------------------------------------------*/
char decode_base(uint8_t base)
{
      switch(base)
      {
      case 0x0 :
	return('A');
      case 0x1 :
	return('C');
      case 0x2 :
	return('G');
      case 0x3 :
	return('T');
      default:
	fprintf(stderr,"print_base: Invalid DNA base : %x\n\n",base);
	exit(0);
      }
}
void print_base(int i, uint8_t base)
{
      switch(base)
      {
      case 0x0 :
	fprintf(stdout,"\t the %d base is 'A'\n",i);
	break;
      case 0x1 :
	fprintf(stdout,"\t the %d base is 'C'\n",i);
	break;
      case 0x2 :
	fprintf(stdout,"\t the %d base is 'G'\n",i);
	break;
      case 0x3 :
	fprintf(stdout,"\t the %d base is 'T'\n",i);
	break;
      default:
	fprintf(stderr,"print_base: Invalid DNA base : %x\n\n",base);
	exit(0);
      }
}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
static uint64_t *power4;

void initialize_power4_array()
{
  int Kmer_size;
  int i;

  Kmer_size = atoi(get_arg_value("Kmer_len"));
  ALLOC(power4,Kmer_size);
  power4[0] = 1;
  for(i = 1; i < Kmer_size; i++) {
    power4[i] = power4[i-1] * 4;
  }
}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
uint64_t *get_power4_array()
{
  return(power4);
}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
uint8_t return_base(uint8_t *v, int pos);

static struct hash_return_value hrv;

struct hash_return_value *hash(uint8_t *v, int pos, int kmer_len, uint32_t HashMaxSize) 
{
  /*
   * This function returns the dividend and the remainder of h when divided by HashMaxSize (a prime number)
   * The remainder is used as a hash value in a hash table and the dividend is used to solve collision issues
   *
   * The k-mer is considered to be a number in base 4. The function first transforms it into a decimal value (h), i.e., 
   * h = v[0] * 4^0 + v[1] * 4^1 + v[2] * 4^2 + ... + v[kmer_len-1] * 4^(kmer_len-1)
   * v contains the K-mer sequence encoded with 2 bits 'A' = 0x0, 'C' = 0x1, 'G' = 0x2, 'T' = 0x3 (0, 1, 2, 3 in decimal)
   * kmer_len is the k-mer length (pos is the k-mer starting position in the read) 
   * The function then computes the dividend and the remainder.
   *
   * To avoid overflow of h for large values of k, the number h can be divided into pieces that are processed independently.
   * For instance, suppose h = h1 + h2 + h3 with
   * h1 = d1 * z + r1
   * h2 = d2 * z + r2
   * h3 = d3 * z + r3
   * Therefore :
   * h = (d1 + d2 + d3) * z + (r1 + r2 + r3). Let (r1 + r2 + r3) = d4 * z + r4. Thus 
   * h = (d1 + d2 + d3 + d4) * z + r4
   *
   */
  uint64_t rr = 0; // r4 in the above equation
  uint64_t dd = 0; // (d1 + d2 + d3 + d4) in the above equation 
  uint64_t ri;     // r1, r2 or r3 in the above equation 
  uint64_t di;     // d1, d2 or d3
  uint64_t hh;     // h1, h2 or h3 (or h if Kmer_len <= 30)
  int i;
  uint8_t base;

  /*
   * For small values of the K-mer it takes less operations to compute h directly then computes the dividend and remainder.
   * The largest integer that can fit in an uint64_t is 2^64 - 1. If the size of the K-mer is 31 the largest integer
   * that we can compute as explained above is 4^32-1 = 2^64-1 which will just fit.
   */
  if(kmer_len > 30 || atoi(get_arg_value("bypass_kmer_len"))) {
    for(i = 0; i < kmer_len; i++) {
      base = return_base(v,pos+i);
      hh = base * power4[i];
      ri = hh % HashMaxSize;
      di = hh / HashMaxSize;
      dd += di + (rr + ri) / HashMaxSize;
      rr = (rr + ri) % HashMaxSize;
    }
  } else {
    hh = 0;
    for(i = 0; i < kmer_len; i++) {
      base = return_base(v,pos+i);
      hh += base * power4[i];
    }
    dd = hh / HashMaxSize;
    rr = hh % HashMaxSize;
  }

  hrv.kmer_tag = dd;
  hrv.hash_val = rr;

  return(&hrv);

}
/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
uint8_t **get_seq_bits();
int *get_seq_lengths();
int get_Nreads();
int get_first_LR();
int get_sorted_len_indx(int i);
void initialize_power4_array();
char **get_seq_titles();

void hash_long_reads(uint32_t HashTableSize)
{

  int ii, jj, lr;
  uint8_t **seq_bits;           // Array containing the compressed sequence of the reads
  int *seq_lengths;             // Vector containing the length of the reads
  int Nreads;                   // Number of reads
  int first_LR;                 // Index of the shortest long read *when the read lengths are sorted by increasing values*
  uint8_t *work;                // Vector containing a read in compressed form
  int Kmer_len;                 // Size of K-mer
  struct hash_return_value *pp; // Pointer to the structure returned by the hash function
  char **seq_titles;
  int nLR;
  struct node *pn;

  /*
   * Get the required data
   */
  seq_bits = get_seq_bits();
  Nreads = get_Nreads();
  first_LR = get_first_LR();
  seq_lengths = get_seq_lengths();
  Kmer_len = atoi(get_arg_value("Kmer_len"));
  seq_titles = get_seq_titles();

  /*
   * Initialize data for the computation of the hash value
   */
  initialize_power4_array();

  /*
   * Hash the reads
   */
  for(ii = first_LR, nLR = 1; ii < Nreads; ii++, nLR++) {   // Note: the first long read has number 1
    lr = get_sorted_len_indx(ii);
    work = seq_bits[lr];
    for(jj = 0; jj <= seq_lengths[lr] - Kmer_len; jj++) {
      pp = hash(work,jj,Kmer_len,HashTableSize);
      if(HashTable[pp->hash_val].LongReadNumber == 0) {
	HashTable[pp->hash_val].LongReadNumber = nLR;    
	HashTable[pp->hash_val].position = jj;              // The position in reads starts at 0
	HashTable[pp->hash_val].kmer_tag = pp->kmer_tag;
      } else {
	ALLOC(pn,1);
	pn->LongReadNumber = nLR;
	pn->position = jj;
	pn->kmer_tag = pp->kmer_tag;
	pn->next = HashTable[pp->hash_val].next;
	HashTable[pp->hash_val].next = pn;
      }
    }
  }

}
/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
int explore_lkd_list(uint32_t hash_val)
{
  struct node *pn;
  int i = 0;

#ifdef DEBUG1
  fprintf(stdout,"Hash value = %d\n",hash_val);
#endif

  for(pn = &HashTable[hash_val]; pn != NULL; pn = pn->next) {
#ifdef DEBUG1
    fprintf(stdout,"\taddress: %p  position=%u  LongReadNumber=%u kmer_tag=%u\n",pn,pn->position,pn->LongReadNumber,pn->kmer_tag);
#endif
    if(pn->LongReadNumber != 0) {
      i++;
    }
  }

  return(i);

}

/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
#define BASE_MASK 0x3

uint8_t return_base(uint8_t *v, int pos)
{
  /*
   * Returns the value of the pos-th nucleotide in a read coded on 2 bits
   * v contains the read sequenced coded on 2 bits
   * A => 0x0
   * C => 0x1
   * G => 0x2
   * T => 0x3
   */

  uint8_t mask;
  uint8_t shift;
  uint8_t base;

  shift = 6 - 2 * (pos % 4);
  mask = BASE_MASK << shift;
  base = (v[pos/4] & mask) >> shift;
  return(base);

}
