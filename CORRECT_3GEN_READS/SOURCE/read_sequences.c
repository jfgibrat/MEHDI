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

void set_seq_bits(uint8_t **value);
void set_seq_lengths(int *value);
void set_seq_titles(char **value);
void set_Nreads(int value);
void dna_bit_encode(uint8_t **m_data, char *dna_str, int dna_len);
char *get_arg_value(char *arg_name);

void read_sequence_file()
{
  /*
   * This function reads the file containing the 3rd generation sequencing technology 'reads' and store them in a compressed form.
   * Each nucleotide {A, C, G, T} is coded with 2 bits. This allows use to pack 4 nucleotides per byte.
   */

  FILE *fp;
  int ii;
  int ll;
  int Nreads;
  char buffer[BUFSIZE];
  uint8_t **seq_bits;
  int *seq_lengths;
  char **seq_titles;

  FOPEN(fp,get_arg_value("SeqInputFile"),"r");

  Nreads = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '>') {
      Nreads++;
    }
  }
  set_Nreads(Nreads);

  rewind(fp);

  ALLOC(seq_lengths,Nreads);
  ALLOC(seq_titles,Nreads);
  ALLOC(seq_bits,Nreads);

  ii = -1;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '>') {
      ii++;
      ll = strlen(buffer);
      ALLOC(seq_titles[ii], ll);
      buffer[ll-1] = '\0';
      strcpy(seq_titles[ii],buffer);
    } else {
      seq_lengths[ii] = strlen(buffer) - 1;
      buffer[seq_lengths[ii]] = '\0';
      dna_bit_encode(&(seq_bits[ii]),buffer,seq_lengths[ii]);
    }
  }

  set_seq_bits(seq_bits);
  set_seq_lengths(seq_lengths);
  set_seq_titles(seq_titles);

}
/************************************************************************************************************
 *                                                                                                          *
 *                                                                                                          *
 ************************************************************************************************************/
 
void dna_bit_encode(uint8_t **m_data, char *dna_str, int dna_len)
{
  /*
   * This function encodes the 4 nucleotides {A, C, G, T} with 2 bits and store the compressed DNA sequence in m_data
   */

  uint8_t shift;
  int i;
  int dna_bytes;

  dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);
  ALLOC((*m_data), dna_bytes);
  memset((*m_data),0,dna_bytes);

  for(i = 0; i < dna_len; i++) {
    shift = 6 - 2 * (i % 4);
    switch(dna_str[i])
      {
      case 'A':
	(*m_data)[i/4] |= 0x0 << shift;
	break;
      case 'C':
	(*m_data)[i/4] |= 0x1 << shift;
	break;
      case 'G':
	(*m_data)[i/4] |= 0x2 << shift;
	break;
      case 'T':
	(*m_data)[i/4] |= 0x3 << shift;
	break;
      default:
	fprintf(stderr,"dna_bit_encode: Invalid DNA base : %c\n\n",dna_str[i]);
	exit(1);
      }
    shift = (shift == 0) ? 6 : shift-2;
  }
  
}
/***********************************************************
 *                                                         *
 *                                                         *
 ***********************************************************/
#define BASE_MASK 0x3

void dna_bit_decode(char **dna_str, uint8_t *m_data, int dna_len, int skip_alloc)
{
  /*
   * This function does the converse of function 'dna_bit_encode'.
   * It recovers, from m_data, the original DNA sequence as a character string
   */

  uint8_t mask;
  uint8_t shift;
  uint8_t base;

  if(skip_alloc) {
    ALLOC((*dna_str), dna_len+1);
  }

  for(int i = 0; i < dna_len; i++) {
    shift = 6 - 2 * (i % 4);
    mask = BASE_MASK << shift;
    // Get the ith DNA base
    base = (m_data[i/4] & mask) >> shift;

    switch(base)
      {
      case 0x0 :
	(*dna_str)[i] = 'A';
	break;
      case 0x1 :
	(*dna_str)[i] = 'C';
	break;
      case 0x2 :
	(*dna_str)[i] = 'G';
	break;
      case 0x3 :
	(*dna_str)[i] = 'T';
	break;
      default:
	fprintf(stderr,"dna_bit_decode: Invalid DNA base : %x\n\n",base);
	exit(1);
      }
  }

  (*dna_str)[dna_len] = '\0';

}
