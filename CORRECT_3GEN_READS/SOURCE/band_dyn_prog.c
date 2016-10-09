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

// Prototypes of functions called by the next function 
char *get_arg_value(char *arg_name);
int get_first_LR();

static char *LRseq, *SRseq;

void initialize_band_dyn_prog()
{
  uint32_t longest_LR, longest_SR;
  int nn;
  int first_LR;

  first_LR = get_first_LR();
  longest_SR = get_sorted_len(first_LR-1);
  nn = (uint32_t) longest_LR * 1.2;
  ALLOC(SRseq, nn);

  longest_LR = get_sorted_len(get_Nreads-1);
  nn = (uint32_t) longest_LR * 1.2;
  ALLOC(LRseq, nn);
}

/************************************************************************************************************
 *                                                                                                          *
 *                                                                                                          *
 ************************************************************************************************************/
// Prototypes of functions called by the next function 
uint8_t **get_seq_bits();
void dna_bit_decode(char **dna_str, uint8_t *m_data, int dna_len, int skip_alloc);
int *get_seq_lengths();

void band_dyn_prog(uint32_t LRnum, uint32_t SRnum)
{
  int band_width;               // The band width used is in fact 2 * band_width + 1
  uint8_t **seq_bits;           // Array containing the compressed sequence of the reads
  int *seq_lengths;

  band_width = atoi(get_arg_value("band_width"));
  seq_bits = get_seq_bits();
  seq_lengths = get_seq_lengths();
  dna_bit_decode(&LRseq, seq_bits[LRnum], seq_lengths[LRnum],1);
  dna_bit_decode(&SRseq, seq_bits[SRnum], seq_lengthsSLRnum],1);

}

