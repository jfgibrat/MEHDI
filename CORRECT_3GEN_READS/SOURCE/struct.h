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

/*
 * Node of the hash table linked lists.
 *
 * Note : if you must change the type of the members in node, 
 * do not forget to modify the sanity_check function accordingly 
 * (in util.c)
 */
struct node		    
{			    
  uint8_t  kmer_tag;        // The K-mer id. (kmer_tag *and the hash value* uniquely identify a given k-mer) 
  uint16_t LongReadNumber;  // The long read number on which the current k-mer matches (numbering starts at 1)
  uint16_t position;	    // The starting position on the long read of the matching k-mer
  struct node *next;	    
};                          

/*
 * Structure returned by the hash function
 */
struct hash_return_value
{
  uint32_t hash_val;        // The remainder of the division by a prime number ==> the hash value 
  uint32_t kmer_tag;        // The dividend  of the division by a prime number ==> the k-mer identifyer
};

/*
 * Node of the long read linked lists
 */
struct node1
{
  uint16_t ShortReadNumber; // Number of a short read whose k-mers map the current long read
  uint16_t Nkmer;           // Numbering of the (non overlapping) k-mer on the short read (starts at 1
  uint16_t position;        // Starting position of the Nkmer-th k-mer on the * current long read*
  struct node1 *next;	    
};
