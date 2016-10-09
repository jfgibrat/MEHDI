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
 * Node of the hash table linked lists
 */
struct node		    
{			    
  uint8_t  kmer_tag;        // If you must change the type of the members in node, 
  uint16_t LongReadNumber;  // do not forget to modify the sanity_check function 
  uint16_t position;	    // accordingly (in util.c)
  struct node *next;	    
};                          

/*
 * Structure returned by the hash function
 */
struct hash_return_value
{
  uint32_t hash_val;
  uint32_t kmer_tag;
};

/*
 * Node of the long read linked lists
 */
struct node1
{
  uint16_t ShortReadNumber;
  uint16_t position;
  uint16_t Nkmer;
  struct node1 *next;	    
};
