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
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "alloc.h"

char *Get_Line(char *buffer, int buffer_size, FILE *fp)
{

/*
 * Read a line in a buffer and *test* whether a buffer overflow occurred.
 */

  char *returned_value;

  returned_value = fgets(buffer, buffer_size, fp);

  if(strlen(buffer) == buffer_size - 1) {
    PRINT_TAG_ERR;
    fprintf(stderr,"Buffer size too small for the line read. Actual size is %d. You must redefine the size correctly\n",buffer_size);
    exit(1);
  }

  return(returned_value);

}
/*==================================================================================================+
 |                                                                                                  |
 |                                                                                                  |
 +==================================================================================================*/
FILE *Get_File_Ptr(char *fname, char *access_mode)
{

/*
 * This function returns a pointer to a file described in string fname.
 * If fname corresponds to an absolute path, i.e., starts with a '/' it is used as such,
 * otherwise the file name is considered to be relative to the directory stored in 
 * environment variable 'PROG_DIR'.
 */

  FILE *fp = NULL;
  char *PROGdir;
  char newFname[BUFSIZE];

  if(access_mode[0] == 'w') {
    /*
     *  Written file names are always taken verbatim
     */
    strcpy(newFname,fname);

  } else {
    if(fname[0] == '/') {
      /*
       * The file is defined by an absolute path.
       */
      strcpy(newFname,fname);
      
    } else if (fname[0] == '.' && fname[1] == '/') {
      /*
       * The file is in the current directory
       */
      strcpy(newFname,fname);

    } else if (fname[0] == '.' && fname[1] == '.') {
      /*
       * The file is given relative to the current directory, e.g., ../foo.txt in the parent directory of the current one.
       */
      strcpy(newFname,fname);

    } else {
      /*
       * The file is defined relative to the directory stored in environment variable 'PROG_DIR'.
       */
      PROGdir = getenv(PROG_DIR);

      if(PROGdir == NULL) {
	PRINT_TAG_ERR;
	fprintf(stderr,"If you do not use an absolute file name path or do not prefix it with ./ or ../ (%s),\n",fname);
	fprintf(stderr,"you must define an environment variable PROG_DIR pointing to the PROG installation directory\n\n");
	exit(1);
      }

      strcpy(newFname,PROGdir);
      strcat(newFname,"/");
      strcat(newFname,fname);

    }
  }

  fp = fopen(newFname,access_mode);

  return(fp);

}
/*******************************************************************************/
/*                                                                             */
/*                      routine iQsort                                         */
/*                                                                             */
/*******************************************************************************/
void swap(int indx[], int i, int j);

void iQsort(int v[], int indx[], int left, int right)
{

/*
 * sort v[left]...v[right] into increasing order. Note that array v itself is left
 * unmodified and only the pointers to the elements of v, indx are really sorted.
 * Do not forget to initialize the indx array before the first call to iQsort, i.e.,
 *      for(i = 0 ; i < N; i++) {indx[i] = i;}
 *      iQsort(v,indx,0,N-1);
 * where N is the size of v
 */

  int i, last;

  if(left >= right)    /* Do nothing if array contains */
    return;            /* fewer than two elements */

  swap(indx, left, (left+right)/2);  /* move partition elements */
  last = left;                       /* to v[0] */

  for(i = left + 1; i <= right; i++)  /* partition */
    if(v[indx[i]] < v[indx[left]])
       swap(indx, ++last, i);

  swap(indx, left, last);    /* restore partition element */

  iQsort(v, indx, left, last-1);
  iQsort(v, indx, last+1, right);

}


/*******************************************************************************/
/*                                                                             */
/*                      routine rQsort                                         */
/*                                                                             */
/*******************************************************************************/
void swap(int indx[], int i, int j);

void rQsort(float v[], int indx[], int left, int right)
{

/*
 * sort v[left]...v[right] into increasing order. Note that array v itself is left
 * unmodified and only the pointers to the elements of v, indx are really sorted.
 *
 */

  int i, last;

  if(left >= right)    /* Do nothing if array contains */
    return;            /* fewer than two elements */

  swap(indx, left, (left+right)/2);  /* move partition elements */
  last = left;                       /* to v[0] */

  for(i = left + 1; i <= right; i++)  /* partition */
    if(v[indx[i]] < v[indx[left]])
       swap(indx, ++last, i);

  swap(indx, left, last);    /* restore partition element */

  rQsort(v, indx, left, last-1);
  rQsort(v, indx, last+1, right);

}

/*******************************************************************************/
/*                                                                             */
/*                               routine swap                                 */
/*                                                                             */
/*******************************************************************************/

void swap(int v[], int i, int j)
{

/*
 * Interchange indx[i] and indx[j]
 *
 */

  int temp;

  temp = v[i];
  v[i] = v[j];
  v[j] = temp;

}
/*******************************************************************************/
/*                                                                             */
/*                               function chomp                                */
/*                                                                             */
/*******************************************************************************/
int chomp(char *str)
{

/*
 * Get rid of the last '\n' character and return the new string length
 */

  int ll;

  ll = strlen(str);

  if(str[ll-1] != '\n') {
    fprintf(stdout,"Function chomp warning: character string (%s) does not end with a newline character\n\n",str);
    return(ll);
  } else {
    str[ll-1] = '\0';
    return(ll-1);
  }

}
/*******************************************************************************/
/*                                                                             */
/*                            function sanity_check                            */
/*                                                                             */
/*******************************************************************************/
void sanity_check(int NberLongReads, int LongestReadLength, int Kmer_len, uint32_t HashMaxSize)
{
  /*
   * This functions verifies that there is no risk of overflow in the hash table.
   *
   * In file hashtable.c, a node of the hash table linked lists is defined as:
   *   struct node		     
   *   {			    
   *     uint8_t  kmer_tag;	    
   *     uint16_t LongReadNumber;  
   *     uint16_t position;	    
   *     struct node *next;	    
   *   };
   *
   *  It is very important to be extremely parcimonious in terms of node size because there are as many nodes
   *  as nucleotides in all the long reads (typically several hundreds of millions). 
   *  Here, the size of the node is 1 + 2 + 2 + 8 = 13 bytes.
   *
   *  The largest value that 'position' can assume is the length of the longest read minus the k-mer length
   *  The largest value that 'longReadNumber' can assume is the total number of long reads
   *  The largest value that 'kmer_tag' can assume is the dividend of the largest integer that can be generated
   *  with the Kmer (see function hash) by the size of the hash table HashMaxSize (usually a very large prime number).
   *  Recall that the largest integers that can fit in uint8_t and uint16_t are respectively 255 and 65535
   */
  int i;
  int nerr = 0;
  uint64_t hh;
  uint64_t ri, di;
  uint64_t dd = 0, rr = 0; // dd dividend, rr remainder of hh / HashMaxSize
  uint64_t power4;
  
  /*
   * Compute the dividend piecewise to avoid overflow of h
   */
  power4 = 1;
  for(i = 0; i < Kmer_len; i++) {
    hh = 3 * power4;
    ri = hh % HashMaxSize;
    di = hh / HashMaxSize;
    dd += di + (rr + ri) / HashMaxSize;
    rr = (rr + ri) % HashMaxSize;
    power4 *= 4;
  }

  if(dd > 255) {
    nerr++;
    fprintf(stderr,"\nThe maximum value that member 'kmer_tag' in structure node can take is 255. here the value is %llu. IT WILL OVERFLOW\n",dd);
    fprintf(stderr,"Upgrade the type of this member from uint8_t to uint16_t in file hashtable.c and recompile the program\n");
  }

  if(LongestReadLength - Kmer_len > 65535) { // Note: the position in the reads starts at 0
    nerr++;
    fprintf(stderr,"\nThe maximum value that member 'position' in structure node can take is 65535. Here the value is %d. IT WILL OVERFLOW\n",LongestReadLength);
    fprintf(stderr,"Upgrade the type of this member from uint16_t to uint32_t in file hashtable.c and recompile the program\n");
  }

  if(NberLongReads > 65535) {
    nerr++;
    fprintf(stderr,"\nThe maximum value that member 'LongReadNumber' in structure node can take is 65535. Here the value is %d. IT WILL OVERFLOW\n",NberLongReads);
    fprintf(stderr,"Upgrade the type of this member from uint16_t to uint32_t in file hashtable.c and recompile the program\n");
  }

  if(nerr != 0) {
    fprintf(stderr,"Of course, after modifying the 'node' member types in hashtable.c, do not forget to modify the 'sanity_check' function (in util.c) accordingly!\n");
    exit(1);
  } 

}
