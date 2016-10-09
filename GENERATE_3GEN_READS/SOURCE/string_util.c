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
#include "alloc.h"

struct tokenNode {
  char *token;
  struct tokenNode *next;
};

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
char *strip_string(char *string)
{
  int len;
  int i;

  /*
   * Remove leading and trailing blanks (' ', '\t', '\r', '\n') of a character string.
   * Note: the initial string is modified since an end of string character ('\0') is put
   * after the last non blank character to get rid of the trailing blanks.
   */

  len = strlen(string);

  for(i = len-1; i >= 0; i--) {
    if(string[i] != ' ' && string[i] != '\t' && string[i] != '\r' && string[i] != '\n') {
      break;
    }
  }

  string[i+1] = '\0';
  len = strlen(string);
  if(len == 0) {
    return(string);
  }

  for(i = 0; i < len; i++) {
    if(string[i] != ' ' && string[i] != '\t' && string[i] != '\r' && string[i] != '\n') {
      break;
    }
  }

  return(string+i);
  
}

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
int token_list(const char *string, char *delim, struct tokenNode **head)
{
/*
 * Return a list of tokens separated by delimiters (specified in delim) that are found in the string.
 * Since the number of tokens is not known in advance these tokens are stored in a linked list.
 * The number of tokens found is returned by the function.
 */

  int Ntok;
  size_t word_len, sep_len;
  struct tokenNode *t = NULL;

  word_len = sep_len = 0;
  Ntok = 0;

  while(*string != '\0') {
    Ntok++;
    word_len = strcspn(string, delim);
    sep_len =   strspn(string + word_len, delim);

    if(Ntok == 1) {
      ALLOC(*head,1);
      t = *head;
    } else {
      ALLOC(t->next,1);
      t = t->next;
    }

    t->token = NULL;
    ALLOC(t->token,word_len + 1);
    strncpy(t->token,string,word_len);
    t->token[word_len] = '\0';
    t->next = NULL;

    string = string + (word_len + sep_len) * sizeof(char);
  }

  return(Ntok);
  
}
/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
void free_token_list(struct tokenNode *head)
{
  struct tokenNode *t, *p = NULL;

  for(t = head; t != NULL; t = t->next) {
    if(p != NULL) {
      FREE(p);
    }
    FREE(t->token);
    p = t;
  }
  FREE(p);
}
