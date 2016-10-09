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
#include <ctype.h>
#include "alloc.h"

struct tokenNode {
  char *token;
  struct tokenNode *next;
};

#define DIM2 5

/*
 * Global variables to this file
 */
static char ***DefaultArguments;
static int maxargs;

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
char *strip_string(char *string);
int   token_list(const char *string, char *delim, struct tokenNode **head);
char *Get_Line(char *buffer, int buffer_size, FILE *fp);
FILE *Get_File_Ptr(char *fname, char *access_mode);
void  free_token_list(struct tokenNode *head);

void initDefaultArguments()
{

/*
 * Read the default arguments on a file, 
 * allocate array DefaultArguments, 
 * fill it with the values read.
 */

  FILE *fp = NULL;
  char buffer[BUFSIZE];
  char filename[BUFSIZE];
  int mxsiz, toklen;
  int i, j;
  int maxcol;
  struct tokenNode *head = NULL, *t = NULL;

/*
 * Open the default file (look first in the current directory then in $PROG_DIR/DATA)
 */

  strcpy(filename,getenv("PWD"));
  strcat(filename,"/DefaultArguments.dat");
  fp = fopen(filename,"r");
  if(fp == NULL) {
    fp = Get_File_Ptr("DATA/DefaultArguments.dat","r");
    if(fp == NULL) {
      fprintf(stderr,"Unable to find file DefaultArguments.dat neither in the current directory nor in $PROG_DIR/DATA\n\n");
      exit(1);
    }
  }



/*
 * Read the file a first time to determine the dimensions of array DefaultArguments.
 */

  maxargs = mxsiz = 0;
  while(Get_Line(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '#' || buffer[0] == '\n') {
      continue;
    }
    maxargs++;
    maxcol = token_list(buffer,"|",&head);
    if(maxcol > DIM2) {
      PRINT_TAG_ERR;
      fprintf(stderr,"The second dimension of array DefaultArguments is set to %d\n",DIM2);
      fprintf(stderr,"If you really need to change this dimension you must modify accordingly a number of functions in file arguments.c\n");
      exit(1);
    }
    for(t = head; t != NULL; t = t->next) {
      toklen = strlen(strip_string(t->token));
      if(mxsiz < toklen) {
	mxsiz =  toklen;
      }
    }
    free_token_list(head);
    head = NULL;
  }
/*
 * Allocate array DefaultArguments
 */

  ALLOC3(DefaultArguments,maxargs,DIM2,mxsiz+1);

/*
 * Read the file a second time to fill array DefaultArguments
 */

  rewind(fp);

  i = 0;
  while(Get_Line(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '#' || buffer[0] == '\n') {
      continue;
    }
    maxcol = token_list(buffer,"|",&head);
    for(t = head, j = 0; t != NULL; t = t->next, j++) {
      strcpy(DefaultArguments[i][j],strip_string(t->token));
    }
    free_token_list(head);
    head = NULL;
    i++;
  }

  FCLOSE(fp);

}
/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
void print_usage(char *progname)
{
/*
 * Print the program help message
 */

  int i;
  int mxs=0;
  char *PROGdir = NULL;
  char buffer[BUFSIZE];
  char fmt[BUFSIZE];

  sprintf(buffer,"\nUsage: %s ",progname);
  for(i = 0; i < maxargs; i++) {
    if(strlen(DefaultArguments[i][0]) > mxs) {
      mxs = strlen(DefaultArguments[i][0]);
    }
    if(strcmp(DefaultArguments[i][3],"compulsory") == 0) {
      strcat(buffer,DefaultArguments[i][1]);
      strcat(buffer," ");
      strcat(buffer,DefaultArguments[i][0]);
      strcat(buffer," ");
    } else {
      strcat(buffer,"[");
      strcat(buffer,DefaultArguments[i][1]);
      strcat(buffer," ");
      strcat(buffer,DefaultArguments[i][0]);
      strcat(buffer,"] ");
    }
  }

  fprintf(stdout,"%s\nwhere\n",buffer);

  sprintf(fmt,"%%-%ds: %%s [%%s]\n",mxs);
  for(i = 0; i < maxargs; i++) {
    fprintf(stdout,fmt,DefaultArguments[i][0],DefaultArguments[i][4],DefaultArguments[i][3]);
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"Note: Names of data files that are created (output files) are always taken verbatim, i.e.,\n");
  fprintf(stdout,"      if you specify 'output.txt' the program will create this file in the working directory.\n\n"); 
  fprintf(stdout,"      Names of data files that are read (input files), if not prefixed with either '/' or './' or '../',\n"); 
  fprintf(stdout,"      will be taken relative to the directory stored in variable %s.\n",PROG_DIR); 
  fprintf(stdout,"      For instance:\n"); 
  fprintf(stdout,"        - if a file is specified as '/home/jones/input.dat' it will be searched for in jones home directory.\n");
  fprintf(stdout,"        - if a file is specified as './input.dat' it will be searched for in the working directory, i.e., the one from which the program has been sent.\n");
  fprintf(stdout,"        - if a file is specified as '../input.dat' it will be searched in the parent working directory.\n");
  fprintf(stdout,"        - otherwise if a file is specified as SUBDIR1/SUBDIR2/input.dat it will be searched for in directory $PROG_DIR/SUBDIR1/SUBDIR2/\n\n");

  if((PROGdir = getenv(PROG_DIR)) != NULL) {
    fprintf(stdout,"Currently, the value of the environment variable %s is: %s\n\n",PROG_DIR,PROGdir);
  } else {
    fprintf(stdout,"The environment variable %s is currently not defined\n\n",PROG_DIR);
  }
}

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
void set_arg_value(char *arg_name, char *new_value);
char *get_arg_value(char *arg_name);
void set_output_levels();
int is_digit(char *str);

void process_line_args(int argc, char *argv[])
{

/*
 * This function processes the arguments read on the command line
 */

  int i, j;
  int found;

/*
 * First read default argument values in $PROG_DIR/DATA/DefaultArguments.
 * This file defines the valid arguments, sets the corresponding tags,
 * provides the default values and gives a brief description of the arguments.
 */

  initDefaultArguments();
  
/*
 * If no argument is provided print the help
 */

  if(argc == 1) {
    print_usage(argv[0]);
    exit(1);
  }

 /*
  * Check that the command line has the correct syntax
  */

  if((argc-1) % 2 != 0) {
    fprintf(stderr,"The syntax of the command line should be: a tag followed by an argument, e.g., -tg1 foo1 -tg2 foo2\n");
    print_usage(argv[0]);
    exit(1);
  }

  for(i = 1; i < argc; i += 2) {
    if(argv[i][0] != '-' || argv[i+1][0] == '-') {
      if(!is_digit(argv[i+1])) {
	fprintf(stderr,"Invalid syntax >>> %s %s <<<\n",argv[i],argv[i+1]);
	fprintf(stderr,"The syntax of the command line should be: a tag followed by an argument, e.g., -tg1 foo1 -tg2 foo2\n");
	print_usage(argv[0]);
	exit(1);
      }
    }
  }

/*
 * Process the command line arguments
 */

  for(i = 1; i < argc; i += 2) {
    found = 0;
    for(j = 0; j < maxargs; j++) {
      if(strcmp(argv[i],DefaultArguments[j][1]) == 0) {
	strcpy(DefaultArguments[j][3],argv[i+1]);
	found = 1;
	break;
      }
    }
    if(!found) {
      fprintf(stderr,"Unknown tag: %s!!\n\n",argv[i]);
      print_usage(argv[0]);
      exit(1);
    }
  }

/*
 * Check that compulsory arguments have been specified
 */

  int nerr = 0;
  for(j = 0; j < maxargs; j++) {
    if(strcmp(DefaultArguments[j][3],"compulsory") == 0) {
      PRINT_TAG_ERR;
      fprintf(stderr,"argument %s is compulsory!\n",DefaultArguments[j][0]);
      nerr++;
    }
  }
  if(nerr != 0) {
    exit(1);
  }
  

/*
 * Set output variables (log and debug)
 */

  set_output_levels();

/*
 * Possibly, process the arguments if required... 
 */

}

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
char *get_arg_value(char *arg_name)
{
  /*
   * Return the default value for argument arg_name
   *
   * Warning! This function returns the argument value as a character string.
   * It is the programmer's responsability to cast it to the appropriate type
   * as needed. For instance if the argument corresponds to a float then this
   * function can be called as foo = atof(get_arg_value(foo_name)) where 
   * variable foo is of type float in the calling function.
   */

  int i;

  for(i = 0; i < maxargs; i++) {
    if(strcmp(arg_name,DefaultArguments[i][0]) == 0) {
      break;
    }
  }

  if(i == maxargs) {
    PRINT_TAG_ERR;
    fprintf(stderr,"Unknown argument name %s\n",arg_name);
    exit(1);
  }

  return(DefaultArguments[i][3]);

}

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
void set_arg_value(char *arg_name, char *new_value)
{
  /*
   * Set the default value to new_value for argument arg_name in DefaultArguments
   */

  int i;

  for(i = 0; i < maxargs; i++) {
    if(strcmp(arg_name,DefaultArguments[i][0]) == 0) {
      strcpy(DefaultArguments[i][3],new_value);
      return;
    }
  }

  if(i == maxargs) {
    PRINT_TAG_ERR;
    fprintf(stderr,"Unknown argument name %s\n",arg_name);
    exit(1);
  }

}

/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
/*
 * Global variables to this part of the file
 */
int prt_log = 0;
int prt_dbg = 0;
FILE *fplog = NULL;
FILE *fpdbg = NULL;

void set_output_levels()
{
  /*
   * Set the print levels for the program (debug and log outputs)
   */

  int i;

  /*
   * Set log output variables
   */

  for(i = 0; i < maxargs; i++) {
    if(strcmp("LOG_LVL",DefaultArguments[i][0]) == 0) {
      if(strcmp(DefaultArguments[i][3],"void") != 0) {
	if(strcmp(DefaultArguments[i][3],"stdout") == 0) {
	  fplog = stdout;
	} else if (strcmp(DefaultArguments[i][3],"stderr") == 0) {
	  fplog = stderr;
	} else {
	  FOPEN(fplog,DefaultArguments[i][3],"w");
	}
	prt_log = 1;
	break;
      } else {
	break;
      }
    }
  }

  if(i == maxargs) {
    PRINT_TAG_ERR;
    fprintf(stderr,"Argument 'LOG_LVL' is required by the function.\n");
    exit(1);
  }

  /*
   * Set debug output variables
   */

  for(i = 0; i < maxargs; i++) {
    if(strcmp("DBG_LVL",DefaultArguments[i][0]) == 0) {
      if(strcmp(DefaultArguments[i][3],"void") != 0) {
	if(strcmp(DefaultArguments[i][3],"stdout") == 0) {
	  fpdbg = stdout;
	} else if (strcmp(DefaultArguments[i][3],"stderr") == 0) {
	  fpdbg = stderr;
	} else {
	  FOPEN(fpdbg,DefaultArguments[i][3],"w");
	}
	prt_dbg = 1;
	break;
      } else {
	break;
      }
    }
  }

  if(i == maxargs) {
    PRINT_TAG_ERR;
    fprintf(stderr,"Argument 'DBG_LVL' is required by the function\n");
    exit(1);
  }

}
/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
int prt_lvl(char *type)
{

  if(strcmp(type,"LOG") == 0) {
    return(prt_log);
  } else if(strcmp(type,"DBG") == 0) {
    return(prt_dbg);
  } else {
    PRINT_TAG_ERR;
    fprintf(stderr,"Invalid argument: type=%s\n",type);
    exit(1);
  }

}
/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
FILE *fp_ptr(char *type)
{

  if(strcmp(type,"LOG") == 0) {
    return(fplog);
  } else if(strcmp(type,"DBG") == 0) {
    return(fpdbg);
  } else {
    PRINT_TAG_ERR;
    fprintf(stderr,"Invalid argument: type=%s\n",type);
    exit(1);
  }

}
/*==============================================================================+
 |                                                                              |
 +==============================================================================*/
int is_digit(char *str)
{
  int i;

  for(i = 0; i < strlen(str); i++) {
    if(str[i] != '.' && str[i] != '-' && ! isdigit(str[i])) {
      return(0);
    }
  }

  return(1);

}
