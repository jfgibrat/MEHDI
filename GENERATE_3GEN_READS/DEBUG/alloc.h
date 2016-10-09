#define BUFSIZE 70000

/*! \def S(x)
 * Stringize operator # : returns a character string corresponding to
 * the name of a variable
 */

#define S(x)      #x

#define EXIT_FAILURE 1

#define PRINT_TAG_OUT\
  fprintf(stdout, "-> Warning <-\n"\
	    "%s: at line %d\n", __FILE__, __LINE__)

#define PRINT_TAG_ERR\
      fprintf(stderr, "=> ERROR <=\n"\
          "%s: at line %d\n", __FILE__, __LINE__)

/*! \def ALLOC3
 * Allocate a three dimensions array
 * It only uses three malloc calls.
 */
#define ALLOC3(ptr, nelem1, nelem2, nelem3) {\
  int alloc3_i;\
  CALLOC2((ptr), (nelem1), (nelem2))\
  ALLOC((ptr)[0][0], (nelem1) * (nelem2) * (nelem3))\
  for (alloc3_i = 1; alloc3_i < ((nelem1) * (nelem2)); alloc3_i++) {\
    (ptr)[0][alloc3_i] = (ptr)[0][alloc3_i - 1] + (nelem3);\
  }\
}

/*! \def CALLOC
 * Same as ALLOC but for CALLOC
 */
#define CALLOC(ptr, nelem) {\
  if ((ptr) != NULL) {\
    PRINT_TAG_OUT;\
    fprintf(stdout, "Allocating (calloc) a non NULL pointer!\n");\
  }\
  (ptr) = calloc((size_t) (nelem), (size_t) sizeof(*ptr));\
  if ((ptr) == NULL) {\
    PRINT_TAG_ERR;\
    fprintf(stderr, "Unable to allocate (calloc) memory (%d bytes) for %s\n",\
        (int) (nelem) * (int) sizeof(*ptr), S(ptr));\
    exit(EXIT_FAILURE);\
  }\
}

/*! \def CALLOC2
 * Same as ALLOC2 but for CALLOC
 */
#define CALLOC2(ptr, nelem1, nelem2) {\
  int calloc2_i;\
  CALLOC((ptr), (nelem1))\
  CALLOC((ptr)[0], (nelem1) * (nelem2))\
  for (calloc2_i = 1; calloc2_i < (nelem1); calloc2_i++) {\
    (ptr)[calloc2_i] = (ptr)[calloc2_i - 1] + (nelem2);\
  }\
}
/*! \def ALLOC
 * malloc with return value check only
 */
#define ALLOC(ptr, nelem) {\
  (ptr) = malloc((size_t) (nelem) * (size_t) sizeof(*ptr));\
  if ((ptr) == NULL) {\
    PRINT_TAG_ERR;\
    fprintf(stderr, "Unable to allocate (malloc) memory (%d bytes) for %s\n",\
        (int) (nelem) * (int) sizeof(*ptr), S(ptr));\
    exit(EXIT_FAILURE);\
  }\
}

/*! \def FOPEN
 * Redefine fopen to include error check
 */

#define FOPEN(fp,fich,mode) { fp = fopen(fich,mode); \
    if(fp == NULL) { fprintf(stderr,"Unable to open file %s\n",fich); exit(1);} }

/*! \def FREE
 * Redefine free to include warning check
 */

#define FREE(p) { if(p == NULL) { fprintf(stderr, \
        "Trying to free a NULL pointer: %s\n", S(p)); } \
        else { free(p); p = NULL; } }

/*! \def FREE2
 * Free an array allocated with CALLOC2
 */
#define FREE2(ptr) {\
  if ((ptr) == NULL) {\
    PRINT_TAG_OUT;\
    fprintf(stdout, "Trying to free a NULL pointer: %s\n", S(ptr));\
  } else {\
    FREE((ptr)[0])\
    FREE(ptr)\
  }\
}

/*! \def FREE3
 * Free an array allocated with ALLOC3
 */
#define FREE3(ptr) {\
  if ((ptr) == NULL) {\
    PRINT_TAG_OUT;\
    fprintf(stdout, "Trying to free a NULL pointer: %s\n", S(ptr));\
  } else {\
    FREE((ptr)[0][0])\
    FREE((ptr)[0])\
    FREE(ptr)\
  }\
}

/*! \def CLOSE
 * Redefine fclose to include warning check
 */

#define FCLOSE(fp) { if(fp == NULL) { fprintf(stderr, \
        "Trying to close a NULL file pointer: %s\n", S(fp)); } \
        else { fclose(fp); fp = NULL; } }

