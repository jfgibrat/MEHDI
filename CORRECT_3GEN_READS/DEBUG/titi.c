#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

int main(int argc, char *argv[])
{

  char str[100];

  strcpy(str,argv[1]);
  fprintf(stdout,"%x\n",str);
  
  return 0;
}

