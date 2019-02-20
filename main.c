#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "allvars.h"
#include "proto.h"


int main(int argc, char *argv[3])
{
  init();
  run();
//  finalize();
  acknowledge();

  return 0;
}

void acknowledge()
{
  T_end = time(NULL);
  printf("---------------------------------------\n");
  printf("Wall clock time lapsed: %.2lf sec\n", difftime(T_end, T_start));
  printf("LIFE IS BEAUTIFUL!!!\n");
}
