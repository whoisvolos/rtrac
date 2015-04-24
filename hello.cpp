#include <omp.h>
#include <stdio.h>
 
int main(void)
{
  int th_id;
  #pragma omp parallel
    th_id = omp_get_thread_num();
    printf("Hello world from thread %i.\n", th_id);
  return 0;
}
