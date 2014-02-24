//#include <omp.h>
//#include <stdio.h>
//#include <stdlib.h>
//int main(int argc, char *argv[]) {
//
//	int nthreads, tid;
//
//	int yend = 30;
//	int xend = 30;
//	int ychange = 1;
//	int xchange = 1;
//	for (int y = 0; y < yend; y += ychange) {
//
//		#pragma omp parallel for
//
//		for (int x = 0; x < xend; x += xchange) {
//			printf("hello.\n");
//		}
//	}
//
////	/* Fork a team of threads giving them their own copies of variables */
////#pragma omp parallel private(nthreads, tid)
////	{
////
////		/* Obtain thread number */
////		tid = omp_get_thread_num();
////		printf("Hello World from thread = %d\n", tid);
////
////		/* Only master thread does this */
////		if (tid == 0)
////		{
////			nthreads = omp_get_num_threads();
////			printf("Number of threads = %d\n", nthreads);
////		}
////
////	} /* All threads join master thread and disband */
//
//}