#include "define_t.h"

void PLMsummarize_cpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/);
void PLMsummarize_gpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/);


int main (int argc, char** argv){

	struct timeval start, end;
	long utime;	
	
	//input data
	int_t numprobes = Y_ROWS;
	int_t numchips = Y_COLS;
	
	double* z_cpu = new double[numprobes * numchips]();
	double* z_gpu = new double[numprobes * numchips]();
	
	double* results_cpu = new double[numchips]();
	double* results_gpu = new double[numchips]();
	
	double* SEresults_cpu = new double[numchips]();
	double* SEresults_gpu = new double[numchips]();

	//init z
	for(int_t i = 0; i < numprobes; i++){
 		for(int_t j = 0; j < numchips; j++){
			z_cpu[i*numchips + j] = rand()%(RANDOM_MAX - RANDOM_MIN + 1) + RANDOM_MIN;	
			z_gpu[i*numchips + j] = z_cpu[i*numchips + j];
		}	
	}

#if 1

#ifdef CPU_COMPUTE
	//CPU computation
	gettimeofday(&start, NULL);

	PLMsummarize_cpu(z_cpu, numprobes, numchips, results_cpu, SEresults_cpu /*double* affinities*/);	
	
	gettimeofday(&end, NULL);
	utime = ((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);
	printf("\nTime PLMsummarize_cpu = %ld us\n",utime);	
#endif

	//GPU computation
	gettimeofday(&start, NULL);
	
	PLMsummarize_gpu(z_gpu, numprobes, numchips, results_gpu, SEresults_gpu /*double* affinities*/);

	gettimeofday(&end, NULL);
	utime = ((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);
	printf("\nTime PLMsummarize_gpu = %ld us\n",utime);	

#ifdef CPU_COMPUTE
	//check for z_err
	double z_err = 0.0;
	for(int_t j=0; j<(numprobes * numchips); j++){
		z_err+=fabs((z_cpu[j] - z_gpu[j])/z_cpu[j]);
	}
	z_err/=(numprobes * numchips);
	printf("\nError z_err = %e\n",z_err);

	//check for results_err
	double results_err = 0.0;
	for(int_t j=0; j<numchips; j++){
		z_err+=fabs((results_cpu[j] - results_gpu[j])/results_cpu[j]);
	}
	results_err/=numchips;
	printf("\nError results_err = %e\n",results_err);

	//check for SEresults_err
	double SEresults_err = 0.0;
	for(int_t j=0; j<numchips; j++){
		SEresults_err+=fabs((SEresults_cpu[j] - SEresults_gpu[j])/SEresults_cpu[j]);
	}
	SEresults_err/=numchips;
	printf("\nError SEresults_err = %e\n",SEresults_err);	
#endif

#else

	if(numchips<2000){
		PLMsummarize_cpu(z_cpu, numprobes, numchips, results_cpu, SEresults_cpu /*double* affinities*/);
	}
	else {
		PLMsummarize_gpu(z_gpu, numprobes, numchips, results_gpu, SEresults_gpu /*double* affinities*/);
	}
	
#endif	
	delete[] z_cpu;		delete[] z_gpu;
	delete[] results_cpu;	delete[] results_gpu;
	delete[] SEresults_cpu;	delete[] SEresults_gpu;

	return 0;
}


