#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "tntblast.h"

#include <stdlib.h>
#include <iostream>

#ifdef _OPENMP
// Under windows, we need to include omp.h to load
// vcomp.dll (which is required for openMP on windows)
#include <omp.h>
#endif // _OPENMP

using namespace std;

// Global variables
int mpi_numtasks;
int mpi_rank;

#ifdef PROFILE
unsigned int num_plus_tm_eval = 0;
unsigned int num_minus_tm_eval = 0;
#endif // PROFILE


int main(int argc, char *argv[])
{	
	int ret_value = EXIT_FAILURE;
	
	#ifdef USE_MPI

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	if(mpi_numtasks < 2){

		#ifdef _OPENMP
		cout << "Running on local machine [" << omp_get_max_threads() << " thread(s)]" << endl;
		#else
		cout << "Running on local machine (1 thread)" << endl;
		#endif // _OPENMP

		ret_value = local_main(argc, argv);
	}
	else{
		if(mpi_rank == 0){ // Master
			ret_value = master(argc, argv);
		}
		else{ // Worker
			ret_value = worker(argc, argv);
		}
	}
	
	MPI_Finalize();
	#else // USE_MPI not defined
	
	#ifdef _OPENMP
	cout << "Running on local machine [" << omp_get_max_threads() << " thread(s)]" << endl;
	#else
	cout << "Running on local machine (1 thread)" << endl;
	#endif // _OPENMP
	
	ret_value = local_main(argc, argv);
	#endif //USE_MPI
	
	return ret_value;
}


