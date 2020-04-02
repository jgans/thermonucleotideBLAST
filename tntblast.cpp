// ThermonucleotideBLAST
// 
// Copyright (c) 2007, Los Alamos National Security, LLC
// All rights reserved.
// 
// Copyright 2007. Los Alamos National Security, LLC. This software was produced under U.S. Government 
// contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
// National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, 
// LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
// If software is modified to produce derivative works, such modified software should be clearly marked, 
// so as not to confuse it with the version available from LANL.
// 
// Additionally, redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 
//      * Redistributions of source code must retain the above copyright notice, this list of conditions 
//        and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
//        and the following disclaimer in the documentation and/or other materials provided with the distribution.
//      * Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, 
//        the U.S. Government, nor the names of its contributors may be used to endorse or promote products 
//        derived from this software without specific prior written permission.
// 
// 
// THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC 
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
		#endif // USE_OPENMP

		return local_main(argc, argv);
	}

	if(mpi_rank == 0){ // Master
		ret_value = master(argc, argv);
	}
	else{ // Worker
		ret_value = worker(argc, argv);
	}

	MPI_Finalize();
	#else // USE_MPI not defined
	
	#ifdef _OPENMP
	cout << "Running on local machine [" << omp_get_num_threads() << " thread(s)]" << endl;
	#else
	cout << "Running on local machine (1 thread)" << endl;
	#endif // USE_OPENMP
	
	ret_value = local_main(argc, argv);
	#endif //USE_MPI
	
	return ret_value;
}


