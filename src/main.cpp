#include "struct.hpp"
#include "mpi_example.hpp"

#include "source.hpp"
#include "slave.hpp"
#include "run.hpp"
#include <iostream>
#include <unistd.h>

int main ( int argc, char *argv[] ){

    MPIIO 	mpiIo(argc, argv);
	int 	slaveNum = mpiIo.getSzProc()-1;


	Choice::parse(argc, argv);

	volatile int a = 0;

	if(mpiIo.isMaster()){
		Choice::print();
	}
	
	run_experiment(Choice::inFileName, Choice::outPath, mpiIo, slaveNum, Choice::method, Choice::budget, Choice::trial);
	
	mpiIo.cleanup();
}
