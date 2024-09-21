#ifndef RUN_EXPERIMENT_HPP
#define RUN_EXPERIMENT_HPP

#include "struct.hpp"
#include "source.hpp"
#include "slave.hpp"
#include <vector>
#include <iostream>
#include <unordered_map>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <limits>


void run_experiment (const char* input, const char* outPath, MPIIO &mpiIo, int slaveNum, Method_t method, int memSize, int repeat, int bufLen=1000, double tolerance=0.2);

double run_experiment_mpi(const char* filename, MPIIO &mpiIo, int slaveNum, const Method_t method, int memSize, int lenBuf, double tolerance, unsigned int seed, std::vector<float> & oLocalCount, double &srcCompCost, double &slaveCompCostMax, double &slaveCompCostSum);

void print_count(const char* outPath, double globalCount, const std::vector<float> &localCount, int id = 0);

#endif 