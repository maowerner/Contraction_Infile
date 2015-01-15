//============================================================================
// Name        : LapHs.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : stochastic LapH code
//============================================================================

#include "GlobalData.h"
#include "QuantumNumbers.h"

int main (int ac, char* av[]) {

  // reading in global parameters from input file
  QuantumNumbers* qns = QuantumNumbers::Instance();
  GlobalData* global_data = GlobalData::Instance();

  global_data->read_parameters(ac, av);
  qns->init_from_infile();

}

