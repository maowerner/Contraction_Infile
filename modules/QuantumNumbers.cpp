#include "QuantumNumbers.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

QuantumNumbers* QuantumNumbers::instance_ = 0;

QuantumNumbers* QuantumNumbers::Instance () {

  if(instance_ == 0) instance_ = new QuantumNumbers;

  return instance_;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// function that compares two pdg structs and checks if the corresponding 
// entries of lookup_corr coincide
bool compare_quantum_numbers_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3) && 
      (in1.gamma == in2.gamma))
    return true;
  else
    return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// function that compares a pdg struct with an operator as defined by the input
// file and checkes if the quantum numbers of pdg are contained in the physical
// situation described by the operator
bool compare_quantum_numbers_of_pdg(const pdg& in1, const Operators& in2){

  if( in1.gamma == in2.gammas &&
      in1.dis3 == in2.dil_vec){
    for(auto mom_sq : in2.mom_vec){
      for(auto mom : mom_sq){
        if(in1.p3 == mom){
          return true;
          //TODO: is that safer or just spam?
          break;
        }
      }
    }
  }

  return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
size_t get_nb_vdaggerv(const std::vector<pdg>& in){

  size_t counter = 0;

  auto it = in.begin();
  while(it != in.end()) {
    auto it2 = it;
    it2++;
    while(it2 != in.end()) {
      if(it->gamma != it2->gamma)
        counter++;
      it2++;
    }
    it++;
  }

  std::cout << counter << std::endl;
  return counter;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void QuantumNumbers::set_index_2pt(const Operators& in1, const Operators& in2) {

  index_2pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      write.index_Q2 = op1.id;
      write.index_Corr = op2.id;

      lookup_2pt.push_back(write);
    }} //loops over sink end here
  }} //loops over source end here

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void QuantumNumbers::set_index_4pt(const Operators& in1, const Operators& in2, 
                               const Operators& in3, const Operators& in4) {

  index_4pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      for(const auto& op3 : lookup_corr){
      if(compare_quantum_numbers_of_pdg(op3, in3)){
        for(const auto& op4 : lookup_corr){
        if(compare_quantum_numbers_of_pdg(op4, in4)){
          write.index_Q2[0]   = op1.id;
          write.index_Corr[0] = op2.id;
          write.index_Q2[1]   = op3.id;
          write.index_Corr[1] = op4.id;

          lookup_4pt.push_back(write);
        }}
      }} //loops over sink end here
    }}
  }} //loops over source end here

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void QuantumNumbers::init_correlators() {

  const Correlator_list correlator_list = global_data->get_correlator_list();
  const std::vector<Operator_list> operator_list = 
    global_data->get_operator_list();
  size_t i = 0;
  // extracting all operators which are used in correlations functions
  std::vector<int> used_operators;
  for(const auto& corr_list : correlator_list)
    used_operators.insert(used_operators.end(), 
                          corr_list.operator_numbers.begin(), 
                          corr_list.operator_numbers.end());
  
  sort(used_operators.begin(), used_operators.end());
  used_operators.erase(std::unique(used_operators.begin(), 
                                   used_operators.end()),
                       used_operators.end());
  // write quantum number in lookup_corr
  for(const auto& op_entry : used_operators){
    for(const auto& individual_operator : operator_list[op_entry]){
      pdg write;
      write.gamma = individual_operator.gammas;
      write.dis3 = individual_operator.dil_vec;
      for(auto mom_sq : individual_operator.mom_vec){
        for(auto mom : mom_sq){
          write.p3 = mom;
          //TODO: creates wrong numbers if operators are doubly counted
          write.id = i++;
          lookup_corr.push_back(write);
        }
      }
    }
  }
  // doubly counted lookup_corr entries are deleted
  auto it = lookup_corr.begin();
  while(it != lookup_corr.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_corr.end()) {
      if(compare_quantum_numbers_of_pdg(*it, *it2))
        lookup_corr.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // Test output for the time being TODO: can be deleted later
  for(const auto& a : lookup_corr){
    std::cout << a.id;
    for(auto b : a.gamma)
      std::cout << "\t" << b;
    for(auto b : a.dis3)
      std::cout << "\t" << b;
    for(auto b : a.p3)
      std::cout << "\t" << b;
    std::cout << std::endl;
  }

//    for(const auto& bla : lookup_corr)
//    std::cout << bla.gamma << std::endl;

//  // nb_op - number of combinations of three-momenta and gamma structures
//  // op    - vector of all three-momenta, three-displacements and gamma 
//  //         structure combinations
  const size_t nb_op = lookup_corr.size();
  std::cout << nb_op << std::endl;
  const size_t nb_VdaggerV = get_nb_vdaggerv(lookup_corr);
//  op_VdaggerV.resize(nb_VdaggerV);
//  const size_t nb_rVdaggerVr = nb_dis*nb_mom;
//  op_rVdaggerVr.resize(nb_rVdaggerVr);
//  set_Corr();
//
//  // nb_op_C2 - number of combinations of absolute values squared of momenta
//  //            and gamma-displacement combinations for 2pt-fct
//  // op_C2    - vector of all combinations for 2pt-fct and vector of 
//  //            op-index-pairs with all corresponding three-vectors and gammas
//  const size_t nb_op_C2 = nb_mom_sq * nb_dg * nb_dg;
//  op_C2.resize(nb_op_C2);
//  set_C2();
//
//  const size_t nb_op_C4 = nb_mom_sq * nb_mom_sq * nb_dg * nb_dg;
//  op_C4.resize(nb_op_C4);
//  set_C4();

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void QuantumNumbers::init_lookup_2pt() {

  const Correlator_list correlator_list = global_data->get_correlator_list();
  const std::vector<Operator_list> operator_list = 
    global_data->get_operator_list();

  // TODO: think about symmetry of switching index_Q2 and index_Corr
  // init lookup_2pt
  for(const auto& corr : correlator_list){

    // must be done for 2pt as well as 4pt function
    if( (corr.type.compare(0,3,"C2+") == 0) || 
        (corr.type.compare(0,5,"C4I2+") == 0) ){
      for(const auto& infile_op_so : operator_list[corr.operator_numbers[0]]){
      for(const auto& infile_op_si : operator_list[corr.operator_numbers[1]]){
        // TODO: give that guy the quark numbers from corr and think about 
        // giving the random numbers as well
        set_index_2pt(infile_op_so, infile_op_si);
      }} //loops over same physical situation end here
    } //case 2pt-function ends here

    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& infile_op_so : operator_list[corr.operator_numbers[2]]){
      for(const auto& infile_op_si : operator_list[corr.operator_numbers[3]]){
        set_index_2pt(infile_op_so, infile_op_si);
      }} //loops over same physical situation end here
    } //case 2pt-function ends here

  } // loop over correlator_list ends here

  // doubly counted lookup_index_C2 entries are deleted in order to not be
  // calculated twice
  auto it = lookup_2pt.begin();
  while(it != lookup_2pt.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_2pt.end()) {
      if( (it->index_Q2 == it2->index_Q2) && 
          (it->index_Corr == it2->index_Corr) )
        lookup_2pt.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // initialize the id's of lookup_2pt
  size_t counter = 0;
  for(auto& op : lookup_2pt){
    op.id = counter++;
  }

  std::cout << std::endl;
  for(const auto& a : lookup_2pt){
    std::cout << a.id << "\t" << a.index_Q2 << "\t" << a.index_Corr<< std::endl;
  }

  // init lookup_2pt_IO
  for(const auto& corr : correlator_list){

    if(corr.type.compare(0,3,"C2+") == 0){
      std::list<size_t> write;
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
        for(auto op : lookup_2pt){
          // TODO: change that to op_Q2[]
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2], 
                                            op1_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr], 
                                            op2_from_list)){

            write.emplace_back(op.id);
          }}
        }
      }}
      lookup_2pt_IO.push_back(write);
    }
  }

  std::cout << std::endl;
  for(const auto& a : lookup_2pt_IO){
    for(const auto& b : a)
      std::cout << b << std::endl;
    std::cout << std::endl;
  }

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void QuantumNumbers::init_lookup_4pt() {

  const Correlator_list correlator_list = global_data->get_correlator_list();
  const std::vector<Operator_list> operator_list = 
    global_data->get_operator_list();

  // initialization of op_C4
  
  for(auto& corr : correlator_list){
    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        set_index_4pt(op1_from_list, op2_from_list, op3_from_list, 
                      op4_from_list);
      }}}} //loops over same physical situation end here
    } //case 2pt-function ends here
  } // loop over correlator_list ends here

  std::cout << std::endl;
  for(auto a : lookup_4pt){
    std::cout << a.index_Q2[0] << "\t" << a.index_Corr[0] << "\t" 
              << a.index_Q2[1] << "\t" << a.index_Corr[1] << std::endl;
  }

}

void QuantumNumbers::init_from_infile() {

  try{
    init_correlators();    
    init_lookup_2pt();
    init_lookup_4pt();

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

