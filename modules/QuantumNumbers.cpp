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
bool compare_mom_dis_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3))
    return true;
  else
    return false;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
static void copy_quantum_numbers(const pdg& in, std::array<int, 6>& out){
  out[0] = in.dis3[0];
  out[1] = in.dis3[1];
  out[2] = in.dis3[2];
  out[3] = in.p3[0];
  out[4] = in.p3[1]; 
  out[5] = in.p3[2];
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
    for(const auto& mom_vec : in2.mom_vec){
      for(auto mom : mom_vec){
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
void QuantumNumbers::set_index_corr(){

  // first number is the operator id, the next three is the displacement vector
  // and the last three the momentum vector
  std::vector<std::array<int, 6> > rvdaggervr_qu_nb;
  std::vector<std::array<int, 6> > vdaggerv_qu_nb;
  size_t counter_rvdvr = 0;
  size_t counter_vdv = 0;
  for(auto& op : lookup_corr){
    std::array<int, 6> write;
    if(op.id != 0){
      copy_quantum_numbers(op, write);
      // ######################################################################
      // check if quantum numbers are already stored in rvdaggervr_qu_nb
      bool is_known_rvdvr = false;
      size_t fast_counter_rvdvr = 0;// this gives the Op id if QN are duplicate
      for(const auto& rvdvr : rvdaggervr_qu_nb){
        if(rvdvr == write){
          is_known_rvdvr = true;
          break;
        }
        fast_counter_rvdvr++;
      }
      if(!is_known_rvdvr){ // setting the unknown quantum numbers
        op.id_rVdaggerVr = counter_rvdvr;
        counter_rvdvr++;
        rvdaggervr_qu_nb.push_back(write);
      }
      else
        op.id_rVdaggerVr = fast_counter_rvdvr;
      // ######################################################################
      // check if quantum numbers are already stored in vdaggerv_qu_nb
      bool is_known_vdv = false;
      size_t fast_counter_vdv = 0;// this gives the Op id if QN are duplicate
      // first check for duplicate quantum numbers
      for(const auto& vdv : vdaggerv_qu_nb){
        if(vdv == write){
          is_known_vdv = true;
          break;
        }
        fast_counter_vdv++;
      }
      if(!is_known_vdv){ // second check for complex conjugate momenta
        fast_counter_vdv = 0;
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;
        for(const auto& vdv : vdaggerv_qu_nb){
          if(vdv == write){
            is_known_vdv = true;
            break;
          }
          fast_counter_vdv++;
        }
        if(!is_known_vdv){
          op.id_VdaggerV = counter_vdv;
          vdaggerv_qu_nb.push_back(write);
          counter_vdv++;
        }
        else{
          op.flag_VdaggerV = -1;
          op.id_VdaggerV = fast_counter_vdv;
        }
      }
      else{
        op.flag_VdaggerV = 1;
        op.id_VdaggerV = fast_counter_vdv;
      }
    }
    else{ // setting the very first entry
      copy_quantum_numbers(op, write);
      rvdaggervr_qu_nb.push_back(write);
      vdaggerv_qu_nb.push_back(write);
      op.id_VdaggerV = counter_vdv;
      op.id_rVdaggerVr = counter_rvdvr;
      counter_rvdvr++;
      counter_vdv++;
    }
  }

  // setting the lookuptables to be able to reconstruct the quantum numbers
  // when computing VdaggerV and rVdaggerVr
  lookup_vdv.resize(vdaggerv_qu_nb.size());
  lookup_rvdvr.resize(rvdaggervr_qu_nb.size());

  size_t index = 0;
  for(auto& op_vdv : lookup_vdv){
    op_vdv.id = index;
    for(const auto& op : lookup_corr){
      if(index == op.id_VdaggerV)
        op_vdv.index = op.id;
    }
    index++;
  }
  index = 0;
  for(auto& op_rvdvr : lookup_rvdvr){
    op_rvdvr.id = index;
    for(const auto& op : lookup_corr){
      if(index == op.id_VdaggerV){
        op_rvdvr.index = op.id;
        if(op.flag_VdaggerV == 1)
          op_rvdvr.adjoint = false;
        else
          op_rvdvr.adjoint = true;
      }
    }
    index++;
  }

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
void QuantumNumbers::init_lookup_corr() {

  const Correlator_list correlator_list = global_data->get_correlator_list();
  const std::vector<Operator_list> operator_list = 
    global_data->get_operator_list();

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
      for(const auto& mom_vec : individual_operator.mom_vec){
        for(auto mom : mom_vec){
          write.p3 = mom;
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
  // sorting lookup_corr for equal momentum and displacement vectors - makes it
  // easier to run over it with auto loops
  std::vector<pdg> dump_write;
  while(lookup_corr.size() != 0){

    it = lookup_corr.begin();
    dump_write.push_back(*it);
    lookup_corr.erase(it);

    auto it2 = dump_write.end()-1;
    while(it != lookup_corr.end()) {
      if(compare_mom_dis_of_pdg(*it, *it2)){
        dump_write.push_back(*it);
        lookup_corr.erase(it);
      }
      else
        it++;
    }
  }
  lookup_corr.swap(dump_write);

  // setting the identification numbers of lookup_corr
  size_t counter = 0;
  for(auto& op : lookup_corr)
    op.id = counter++;

  // final setting lookuptables for vdaggerv and so on
  set_index_corr();

  for(auto a : lookup_corr){
    std::cout << a.id;
    for(auto b : a.gamma)
    std::cout << " " << b;
    for(auto b : a.dis3)
      std::cout << " " << b;
    for(auto b : a.p3)
      std::cout << " " << b;
    std::cout << std::endl;
  }

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
    init_lookup_corr();    
    init_lookup_2pt();
    init_lookup_4pt();

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

