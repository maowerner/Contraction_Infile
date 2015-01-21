#include "global_data.h"
#include "global_data_utils.h"

namespace gdu = ::global_data_utils;

namespace {

// need some functions from namespace global_data_utils
using gdu::compare_quantum_numbers_of_pdg;
using gdu::compare_mom_dis_of_pdg;
using gdu::set_index_corr;
using gdu::set_index_2pt;
using gdu::set_index_4pt;

void init_lookup_corr(const Correlator_list& correlator_list, 
                      const std::vector<Operator_list>& operator_list,
                      vec_pdg_Corr& lookup_corr, vec_pd_VdaggerV& lookup_vdv,
                      vec_pd_rVdaggerVr& lookup_rvdvr) {

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
  set_index_corr(lookup_corr, lookup_vdv, lookup_rvdvr);

  for(auto a : lookup_corr){
    std::cout << a.id;
    for(auto b : a.gamma)
    std::cout << " " << b;
    for(auto b : a.dis3)
      std::cout << " " << b;
    for(auto b : a.p3)
      std::cout << " " << b;
    std::cout << "\t" << a.id_vdv;
    std::cout << " " << a.first_vdv;
    std::cout << " " << a.negative_momentum;
    std::cout << " " << a.id_rvdvr;
    std::cout << std::endl;
  }

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_2pt(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     const vec_pdg_Corr& lookup_corr, vec_index_2pt& lookup_2pt,
                     std::vector<indexlist_1>& lookup_2pt_IO,
                     std::vector<indexlist_2>& lookup_4pt_1_IO, 
                     std::vector<indexlist_2>& lookup_4pt_2_IO) {

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
        set_index_2pt(infile_op_so, infile_op_si, lookup_corr, lookup_2pt);
      }} //loops over same physical situation end here
    } //case 2pt-function ends here

    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& infile_op_so : operator_list[corr.operator_numbers[2]]){
      for(const auto& infile_op_si : operator_list[corr.operator_numbers[3]]){
        set_index_2pt(infile_op_so, infile_op_si, lookup_corr, lookup_2pt);
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
          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
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

  // init lookup_4pt_IO
  for(const auto& corr : correlator_list){

    if(corr.type.compare(0,5,"C4I2+") == 0){
      std::list<std::pair<size_t, size_t> > write;
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        for(auto op : lookup_2pt){
          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2], 
                                            op1_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr], 
                                            op2_from_list)){
          for(auto op2 : lookup_2pt){
            if(compare_quantum_numbers_of_pdg(lookup_corr[op2.index_Q2], 
                                              op3_from_list)){
            if(compare_quantum_numbers_of_pdg(lookup_corr[op2.index_Corr], 
                                              op4_from_list)){

              write.emplace_back(std::pair<size_t, size_t>(op.id, op2.id));
            }}
          }
          }}
        }
      }}}}
      lookup_4pt_1_IO.push_back(write);
      lookup_4pt_2_IO.push_back(write);
    }
  }

  std::cout << std::endl;
  for(const auto& a : lookup_4pt_1_IO){
    for(const auto& b : a)
      std::cout << b.first << "\t" << b.second << std::endl;
    std::cout << std::endl;
  }

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void init_lookup_4pt(const Correlator_list& correlator_list, 
                     const std::vector<Operator_list>& operator_list,
                     const vec_pdg_Corr& lookup_corr, vec_index_4pt& lookup_4pt,
                     std::vector<indexlist_1>& lookup_4pt_3_IO) {

  // initialization of op_C4
  
  for(auto& corr : correlator_list){
    // must be down in addition for 4pt function
    if(corr.type.compare(0,5,"C4I2+") == 0){
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        set_index_4pt(op1_from_list, op2_from_list, op3_from_list, 
                      op4_from_list, lookup_corr, lookup_4pt);
      }}}} //loops over same physical situation end here
    } //case 4pt-function ends here
  } // loop over correlator_list ends here

  // doubly counted lookup_index_C4 entries are deleted in order to not be
  // calculated twice
  auto it = lookup_4pt.begin();
  while(it != lookup_4pt.end()) {
    auto it2 = it;
    it2++;
    while(it2 != lookup_4pt.end()) {
      if( (it->index_Q2[0] == it2->index_Q2[0]) && 
          (it->index_Corr[0] == it2->index_Corr[0]) &&
          (it->index_Q2[1] == it2->index_Q2[1]) && 
          (it->index_Corr[1] == it2->index_Corr[1]) )
        lookup_4pt.erase(it2);
      else
        it2++;
    }
    it++;
  }

  // initialize the id's of lookup_2pt
  size_t counter = 0;
  for(auto& op : lookup_4pt){
    op.id = counter++;
  }

  std::cout << std::endl;
  for(auto a : lookup_4pt){
    std::cout << a.id << "\t" << a.index_Q2[0] << "\t" << a.index_Corr[0] 
              << "\t" << a.index_Q2[1] << "\t" << a.index_Corr[1] << std::endl;
  }

  // init lookup_4pt_IO
  for(const auto& corr : correlator_list){

    if(corr.type.compare(0,5,"C4I2+") == 0){
      std::list<size_t> write;
      for(const auto& op1_from_list : operator_list[corr.operator_numbers[0]]){
      for(const auto& op2_from_list : operator_list[corr.operator_numbers[1]]){
      for(const auto& op3_from_list : operator_list[corr.operator_numbers[2]]){
      for(const auto& op4_from_list : operator_list[corr.operator_numbers[3]]){
        for(auto op : lookup_4pt){
          // TODO: change lookup_corr[op.index_Q2] to lookup_Q2[op.index_Q2]
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2[0]], 
                                            op1_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr[0]], 
                                            op2_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Q2[1]], 
                                            op3_from_list)){
          if(compare_quantum_numbers_of_pdg(lookup_corr[op.index_Corr[1]], 
                                            op4_from_list)){
            write.emplace_back(op.id);
          }}}}
        }
      }}}}
      lookup_4pt_3_IO.push_back(write);
    }
  }

  std::cout << std::endl;
  for(const auto& a : lookup_4pt_3_IO){
    for(const auto& b : a)
      std::cout << b << std::endl;
    std::cout << std::endl;
  }

}


} // end of unnamed namespace

void GlobalData::init_lookup_tables() {

  try{
    init_lookup_corr(correlator_list, operator_list, lookup_corr, lookup_vdv, 
                     lookup_rvdvr);    
    init_lookup_2pt(correlator_list, operator_list, lookup_corr, lookup_2pt,
                    lookup_2pt_IO, lookup_4pt_1_IO, lookup_4pt_2_IO);
    init_lookup_4pt(correlator_list, operator_list, lookup_corr, lookup_4pt,
                    lookup_4pt_3_IO);

  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}

