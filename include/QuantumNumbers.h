#ifndef QUANTUMNUMBERS_H_
#define QUANTUMNUMBERS_H_

#include <array>
#include <vector>

#include "GlobalData.h"
#include "Operators.h"
#include "typedefs.h"

class QuantumNumbers {

  private:

    static QuantumNumbers* instance_;

    void set_index_2pt(const Operators& in1, const Operators& in2);
    void set_index_4pt(const Operators& in1, const Operators& in2, const Operators& in3, 
                const Operators& in4);

    void init_correlators();
    void init_lookup_2pt();
    void init_lookup_4pt();

    vec_pdg_Corr lookup_corr;
    vec_index_2pt lookup_2pt;
    vec_index_4pt lookup_4pt;
    std::vector<std::list<size_t> > lookup_2pt_IO;
    std::vector<std::list<size_t> > lookup_4pt_IO;
    vec_pd_VdaggerV lookup_vdv;
    vec_pd_rVdaggerVr lookup_rvdvr;

  public:
    static QuantumNumbers* Instance ();
  
    void init_from_infile();

    inline const vec_pdg_Corr& get_lookup_corr() {
      return lookup_corr;
    }
    inline const vec_index_2pt& get_lookup_2pt_trace() {
      return lookup_2pt;
    }
    inline const vec_index_4pt& get_lookup_4pt_trace() {
      return lookup_4pt;
    }
    inline const std::vector<std::list<size_t> >& get_lookup_2pt_IO() {
      return lookup_2pt_IO;
    }
    inline const std::vector<std::list<size_t> >& get_lookup_4pt_IO() {
      return lookup_4pt_IO;
    }
//    inline const indexlist_2& get_rnd_vec_C2() {
//      return rnd_vec_C2;
//    }
//    inline const indexlist_4& get_rnd_vec_C4() {
//      return rnd_vec_C4;
//    }
//    inline const size_t get_index_of_unity() {
//      return index_of_unity;
//    }
    inline const vec_pd_VdaggerV get_lookup_vdaggerv() {
      return lookup_vdv;
    }
    inline const vec_pd_rVdaggerVr get_lookup_rvdaggervr() {
      return lookup_rvdvr;
    }

  protected:
    QuantumNumbers () {
    }
    QuantumNumbers (const QuantumNumbers& other) {
    }
    virtual ~QuantumNumbers () {
    }

};

#endif // QUANTUMNUMBERS_H_
