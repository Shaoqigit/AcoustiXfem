#pragma once
#include "xData.h"
#include "xMesh.h"
#include "xIntegrationRule.h"
#include "xParseData.h"
#include "xLevelSet.h"
#include "xGeomElem.h"
#include "xField.h"
#include "mPoint.h"
#include "mEdge.h"
#include "xAttachableGP.h"
#include "xAssembler.h"
#include "xSimpleGeometry.h"
#include <memory>

#include "xPhysSurfByTagging.h"
#include "xNitscheFile.h"
using namespace Trellis_Util;
using namespace xfem;
using namespace AOMD;
using namespace xtensor;


class helmholtzNd
{
public:
    //constructor
    helmholtzNd (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos);
    //function for constructing the FEM form
    void TreatmentOfFormulation (xParseData &parsedinfos) ;

    template <class FIELD , class ITER >
    void TreatmentOfEssEnv(const FIELD& fct, ITER itenv, ITER endenv,
                           const std::vector<string>& vars, xfem::xMesh* mesh);


    template <class ASSEMBLER, class FIELD >
    void  TreatmentOfMixedBCs(xfem::xData &data, ASSEMBLER &assembler, const xfem::xIntegrationRule& integration_rule,
                              const FIELD& press_l, double omega);

    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnv(const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                           xfem::xEval<std::complex<double> > &exact_velocity, 
                           std::complex<double> minusjrhock, std::complex<double> minusjomegaexp, double beta);
    
    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnv(const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                           xfem::xEval<xtensor::xVector<std::complex<double> > > &exact_velocity,
                           std::complex<double> minusjrhock, std::complex<double> minusjomegaexp, double beta);



    void setLevelSet(xfem::xLevelSet &lset_){
        if(lset) delete lset;
        lset = &lset_;
    }

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    //Parallel stuff
    string pid;
    int proc_id,nb_proc;

    // Xfiles members
    xfem::xData &data;

    //unique_ptr are used in order to avoid to bother with memory release
    //not for lset as it crashes while releasing the pointer
    xfem::xLevelSet *lset{nullptr};
    std::shared_ptr<xcut::xPhysSurfByTagging> inclusion;

};



template<typename T>
class xEvalInvDensity : public xfem::xEval<T> {

public :
    xEvalInvDensity(xfem::xEval<T> &ev_inv_dens_) : ev_inv_dens(ev_inv_dens_) {}
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, T& res) const{
        ev_inv_dens(geo_appro, geo_appro, res);

    //    cout<<"inv_dens "<<res<<endl;
        return;
    }

private:
    xfem::xEval<T> &ev_inv_dens;
};


template<typename T>
class xEvalJump : public xfem::xEval<T> {

public :
    xEvalJump(xfem::xEval<T> &ev_field_, xfem::xField<std::complex<double>> *field_=nullptr) : ev_field(ev_field_), 
    field(field_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag")) {}

    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, T& res) const{

         //Positive side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        auto resPos = res * 0.;
        ev_field(geo_appro, geo_integ, resPos);
        // cout<<"eval_jump+ : "<<resPos<<endl;


         //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        auto resNeg = res * 0.;
        if(field) field->resetStorage(field->beginFcts(geo_appro->getEntity()), field->endFcts(geo_appro->getEntity()), geo_appro, geo_integ);
        ev_field(geo_appro, geo_integ, resNeg);
        // cout<<"eval_jump- : "<<resNeg<<endl;

         //Jump
        res = resPos - resNeg;
        // res = 1;

         //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);


        // cout<<"eval_jump : "<<res<<endl;
        return;
     }

private:
    xfem::xEval<T> &ev_field;
    xfem::xField<T> *field;
    const unsigned int side_tag;
};


// evaluator for stablization parameter as a mater

// template <typename T>
// class  xEvalStabParam : public xEval<T> {

// using valtype = std::complex<double>;
// public:
// typedef typename xEval<T>::result_type result_type;
// xEvalStabParam(xMesh& mInterf_, xfem::xField<valtype> &field_, xfem::xLevelSet& lst, xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, const char * side_tag_name = "side_tag"):
// mInterf(mInterf_), eval_normal(lst), classNeg(classNeg_), classPos(classPos_), field(field_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {
//     Nitsche_param_XFEM getgamma(field, eval_normal, classNeg, classPos);
// }
// void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& gamma) const {
    
//     Nitsche_param_XFEM getgamma(field, eval_normal, classNeg, classPos);
//     // AOMD::mEntity* e = geo_appro->getEntity();
//     // getgamma.ComputeParameter(mInterf);
//     // gamma = getgamma.GetParameter(e);

//     return;
// }

// private:
// Nitsche_param_XFEM &getgamma
// xMesh& mInterf;
// xEvalNormalFronLevelSet eval_normal;
// xfem::xClassifyOn &classNeg;
// xfem::xClassifyOn &classPos;
// xfem::xField<valtype> &field;
// unsigned int side_tag;

// };


class xEvalPointSource: public xfem::xEval<std::complex<double> > {
    public : 
    xEvalPointSource(double omega_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;
private:
    double omega;
    xtensor::xPoint o;
};

void lineExporterFE(xMesh &m, xtensor::xPoint p1, xtensor::xPoint p2, xEval<double> &eval, int nbPoints, string filename);




