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
#include "xSimpleGeometry.h"
using namespace Trellis_Util;
using namespace xfem;
using namespace AOMD;
using namespace xtensor;
using namespace xgeom;


class helmholtzNdbiot
{
public:
    //constructor
    helmholtzNdbiot (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos);
    //function for constructing the FEM form
    void TreatmentOfFormulation (xParseData &parsedinfos, xTranslate s_f) ;  // xUnion for two seats problems

    template < typename FIELD >
    void TreatmentOfEssEnv   (const FIELD& listFunctionSpace, xfem::xData &data);

    template <class ASSEMBLER, class FIELD >
    void  TreatmentOfMixedBCs(xfem::xData &data, ASSEMBLER &assembler, const xfem::xIntegrationRule& integration_rule,
                              const FIELD& press_l, double omega);

    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnvAcoustic(xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                            xfem::xEval<xtensor::xVector<std::complex<double> > > &exact_velocity, 
                            double omega_sq, double beta);

    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnvAcoustic(xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                            xfem::xEval<std::complex<double> > &exact_velocity, 
                            double omega_sq);

    
    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnvPEMFluid(const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule, 
                           xfem::xEval<xtensor::xVector<std::complex<double> > > &disp_total);

    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnvPEMSolid(const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule, 
                           xfem::xEval<xtensor::xTensor2<std::complex<double> > > &stress);


    void setLevelSet(xfem::xLevelSet &lset_){lset.reset(&lset_);}

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    //Parallel stuff
    string pid;
    int proc_id,nb_proc;

    // Xfiles members
    xfem::xData data;

    //unique_ptr are used in order to avoid to bother with memory release
    std::unique_ptr<xfem::xLevelSet> lset;
    std::shared_ptr<xcut::xPhysSurfByTagging> inclusion;

};

template<typename T>
class xEvalInvDensity : public xfem::xEval<T> {

public :
    xEvalInvDensity(xfem::xEval<T> &ev_inv_dens_) : ev_inv_dens(ev_inv_dens_) {}
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, T& res) const {

        throw;
        cout<<"Hello !\n";

        ev_inv_dens(geo_appro, geo_integ, res);

       auto xyz= geo_integ->getXYZ();
       cout<<"Pos "<<xyz(0)<<" "<<xyz(1)<<" inv_dens "<<res<<endl;
       //cout<<"inv_dens "<<res<<endl;
        return;
    }

private:
    xfem::xEval<T> &ev_inv_dens;
};


class xEvalInvDensity__ : public xfem::xEval<std::complex<double> > {

public :
    xEvalInvDensity__(xfem::xEval<std::complex<double> > &ev_inv_dens_) : ev_inv_dens(ev_inv_dens_) {}
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double> & res) const {

        //throw;
        cout<<"Hello !\n";

        ev_inv_dens(geo_appro, geo_integ, res);

       auto xyz= geo_integ->getXYZ();
       cout<<"Pos "<<xyz(0)<<" "<<xyz(1)<<" inv_dens "<<res<<endl;
       //cout<<"inv_dens "<<res<<endl;
        return;
    }

private:
    xfem::xEval<std::complex<double> > &ev_inv_dens;
};




template<typename T>
class xEvalJumpB : public xfem::xEval<T> {

public :
    xEvalJumpB(xfem::xEval<T> &ev_field_, xfem::xField<std::complex<double>> *field_=nullptr) : ev_field(ev_field_), 
    field(field_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag")) {}

    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, T& res) const{

         //Positive side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        auto resPos = res * 0.;
        ev_field(geo_appro, geo_integ, resPos);
        // cout<<"eval_jump+ : "<<resPos<<endl;


         //Negative side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        auto resNeg = res * 0.;
        if(field) field->resetStorage(field->beginFcts(geo_appro->getEntity()), field->endFcts(geo_appro->getEntity()), geo_appro, geo_integ);
        ev_field(geo_appro, geo_integ, resNeg);
        // cout<<"eval_jump- : "<<resNeg<<endl;

         //Jump
        res = resNeg - resPos;
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


// class xUpperCreatorOnMeshC {
// public:
//   xUpperCreatorOnMeshC(xMesh *mC_) {mC=mC_;}
//    AOMD::mEntity* operator()(AOMD::mEntity* e);

//   xMesh *mC;
// };




void lineExporterbt(xfem::xMesh &m, xtensor::xPoint p1, xtensor::xPoint p2, xfem::xEval<double> &eval, int nbPoints, string filename);



