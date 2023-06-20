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
#include "xOperators.h"
#include "hoGeometricalModel_hoOctree.h"
#include <memory>

#include "xPhysSurfByTagging.h"
#include "nitscheOperators.h"
#include "xMaterialSensitivity.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace Trellis_Util;
using namespace xfem;
using namespace AOMD;
using namespace xtensor;
using namespace Eigen;
using namespace xgeom;


class helmholtzNdPSG
{
public:
    //constructor
    helmholtzNdPSG (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos);
    //function for constructing the FEM form
    void TreatmentOfFormulation (xParseData &parsedinfos, xMinus m2) ;

    template <class FIELD , class ITER >
    void TreatmentOfEssEnv(const FIELD& fct, ITER itenv, ITER endenv,
                           const std::vector<string>& vars, xfem::xMesh* mesh);


    template <class ASSEMBLER, class FIELD >
    void  TreatmentOfMixedBCs(xfem::xData &data, ASSEMBLER &assembler, const xfem::xIntegrationRule& integration_rule,
                              const FIELD& press_l, double omega);
    
    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnv(xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                           xfem::xEval<xtensor::xVector<std::complex<double> > > &exact_velocity);

    template < typename ASSEMBLER, typename FIELD >
    void TreatmentOfNatEnv(xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                           xfem::xEval<xtensor::xVector<std::complex<double> > > &exact_velocity, double omega_sq);


    void setLevelSet(xfem::xLevelSet &lset_){
        if(lset) delete lset;
        lset = &lset_;
        // void setLevelSet(xfem::xLevelSet &lset_){lset.reset(&lset_);
    }

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    //Parallel stuff
    string pid;
    int proc_id,nb_proc;

    // Xfiles members
    xfem::xData &data;
    // hoGeometricalModel_hoOctree &geoModel;

    //unique_ptr are used in order to avoid to bother with memory release
    //not for lset as it crashes while releasing the pointer
    xfem::xLevelSet *lset{nullptr};
    std::shared_ptr<xcut::xPhysSurfByTagging> inclusion;

    // std::unique_ptr<xfem::xLevelSet> lset;
    // std::shared_ptr<xcut::xPhysSurfByTagging> inclusion;

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




// double ComputeParameter(AOMD::mEntity* e, mEntity* ed, xfem::xField<std::complex<double>> &field,
//                         const xfem::xEval<xtensor::xVector<> >& eval_normal_, xfem::xClassifyOn &classNeg, xfem::xClassifyOn &classPos)
// {

//     using valtype=std::complex<double>; 
//     xUniformMaterialSensitivity<valtype > fluid_inv_density("inv_density");
//     xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
//             xEval<valtype>,
//             xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > B_bilin(fluid_inv_density);  // result type is valtype
    

//     double gamma = 0.5;
//     xComputeNormalVelocity opNormalvelocity(eval_normal_, fluid_inv_density, classNeg, classPos, 1.);
//     xGradMeanOperatorWithParamInterConj<xComputeNormalVelocity> opMeanNormalVelocity_cong(gamma, opNormalvelocity, true);
//     xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(gamma, opNormalvelocity);
//     xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
//             xGradMeanOperatorWithParamInterConj<xComputeNormalVelocity>, valtype > A_bilin(opMeanNormalVelocity, opMeanNormalVelocity_cong); 


//     xIntegrationRulePartition integration_rule_neg(xMesh::getPartition, xAccept("material_neg"), 20);
//     xIntegrationRulePartition integration_rule_pos(xMesh::getPartition, xAccept("material_pos"), 20);
//     xIntegrationRuleBasic integration_rule_(20);
//     xIntegrateFormCommand commandA(&A_bilin);
//     AOMD::mEntity* e_integ = e;
//     AOMD::mEntity* e_appro = e_integ;
//     xIntegrateFormCommand commandB(&B_bilin);
//     AOMD::mEntity* ed_integ = ed;
//     AOMD::mEntity* ed_appro = e_appro;
//     xfem::femFcts_t trialfcts(field.beginFcts(e_appro), field.endFcts(e_appro));
//     if (!trialfcts.empty())
//          {
//         //    negtive side
//             xGeomElem geo_approA(e_appro);
//             A_bilin.init(&geo_approA, &trialfcts);
//             integration_rule_.accept(commandA, ed);
//             A_bilin.finalize();
//             auto A = A_bilin.getMatrix();

//             xGeomElem geo_approB(e_appro);
//             B_bilin.init(&geo_approB, &trialfcts);
//             integration_rule_.accept(commandB, e_appro);
//             B_bilin.finalize();
//             auto B = B_bilin.getMatrix();

//             Map<Matrix<valtype, Dynamic, Dynamic, RowMajor>> MapA(A.ReturnValPointer(), A.getNbRow(), A.getNbCol()); 
//             Map<Matrix<valtype, Dynamic, Dynamic, RowMajor>> MapB(B.ReturnValPointer(), B.getNbRow(), B.getNbCol()); 

            
//             MatrixXd absA = MapA.real();
//             MatrixXd absB = MapB.real();
//             // cout<<"matA:"<<rA<<endl;
//             // cout<<"matB:"<<rB<<endl;
//             std::ofstream out("matA_eigen.txt",ios::app);
//             out<<"matA:"<<absA<<endl;
//             out.close();
//             std::ofstream outB("matB_eigen.txt",ios::app);
//             outB<<"matB:"<<absB<<endl;
//             outB.close();
//             GeneralizedEigenSolver<MatrixXd> ges;
//             ges.compute(absA, absB);
//             VectorXcd  Lambda_set = ges.eigenvalues();
//             // double Lambda = Lambda_set.maxCoeff();
//             // cout<<"eigen value"<<Lambda<<endl;
        
//          }


//     return 0.;

// };


void lineExporter(xfem::xMesh &m, Trellis_Util::mPoint p1, Trellis_Util::mPoint p2, xfem::xEval<double> &eval, int nbPoints);