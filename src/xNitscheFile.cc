#pragma once

#include "xNitscheFile.h"
#include "xEval.h"
#include "xForm.h"
#include "xOperators.h"
#include "xApproxFctPtr.h"
#include "nitscheOperators.h"
#include "xMaterialSensitivity.h"

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "xGenericSparseMatrixTraitPolicy.h"
#include "Eigen/Eigenvalues"


using namespace Eigen;


// #include "xGenericSparseMatrix.h"
using namespace xfem;
using namespace AOMD;

// #define OCTREE


// take a xEval as enter for the class
Nitsche_param::Nitsche_param(xfem::xMesh* meshC_, xfem::xField<valtype> &field_,  xfem::xLevelSet& lst_, xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, const char * side_tag_name):
 uppercreatoronmc(meshC_), eval_normal(lst_), classNeg(classNeg_), classPos(classPos_), field(field_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {};

double Nitsche_param::GetStable(AOMD::mEntity* e)
{
    double stable = *datamanger_beta.getData(*e);
    return stable;
}

double Nitsche_param::GetWeighting(AOMD::mEntity* e)
{
    double gamma_1 = *datamanger_gamma.getData(*e);
    return gamma_1;
}


void Nitsche_param::ComputeParameter(xMesh& mInterf, valtype inv_rho) {

    for (mEntity* ed: mInterf.range(mInterf.dim())) {
        // xPartition partition;
        // xMesh::getPartition(e, partition);
#ifdef OCTREE
        mEntity* e = uppercreatoronmc(ed);
        double k_2 = 1./1.213;// air density
        double k_1 = std::abs(inv_rho);
        // double k_1 = 1./1.213;// neg density
        // double k_2 = std::abs(inv_rho); // pos density
#else
        mEntity* e = uppercreator(ed);
        double k_1 = 1./1.213;// neg density
        double k_2 = std::abs(inv_rho); // pos density
        // double k_2 = 1./1.213;// air density
        // double k_1 = std::abs(inv_rho);
#endif
        if (e) {
            double *parama;
            parama = ComputeParameter(e, ed);
            double C_1 = *(parama);
            double C_2 = *(parama+1);

            // cout<<"mdiscr constant C1: "<<C_1<<endl;
            // cout<<"mdiscr constant C2: "<<C_2<<endl;
            double C = 1. / (1./(C_1*k_1)+1./(C_2*k_2));
            double gamma_1 = (1./(C_1*k_1))*C;
            // cout<<"eigen value final: "<<16.*C<<endl;

            datamanger_beta.setData(*e) = C;
            datamanger_gamma.setData(*e) = gamma_1;
            // cout<<"mean weighting negative : "<<gamma_1<<endl;
            
            } 
    }
    return;

}
double * Nitsche_param::ComputeParameter(AOMD::mEntity* e, mEntity* ed)
{
    xFormBilinearWithoutLaw<xGradOperator<xtensor::xCastToVectorComplex<double>>,
                xGradOperator<xtensor::xCastToVectorComplex<double>>, valtype> B_bilin;
    

    xComputeNormalFlux opNormalflux(eval_normal, classNeg, classPos);
    xGradOperatorWithParam<xComputeNormalFlux> gradflux(opNormalflux);
    xFormBilinearWithoutLaw<xGradOperatorWithParam<xComputeNormalFlux>,
                xGradOperatorWithParam<xComputeNormalFlux>, valtype> A_bilin(gradflux);
    int degree = 3;
    xIntegrationRulePartition integration_rule_neg(xMesh::getPartition, xAccept("material_neg"), 2*(degree+1));
    xIntegrationRulePartition integration_rule_pos(xMesh::getPartition, xAccept("material_pos"), 2*(degree+1));
    xIntegrationRuleBasic integration_rule_(4*(degree+1));
    xIntegrateFormCommand commandA(&A_bilin);
    AOMD::mEntity* e_integ = e;
    AOMD::mEntity* e_appro = e_integ;
    xIntegrateFormCommand commandB(&B_bilin);
    AOMD::mEntity* ed_integ = ed;
    AOMD::mEntity* ed_appro = e_appro;
    static double C[2];
    xfem::femFcts_t trialfcts(field.beginFcts(e_appro), field.endFcts(e_appro));
    if (!trialfcts.empty())
         {

        //    negtive side
            ed -> attachInt(side_tag, -1);
            e_appro -> attachInt(side_tag, -1);
            xGeomElem geo_approA(e_appro);
            A_bilin.init(&geo_approA, &trialfcts);
            integration_rule_.accept(commandA, ed);
            A_bilin.finalize();
            auto An = A_bilin.getMatrix();

            e->attachInt(side_tag, -1);        
            xGeomElem geo_approB(e_appro);
            B_bilin.init(&geo_approB, &trialfcts);
            integration_rule_neg.accept(commandB, e_appro);
            B_bilin.finalize();
            auto Bn = B_bilin.getMatrix();

            Map<Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>> MapAn(An.ReturnValPointer(), An.getNbRow(), An.getNbCol()); 
            Map<Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>> MapBn(Bn.ReturnValPointer(), Bn.getNbRow(), Bn.getNbCol()); 

            MatrixXd nA = MapAn.real();
            MatrixXd nB = MapBn.real();

            MatrixXd T(nB.rows(), nB.cols()); // transfer matrices
            T.topLeftCorner(nB.rows()/2, nB.cols()/2) = 0.5 * MatrixXd::Identity(nB.rows()/2, nB.cols()/2);
            T.topRightCorner(nB.rows()/2, nB.cols()/2) = 0.5 * MatrixXd::Identity(nB.rows()/2, nB.cols()/2);
            T.bottomLeftCorner(nB.rows()/2, nB.cols()/2) = 0.5 * MatrixXd::Identity(nB.rows()/2, nB.cols()/2);
            T.bottomRightCorner(nB.rows()/2, nB.cols()/2) = -0.5 * MatrixXd::Identity(nB.rows()/2, nB.cols()/2);
            // cout<<"Transfer matrices:\n"<<T<<endl;

// Transfer to Hansbo's basis
            nB = T * nB * T;
            nA = T * nA * T;
            // cout<<"matA:\n"<<nA<<endl;
            // cout<<"matB:\n"<<nB<<endl;
            
            nB = nB.bottomRightCorner(nB.rows()/2, nB.cols()/2);
            nA = nA.bottomRightCorner(nA.rows()/2, nA.cols()/2);

            // FullPivLU<MatrixXd> lunB(nB);
            // MatrixXd nB_null_space = lunB.kernel();
            // // cout<<"null space nB: \n"<<nB_null_space<<endl;
            // // cout<<"matB 1:"<<nB<<endl;
            // EigenSolver<MatrixXd> evnB(nB);
            // cout << "The eigenvectors of the matrix nB are:\n" << endl << evnB.eigenvectors() << endl;

            MatrixXd B_null_space(nB.rows(), 1);
            B_null_space.setOnes();
 
            MatrixXd a(1, nA.cols());
            MatrixXd b(1, nB.cols());
            MatrixXd A_til;
            MatrixXd B_til;

            // cout<<"deflation inverse: "<<endl;
            a.row(0) = nA.row(0);
            b.row(0) = nB.row(0);
            A_til = nA - B_null_space*a;
            B_til = nB - B_null_space*b;

            // cout<<"matB_til:\n"<<B_til<<endl;

            MatrixXd A_hat(A_til.rows()-1, A_til.cols()-1);
            A_hat = A_til.block(1, 1, A_til.rows()-1, A_til.cols()-1);
            MatrixXd B_hat(B_til.rows()-1, B_til.cols()-1);
            B_hat = B_til.block(1, 1, B_til.rows()-1, B_til.cols()-1);
            // cout<<"matB_hat:\n"<<B_hat<<endl;

            GeneralizedEigenSolver<MatrixXd> ges;
            ges.compute(A_hat, B_hat);
            VectorXd  Lambda_Neg_set = ges.eigenvalues().real();
            double C_neg = Lambda_Neg_set.maxCoeff();
            // cout<<"eigen value"<<C_neg<<endl;
//  positive side

            ed -> attachInt(side_tag, 1);
            e_appro -> attachInt(side_tag, 1);
            xGeomElem geo_approAp(e_appro);
            A_bilin.init(&geo_approAp, &trialfcts);
            integration_rule_.accept(commandA, ed);
            A_bilin.finalize();
            auto Ap = A_bilin.getMatrix();

            e->attachInt(side_tag, 1);        
            xGeomElem geo_approBp(e_appro);
            B_bilin.init(&geo_approBp, &trialfcts);
            integration_rule_pos.accept(commandB, e_appro);
            B_bilin.finalize();
            auto Bp = B_bilin.getMatrix();

            Map<Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>> MapAp(Ap.ReturnValPointer(), Ap.getNbRow(), Ap.getNbCol()); 
            Map<Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>> MapBp(Bp.ReturnValPointer(), Bp.getNbRow(), Bp.getNbCol()); 

            MatrixXd pA = MapAp.real();
            MatrixXd pB = MapBp.real();
            // cout<<"matA:"<<rA<<endl;
            
            pB = T * pB * T;
            pA = T * pA * T;
            // cout<<"matA in positive:\n"<<pA<<endl;
            // cout<<"matB in positive:\n"<<pB<<endl;

            pB = pB.topLeftCorner(pB.rows()/2, pB.cols()/2);
            pA = pA.topLeftCorner(pA.rows()/2, pA.cols()/2);

            a.setZero();
            b.setZero();
            A_til.setZero();
            B_til.setZero();

            // cout<<"deflation inverse: "<<endl;
            a.row(0) = pA.row(0);
            b.row(0) = pB.row(0);
            A_til = pA - B_null_space*a;
            B_til = pB - B_null_space*b;

            // cout<<"matB_til positive:\n"<<B_til<<endl;

            A_hat.setZero();
            B_hat.setZero();
            A_hat = A_til.block(1, 1, A_til.rows()-1, A_til.cols()-1);
            B_hat = B_til.block(1, 1, B_til.rows()-1, B_til.cols()-1);

            GeneralizedEigenSolver<MatrixXd> gesp;
            gesp.compute(A_hat, B_hat);
            VectorXd  Lambda_pos_set = gesp.eigenvalues().real();
            double C_pos = Lambda_pos_set.maxCoeff();
            // cout<<"eigen value pos"<<C_pos<<endl;

            C[0] = C_neg;
            C[1] = C_pos;


        
         }
    
    return C;

}


Nitsche_param_XFEM::Nitsche_param_XFEM(xfem::xMesh* meshC_, xfem::xField<valtype> &field_, xfem::xLevelSet& lst, xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, const char * side_tag_name):
 uppercreatoronmc(meshC_), eval_normal(lst), classNeg(classNeg_), classPos(classPos_), field(field_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {};

double Nitsche_param_XFEM::GetParameter(AOMD::mEntity* e)
{

    double gamma = *datamanger.getData(*e);
    return gamma;
}


void Nitsche_param_XFEM::ComputeParameter(xMesh& mInterf) {
    
    // double gamma = 0.;

    for (mEntity* ed: mInterf.range(mInterf.dim())) {
        // xPartition partition;
        // xMesh::getPartition(e, partition);
#ifdef OCTREE
        mEntity* e = uppercreatoronmc(ed);
#else
        mEntity* e = uppercreator(ed);
#endif
        // ed->print();
        // cout<<"ed: "<<e<<endl;
        if (e) {
            double parama = ComputeParameter(e, ed);
            datamanger.setData(*e) = parama;
            //  if (gamma < parama)
            //     {
            //         std::swap(gamma, parama);
            //     }
            
            }
    }
    
    return;

}

double Nitsche_param_XFEM::ComputeParameter(AOMD::mEntity* e, mEntity* ed)
{
    xUniformMaterialSensitivity<valtype > fluid_inv_density("abs_inv_density");
    xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
            xEval<valtype>,
            xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > B_bilin(fluid_inv_density);  // result type is valtype

    using valtype = std::complex<double>;
    double gamma = 0.5;
    xComputeNormalVelocity opNormalvelocity(eval_normal, fluid_inv_density, classNeg, classPos, 1.);
    xGradMeanOperatorWithParamInterConj<xComputeNormalVelocity> opMeanNormalVelocity_cong(gamma, opNormalvelocity, true);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(gamma, opNormalvelocity);
    // xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
    //         xGradMeanOperatorWithParamInterConj<xComputeNormalVelocity>, valtype > A_bilin(opMeanNormalVelocity, opMeanNormalVelocity_cong); 
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > A_bilin(opMeanNormalVelocity, opMeanNormalVelocity); 

    int degree = 3;
    xIntegrationRuleBasic integration_rule_(4*(degree+1));
    xIntegrationRulePartition integration_rule_part(xMesh::getPartition,  4*(degree+1));
    xIntegrateFormCommand commandA(&A_bilin);
    AOMD::mEntity* e_integ = e;
    AOMD::mEntity* e_appro = e_integ;
    xIntegrateFormCommand commandB(&B_bilin);
    AOMD::mEntity* ed_integ = ed;
    AOMD::mEntity* ed_appro = e_appro;
    // e_appro->print();
    double gamma_max;
    xfem::femFcts_t trialfcts(field.beginFcts(e_appro), field.endFcts(e_appro));
    if (!trialfcts.empty())
         {
            xGeomElem geo_approA(e_appro);
            A_bilin.init(&geo_approA, &trialfcts);
            integration_rule_.accept(commandA, ed);
            A_bilin.finalize();
            auto A = A_bilin.getMatrix();

            xGeomElem geo_approB(e_appro);
            B_bilin.init(&geo_approB, &trialfcts);
            integration_rule_part.accept(commandB, e_appro);
            B_bilin.finalize();
            auto B = B_bilin.getMatrix();

            Map<Matrix<valtype, Dynamic, Dynamic, RowMajor>> MapA(A.ReturnValPointer(), A.getNbRow(), A.getNbCol()); 
            Map<Matrix<valtype, Dynamic, Dynamic, RowMajor>> MapB(B.ReturnValPointer(), B.getNbRow(), B.getNbCol()); 
            
            MatrixXcd realA = MapA;
            MatrixXcd realB = MapB;

            MatrixXd absA = MapA.real();
            MatrixXd absB = MapB.real(); 

            // FullPivLU<MatrixXd> luB(absB);
            // MatrixXd B_null_space = luB.kernel();
            MatrixXd B_null_space(absB.rows(), 2);
            B_null_space.setZero();
            B_null_space.topLeftCorner(absB.cols()/2, 1) = VectorXd::Ones(absB.cols()/2);
            B_null_space.bottomRightCorner(absB.cols()/2, 1) = VectorXd::Ones(absB.cols()/2);
            // MatrixXd B_null_space(absB.cols(), 2);
            // B_null_space(B_null_space.rows(),0) = 1.;
            // cout<<"B: \n"<<absB<<endl;
            EigenSolver<MatrixXd> ev(absB);
            // cout << "The eigenvalues of the matrix are:" << endl << ev.eigenvalues() << endl;
            // cout << "The eigenvectors of the matrix are:" << endl << ev.eigenvectors() << endl;
            // cout<<"null space A: \n"<<A_null_space<<endl;
            // cout<<"null space B: \n"<<B_null_space<<endl;
            // std::ofstream outB("matB_eigen.txt",ios::app);
            // outB<<"matB:"<<absB<<endl;
            // outB.close();
            // cout<<"Size of null space B: \n"<<B_null_space.rows()<<"X"<<B_null_space.cols()<<endl;

            MatrixXd a(2, absA.cols());
            MatrixXd b(2, absB.cols());
            MatrixXd A_til;
            MatrixXd B_til;

            // cout<<"deflation inverse: "<<endl;
            a.row(1) = absA.row(absA.rows()-1);
            a.row(0) = absA.row(0);
            b.row(1) = absB.row(absB.rows()-1);
            b.row(0) = absB.row(0);
            A_til = absA - B_null_space*a;
            B_til = absB - B_null_space*b;
  
            // cout<<"A_til: \n"<<A_til<<endl;
            // cout<<"B_til: \n"<<B_til<<endl;
            MatrixXd A_hat(A_til.rows()-2, A_til.cols()-2);
            A_hat = A_til.block(1, 1, A_til.rows()-2, A_til.cols()-2);
            MatrixXd B_hat(B_til.rows()-2, B_til.cols()-2);
            B_hat = B_til.block(1, 1, B_til.rows()-2, B_til.cols()-2);
            // cout<<"A_hat: \n"<<A_hat<<endl;
            // cout<<"B_hat: \n"<<B_hat<<endl;
            // FullPivLU<MatrixXd> luB(B_hat);
            // MatrixXd B_hat_null_space = luB.kernel();
            // cout<<"null space B hat: \n"<<B_hat_null_space<<endl;
            GeneralizedEigenSolver<MatrixXd> ges;
            ges.compute(A_hat, B_hat);
            VectorXd  Lambda_set = ges.eigenvalues().real();
            double Lambda = Lambda_set.maxCoeff();
            // cout<<"lambdas :"<<Lambda<<endl;
            gamma_max = Lambda;
         }
    

    

    return gamma_max;

}

