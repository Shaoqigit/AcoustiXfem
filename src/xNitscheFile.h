#pragma once

#include "xEval.h"
#include "xMesh.h"
#include "xDomain.h"
#include "xLevelSet.h"
#include "xField.h"
// #include "nitscheOperators.h"
#include "hoGeometricalModel_hoOctree.h"
#include "xEntityToEntity.h"
#include "xApproxFctPtr.h"
#include "xEntityFilter.h"
#include "xEvalStorage.h"
#include "xValKey.h"
#include "xVector.h"
// #include "xVector.h"

// #include "xGenericSparseMatrix.h"
using namespace xfem;
using namespace AOMD;



class xUpperCreatorOnMeshC {
public:
  xUpperCreatorOnMeshC(xMesh *mC_) {mC=mC_;}
   AOMD::mEntity* operator()(AOMD::mEntity* e);

  xMesh *mC;
};

class xEvalNormalFronLevelSet : public xfem::xEval<xtensor::xVector<double> > {
// In order to receve the const xfem::xLevelSet, change the input type of ls with type const
public:
    xEvalNormalFronLevelSet(xfem::xLevelSet &ls_): ls(ls_), crea(new xPartitionCreatorInterface) {};
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<double>& res) const{

        // Gradient is constant on every element, and computed on approx mesh
        // cout<<ls.isDefinedOnElement(geo_appro->getEntity())<<endl;
        mEntity* e_integ = geo_integ->getEntity();
        mEntity* e_appro = geo_appro->getEntity();
        mEntity * e_integ_ls = nullptr;

        if (!ls.isDefinedOnElement(e_appro)){
            if (!ls.isDefinedOnElement(e_integ)){
                e_integ_ls = (*crea)(e_integ);
            }else{
                e_integ_ls = e_integ;
            }
        }else{
            e_integ_ls = e_appro;
        }
        if(!e_integ_ls) throw;
        // cout<<ls.isDefinedOnElement(e_integ_ls)<<endl;

        res = ls.getGrad(e_integ_ls);

        // res = ls.getGrad(geo_integ->getEntity());
        // cout<<res<<endl;
        //Warning !!!
        //normValue returns the norm of the vector (before normalization), and NORMALIZE IT IN PLACE !
        //we do it just in case the levelset was not properly reinitialized.

        res.normValue();

    }

private:
    xfem::xLevelSet &ls;
    xPartitionCreatorInterface *crea;

};


class xComputeNormalFlux  : public std::unary_function<xtensor::xVector<>, std::complex<double> >
{
public:
    xComputeNormalFlux(const xfem::xEval<xtensor::xVector<> >& _eval_vec,
                           xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_):
        eval_vec(_eval_vec), 
        //zone_tag(xfem::xDomainStringManager::getDomainId("zone")),
        side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag")),
        classNeg(classNeg_), classPos(classPos_){};

    // Set geom is called before applying the operator on the test function
    void setGeom(xfem::xGeomElem* geo_appro, xfem::xGeomElem* geo_integ){
        eval_vec(geo_appro, geo_integ, vec);


        //Tag the right material on geo_integ depending on the side...
        int side = geo_integ->getEntity()->getAttachedInt(side_tag);
        if(!side) throw;//side_tag should have been set in xGradMeanOperatorWithParam
        else if(side > 0) classPos(geo_integ->getEntity());
        else classNeg(geo_integ->getEntity());


        //Clean material classification
        unClassifyMaterial(geo_integ->getEntity());
    }

    std::complex<double> operator()(xtensor::xVector<>& in) const{

        // cout<<"in "<<in<<endl;
        // cout<<"normale "<<vec<<endl;
        // cout<<"invDens "<<inv_density<<endl;
        // cout<<"res "<<(in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2)) * inv_density<<endl;
        
        // cout<<"jump for Biot"<<endl;
        return (in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2));

        
    }
private:
    const xfem::xEval<xtensor::xVector<> >& eval_vec;
    xtensor::xVector<> vec;

    //const unsigned int zone_tag, side_tag;
    const unsigned int side_tag;
    xfem::xClassifyOn &classNeg, &classPos;

    //To unclassify
    void unClassifyMaterial(AOMD::mEntity* e)
    {
        //e->deleteData(zone_tag);
        xDomain::del(*e);
    }
};


class Nitsche_param{

using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD<double>;


using valtype = std::complex<double>;
const valtype j{0.,1.};
public:

Nitsche_param(xfem::xMesh* meshC_, xfem::xField<valtype> &field_, xfem::xLevelSet& lst_, xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, const char * side_tag_name = "side_tag");
double GetStable(AOMD::mEntity* e);
double GetWeighting(AOMD::mEntity* e);
void ComputeParameter(xfem::xMesh& m, valtype inv_rho);
double * ComputeParameter(AOMD::mEntity* e, AOMD::mEntity* ed);


private:
datamanager_t datamanger_gamma;
datamanager_t datamanger_beta;
xEvalNormalFronLevelSet eval_normal;
xfem::xClassifyOn &classNeg;
xfem::xClassifyOn &classPos;
xfem::xField<valtype> &field;
xfem::xMesh* meshC;
xUpperCreatorOnMeshC uppercreatoronmc;
xUpperCreator uppercreator;
unsigned int side_tag;

};

class Nitsche_param_XFEM{

using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD<double>;


using valtype = std::complex<double>;
const valtype j{0.,1.};
public:

Nitsche_param_XFEM(xfem::xMesh* meshC_, xfem::xField<valtype> &field_, xfem::xLevelSet& lst, xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, const char * side_tag_name = "side_tag");
double GetParameter(AOMD::mEntity* e);
void ComputeParameter(xfem::xMesh& m);
double ComputeParameter(AOMD::mEntity* e, AOMD::mEntity* ed);


private:
datamanager_t datamanger;
xEvalNormalFronLevelSet eval_normal;
xfem::xClassifyOn &classNeg;
xfem::xClassifyOn &classPos;
xfem::xField<valtype> &field;
xfem::xMesh* meshC;
xUpperCreatorOnMeshC uppercreatoronmc;
xUpperCreator uppercreator;
std::complex<double> inv_density;
unsigned int side_tag;
// const xEval<double> &eval;

};

template <typename T>
class  xEvalStabParam : public xEval<std::complex<double>> {

using valtype = std::complex<double>;
public:
// typedef typename xEval<T>::result_type result_type;
xEvalStabParam(Nitsche_param_XFEM& nitsche_, xfem::xMesh& mInterf_, valtype alpha_, xfem::xLevelSet &ls_):
nitsche(nitsche_), mInterf(mInterf_), alpha(alpha_), ls(ls_) {
    // nitsche.ComputeParameter(mInterf);
}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double> & lambda) const {
    
    // Nitsche_param_XFEM getgamma(field, eval_normal, classNeg, classPos);
    
    AOMD::mEntity* ee = geo_appro->getEntity();
    AOMD::mEntity* ei = geo_integ->getEntity();
    // if (!ls.isDefinedOnElement(ei))
    // {
    double gamma = 16. * nitsche.GetParameter(ee);
    lambda = 1./(alpha+1./gamma);
    // cout<<"runing nitsche calculation: gamma: "<<gamma<<endl;
    // }else
    // { lambda = 0.;}

// to be verified by lg, and improve the code ? no necessary for line 83 ?
    }


private:
Nitsche_param_XFEM &nitsche;
xfem::xMesh& mInterf;
xfem::xLevelSet &ls;
valtype alpha;

};


template <typename T>
class  xEvalGamma : public xEval<double> {

using valtype = std::complex<double>;

public:
// typedef typename xEval<T>::result_type result_type;
xEvalGamma(Nitsche_param_XFEM& nitsche_, xfem::xMesh& mInterf_, valtype alpha_, xfem::xLevelSet &ls_):
nitsche(nitsche_), mInterf(mInterf_), alpha(alpha_), ls(ls_) {
    // nitsche.ComputeParameter(mInterf);
}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, double & gamma) const {
    
    AOMD::mEntity* ee = geo_appro->getEntity();
    AOMD::mEntity* ei = geo_integ->getEntity();
    
    if (!ls.isDefinedOnElement(ei))
    {
    gamma = 16. * nitsche.GetParameter(ee);
    }
    else
    {gamma = 0.;
    }
    // cout<<"runing nitsche calculation: gamma: "<<gamma<<endl;
    }


private:
Nitsche_param_XFEM &nitsche;
xfem::xMesh& mInterf;
xfem::xLevelSet &ls;
valtype alpha;

};

template <typename T>
class  xEvalStabParamJiang : public xEval<std::complex<double>> {

using valtype = std::complex<double>;
public:
// typedef typename xEval<T>::result_type result_type;
xEvalStabParamJiang(Nitsche_param& nitsche_, xfem::xMesh& mInterf_, valtype alpha_, xfem::xLevelSet &ls_, valtype inv_rho_, double beta_f_):
nitsche(nitsche_), mInterf(mInterf_), alpha(alpha_), ls(ls_), inv_rho(inv_rho_), beta_f(beta_f_) {
    // nitsche.ComputeParameter(mInterf, inv_rho);
}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double> & lambda) const {
    
    // Nitsche_param_XFEM getgamma(field, eval_normal, classNeg, classPos);
    
    AOMD::mEntity* ee = geo_appro->getEntity();
    AOMD::mEntity* ei = geo_integ->getEntity();
    // if (!ls.isDefinedOnElement(ei))
    // {
    double beta = 16. * nitsche.GetStable(ee);
    // double beta = beta_f;
    lambda = 1./(alpha+1./beta);
    // cout<<"runing nitsche calculation: gamma: "<<beta<<endl;
    // }else
    // { lambda = 0.;}

// to be verified by lg, and improve the code ? no necessary for line 83 ?
    }


private:
Nitsche_param &nitsche;
xfem::xMesh& mInterf;
xfem::xLevelSet &ls;
valtype alpha;
double beta_f;
valtype inv_rho;
};


template <typename T>
class  xEvalBetaJiang : public xEval<double> {

using valtype = std::complex<double>;

public:
// typedef typename xEval<T>::result_type result_type;
xEvalBetaJiang(Nitsche_param& nitsche_, xfem::xMesh& mInterf_, valtype alpha_, xfem::xLevelSet &ls_, valtype inv_rho_):
nitsche(nitsche_), mInterf(mInterf_), alpha(alpha_), ls(ls_), inv_rho(inv_rho_) {
    // nitsche.ComputeParameter(mInterf, inv_rho);
}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, double & beta) const {
    
    AOMD::mEntity* ee = geo_appro->getEntity();
    AOMD::mEntity* ei = geo_integ->getEntity();
    
    if (!ls.isDefinedOnElement(ei))
    {
    beta = 16. * nitsche.GetStable(ee);
    }
    else
    {beta = 0.;
    }
    // cout<<"runing nitsche calculation: gamma: "<<gamma<<endl;
    }


private:
Nitsche_param &nitsche;
xfem::xMesh& mInterf;
xfem::xLevelSet &ls;
valtype alpha;
valtype inv_rho;

};


