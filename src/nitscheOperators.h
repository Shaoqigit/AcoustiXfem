// #ifndef _OPERATORS_H
// #define _OPERATORS_H
#pragma once

#include "xLevelSet.h"
#include "xEval.h"
#include "xDomain.h"
#include "xMesh.h"

#include <vector>
#include "xApproxFctPtr.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xApproxFunction.h"
#include "hoGeometricalModel_hoOctree.h"
#include "xComplexUtils.h"
#include "xNitscheFile.h"

using namespace xfem;
using namespace AOMD;



class xComputeNormalVelocity  : public std::unary_function<xtensor::xVector<>, std::complex<double> >
{
public:
    xComputeNormalVelocity(const xfem::xEval<xtensor::xVector<> >& _eval_vec, const xfem::xEval<std::complex<double> >& _eval_inv_density,
                           xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_, double indicator_):
        eval_vec(_eval_vec), eval_inv_density(_eval_inv_density), 
        //zone_tag(xfem::xDomainStringManager::getDomainId("zone")),
        side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag")),
        classNeg(classNeg_), classPos(classPos_), indicator(indicator_){};

    // Set geom is called before applying the operator on the test function
    void setGeom(xfem::xGeomElem* geo_appro, xfem::xGeomElem* geo_integ){
        eval_vec(geo_appro, geo_integ, vec);


        //Tag the right material on geo_integ depending on the side...
        int side = geo_integ->getEntity()->getAttachedInt(side_tag);
        if(!side) throw;//side_tag should have been set in xGradMeanOperatorWithParam
        else if(side > 0) classPos(geo_integ->getEntity());
        else classNeg(geo_integ->getEntity());

        eval_inv_density(geo_appro, geo_integ, inv_density);

        //Clean material classification
        unClassifyMaterial(geo_integ->getEntity());
    }

    std::complex<double> operator()(xtensor::xVector<>& in) const{

        // cout<<"in "<<in<<endl;
        // cout<<"normale "<<vec<<endl;
        // cout<<"invDens "<<inv_density<<endl;
        // cout<<"res "<<(in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2)) * inv_density<<endl;
        
    if (indicator == 0.){
        // cout<<"jump for Biot"<<endl;
        return (in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2)) / 1.213;
        
    }
    else {
        // cout<<"jump for Biot"<<endl;
        return ((in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2))) * inv_density;
    }
        
    }
private:
    const xfem::xEval<xtensor::xVector<> >& eval_vec;
    const xfem::xEval<std::complex<double> >& eval_inv_density;
    xtensor::xVector<> vec;
    std::complex<double> inv_density;
    double indicator;

    //const unsigned int zone_tag, side_tag;
    const unsigned int side_tag;
    xfem::xClassifyOn &classNeg, &classPos;

    //To unclassify
    void unClassifyMaterial(AOMD::mEntity* e)
    {
        //e->deleteData(zone_tag);
        xfem::xDomain::del(*e);
    }
};


class xComputeNormalDisplacement  : public std::unary_function<xtensor::xVector<>, std::complex<double> >
{
public:
    xComputeNormalDisplacement(const xfem::xEval<xtensor::xVector<> >& _eval_vec, const xfem::xEval<std::complex<double> >& _eval_gamma_til,
                           xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_):
        eval_vec(_eval_vec), eval_gamma_til(_eval_gamma_til), 
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

        // if(side > 0) xDomain::set(*(geo_integ->getEntity())) = 1;
        // else xDomain::set(*(geo_integ->getEntity())) = 0.;

        eval_gamma_til(geo_appro, geo_integ, gamma_til);

        //Clean material classification
        unClassifyMaterial(geo_integ->getEntity());
    }

    std::complex<double> operator()(xtensor::xVector<>& in) const{

        // cout<<"in "<<in<<endl;
        // cout<<"normale "<<vec<<endl;
        // cout<<"invDens "<<inv_density<<endl;
        // cout<<"res "<<(in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2)) * inv_density<<endl;
        
  
        // cout<<"jump for Biot"<<endl;
        // minus in front of the equation to inverse the levelset normal
        return (-(in(0)*vec(0)+in(1)*vec(1)+in(2)*vec(2))) * gamma_til;

        
    }
private:
    const xfem::xEval<xtensor::xVector<> >& eval_vec;
    const xfem::xEval<std::complex<double> >& eval_gamma_til;
    xtensor::xVector<> vec;
    std::complex<double> gamma_til;

    //const unsigned int zone_tag, side_tag;
    const unsigned int side_tag;
    xfem::xClassifyOn &classNeg, &classPos;

    //To unclassify
    void unClassifyMaterial(AOMD::mEntity* e)
    {
        //e->deleteData(zone_tag);
        xfem::xDomain::del(*e);
    }
};

#if 1
class xComputeNormalstress  : public std::unary_function<xtensor::xTensor2< > , xtensor::xVector<std::complex<double> > >
{
public:
    xComputeNormalstress(const xfem::xEval<xtensor::xVector<double> >& _eval_norm, 
    const xfem::xEval<xtensor::xTensor4<std::complex<double> > >& _eval_hook,
                           xfem::xClassifyOn &classNeg_, xfem::xClassifyOn &classPos_):
        eval_norm(_eval_norm), eval_hook(_eval_hook), 
        //zone_tag(xfem::xDomainStringManager::getDomainId("zone")),
        side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag")),
        classNeg(classNeg_), classPos(classPos_){};

    // Set geom is called before applying the operator on the test function
    void setGeom(xfem::xGeomElem* geo_appro, xfem::xGeomElem* geo_integ){
        eval_norm(geo_appro, geo_integ, normale);


        //Tag the right material on geo_integ depending on the side...
        int side = geo_integ->getEntity()->getAttachedInt(side_tag);
        if(!side) throw;//side_tag should have been set in xGradMeanOperatorWithParam
        else if(side > 0) classPos(geo_integ->getEntity());
        else classNeg(geo_integ->getEntity());

        eval_hook(geo_appro, geo_integ, hook);

        //Clean material classification
        unClassifyMaterial(geo_integ->getEntity());
    }

     xtensor::xVector<std::complex<double> >  operator()(xtensor::xTensor2< > & in) const{
        // cout<<"in "<<in<<endl;
        return (hook * castT2(in)) * castV(normale);
        // cout<<"res "<<(hook * ) * castV(normale)<<endl;
    }
private:
    const xfem::xEval<xtensor::xVector<double> >& eval_norm;
    const xfem::xEval<xtensor::xTensor4< std::complex<double > > >& eval_hook;
    xtensor::xVector<double> normale;
    xtensor::xTensor4< std::complex<double> > hook;
    xtensor::xCastToVectorComplex<double> castV;
    xtensor::xCastToTensor2Complex<double> castT2;
    //const unsigned int zone_tag, side_tag;
    const unsigned int side_tag;
    xfem::xClassifyOn &classNeg, &classPos;

    //To unclassify
    void unClassifyMaterial(AOMD::mEntity* e)
    {
        //e->deleteData(zone_tag);
        xfem::xDomain::del(*e);
    }
};
#endif


template <class UnaryOperator>
class xValImpJumpOperatorinter
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   xValImpJumpOperatorinter(std::complex<double> imp1_, std::complex<double> imp2_, const char* side_tag_name = "side_tag")
       : imp1(imp1_), imp2(imp2_), funct(), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name))
   {
   }
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      const int rs = vals.size() + f->size();
      std::vector<result_type> valsNeg = vals;
      vals.reserve(rs);
      valsNeg.reserve(rs);

      // Positive side
      geo_appro->getEntity()->attachInt(side_tag, 1);
      geo_integ->getEntity()->attachInt(side_tag, 1);
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      // Negative side
      geo_appro->getEntity()->attachInt(side_tag, -1);
      geo_integ->getEntity()->attachInt(side_tag, -1);
      xField<>::getFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

      // Untag side:
      geo_appro->getEntity()->deleteData(side_tag);
      geo_integ->getEntity()->deleteData(side_tag);

      // Reset Storage before leaving (test)
      xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

      // Compute jump :
      // std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::minus<result_type>()); //Bug ?
      for (unsigned int i = 0; i < vals.size(); ++i)
      {
         vals[i] = imp1*valsNeg[i] - imp2*vals[i];
      }
   }

  private:
   unsigned int side_tag;
   std::complex<double> imp1, imp2;
};


template <class UnaryOperator>
class xGradMeanOperatorWithParamInter {
public:
    typedef typename UnaryOperator::result_type   result_type;
public:
    //GradMeanOperatorWithParamInter(double gamma_, UnaryOperator& f, const bool conjugate_=false, const char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), conjugate(conjugate_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){}
    xGradMeanOperatorWithParamInter(double gamma_, UnaryOperator& f, const char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){}
    UnaryOperator& funct;
    void eval(xfem::femFcts_t* f, xfem::xGeomElem* geo_appro,
              xfem::xGeomElem* geo_integ, std::vector<result_type>& vals) const
    {
        const int rs = vals.size()+f->size();
        std::vector<result_type> valsNeg = vals;
        vals.reserve(rs);
        valsNeg.reserve(rs);

        //Positive side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

        //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct,true);

        //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);

        //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

        //Compute mean :
        // cout<<"compute mean gradvals"<<endl;
        //std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

        //A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE VA PAS DANS LE TRANSFORM....
        //     std::transform(vals.begin(),vals.end(),vals.begin(),bind2nd(std::multiplies<result_type>(),0.5));
        //for(int i =0; i<vals.size();++i) vals[i]*=0.5;
        double gamma_2 = 1.-gamma;
        for(unsigned int i =0; i<vals.size();++i) vals[i]=vals[i]*gamma+valsNeg[i]*gamma_2;
        // for(unsigned int i =0; i<vals.size();++i) vals[i]=vals[i]*gamma_2+valsNeg[i]*gamma;

    }

private:
    unsigned int side_tag;
    double gamma;

};

// weighting mean flux from jiang 2016
template <class UnaryOperator>
class xGradWeightMeanOperatorWithParamInter {
public:
    typedef typename UnaryOperator::result_type   result_type;
public:
    //GradMeanOperatorWithParamInter(double gamma_, UnaryOperator& f, const bool conjugate_=false, const char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), conjugate(conjugate_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){}
    xGradWeightMeanOperatorWithParamInter(xfem::xMesh& mInterf_, Nitsche_param& nitsche_, std::complex<double> inv_rho_, UnaryOperator& f, const char * side_tag_name = "side_tag") : mInterf(mInterf_), nitsche(nitsche_), inv_rho(inv_rho_), funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){
        // nitsche.ComputeParameter(mInterf, inv_rho);
// if(!nitsche.isInitialized()) throw;
    }
    UnaryOperator& funct;
    void eval(xfem::femFcts_t* f, xfem::xGeomElem* geo_appro,
              xfem::xGeomElem* geo_integ, std::vector<result_type>& vals) const
    {
        const int rs = vals.size()+f->size();
        std::vector<result_type> valsNeg = vals;
        vals.reserve(rs);
        valsNeg.reserve(rs);

        //Positive side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

        //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct,true);
        mEntity* ee = geo_appro->getEntity();
        double gamma_neg = nitsche.GetWeighting(ee);
        //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);

        //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

        //Compute mean :
        // cout<<"compute mean gradvals"<<endl;
        //std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

        //A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE VA PAS DANS LE TRANSFORM....
        //     std::transform(vals.begin(),vals.end(),vals.begin(),bind2nd(std::multiplies<result_type>(),0.5));
        //for(int i =0; i<vals.size();++i) vals[i]*=0.5;
        double gamma_pos = 1 - gamma_neg;
        // cout<<"weighting coefficient: "<<gamma_neg<<endl;
        for(unsigned int i =0; i<vals.size();++i) vals[i]=vals[i]*gamma_pos+valsNeg[i]*gamma_neg;
        // for(unsigned int i =0; i<vals.size();++i) vals[i]=vals[i]*gamma_neg+valsNeg[i]*gamma_pos;

    }
private:
    Nitsche_param& nitsche;
    xfem::xMesh& mInterf;
    unsigned int side_tag;
    double gamma;
    std::complex<double> inv_rho;

};


template <class UnaryOperator>
class xGradMeanOperatorWithParamInterConj {
public:
    typedef typename UnaryOperator::result_type   result_type;
public:
    //GradMeanOperatorWithParamInter(double gamma_, UnaryOperator& f, const bool conjugate_=false, const char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), conjugate(conjugate_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){}
    xGradMeanOperatorWithParamInterConj(double gamma_, UnaryOperator& f, const bool conjugate_=false, char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), conjugate(conjugate_), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)){}
    UnaryOperator& funct;
    void eval(xfem::femFcts_t* f, xfem::xGeomElem* geo_appro,
              xfem::xGeomElem* geo_integ, std::vector<result_type>& vals) const
    {
        const int rs = vals.size()+f->size();
        std::vector<result_type> valsNeg = vals;
        vals.reserve(rs);
        valsNeg.reserve(rs);

        //Positive side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

        //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getGradFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct,true);

        //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);

        //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

        //Compute mean :
        // cout<<"compute mean gradvals"<<endl;
        //std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

        //A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE VA PAS DANS LE TRANSFORM....
        //     std::transform(vals.begin(),vals.end(),vals.begin(),bind2nd(std::multiplies<result_type>(),0.5));
        //for(int i =0; i<vals.size();++i) vals[i]*=0.5;
        double gamma_2 = 1.-gamma;
        for(unsigned int i =0; i<vals.size();++i) vals[i]=vals[i]*gamma+valsNeg[i]*gamma_2;

        if (conjugate)
        {
            for(unsigned int i =0; i<vals.size();++i) vals[i]=std::conj(vals[i]);

        }
    }

private:
    unsigned int side_tag;
    double gamma;
    const bool conjugate;
};

template <class UnaryOperator>
class xMeanValOperatorWithParamInter {
public:
    typedef typename UnaryOperator::result_type   result_type;
public:
    xMeanValOperatorWithParamInter(double gamma_, UnaryOperator& f, const char * side_tag_name = "side_tag") : gamma(gamma_), funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
    UnaryOperator& funct;
    void eval(xfem::femFcts_t* f, xfem::xGeomElem* geo_appro,
              xfem::xGeomElem* geo_integ, std::vector<result_type>& vals) const
    {
        const int rs = vals.size()+f->size();
        std::vector<result_type> valsNeg = vals;
        vals.reserve(rs);
        valsNeg.reserve(rs);

        //Positive side
        geo_appro->getEntity()->attachInt(side_tag, 1);
        geo_integ->getEntity()->attachInt(side_tag, 1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

        //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);
        funct.setGeom(geo_appro, geo_integ);
        xfem::xField<>::getFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

        //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);

        
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);
        // cout<<"compute mean vals"<<endl;
        //Compute mean :
        //std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

        //A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE VA PAS DANS LE TRANSFORM....
        //     std::transform(vals.begin(),vals.end(),vals.begin(),bind2nd(std::multiplies<result_type>(),0.5));
        //for(int i =0; i<vals.size();++i) vals[i]*=0.5;
        double gamma_2 = 1.-gamma;
        for(unsigned int i =0; i<vals.size();++i)  { vals[i]=vals[i]*gamma+valsNeg[i]*gamma_2;}
    }

private:
    unsigned int side_tag;
    double gamma;
};

template <class UnaryOperator>
class xValOperatorWithParamSolidPart
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xValOperatorWithParamSolidPart(UnaryOperator& f, const char * side_tag_name = "side_tag") : funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
       //Negative side
        geo_appro->getEntity()->attachInt(side_tag, -1);
        geo_integ->getEntity()->attachInt(side_tag, -1);

      funct.setGeom(geo_appro, geo_integ);
      vals.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct, true);
    //   cout<<"compute neg side vals"<<endl;
      //Untag side:
        geo_appro->getEntity()->deleteData(side_tag);
        geo_integ->getEntity()->deleteData(side_tag);
    //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);
   }
   private:
   unsigned int side_tag;
};

template <class UnaryOperator>
class xValOperatorNeg
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   typedef typename UnaryOperator::argument_type argument_type;
   xValOperatorNeg(const char * side_tag_name = "side_tag") : funct(), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
    //    cout<<"eval solid part displacement"<<endl;
    //Negative side
    geo_appro->getEntity()->attachInt(side_tag, -1);
    geo_integ->getEntity()->attachInt(side_tag, -1);

    vals.reserve(vals.size() + f->size());
    xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct, true);;  // FF are assumed double valued

    //Untag side:
    geo_appro->getEntity()->deleteData(side_tag);
    geo_integ->getEntity()->deleteData(side_tag);
    //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);
   }

   private:
   unsigned int side_tag;
};

template <class UnaryOperator>
class xValOperatorPos
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   typedef typename UnaryOperator::argument_type argument_type;
   xValOperatorPos(const char * side_tag_name = "side_tag") : funct(), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {

    //Positive side
    geo_appro->getEntity()->attachInt(side_tag, 1);
    geo_integ->getEntity()->attachInt(side_tag, 1);

      vals.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);  // FF are assumed double valued

    //Untag side:
    geo_appro->getEntity()->deleteData(side_tag);
    geo_integ->getEntity()->deleteData(side_tag);
    //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);
   }

   private:
   unsigned int side_tag;
};

template <class UnaryOperator>
class xGradOperatorWithParamNeg
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xGradOperatorWithParamNeg(UnaryOperator& f, const char * side_tag_name = "side_tag") : funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {

    //Negative side
    geo_appro->getEntity()->attachInt(side_tag, -1);
    geo_integ->getEntity()->attachInt(side_tag, -1);

      funct.setGeom(geo_appro, geo_integ);
      vals.reserve(vals.size() + f->size());
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

    //Untag side:
    geo_appro->getEntity()->deleteData(side_tag);
    geo_integ->getEntity()->deleteData(side_tag);
    //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

   }

   private:
   unsigned int side_tag;
};


template <class UnaryOperator>
class xGradOperatorWithParamPos
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xGradOperatorWithParamPos(UnaryOperator& f, const char * side_tag_name = "side_tag") : funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name)) {}
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {

    //Positive side
    geo_appro->getEntity()->attachInt(side_tag, 1);
    geo_integ->getEntity()->attachInt(side_tag, 1);

      funct.setGeom(geo_appro, geo_integ);
      vals.reserve(vals.size() + f->size());
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

    //Untag side:
    geo_appro->getEntity()->deleteData(side_tag);
    geo_integ->getEntity()->deleteData(side_tag);
    //Reset Storage before leaving (test)
        xfem::xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

   }

   private:
   unsigned int side_tag;
};
