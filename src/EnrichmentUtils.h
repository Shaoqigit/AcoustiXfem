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
using namespace Trellis_Util;
using namespace xfem;
using namespace AOMD;
using namespace xtensor;


void lineExporter(xfem::xMesh &m, Trellis_Util::mPoint p1, Trellis_Util::mPoint p2, xfem::xEval<double> &eval, int nbPoints, string filename);

struct xAcceptSupportHasSubmesh{
    xAcceptSupportHasSubmesh() {};

    bool operator()(AOMD::mEntity* e)
    {
//        if(e->getLevel() == 2) cout<<"size2 "<<e->size(2)<<endl;

        for(int isupp=0; isupp != e->size(2); ++isupp){
            mEntity *esupp = e->get(2,isupp);


            xPartition part;
            xMesh::getPartition(esupp,part);

            if(part.size()>1) return true;

        }

        if(e->getLevel() == 2){
            xPartition part;
            xMesh::getPartition(e,part);
            if(part.size()>1) return true;
        }

        return false;

    }

//    xLevelSet &ls;
};

class xScalarFunctionDerivDiscAnalytical : public xApproxFunction {
public:
  xScalarFunctionDerivDiscAnalytical(xPoint center_, double rad_);
  void  getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& )const;
  void  getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xVector<>&)const;
  void resetStorage(){storage.reset();}
private:
  xPoint center;
  const double rad;
  mutable xEvalStorage storage;
};


struct xAcceptInTestCaseCDG1
{
  xAcceptInTestCaseCDG1(xPoint center_=xPoint(0.,0.,0.),double rad=1.) {r=rad; center=center_;}
  bool operator()(AOMD::mEntity* e)
  {
    std::vector<xPoint> nds;
    if (e->getLevel()>0)
    {
//      for(int i=0;i<e->size(0);i++)
//      {
//        nds.push_back(((mVertex*)e->get(0,i))->point());
//      }
      xGeomElem gg(e);
      nds.push_back(gg.getCDGxyz());
    }
    else
    {
      nds.push_back(((mVertex*)e)->point());
    }
    bool ok=false;
    for( auto it=nds.begin();it!=nds.end();++it)
    {
//      xPoint P=(*it)->point();
      xPoint P=(*it);
      P(0)-=center(0);
      P(1)-=center(1);
      //      cout << P ;
      double rr=sqrt(P(0)*P(0)+P(1)*P(1));
      //    cout << " " << rr << endl;
      if (rr<r)
      {
        ok=true;break;
      }
    }
    return ok;
  }

  double r;
  xtensor::xPoint center;
};

struct xAcceptCyl
{
  xAcceptCyl(double rad) {r=rad;}
  bool operator()(AOMD::mEntity* e)
  {
    std::vector<mVertex *> nds;
    if (e->getLevel()>0)
    {
      for(int i=0;i<e->size(0);i++)
      {
  nds.push_back((mVertex*)e->get(0,i));
      }
    }
    else
    {
      nds.push_back((mVertex*)e);
    }
    bool ok=false;
    for( std::vector<mVertex *>::iterator it=nds.begin();it!=nds.end();++it)
    {
      mPoint P=(*it)->point();
//      cout << P ;
      double rr=sqrt(P(0)*P(0)+P(1)*P(1));
//    cout << " " << rr << endl;
      if (rr<r)
      {
  ok=true;break;
      }
    }
    return ok;
  }
  double r;
};


class xScalarFunctionDiscXFEMPSG : public xApproxFunction
{
  public:
   xScalarFunctionDiscXFEMPSG(const xLevelSet& ls_);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void resetStorage() override { storage.reset(); }

  private:
   const xLevelSet& ls;
   mutable xEvalStorage storage;
};