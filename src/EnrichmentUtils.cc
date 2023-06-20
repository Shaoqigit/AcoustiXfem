#pragma once

#include "xMesh.h"
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
#include "EnrichmentUtils.h"
using namespace Trellis_Util;
using namespace xfem;
using namespace AOMD;
using namespace xtensor;





void lineExporter(xMesh &m, xtensor::xPoint p1, xtensor::xPoint p2, xEval<double> &eval, int nbPoints, string filename){

  xtensor::xVector<> v12(p1,p2);
  double length = v12.normValue();
  double step = length / (nbPoints-1);


  //Evaluation points
  std::vector<xtensor::xPoint> points;
  points.reserve(nbPoints);
  //   std::vector<double> curvAbs(nbPoints);
  points.push_back(p1);
  //   curvAbs.push_back(0);

  for(int iPt=1; iPt < nbPoints; ++iPt){
      points.push_back(xtensor::xPoint(p1(0) + iPt * step * v12(0) , p1(1) + iPt * step * v12(1), p1(2) + iPt * step * v12(2)));
      //     curvAbs.push_back(iPt * step / length);
    }

  ofstream file(filename, ofstream::out);

  double curvAbs{0.};
  for(xtensor::xPoint pt: points){
      std::set<mEntity *> elts;
      m.locateElement(pt);


      if(!elts.empty()){
          xGeomElem geo(*elts.begin());
          geo.setUVWForXYZ(pt);
          double value;
          eval(&geo,&geo,value);

          file << curvAbs << " " << pt(0)<< " "<<  pt(1)<< " "<<  pt(2)<< " "<< value<<endl;

        }

      curvAbs += step / length;
    }

  file.close();
}

xScalarFunctionDerivDiscAnalytical::xScalarFunctionDerivDiscAnalytical(xPoint center_, double rad_) : center(center_), rad(rad_) {}

void xScalarFunctionDerivDiscAnalytical::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
  if (storage.exist(geo_appro, geo_integ, res)) return;
  xPoint xyz = geo_integ->getXYZ();
  xVector<> pos(center,xyz);
  res = abs(pos.normValue()-rad);
//  res = xyz(0)*xyz(0)*xyz(0)*xyz(0);
  storage.store(geo_appro, geo_integ, res);
}

void xScalarFunctionDerivDiscAnalytical::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xVector<>& res) const
{
    xPoint xyz = geo_integ->getXYZ();
    xVector<> pos(center,xyz);
    res = pos.norm();

     xVector<> pos2(center,xyz);
    if(pos2.normValue() < rad) res*=-1;
//  res[0]=0.0; res[1]=0.0; res[2]=0.0;
}

xScalarFunctionDiscXFEMPSG::xScalarFunctionDiscXFEMPSG(const xLevelSet& ls_) : ls(ls_) {}

void xScalarFunctionDiscXFEMPSG::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   res = (double)ls.side_of_Octree(geo_appro, geo_integ);
   storage.store(geo_appro, geo_integ, res);
}

void xScalarFunctionDiscXFEMPSG::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   res[0] = 0.0;
   res[1] = 0.0;
   res[2] = 0.0;
}



