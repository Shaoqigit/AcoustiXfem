/*
    xfem : C++ Finite Element Library
    Copyright (c) 2001-2003 Nicolas MOËS and Jean-François REMACLE
    Contact : <nicolas.moes@ec-nantes.fr>
    This file is part of the xfem library.
    See the NOTICE & LICENSE files for conditions.
 */
//#define EIGEN_RUNTIME_NO_MALLOC
//#define EIGEN_USE_BLAS
#include <fstream>
#include <cassert>

#include "main.h"
#include "hoPhysicalModel.h"
#include "hoParser.h"
#include "xEval.h"
#include "hoExactSolutions.h"
#include "xAlgorithm.h"
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
#include "hoGeometricalModel_hoOctree.h"

#include <iomanip>


using namespace xtool;
using namespace xfem;
using namespace xlinalg;
using namespace std;
using namespace xtensor;
using namespace xexport;
using namespace Trellis_Util;

//#include "boost/program_options.hpp"



//#define XOCTREE_GEO_INTERFACE 0
#define HOOCTREE_GEO_INTERFACE 1

// sphere
double ls_circle(const double& x, const double& y, const double& z)
{
      return - sqrt((x-00)*(x-00)+(y-0.)*(y-0.)) + 0.3;//0.69
//      return - sqrt((x-100)*(x-100)+(y-100.)*(y-100.)) + 0.3;
//  return - sqrt((x-1)*(x-1)+(y-0.5)*(y-0.5)) + 0.3;
//  return - sqrt((x-0.67)*(x-0.67)+(y-0.5)*(y-0.5)) + 0.3;
//  return - sqrt((x-0.)*(x-0.)+(y-1.)*(y-1.)) + 0.3;
}

int main(int argc, char *argv[])
{
    xMPIEnv::init(argc,argv);
    
    /*    
  //Initialise command-line parser
  namespace po = boost::program_options;
  po::options_description description("testPhaseField Usage");
  description.add_options()
      ("help,h", "Display this help message")
      ("data,d", po::value<string>(), "main.dat (BCs & mesh definition file)")
      ("param,p", po::value<string>(), "param.dat (physical and numerical definitions)");


  //Read command-line
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
  po::notify(vm);

  if(vm.count("help")){
      std::cout << description;
      return 0;
      }*/
    
    string pname = "data/main.dat";
    string iname = "data/param.dat";

    //xParseData parsedinfos = getInfos(iname);

    std::array<double, 3> box_inf{-1,-1,0};
    std::array<double, 3> box_sup{1,1,0};
    std::array<int, 3> per{0,0,0};
    const int Dim = 2;

    string in{"matrix"};
    string out{"out"};

    if(Dim > 2){
        box_inf[2] = -1;
        box_sup[2] = 1;
    }


    xIntegrationRuleBasic integ_basic(2);
    xIntegrationRulePartition integ_part(xMesh::getPartition, 2);
    xIntegrationRuleBasic integ_basicIn(xAccept(in), 2);
    xIntegrationRulePartition integ_partIn(xMesh::getPartition, xAccept(in), 2);
    xIntegrationRulePartition integ_partOut(xMesh::getPartition, xAccept(out), 2);

    xEvalConstant<double> one(1.);
    xExportGmshAsciiSort pexport;




#ifdef HOOCTREE_GEO_INTERFACE
    // Geometry Class
    hoGeometricalModel_hoOctree geoModel(Dim);

    std::string meshFile = "data/square.msh";
    xMesh loadedMesh(meshFile);
    const int nLevels = 3;
    auto lset_ = [](std::array<double,3> xyz) -> double { return -(xyz[0]*xyz[0] + xyz[1]*xyz[1] - 0.3*0.3);};

//    geoModel.createMeshG(loadedMesh, nLevels);
    geoModel.createMeshGFromLs(loadedMesh, lset_, nLevels);

    geoModel.createLsG(lset_);
    geoModel.createPhysicalSurface(in, out);
    geoModel.classifyElementsC();
    geoModel.preprocessFilters();



    Export(one, pexport, "PART", integ_part, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));

    Export(one, pexport, "PARTIN", integ_partIn, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));
    Export(one, pexport, "PARTOUT", integ_partOut, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));
    Export(one, pexport, "PARTBIN", integ_basicIn, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));

    xRegion all(geoModel.getMeshC());
    xFilteredRegion < xIter, xAccept >  matrix(all.begin(), all.end(), xAccept(in));
    Export(one, pexport, "MATRIX", integ_basic, matrix.begin(), matrix.end());



    // Retag mesh boundaries
//    for(xIter it = geoModel.getMeshC()->begin(1); it != geoModel.getMeshC()->end(1); ++it){
//        mEntity *ed = *it;
//        auto centroid = ed->getCentroid();
//        if(abs(centroid(1) + 1) < 1.e-13) ed->classify(geoModel.getMeshC()->getMesh().getGEntity(1, 1)); // bottom
//        if(abs(centroid(1) - 1) < 1.e-13) ed->classify(geoModel.getMeshC()->getMesh().getGEntity(3, 1)); // top
//        if(abs(centroid(0) + 1) < 1.e-13) ed->classify(geoModel.getMeshC()->getMesh().getGEntity(2, 1)); // left
//        if(abs(centroid(0) - 1) < 1.e-13) ed->classify(geoModel.getMeshC()->getMesh().getGEntity(4, 1)); // right
//    }

//    for(xIter it = geoModel.getMeshG()->begin(1); it != geoModel.getMeshG()->end(1); ++it){
//        mEntity *ed = *it;
//        auto centroid = ed->getCentroid();
//        if(abs(centroid(1) + 1) < 1.e-13) ed->classify(geoModel.getMeshG()->getMesh().getGEntity(1, 1)); // bottom
//        if(abs(centroid(1) - 1) < 1.e-13) ed->classify(geoModel.getMeshG()->getMesh().getGEntity(3, 1)); // top
//        if(abs(centroid(0) + 1) < 1.e-13) ed->classify(geoModel.getMeshG()->getMesh().getGEntity(2, 1)); // left
//        if(abs(centroid(0) - 1) < 1.e-13) ed->classify(geoModel.getMeshG()->getMesh().getGEntity(4, 1)); // right
//    }

     AOMD_Util::Instance()->ex_port("meshGGG.msh", &geoModel.getMeshG()->getMesh());
#endif


    // wait every destruction done
    MPI_Barrier(MPI_COMM_WORLD);

    return xMPIEnv::finalize();

}

