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

#include "helmholtzPSG.h"

// #include "hoParser.h"
#include "heParser.h"
#include "xEval.h"

#include "xAlgorithm.h"
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
// #include "hoGeometricalModel_hoOctree.h"
#include "xMPIEnv.h"
#include "xSimpleGeometry.h"
#include <iomanip>


using namespace xtool;
using namespace xfem;
using namespace xlinalg;
using namespace std;
using namespace xtensor;
using namespace xexport;
using namespace Trellis_Util;
using namespace xgeom;

#define HOOCTREE_GEO_INTERFACE 1


int main(int argc, char *argv[])
{
    xMPIEnv::init(argc,argv);
    
    xData data;
    std::string dataname(argv[1]);
    data.ReadInfo(dataname.c_str());

    std::string paramname;
    if(argc>2) paramname = argv[2];

    xParseData parsedinfos = getInfos(paramname);



    data.ReadMesh();

    // xRegion all(data.getMesh());
    // xLevelSet ls(all, xCylinder(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(0., 0., 1.), 0.1));

    helmholtzNdPSG helmholtzNd(data, MPI_COMM_WORLD, parsedinfos);
    // helmholtzNd.setLevelSet(ls);
    helmholtzNd.TreatmentOfFormulation(parsedinfos);

    // wait every destruction done
    MPI_Barrier(MPI_COMM_WORLD);

    return xMPIEnv::finalize();

}

