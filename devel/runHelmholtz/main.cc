#include "helmholtz.h"
#include "xData.h"
#include "heParser.h"

#include <string>

// xTool
#include "xMPIEnv.h"

#include "xSimpleGeometry.h"

using namespace xfem;
using namespace xtool;
using namespace xgeom;



int main(int argc, char *argv[])
{
    xMPIEnv::init(argc,argv);

    xData data;
    std::string dataname(argv[1]);
    data.ReadInfo(dataname.c_str());

    std::string paramname;
    if(argc>2) paramname = argv[2];

    xParseData parsedinfos = getInfos(paramname);

    //data.ReadInfo("data/main_test.dat");
    data.ReadMesh();


//define interface
    // xMesh *pmesh = data.mesh
    // xRegion zone1(data.mesh, "21");
    xRegion all(data.getMesh());
    // 
    xLevelSet ls(all, xPlane(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(1.,0.,0.)));
    // xLevelSet ls(all, xSphere(xtensor::xPoint(0.07,0.025,0.025), 0.018));
    // xLevelSet ls(all, xCylinder(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(0., 0., 1.), 0.1));
    // xLevelSet ls(all, xCompl(xCylinder(xtensor::xPoint, xtensor::xVector<>(0., 0., 1.), 0.2)));

    helmholtzNd helmholtz(data, MPI_COMM_WORLD, parsedinfos);
    helmholtz.setLevelSet(ls);
    helmholtz.TreatmentOfFormulation(parsedinfos);

    // wait every destruction done
    MPI_Barrier(MPI_COMM_WORLD);

    return xMPIEnv::finalize();

}
