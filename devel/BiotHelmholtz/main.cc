#include "helmholtz_biot.h"
#include "xData.h"
#include "heParser.h"

#include <string>

// xTool
#include "xMPIEnv.h"

#include "xSimpleGeometry.h"

using namespace xfem;
using namespace xtool;
using namespace xgeom;


#include "xExportGmsh.h"
#include "xExportAlgorithm.h"


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
    xRegion all(data.getMesh());

    // hoGeometricalModel_hoOctree geoModel(2);

//define interface of car seat
#if 1
    // xPlane created by a point, xVector for normal direction of levelset
    // xPlane p1(xtensor::xPoint(1.1625,-0.4625,0.), xtensor::xVector<>(0.1071, 0.254,0.));
    xPlane p0(xtensor::xPoint(0.9085 , -0.3554,0.), xtensor::xVector<>(-0.1504, 0.2944000000000001,0.));  //30 28
    xPlane p1(xtensor::xPoint(1.1625,-0.4625,0.), xtensor::xVector<>(-0.1071, -0.254,0.));  //29 28
    xPlane p2(xtensor::xPoint(1.1625,-0.4625,0.), xtensor::xVector<>(0.25750000000000006, -0.04039999999999999,0.));  //29 30

    xPlane p00(xtensor::xPoint(0.725 , -0.3375,0.), xtensor::xVector<>(-0.01789999999999997, -0.1835,0.));  //27 28
    xPlane p01(xtensor::xPoint(0.725 , -0.3375,0.), xtensor::xVector<>(0.01789999999999997, 0.1835,0.));  //28 30
    xPlane p02(xtensor::xPoint(1.2029 , -0.205,0.), xtensor::xVector<>(0.205, -0.17209999999999992,0.));  //30 31
    xPlane p03(xtensor::xPoint(0.9085 , -0.3554,0.), xtensor::xVector<>(0.1504, -0.2944000000000001,0.));  //30 28 -

    // xPlane p2(xtensor::xPoint(1.4125 , -0.00625,0.), xtensor::xVector<>(0.0375, 0.03749999999999987,0.));  //315 316
    xPlane p3(xtensor::xPoint(1.375,0.0,0.), xtensor::xVector<>(-0.00625, -0.03750000000000009,0.));  //31
    xPlane p4(xtensor::xPoint(1.45 , -0.075,0.), xtensor::xVector<>(0.03125, -0.0,0.));  //33
    xPlane p5(xtensor::xPoint(1.275 , -0.475,0.), xtensor::xVector<>(0.39999999999999997, -0.17500000000000004,0.));  //31
    xPlane p6(xtensor::xPoint(1.2 , -0.5625,0.), xtensor::xVector<>(0.08750000000000002, -0.07499999999999996,0.));  //35
    xPlane p7(xtensor::xPoint(1.0875 , -0.5375,0.), xtensor::xVector<>(-0.025000000000000022, -0.11250000000000004,0.));  //36
    xPlane p8(xtensor::xPoint(0.7625 , -0.5375,0.), xtensor::xVector<>(-0.025000000000000022, -0.4375,0.));  //37 35
    xPlane p9(xtensor::xPoint(0.675 , -0.4,0.), xtensor::xVector<>(-0.13749999999999996, -0.08749999999999991,0.));  //38

    xPlane p10(xtensor::xPoint(0.725 , -0.3375,0.), xtensor::xVector<>(-0.0625, 0.04999999999999993,0.));  //27 
    xPlane p11(xtensor::xPoint(1.45 , -0.04375,0.), xtensor::xVector<>(0.03125, 0.0,0.));  //316
    xPlane p12(xtensor::xPoint(1.4125 , -0.00625,0.), xtensor::xVector<>(0.0375, 0.03749999999999987,0.));  //315
    xPlane p13(xtensor::xPoint(1.375 , 0.0, 0.), xtensor::xVector<>(0.00625, 0.03750000000000009,0.));  //31
    xPlane p14(xtensor::xPoint(0.725 , -0.3375, 0.), xtensor::xVector<>(-0.3375, 0.65,0.));  //27 31


    xUnion u1(p2, p3);
    xInter i1(p5, p6);
    xInter i2(i1, p8);
    xInter i3(i2, p9);
    xInter i4(i3, p10);
    xInter i5(i4, p11);
    xInter i6(i5, p12);
    xInter i7(i6, p13);
    xInter i8(i7, p14);

    xInter i9(p1,p2);
    xInter i10(i9, p0);
    xMinus m1(i8,i10);

    xInter i11(p00, p03);
    xInter i12(i11, p02);
    xInter i13(i12,p14);

    xMinus m2(m1, i13);
    // translate the car seat profile to 0.3, 0.03
    // two seat
    xTranslate t0(m2, xtensor::xVector<>(-0, 0.0));
    xTranslate t1(m2, xtensor::xVector<>(-0.3, 0.05));
    xTranslate t2(m2, xtensor::xVector<>(+0.5, 0.));
    xUnion s_f(t1, t2);
    xTranslate t3(s_f, xtensor::xVector<>(-0, 0.0));

    // xLevelSet ls(all, m2);
    // xLevelSet ls(all, t1);
    // two seats conf for article 2
    xLevelSet ls(all, t3);
#endif

#if 0
  xLevelSet ls(all, xPlane(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(-1.,0.,0.)));
//   xLevelSet ls(all, xPlane(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(-1.0,1.0,0.))); //45 °surface
//   xLevelSet ls(all, xSphere(xtensor::xPoint(0.08,0.025,0.025), 0.05));
    // xLevelSet ls(all, xCylinder(xtensor::xPoint(0.25,0.1,0.), xtensor::xVector<>(0., 0., 1.), 0.26925824));
#endif

    // xCylinder created by a point, xVector for normal direction of levelset
    // xCylinder cyl1(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(0., 0., 1.), 0.15);
    // xCylinder cyl2(xtensor::xPoint(-0.,0.,0.), xtensor::xVector<>(0., 0., 1.), 0.1);
    // xMinus u_cyl(cyl1, cyl2);
    // xLevelSet ls(all, u_cyl);
    // xLevelSet ls(all, xCylinder(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(0., 0., 1.), 0.1));
    // xLevelSet ls(all, xCompl(xCylinder(xtensor::xPoint(0.,0.,0.), xtensor::xVector<>(-1.,0.,0.))));
    helmholtzNdbiot helmholtzbiot(data, MPI_COMM_WORLD, parsedinfos);
    helmholtzbiot.setLevelSet(ls);
    helmholtzbiot.TreatmentOfFormulation(parsedinfos, t3);


    // wait every destruction done
    MPI_Barrier(MPI_COMM_WORLD);

    return xMPIEnv::finalize();

}
