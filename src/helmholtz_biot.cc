#include "helmholtz_biot.h"
#include <iomanip>
#include <cmath>
#include <time.h>
// Export
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"


//xfem::distmesh
#include "xLoadBalanceTools.h"
#include "xSpaceOLD.h"

// Xfiles
#include "xSpacePolynomial.h"
#include "xValueCreators.h"
#include "xField.h"
#include "xAlgorithm.h"
#include "xDistStateOfValue.h"
#include "xLinearSystemSolverMumps.h"
#include "xGenericOperations.h"
#include "xForm.h"
#include "xMaterialSensitivity.h"
#include "xPhysSurf.h"
#include "ParUtil.h"

// From helmholtzNd
#include "exactSolutionsGenericInterface.h"
#include "AnalyticalMaterial.h"
#include "xComplexUtils.h"
#include "exactSolutions.h"
#include "BiotexactSolutions.h"
#include "BiotCoupexactSolutions.h"
#include "CylinderExactSolutions.h"
#include "acousticMaterial.h"
#include "heParser.h"
#include "nitscheOperators.h"
#include "EnrichmentUtils.h"
#include "xNitscheFile.h"

#include "xSimpleGeometry.h"
#include "hoGeometricalModel_hoOctree.h"

// #define CHECK_WITH_EIGEN 1
#define GENERIC_INTERFACE
// #define XFEM_CONTINUE
#define FEM_STANDARD
// #define XFEM_DISCONTINUE
#define AIRBIOT
// #define BIOTBIOT
// #define CYL
// #define THREE_DIM

// #define CAR
// #define USE_NITSCHE_CONTINUE
// #define NURBS
// #define OCTREE


// from NEXFEM
#include "nurbsUtilities.h"
#include "sisl.h"
#include "xSpacePolynomialOLD.h"
#include "orthoEnrichFunction2.h"

// matrix manipulation
// #ifdef CHECK_WITH_EIGEN
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/QR" 
#include "unsupported/Eigen/MatrixFunctions"
#include "unsupported/Eigen/SparseExtra"
#include "xGenericSparseMatrixTraitPolicy.h"
// #endif

using namespace std;
using namespace xfem;
using namespace xtool;
using namespace xtensor;
using namespace AOMD;
using namespace xlinalg;
using namespace xcut;       
using namespace Eigen;


void lineExporterbt(xMesh &m, xMesh &mG, xtensor::xPoint p1, xtensor::xPoint p2, xEval<double> &eval, int nbPoints, string filename){

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
      //std::set<mEntity *> elts;
      auto elt_uvw_list = m.locateElement(pt);

      //std::set<mEntity *> eltsG;
      auto elt_uvw_listG = mG.locateElement(pt);



      if(!elt_uvw_list.empty()){
          auto elt_uvw = elt_uvw_list.front();
          AOMD::mEntity* pe = const_cast<AOMD::mEntity*>(elt_uvw.first);

            auto elt_uvwG = elt_uvw_listG.front();
          AOMD::mEntity* peG = const_cast<AOMD::mEntity*>(elt_uvwG.first);

          xGeomElem geo(pe);
          //geo.setUVWForXYZ(elt_uvw.second);
            geo.setUVW(elt_uvw.second);
         xGeomElem geoG(peG);
          //geoG.setUVWForXYZ(elt_uvwG.second);
            geoG.setUVW(elt_uvwG.second);

          double value;
          eval(&geo,&geoG,value);

          file << curvAbs << " " << pt(0)<< " "<<  pt(1)<< " "<<  pt(2)<< " "<< value<<endl;

        }

      curvAbs += step / length;
    }

  file.close();
}

double PointExporterbt(xMesh &m, xMesh &mG, xtensor::xPoint p, xEval<double> &eval){


      std::set<mEntity *> elts;
      auto elt_uvw_list = m.locateElement(p);
      auto elt_uvw_listG = mG.locateElement(p);


      if(!elt_uvw_list.empty()){
        auto elt_uvw = elt_uvw_list.front();
        AOMD::mEntity* pe = const_cast<AOMD::mEntity*>(elt_uvw.first);
        auto elt_uvwG = elt_uvw_listG.front();
        AOMD::mEntity* peG = const_cast<AOMD::mEntity*>(elt_uvwG.first);

          xGeomElem geo(pe);
          geo.setUVW(elt_uvw.second);
          xGeomElem geoG(peG);
          geoG.setUVW(elt_uvwG.second);
          double value;
          eval(&geo,&geoG,value);
        return value;
        }
        

}


helmholtzNdbiot::helmholtzNdbiot (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos) : data(data_) {

    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);
    // read materials
    //Give the pulsation to the porous material
    xfem::xAcousticEquivalentFluid::setOmega(parsedinfos.getDouble("Frequency_Hz"));
    xfem::xAcousticPEMBiot::setOmega(parsedinfos.getDouble("Frequency_Hz"));
    xfem::xAcousticLimpFluid::setOmega(parsedinfos.getDouble("Frequency_Hz"));
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_fluid_porous", xfem::xMaterialCreator < xfem::xAcousticEquivalentFluid >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_air", xfem::xMaterialCreator < xfem::xAcousticAir >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_limp", xfem::xMaterialCreator < xfem::xAcousticLimpFluid >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_biot", xfem::xMaterialCreator< xfem::xAcousticPEMBiot>() );
    xexport::xExportGmshAsciiSort pexport;

    xfem::xMesh* mesh=data.getMesh();
    AOMD::mMesh* mmesh=&mesh->getMesh();

    //Redirection of all output to individual files
#if 0
    pid = std::to_string(proc_id);
    string fo = "proc_"+pid+"_output.txt";
    FILE * ok = freopen (fo.c_str(),"w",stdout);
    if (!ok) {std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<< std::endl; throw; }
#endif

    // set  partition with parmetis
    xinterface::aomd::xAttachedDataManagerAOMD < int > target_proc = xmeshtool::setParmetisPartition(*mmesh, mesh->getPartitionManager(),  nb_proc  );
    xmeshtool::migrateEntities(*mmesh, mesh->getPartitionManager(), target_proc);

    xmeshtool::exportGMSH_V2_Dist(*mmesh,mesh->getPartitionManager(), "partmesh");
    cout<<"starting reading zone"<<endl;

    //Read zones
    data.ReadZones();
    cout<<"reading zone"<<endl;

}

void helmholtzNdbiot::TreatmentOfFormulation (xParseData &parsedinfos, xTranslate s_f) {

    const int Dim = 2;
#ifdef OCTREE

    cout<<"Octree mesh construct" <<endl;
    xfem::xMesh* meshO=data.getMesh(); // read the original mesh from data
    // AOMD::mMesh* mmeshO=&meshO->getMesh();

    
    hoGeometricalModel_hoOctree geoModel(2);

    xIntegrationRuleBasic integ_basic(2);
    xIntegrationRulePartition integ_part(xMesh::getPartition, 2);
    xIntegrationRuleBasic integ_basicPos(xAccept("material_pos"), 2);
    xIntegrationRulePartition integ_partPos(xMesh::getPartition, xAccept("material_pos"), 2);
    xIntegrationRulePartition integ_partNeg(xMesh::getPartition, xAccept("material_neg"), 2);

    const int nLevels = parsedinfos.getInt("levelG");;
    // auto lset_ = [](std::array<double,3> xyz) -> double { return -sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 0.1;};
    // auto lset_ = [](std::array<double,3> xyz) -> double { return -(xyz[0]*xyz[0] + xyz[1]*xyz[1] - 0.1*0.1);};
    // cylinder surface in the tube
    cout<<"tube cylinder case"<<endl;
    auto lset_ = [](std::array<double,3> xyz) -> double { return sqrt((xyz[0]-0.2)*(xyz[0]-0.2) + (xyz[1]-0.1)*(xyz[1]-0.1)) - sqrt(0.02);};
    // plane wave scattering
    // auto lset_ = [](std::array<double,3> xyz) -> double { return sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) - 0.12;};
    // auto lset_ = [](std::array<double,3> xyz) -> double { return sqrt(xyz[0]*xyz[0] + (xyz[1]-0.25)*(xyz[1]-0.25)) - 0.12;};
    // car
    // auto lset_ = [&s_f] (std::array<double,3> xyz) -> double { xPoint p(xyz[0], xyz[1], xyz[2]); return s_f(p);};
    geoModel.createMeshGFromLs(*meshO, lset_, nLevels);
    geoModel.createLsG(lset_);
    // geoModel.createPhysicalSurface("material_pos", "material_neg");
    geoModel.createPhysicalSurface("material_neg", "material_pos");
    geoModel.classifyElementsC();
    geoModel.preprocessFilters();

{
    xEvalConstant<double> one(1.);
    xexport::xExportGmshAsciiSort pexport;

    Export(one, pexport, "PARTIN", integ_partPos, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));
    Export(one, pexport, "PARTOUT", integ_partNeg, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));
    // xEvalConstant<double> one_(1.);
    // xexport::xExportGmshMPIIO pexportG;
    // Export(lset_, pexportG,"LSET");
}
// throw;

    std::cout<<"==============finish the octree mesh===================="<<endl;

    xfem::xMesh* mesh=geoModel.getMeshC(); // set mesh as MeshC
    xfem::xMesh* meshG = geoModel.getMeshG();
    xRegion all(mesh); // define the region

#endif
#ifndef OCTREE
    xfem::xMesh* mesh=data.getMesh();
    AOMD::mMesh* mmesh=&mesh->getMesh();
    xRegion all(mesh);
#endif
    // int Dim = mesh->dim();
    cout << "Start creating the function space" << endl;


    const int degree = parsedinfos.getInt("order_approx_min");
    // int  = parsedinfos.getInt("odebug")
    double freq = parsedinfos.getDouble("Frequency_Hz");
    double omega = 2.*M_PI*freq;
    double beta = parsedinfos.getDouble("H_wave_angle");  // degree
    double wave_angle = beta*M_PI/180.;  // radians
    double sigmaF = parsedinfos.getDouble("Resistivity");
    double d = parsedinfos.getDouble("Mat_Thickness");
    double alpha_ = parsedinfos.getDouble("Nitsche_fluid");
    double alpha_2 = parsedinfos.getDouble("Nitsche_solid");
    double gamma = parsedinfos.getDouble("Average_gamma");
    double h = parsedinfos.getDouble("Element_size");
    std::string coupling_type = parsedinfos.getString("Coupling_type");
    double radius = parsedinfos.getDouble("Cylinder_radius");

// #ifdef CAR
// for (double freq = 1.; freq <=1000; freq = freq + 10.) {
//     double omega = 2*M_PI*freq;
// #endif

// Generailised transfer matrix for limp model
#ifdef GENERIC_INTERFACE
    double phi_F = parsedinfos.getDouble("Porosity");
    double alpha_F = parsedinfos.getDouble("Tortuosity");
    double lambda_prime_F = parsedinfos.getDouble("Thermal_length");
    double lambda_F = parsedinfos.getDouble("Viscous_length");

    double rho_1_F = 809.;
    double nu_F = 0.3;
    double E_F = 260e6;
    double eta_F = 0.5;

    Analytical_MatLimp mat_F(omega, wave_angle, phi_F, sigmaF, alpha_F, lambda_prime_F, lambda_F, rho_1_F, nu_F, E_F, eta_F);
    valtype Z = mat_F.rho_limp_til * mat_F.c_eq_til;
    valtype k_a = mat_F.k_a;
    valtype k_ay = k_a*sin(wave_angle);
    valtype k_eqx = sqrt(pow(mat_F.rho_limp_til,2)-pow(k_ay,2));
    std::cout<<"Analytical rho limp til: "<<mat_F.rho_limp_til<<std::endl;

    Analytical_MatLimp mat_F2(omega, wave_angle, 0.72, 87e3, 1.02, 480e-6, 480e-6, 171, 0., 50e3, 0.5);
    valtype Z2 = mat_F2.rho_limp_til * mat_F.c_eq_til;
    // valtype k_a = mat_F.k_a;
    // valtype k_ay = k_a*sin(wave_angle);
    valtype k_eqx2 = sqrt(pow(mat_F2.rho_limp_til,2)-pow(k_ay,2));
    
    cout<<"k_ay: "<<k_ay<<endl;
    cout<<"k_eqx: "<<k_eqx<<endl;
    cout<<"rho_eqx: "<<mat_F.rho_limp_til<<endl;
    double  d1 = 0.0005628;
    Matrix2cd M1;
    M1(0,0) = cos(k_eqx * d1); 
    M1(0,1) = -1. * pow(omega,2) * mat_F.rho_limp_til/k_eqx * sin(k_eqx * d1); 
    M1(1,0) = k_eqx/(pow(omega,2) * mat_F.rho_limp_til) * sin(k_eqx * d1);
    M1(1,1) = cos(k_eqx * d1); 

    Matrix2cd M2;
    // double d2 = 0.0006/sqrt(2.);
    double d2 = 0.0004236;
    M2(0,0) = cos(k_eqx2 * d2); 
    M2(0,1) = -1. * pow(omega,2) * mat_F2.rho_limp_til/k_eqx2 * sin(k_eqx2 * d2); 
    M2(1,0) = k_eqx2/(pow(omega,2) * mat_F2.rho_limp_til) * sin(k_eqx2 * d2);
    M2(1,1) = cos(k_eqx2 * d2); 

    // M1(0,0) = 9.99184026e-01+7.89194092e-04*j; 
    // M1(0,1) = -1.58702392e+07+1.09508537e+07*j; 
    // M1(1,0) = 1.16113455e-10-1.92535727e-11*j;
    // M1(1,1) = 9.99184026e-01+7.89194092e-04*j; 

    Matrix2cd M;
    M = M2*M1*M2;
    // M = M1;
    cout<<"Film wave number: "<<k_eqx<<endl;
    valtype M11 = M(0,0); cout<<"Coeff M_11: "<<std::setprecision(5)<<M11<<endl;
    valtype M12 = M(0,1); 
    // valtype M12 = 0.;
    cout<<"Coeff M_12: "<<std::setprecision(5)<<M12<<endl;
    valtype M21 = M(1,0); cout<<"Coeff M_21: "<<std::setprecision(5)<<M21<<endl;
    valtype M22 = M(1,1); cout<<"Coeff M_22: "<<std::setprecision(5)<<M22<<endl;
    cout<<"Pressure jump: "<<-sigmaF *d<<endl;
#endif

#ifdef NURBS
            //Create interface:
            SISLCurve *pc=nullptr;

            int kdim=2, kn, kk;
            //  double aepsge=0.000001; /* geometric tolerance */
            double aepsge=1.e-12;
            double aepsco=1.e-10; /* computational tolerance */


            const double radiusCircle=radius;
            xtensor::xPoint centerCircle(0.,0.,0.);

            //Spline cercle
            double pi = 4.*atan(1.);
            double startpt[3]={centerCircle(0)+radiusCircle*cos(-0.*pi),centerCircle(1)+radiusCircle*sin(-0.*pi),0. };
            double cepsge=1.e-12;
            double angle=4*pi;
            double centrept[3]={centerCircle(0),centerCircle(1),0.};
            double axis[3]={0.,0.,1.};
            int stat;

            if(pc) freeCurve(pc);
            s1303(startpt, cepsge, angle, centrept, axis, kdim, &pc, &stat);

            std::vector<double> cZeroKnots;

            //Arc de cercle -> si R=1.0, surface exacte=3.21460183660255
            kn=3;//number nodes
            kk=3;//order
            double knots[8],coeffs[9];

            //Knot vector
            knots[0]=0.; knots[1]=0.; knots[2]=0.;
            knots[3]=3.; knots[4]=3.;
            knots[5]=6.; knots[6]=6.; knots[7]=6.;


            coeffs[0]=centerCircle(0)+radiusCircle; coeffs[1]=centerCircle(1); coeffs[2]=1.;
            coeffs[3]=centerCircle(0)+radiusCircle; coeffs[4]=centerCircle(1)+radiusCircle; coeffs[5]=1./sqrt(2.);
            coeffs[6]=centerCircle(0); coeffs[7]=centerCircle(1)+radiusCircle; coeffs[8]=1.;
            
            // int stat;
            // int kind=2;
            // //Create Spline curve
            // if(pc) freeCurve(pc);
            // pc = newCurve(kn,kk,knots,coeffs,kind,kdim,1);
            if(!pc) cout<<"courbe pas construite\n";

            double startpar, endpar;
            s1363(pc,&startpar, &endpar, &stat );
            cout<<"Start Par "<<startpar<<" End Par "<<endpar<<endl;
 
// throw;
            xfem::xMesh surface(MPI_COMM_WORLD,2);
            AOMD::mMesh &surfaceAOMD = surface.getMesh();
            // On evalue la position du ou des pts:
            for(int ii=0;ii<500+1;++ii){
                int stat_;
                int computeDerivative=0;
                int noeudADroite2=0;
                double derivees2[1*kdim];//si on veut seulement la pos.
                double intersParam_[1]={(endpar-startpar)/500*ii};
                s1227(pc,computeDerivative, intersParam_[0], &noeudADroite2, derivees2, &stat_);
                //    cout<<derivees2[0]<<" "<<derivees2[1]<<endl;
                surfaceAOMD.createVertex(derivees2[0],derivees2[1],0.,surfaceAOMD.getGEntity(0,0));

            }

            AOMD_Util::Instance()->ex_port("surface.msh",&surfaceAOMD);

              
            xfem::xMesh slice(MPI_COMM_WORLD,199);
            AOMD::mMesh *sliceAOMD = &slice.getMesh();
            //  int degpTheta=2*degp; //polynome de degre degp+1 selon theta -> grad en degp
            //  int degpLambda=2*(degp+3); //polynome de degre degp+3 (cubic spline) selon theta -> grad en degp+2 -> degp+3 car nurbs
            //  int degpClassical=2*(degp-1); //polynome de degre degp -> grad en degp-1


            int degpTheta=100*degree; //polynome de degre degp+1 selon theta -> grad en degp
            int degpLambda=100*(degree+3); //polynome de degre degp+3 (cubic spline) selon theta -> grad en degp+2 -> degp+3 car nurbs
            int degpClassical=100*(degree-1); //polynome de degre degp -> grad en degp-1

            degpTheta=100;
            degpLambda=100;
            degpClassical=100;


            const int maxLevel = 5;
            subdivideAndAttachNurbsEnhancedGaussPoints(data.getMesh(), pc, degpTheta, degpLambda, degpClassical,&slice,&cZeroKnots, maxLevel);
            AOMD_Util::Instance()->ex_port("slice.msh",sliceAOMD);



            xIntegrationRulePartition integration_ruleC(xMesh::getPartition, 4*degree);
            xIntegrationRuleBasic integration_ruleB(2*degree);
            xexport::xExportGmshAsciiSort pexportA;
            xEvalConstant<double > one(1.);
            Export(one, pexportA, "ONE", integration_ruleC, all.begin(), all.end());
            Export(one, pexportA, "ONEB", integration_ruleB, all.begin(), all.end());

            double result=0.;
            xIntegrateEvalCommand<xEvalConstant<double> > integOne(one, result);

            xAcceptAll acc;
            xRegion matrix(data.getMesh());
            Export(one, pexportA, "AONE", integration_ruleC, matrix.begin(), matrix.end());
            xIntegrationRulePartition integration_rule_envT(xMesh::getPartition,  50);
            xIntegrationRulePartition integration_ruleTest(xMesh::getPartition,  30);

            
            xIntegrationRuleStoredPartition integration_ruleSB(xMesh::getPartition, 5*(degree+1));



            cout<<"construct the spline"<<endl;
#endif

    cout<<"Approximation order is "<<degree<<endl;

// Function for exact solution
#ifdef AIRBIOT
cout<<"air-Biot coupling starting..."<<endl;
    xEvalExactBiotUs Biot_exact_disp(omega, wave_angle, sigmaF, d);
    xEvalExactBiotTraction Biot_exact_stress(omega, wave_angle, sigmaF, d);
    xEvalExactBiotUtotal Biot_disp_total(omega, wave_angle, sigmaF, d);
    xEvalExactBiotPressure Biot_exact_press(omega, wave_angle, sigmaF, d);
    xEvalExactBiotVelocity Exact_velocity(omega, wave_angle, sigmaF, d);
#endif
#ifdef BIOTBIOT
cout<<"Biot-Biot coupling starting..."<<endl;
    xEvalExactCoupBiotUs Coup_Biot_exact_disp(omega, wave_angle, sigmaF, d);
    xEvalExactCoupBiotTraction Coup_Biot_exact_stress(omega, wave_angle, sigmaF, d);
    xEvalExactCoupBiotUtotal Coup_Biot_disp_total(omega, wave_angle, sigmaF, d);
    xEvalExactCoupBiotPressure Coup_Biot_exact_press(omega, wave_angle, sigmaF, d);
#endif
#ifdef CYL
cout<<"Air-Biot scattering problem starting..."<<endl;
    xEvalExactCylinderPressure Cylinder_pressure(omega, radius, sigmaF, d);
    xEvalExactCylinderUs Cylinder_displacement(omega, radius, sigmaF, d);
    xEvalExactCylinderUtotal Cylinder_displacement_total(omega, radius, sigmaF, d);
    xEvalExactCylinderVelocity  Cylinder_velocity(omega, radius, sigmaF, d);
    xEvalExactCylinderTraction Cylinder_stress(omega, radius, sigmaF, d);
#endif



// create the polynomial space for different fields
    xSpacePolynomialBernstein dispx("DISPLACEMENT_X", xSpace::VECTOR_X, degree);
    xSpacePolynomialBernstein dispy("DISPLACEMENT_Y", xSpace::VECTOR_Y, degree);
#ifdef THREE_DIM
    xSpacePolynomialBernstein dispz("DISPLACEMENT_Z", xSpace::VECTOR_Z, degree);
    xSpaceComposite disp(dispx, dispy, dispz);
    // OLD::xSpacePolynomialOLDLagrange lagEx("DISPLACEMENT_X", xSpace::VECTOR_X, degree);
    // OLD::xSpacePolynomialOLDLagrange lagEy("DISPLACEMENT_Y", xSpace::VECTOR_Y, degree);
#else
    xSpaceComposite disp(dispx, dispy);
#endif
    // xSpaceComposite  lagrangeE(lagEx, lagEy);

    xSpacePolynomialBernstein press("PRESSURE", xSpace::SCALAR, degree);

    xfem::xValueManagerDist<double> double_manager;
    xfem::xValueManagerDist<valtype> complex_manager;
    xField<valtype> press_l(&complex_manager, press);

#ifndef OCTREE
    xClassifyOn classNeg("material_neg");
    xClassifyOn classPos("material_pos");
    for(xIter ita = all.begin(1); ita != all.end(1); ++ita){
         classPos(*ita);
     }
#endif

#ifdef NURBS
 for(auto it=matrix.begin() ; it!= matrix.end(); ++it){
                mEntity *elem = *it;
                xGeomElem geo(elem);
                xPartition partition;
                data.getMesh()->getPartition(elem,partition);


                for(xPartition::iterator itp = partition.begin() ; itp!= partition.end(); ++itp){
                    mEntity *elem_sub = *itp;
                    xGeomElem geo_sub(elem_sub);
                    xPoint cdg_sub = geo_sub.getCDGxyz();
                    geo.setUVWForXYZ(cdg_sub);
                    xPoint cdg = geo.getUVW();
                    double val=0.;
//                    lsHole.getVal(elem,cdg, val);
                    val = sqrt((cdg_sub(0)-0.)*(cdg_sub(0)-0.) + (cdg_sub(1)-0.)*(cdg_sub(1)-0.))-radius;

                    if(val>0.) classPos(elem_sub);
                    else{
                        classNeg(elem_sub);
                     classNeg(elem);
                    }

                    // cout<<"CLASS\n";
                }

                // cout<<"----\n";

            }

            // xIntegrationRuleStoredPartition integration_ruleSBM(xMesh::getPartition, xAccept("material_pos"), 10*degree);
//              Export(one, pexportA, "MATRIX", integration_ruleSBM, all.begin(), all.end());

            xIntegrationRuleStoredPartition integration_ruleSBI(xMesh::getPartition, xAccept("material_neg"), 10*degree);
            pexportA.setNbSplit(2*degree);
            Export(one, pexportA, "INCLUSION", integration_ruleSBI, all.begin(), all.end());
//             cout<<"classify\n";

        // xcut::xPhysSurfParameter param_surf(classNeg, classPos);
        // inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));
        // xSpaceFiltered::filter_t filter_C = xAcceptCyl(0.2);
        // xFilteredRegion<xIter, xSpaceFiltered::filter_t>  porous(all.begin(), all.end(), filter_C);




        xFilteredRegion<xIter, xAccept> porous(all.begin(), all.end(), xAccept("material_neg"));                
// {
//         xEvalConstant<double> one_(1.);
//         xexport::xExportGmshMPIIO pexport;
//         Export(one_, pexportA, "POROUS2", integration_ruleSBI, porous.begin(), porous.end());
//         Export(one_, pexportA, "POROUS2B", integration_ruleB, porous.begin(), porous.end());
// //         Export(*lset, pexport,"LSET");
// }
        // xSpaceFiltered::filter_t filter_d(std::bind(&xcut::xPhysSurfByTagging::supportCoversIn, inclusion, std::placeholders::_1));
        
        // xField<valtype> disp_l(&complex_manager,  xSpaceFiltered(disp, filter_d));
        xField<valtype> disp_l(&complex_manager, disp);
// throw;
#endif

#ifdef FEM_STANDARD
    int nb_air = 110;
    int nb_porous = 114;
    // xClassRegion porous(mesh, 122, 2);
    // xClassRegion porous(mesh, 105, 2);
    // xClassRegion porous(mesh, 114, 2);
    xClassRegion porous(mesh, nb_porous, 2);
    // xClassRegion porous_F1(mesh, 111, 2);
    // xClassRegion porous_F2(mesh, 112, 2);
    // xClassRegion porous_F3(mesh, 113, 2);
    xField<valtype> disp_l(&complex_manager, disp);
    // xValueCreator < xSingleValue<valtype> >  creator;
    // DeclareInterpolation<xField<valtype> >(disp_l, creator, porous.begin(), porous.end());
#endif

#ifdef XFEM_CONTINUE

        // direction of levelset, normal vector
    ////////////////////////////////////////////////////////
    //                         //                         //
    //         +       <-------//             -           //
    //                         //                         //
    ////////////////////////////////////////////////////////

        //--- Create physical surface object
        xcut::xPhysSurfParameter param_surf(classNeg, classPos);
        inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));
#if defined(AIRBIOT)||defined(CYL)
        xFilteredRegion < xIter, xAccept >  porous(all.begin(), all.end(), xAccept("material_neg"));
        xSpaceFiltered::filter_t filter_d(std::bind(&xcut::xPhysSurfByTagging::supportCoversIn, inclusion, std::placeholders::_1));
        xField<valtype> disp_l(&complex_manager,  xSpaceFiltered(disp, filter_d));
                
{
        xEvalConstant<double> one_(1.);
        xIntegrationRulePartition integr(xMesh::getPartition,xAccept("material_neg"),1);
        xexport::xExportGmshMPIIO pexport;
        pexport.setNbSplit(2*degree);
        Export(one_, pexport, "POROUS2", integr, porous.begin(), porous.end());
        Export(*lset, pexport,"LSET");
}
// throw;
#endif
        // xField<valtype> disp_l(&complex_manager, disp);

        //--- Setup enriched space for pressure field
        //a/ keyword to create the enriched keys
        xValKeyExtend key_modifier("_MAT_ENRICH");
        //b*/ enrichment function : Ridge
        xScalarFunctionDerivDiscXFEM enrichment_function(inclusion->getLevelSet());
        //c*/ filter for ridge function "supportCutStrictlyEltWise", do not cut at bord of element
        xSpaceFiltered::filter_t filter(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictlyEltWise, inclusion, std::placeholders::_1));
        //d/ create enriched space (but ALL the dofs are enriched...)
        xSpaceXFEM space_full(press, enrichment_function, key_modifier);
        //e/ filter the space in order to keep only the relevant nodes
        xSpaceFiltered enriched(space_full, filter);
        //f/ add the enriched space to the approximation
        press_l.insert(enriched);

#ifdef BIOTBIOT
        // xScalarFunctionDiscXFEM enrichment_function_discont(inclusion->getLevelSet());
        xSpaceXFEM space_full_disp(disp, enrichment_function, key_modifier);
        xSpaceFiltered enriched_disp(space_full_disp, filter);
        disp_l.insert(enriched_disp);
#endif
    
    cout<<"End XFEM continue enrichment "<<endl;
#endif

#ifdef XFEM_DISCONTINUE
cout<<"start XFEM discontinue enrichment "<<endl;
        //--- Create physical surface object

#ifndef OCTREE
        xcut::xPhysSurfParameter param_surf(classNeg, classPos);
        inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));
#endif
        // coupling air-Biot
#if defined(AIRBIOT)||defined(CYL)
        xFilteredRegion < xIter, xAccept >  porous(all.begin(), all.end(), xAccept("material_neg"));
        // xSpaceFiltered::filter_t filter_d(std::bind(&xcut::xPhysSurfByTagging::supportCoversIn, inclusion, std::placeholders::_1));

// {
//         xEvalConstant<double> one_(1.);
//         xIntegrationRulePartition integr(xMesh::getPartition,xAccept("material_neg"),1);
//         xexport::xExportGmshMPIIO pexport;
//         pexport.setNbSplit(2*degree);
//         Export(one_, pexport, "POROUS2", integr, porous.begin(), porous.end());
//         Export(*lset, pexport,"LSET");
// }
// throw;

#ifdef OCTREE
    auto filter_d = bind1st(mem_fun(&hoGeometricalModel::filterSupportCoversInC), &geoModel);
    // filter the porous field for the displacement with coverIn
    {
    xEvalConstant<double> one(1.);
    xexport::xExportGmshAsciiSort pexport;
    xIntegrationRuleBasic integ_basicCover_(filter_d, 4);
    Export(one, pexport, "COVER___S", integ_basicCover_, geoModel.getMeshC()->begin(2), geoModel.getMeshC()->end(2));
    xIntegrationRulePartition integration_porous(xMesh::getPartition,xAccept("material_neg"), 2*(degree+1));
    Export(one, pexport, "POROUS2", integration_porous, porous.begin(), porous.end());
    }
#endif
// throw;
        xField<valtype> disp_l(&complex_manager,  xSpaceFiltered(disp, filter_d));
        // xField<valtype> disp_l(&complex_manager, disp);
#else
        // xField<valtype> disp_l(&complex_manager, disp);
#endif
        //--- Setup enriched space
        //a/ keyword to create the enriched keys
        xValKeyExtend key_modifier("_MAT_ENRICH");
        //b/ enrichment function : Heaviside
        // xScalarFunctionDiscXFEM enrichment_function(inclusion->getLevelSet());
        
        //c/ filter to enrich the necessary dofs only (supportCutStrictlyEltWise is defined in xPhysSurfByTagging)
        // xSpaceFiltered::filter_t filter(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictly, inclusion, std::placeholders::_1));
    
#ifdef OCTREE
    //b/ enrichment function : Heaviside
    xScalarFunctionDiscXFEM_ enrichment_function(*geoModel.getLsG());
    //c/ filter to enrich the necessary dofs only (supportCutStrictlyEltWise is defined in xPhysSurfByTagging)
    auto filter = bind1st(mem_fun(&hoGeometricalModel_hoOctree::filterSupportCutStrictly), &geoModel);
#endif 

    {
    // xEvalConstant<double> one(1.);
    // xexport::xExportGmshAsciiSort pexport;
    // xIntegrationRuleBasic integ_basicCut_(filter, 4);
    // Export(one, pexport, "CUT___S", integ_basicCut_, geoModel.getMeshC()->begin(2), geoModel.getMeshC()->end(2));
    // Export(one, pexport, "CUT___E", integ_basicCut_, geoModel.getMeshC()->begin(1), geoModel.getMeshC()->end(1));
    // Export(one, pexport, "CUT___V", integ_basicCut_, geoModel.getMeshC()->begin(0), geoModel.getMeshC()->end(0));
    // xEvalApproxFunction<double> evalappro(enrichment_function);
    // Export(evalappro, pexport, "Heaviside", integ_part, geoModel.getMeshG()->begin(2), geoModel.getMeshG()->end(2)); 

    }
        //d/ create enriched space (but ALL the dofs are enriched...)
        xSpaceXFEM space_full(press, enrichment_function, key_modifier);
        //e/ filter the space in order to keep only the relevant nodes
        xSpaceFiltered enriched(space_full, filter);
        //f/ add the enriched space to the approximation
        press_l.insert(enriched);

#ifdef BIOTBIOT
        xValKeyExtend key_modifier_disp("_RIDGE_ENRICH");
        xScalarFunctionDerivDiscXFEM ridge_function(inclusion->getLevelSet());
        xSpaceFiltered::filter_t filter_disp(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictlyEltWise, inclusion, std::placeholders::_1));
        xSpaceXFEM space_full_disp(disp, ridge_function, key_modifier_disp);
        xSpaceFiltered enriched_disp(space_full_disp, filter_disp);
        disp_l.insert(enriched_disp);
        cout<<"End XFEM discontiuity enrichment "<<endl;
#endif

#endif
// throw;

#ifdef NURBS
    xSpaceFiltered::filter_t filter_E = xAcceptSupportHasSubmesh();
    xSpaceFiltered::filter_t filter(acc);
    OLD::xSpacePolynomialOLDLagrange pressE("PRESSURE", xSpace::SCALAR, degree);
    
    xValKeyExtend key_modifier_E("_Enrichment");
    xScalarFunctionDerivDiscAnalytical enrichment_function_E(xPoint(0.,0.,0), radius);  

    xSpacePolynomialLagrange lagScalar("FCTQ", xSpace::SCALAR, 1);
    xField<> alphaField(&double_manager, lagScalar);
    createDomainRamp cramp(all.begin(),all.end(), data.getMesh(), filter_E, &double_manager, alphaField);

    xEvalField<xIdentity<double> > eval_ramp(alphaField);
    xEvalGradField<xIdentity<xVector<>> > eval_Gramp(alphaField);


    correctedEnrichmentFunction enrichment_function_C(enrichment_function_E, alphaField);//Corrected function
    xSpaceFiltered::filter_t filterC=bind1st(mem_fun(&createDomainRamp::isInDomain), &cramp);

    //OK
    xSpaceFiltered enriched_tip_C(xSpaceXFEMShiftedCorrectedOLD(pressE, enrichment_function_E,key_modifier_E, &eval_ramp, &eval_Gramp), filterC);
    press_l.insert(enriched_tip_C);
cout<<"End NURBs discontiuity enrichment "<<endl;
#endif

    xValueCreator < xSingleValue<valtype> >  creator;
#if defined(AIRBIOT)||defined(CYL)
    DeclareInterpolation<xField<valtype> >(disp_l, creator, porous.begin(), porous.end());
#else
    DeclareInterpolation<xField<valtype> >(disp_l, creator, all.begin(), all.end()); 
#endif
    DeclareInterpolation<xField<valtype> >(press_l, creator, all.begin(), all.end());

    TreatmentOfEssEnv(press_l, data);
    TreatmentOfEssEnv(disp_l, data);

    complex_manager.PrintForDebug("cplx_dcl_"+pid+".dbg");
    // In case some one use this test with Dirichlet BC we need to uniformize fixed values
    finalizeFixedUniform(mesh->getPartitionManager(),complex_manager);
    complex_manager.PrintForDebug("cplx_ess_"+pid+".dbg");

    xDistStateDofCreatorUniform<xMesh::partman_t, xTrue, valtype> snh(mesh->getPartitionManager(),complex_manager, "dofs");

#if defined(AIRBIOT)||defined(CYL)
    DeclareState<xField<valtype> >(disp_l, snh, porous.begin(), porous.end()); 
#else
    DeclareState<xField<valtype> >(disp_l, snh, all.begin(), all.end());
#endif
    DeclareState<xField<valtype> >(press_l, snh, all.begin(), all.end());

    std::shared_ptr < xlinalg::xDistIndex > dist_index = finalizeDofUniform(snh);
    cout << "End creating the function space" << endl;

    cout << "Start Building the Linear System" << endl;
    complex_manager.PrintForDebug("cplx_sym_"+pid+".dbg");

    xIntegrationRulePartition integration_rule_env(xMesh::getPartition, 4*degree);
       //2*degree is the degree to integrate
    // It is necessay to instance different partition integration for levelset
 
#if defined(XFEM_CONTINUE)||defined(XFEM_DISCONTINUE)
    xIntegrationRulePartition integration_rule(xMesh::getPartition, 2*(degree+1)); 
    xIntegrationRulePartition integration_rule_neg(xMesh::getPartition,xAccept("material_neg"), 2*(degree+1));
    xIntegrationRulePartition integration_rule_pos(xMesh::getPartition,xAccept("material_pos"), 2*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc_neg(xMesh::getPartition, xAccept("material_neg"), 4*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc_pos(xMesh::getPartition, xAccept("material_pos"), 4*(degree+1));
#endif

#ifdef NURBS
    xIntegrationRuleStoredPartition integration_rule(xMesh::getPartition, 10*degree);
    xIntegrationRuleStoredPartition integration_rule_neg(xMesh::getPartition, xAccept("material_neg"), 10*degree);
    xIntegrationRuleStoredPartition integration_rule_pos(xMesh::getPartition, xAccept("material_pos"), 10*degree);
    xIntegrationRuleStoredPartition integration_bc_pos(xMesh::getPartition, 50);
#endif

 #ifdef FEM_STANDARD   
    xIntegrationRulePartition integration_rule(xMesh::getPartition, 2*(degree+1)); 
    xIntegrationRulePartition integration_rule_neg(xMesh::getPartition, xAccept(nb_porous), 2*(degree+1));
    xIntegrationRulePartition integration_rule_pos(xMesh::getPartition, xAccept(nb_air), 2*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc_neg(xMesh::getPartition, xAccept(nb_porous), 2*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc_pos(xMesh::getPartition, xAccept(nb_air), 2*(degree+1));
    // xIntegrationRulePartition integration_bc_pos(xMesh::getPartition, 2*(degree+1)); 
    // There is problem on integration when using FEM standard for multilayer problem
#endif

    xlinalg::xDistVector<valtype> b(*dist_index);
    xlinalg::xDistVector<valtype> sol(*dist_index);
    xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),1);
    // xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),false);
#if defined(AIRBIOT)||defined(CYL)
        AssembleGraph(graphsym, disp_l, porous.begin(), porous.end() );
        AssembleGraph(graphsym,press_l, all.begin(), all.end() );
        AssembleGraph(graphsym, press_l, disp_l, porous.begin(), porous.end() );
        // AssembleGraph(graphsym, disp_l, press_l, porous.begin(), porous.end() ); // unsymmetry
#endif
    
#ifdef BIOTBIOT
        AssembleGraph(graphsym,disp_l, all.begin(), all.end() );
        AssembleGraph(graphsym,press_l, all.begin(), all.end() );
        AssembleGraph(graphsym,disp_l, press_l, all.begin(), all.end() );
#endif
    //std::ofstream graph("graph.txt");
    //graphsym.printCOO(graph);
    
#ifdef USE_ZMUMPS
    typedef xlinalg::xLinearSystemSolverMumpsDistributed < valtype > solver_t_;
    typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
    // typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
#endif
    Matrix_SDP_t_ A(graphsym);
    xAssemblerBasic < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler(A, b);
    xAssemblerBasicAndTranspose < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler_transpose(A, b);

 TreatmentOfNatEnvAcoustic(mesh, press_l, assembler, integration_rule_env, Exact_velocity, omega*omega, 0.);
#ifdef AIRBIOT
    // Natural boundary conditions for the fluid part
    TreatmentOfNatEnvAcoustic(mesh, press_l, assembler, integration_bc_pos, Exact_velocity, omega*omega, 0.);
    TreatmentOfNatEnvPEMFluid(press_l, assembler, integration_bc_neg, Biot_disp_total);
    TreatmentOfNatEnvPEMSolid(disp_l, assembler, integration_bc_neg, Biot_exact_stress);
#endif
// point source
            // xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
            //                       xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);
// cout<<"Initialize the point source excitation..."<<endl;
//     xtensor::xPoint position_source(-0.1, 0.02, 0);
//     xEvalPointSourceCplx point_source(omega/340, position_source);
//     xFormLinearWithLoad < xValOperator<xtensor::xCastToComplex<double> >,
//                                   xEval<valtype>, valtype > lin_source(point_source);
//     Assemble(lin_source, assembler, integration_rule, press_l, all.begin(), all.end());

#ifdef CYL
    TreatmentOfNatEnvAcoustic(mesh, press_l, assembler, integration_bc_pos, Cylinder_velocity, omega*omega, 0.);
#endif

#ifdef BIOTBIOT
    TreatmentOfNatEnvPEMFluid(press_l, assembler, integration_rule_env, Coup_Biot_disp_total);
    // Natural boundary conditions for the solid part
    TreatmentOfNatEnvPEMSolid(disp_l, assembler, integration_rule_env, Coup_Biot_exact_stress);
#endif
    //    std::ofstream vecb_("vecb_.txt");//Export per-process vectors...
    //    b.Printdata(vecb_);

    cout << "Start constructing the weak form" << endl;
    
    //Biot formulation solid phase first term
    xUniformMaterialSensitivity<xtensor::xTensor4<valtype> > hooke("strain");
    xFormBilinearWithLaw<xGradOperator<xtensor::xCastToTensor2Complex<double> >, 
                        xUniformMaterialSensitivity<xtensor::xTensor4<valtype> >, 
                        xGradOperator<xtensor::xCastToTensor2Complex<double> >, valtype > diffusive_s(hooke);
// Biot formulation solid mass term
    xUniformMaterialSensitivity<valtype > rho_til("rho_til");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToVectorComplex<double> >,
                        xEval<valtype>,
                        xValOperator<xtensor::xCastToVectorComplex<double> >, valtype > mass_s(rho_til);
// Biot formulation coupling term
    xUniformMaterialSensitivity<valtype > gamma_til("gamma_til");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToVectorComplex<double> >,
                        xEval<valtype>,
                        xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > coupled_s(gamma_til);

// Biot formulation fluid stiff term
    xUniformMaterialSensitivity<valtype > fluid_inv_density("inv_density");
    xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
            xEval<valtype>,
            xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > diffusive_f(fluid_inv_density);  // result type is valtype

    // Biot fluid second mass term -1/K_eq_til
    xUniformMaterialSensitivity<valtype> inv_K_til("inv_bulk");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xValOperator<xtensor::xCastToComplex<double> > ,valtype > mass_f(inv_K_til);

    // xUniformMaterialSensitivity<valtype > gamma_til("gamma_til");        
    xFormBilinearWithLaw<xGradOperator<xtensor::xCastToVectorComplex<double> >,
                        xEval<valtype>,
                        xValOperator<xtensor::xCastToVectorComplex<double> >, valtype > coupled_f(gamma_til);
    // assembly the bilinear
#ifdef BIOTBIOT
    Assemble(diffusive_s, assembler, integration_rule, disp_l, disp_l, all.begin(), all.end());
    assembler.setCoeff(-omega*omega);
    Assemble(mass_s, assembler, integration_rule, disp_l, disp_l, all.begin(), all.end());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    Assemble(coupled_s, assembler, integration_rule, disp_l, press_l, all.begin(), all.end()); // unsymmetry
    assembler.setCoeff(1.);
#endif

#if defined(AIRBIOT)||defined(CYL)
    Assemble(diffusive_s, assembler, integration_rule_neg, disp_l, disp_l, porous.begin(), porous.end());
    
    assembler.setCoeff(-omega*omega);
    Assemble(mass_s, assembler, integration_rule_neg, disp_l, disp_l, porous.begin(), porous.end());
    assembler.setCoeff(1.);
    
    assembler.setCoeff(-1.);
    Assemble(coupled_s, assembler, integration_rule_neg, disp_l, press_l, porous.begin(), porous.end());
    assembler.setCoeff(1.);

    assembler.setCoeff(-1.);
    // Assemble(coupled_f, assembler, integration_rule_neg, press_l, disp_l, porous.begin(), porous.end()); // unsymmetry
    assembler.setCoeff(1.);
#endif
#ifdef FEM_STANDARD
// FILM 1
//     Assemble(diffusive_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F1.begin(), porous_F1.end());
    
//     assembler.setCoeff(-omega*omega);
//     Assemble(mass_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F1.begin(), porous_F1.end());
//     assembler.setCoeff(1.);
    
//     assembler.setCoeff(-1.);
//     Assemble(coupled_s, assembler, integration_rule_neg, disp_l, press_l, porous_F1.begin(), porous_F1.end());
//     assembler.setCoeff(1.);

//     assembler.setCoeff(-1.);
//     Assemble(coupled_f, assembler, integration_rule_neg, press_l, disp_l, porous_F1.begin(), porous_F1.end()); 
// // FILM 2
//     Assemble(diffusive_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F2.begin(), porous_F2.end());
    
//     assembler.setCoeff(-omega*omega);
//     Assemble(mass_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F2.begin(), porous_F2.end());
//     assembler.setCoeff(1.);
    
//     assembler.setCoeff(-1.);
//     Assemble(coupled_s, assembler, integration_rule_neg, disp_l, press_l, porous_F2.begin(), porous_F2.end());
//     assembler.setCoeff(1.);

//     assembler.setCoeff(-1.);
//     Assemble(coupled_f, assembler, integration_rule_neg, press_l, disp_l, porous_F2.begin(), porous_F2.end()); 
// // FILM 3
//     Assemble(diffusive_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F3.begin(), porous_F3.end());
    
//     assembler.setCoeff(-omega*omega);
//     Assemble(mass_s, assembler, integration_rule_neg, disp_l, disp_l, porous_F3.begin(), porous_F3.end());
//     assembler.setCoeff(1.);
    
//     assembler.setCoeff(-1.);
//     Assemble(coupled_s, assembler, integration_rule_neg, disp_l, press_l, porous_F3.begin(), porous_F3.end());
//     assembler.setCoeff(1.);

//     assembler.setCoeff(-1.);
//     Assemble(coupled_f, assembler, integration_rule_neg, press_l, disp_l, porous_F3.begin(), porous_F3.end()); 

#endif


    assembler.setCoeff(1./(omega*omega));
    // Assemble(diffusive_f, assembler, integration_ruleSB, press_l, press_l, all.begin(), all.end());
    Assemble(diffusive_f, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    // Assemble(mass_f, assembler, integration_ruleSB, press_l, press_l, all.begin(), all.end());
    Assemble(mass_f, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    // Assemble(coupled_f, assembler, integration_rule_neg, press_l, disp_l, porous.begin(), porous.end());
    assembler.setCoeff(1.);   

    cout << "end the weak form" << endl; 
    // interface forms


#ifdef XFEM_DISCONTINUE
    // xMesh* interface = inclusion->getMesh_bnd();
    // xEvalNormalFronLevelSet eval_normal_(*(lset.get()));
#ifdef OCTREE
    xMesh* interface = geoModel.getBnd();
    xEvalNormalFronLevelSet eval_normal_(*geoModel.getLsG());
    xClassifyOn classPos(geoModel.getInString());
    xClassifyOn classNeg(geoModel.getOutString());
#endif
    // eigen value problem

    xComputeNormalVelocity opNormalvelocity(eval_normal_, fluid_inv_density, classNeg, classPos, omega);
    xComputeNormalVelocity opNormalvelocity2(eval_normal_, fluid_inv_density, classNeg, classPos, 0.);
    xComputeNormalDisplacement opNormaldisp(eval_normal_, gamma_til, classNeg, classPos);
    xMeanValOperatorWithParamInter<xComputeNormalDisplacement> opMeandispnormal(gamma, opNormaldisp);
    xValOperatorWithParamSolidPart<xComputeNormalDisplacement> opdispnormal(opNormaldisp);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(gamma, opNormalvelocity);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity2(1., opNormalvelocity2);

    xComputeNormalstress opNormalstress(eval_normal_, hooke, classNeg, classPos);
    xGradMeanOperatorWithParamInter<xComputeNormalstress> opMeanNormalstress(0.5, opNormalstress);


// lambda * [q][p]-delta p
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > penalty_press;
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToVectorComplex<double> >,
            xValJumpOperator<xtensor::xCastToVectorComplex<double> >, valtype > penalty_disp;
    
    
    xFormBilinearWithoutLaw<xMeanValOperatorWithParamInter<xComputeNormalDisplacement>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > disp_jumpp(opMeandispnormal);      
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xMeanValOperatorWithParamInter<xComputeNormalDisplacement>, valtype > jumpp_disp(opMeandispnormal);     
  
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xMeanValOperatorWithParamInter<xComputeNormalDisplacement>, valtype > averagep_disp(opMeanNormalVelocity, opMeandispnormal); 
    xFormBilinearWithoutLaw<xMeanValOperatorWithParamInter<xComputeNormalDisplacement>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > disp_averagep(opMeandispnormal, opMeanNormalVelocity); 
    xFormBilinearWithoutLaw<xMeanValOperatorWithParamInter<xComputeNormalDisplacement>,
            xMeanValOperatorWithParamInter<xComputeNormalDisplacement>, valtype > disp_disp(opMeandispnormal); 

    // Nitsche with total displacement
 
// [[q]]<p/n>
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > jump_average_press(opMeanNormalVelocity);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > average_jump_press(opMeanNormalVelocity);
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > jump_average_press2(opMeanNormalVelocity2);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > average_jump_press2(opMeanNormalVelocity2);

// <1/rho*q/n><1/rho*p/n>
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average_press(opMeanNormalVelocity); 
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average_press_1(opMeanNormalVelocity, opMeanNormalVelocity2); 
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average_press_2(opMeanNormalVelocity2, opMeanNormalVelocity); 
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average_press_3(opMeanNormalVelocity2, opMeanNormalVelocity2); 
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToVectorComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalstress>, valtype > jump_average_disp(opMeanNormalstress);
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalstress>,
            xValJumpOperator<xtensor::xCastToVectorComplex<double> >, valtype > average_jump_disp(opMeanNormalstress);


// interface coefficient
    // valtype alpha = -sigmaF*d*omega*j;
    // valtype alpha = sigmaF*d/(j*omega);
    valtype alpha = j*omega*sigmaF*d;
    valtype K = j*sigmaF*d/(omega);
    valtype lambda = 1./(-alpha+1./alpha_);
    valtype epsilon = alpha_/(omega*omega);
    // cout<<"Nitsche para :"<<alpha_<<endl;
    cout<<"epsilon: "<<(lambda)<<endl;
    cout<<"stabilization"<<alpha_<<endl;
    cout<<"pressure jump: "<<-1./(alpha)<<endl;

// Nitsche perfect interface
if (alpha_ ==0.) {
#ifndef GENERIC_INTERFACE
    assembler.setCoeff(-1./(alpha));
    // Assemble(penalty_press, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    Assemble(penalty_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
#else
    cout<<"Generalized interface direct form"<<endl;
    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
        xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > qnpn;
    
    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
        xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > qnpp;

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
        xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > qppn;

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
        xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > qppp;


    cout<<"Generalized Fluid interface"<<endl;
    valtype F11 = M22/M12; cout<<"F11: "<<F11*j*omega<<endl;
    valtype F12 = (M12*M21-M11*M22)/M12; cout<<"F12: "<<F12*j*omega<<endl;
    valtype F21 = 1./M12; cout<<"F21: "<<F21*j*omega<<endl;
    valtype F22 = -M11/M12; cout<<"F22: "<<F22*j*omega<<endl;

    assembler.setCoeff(F22);
    // Assemble(qnpn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    Assemble(qnpn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(F21);
    // Assemble(qnpp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    Assemble(qnpp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(-F12);
    // Assemble(qppn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    Assemble(qppn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(-F11);
    // Assemble(qppp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    Assemble(qppp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    // if film is modelled by Biot
#if 0
    cout<<"Generalized Biot interface"<<endl;
    Analytical_MatBiot mat_F(omega, beta, 0.04, 775e3, 1.15, 230e-6, 230e-6, 809, 0.3, 260e6, 0.5);
    valtype k_a = mat_F.k_a;
    valtype ky = k_a*sin(wave_angle);
    valtype k_eqx = sqrt(pow(mat_F.k_eq_til,2)-pow(ky,2));

    MatrixXcd TM_(6,6);
    TM_ <<  0., 0., 0., j*ky*mat_F.A_hat/mat_F.P_hat, j*ky*mat_F.gamma_til, -(pow(mat_F.A_hat,2)-pow(mat_F.P_hat,2))/mat_F.P_hat*pow(ky,2)-mat_F.rho_til*pow(omega,2),
        0., 0., 0., 1./mat_F.P_hat, 0., j*ky*mat_F.A_hat/mat_F.P_hat,
        0., 0., 0., 0., -1./mat_F.K_eq_til+pow(ky,2)/(mat_F.rho_eq_til*pow(omega,2)), -j*ky*mat_F.gamma_til,
        j*ky, -mat_F.rho_s_til*pow(omega,2), -mat_F.rho_eq_til*mat_F.gamma_til*pow(omega,2), 0., 0., 0.,
        0., mat_F.rho_eq_til*mat_F.gamma_til*pow(omega,2), mat_F.rho_eq_til*pow(omega,2), 0., 0., 0.,
        1./mat_F.N, j*ky, 0., 0., 0., 0.;

    Analytical_MatBiot mat_F2(omega, beta, 0.72, 87e3, 1.02, 480e-6, 480e-6, 171, 0., 50e3, 0.5);
    valtype k_eqx2 = sqrt(pow(mat_F2.k_eq_til,2)-pow(ky,2));

    MatrixXcd TM_2(6,6);
    TM_2 <<  0., 0., 0., j*ky*mat_F2.A_hat/mat_F2.P_hat, j*ky*mat_F2.gamma_til, -(pow(mat_F2.A_hat,2)-pow(mat_F2.P_hat,2))/mat_F2.P_hat*pow(ky,2)-mat_F2.rho_til*pow(omega,2),
        0., 0., 0., 1./mat_F2.P_hat, 0., j*ky*mat_F2.A_hat/mat_F2.P_hat,
        0., 0., 0., 0., -1./mat_F2.K_eq_til+pow(ky,2)/(mat_F2.rho_eq_til*pow(omega,2)), -j*ky*mat_F2.gamma_til,
        j*ky, -mat_F2.rho_s_til*pow(omega,2), -mat_F2.rho_eq_til*mat_F2.gamma_til*pow(omega,2), 0., 0., 0.,
        0., mat_F2.rho_eq_til*mat_F2.gamma_til*pow(omega,2), mat_F2.rho_eq_til*pow(omega,2), 0., 0., 0.,
        1./mat_F2.N, j*ky, 0., 0., 0., 0.;

    MatrixXcd I(6,6);
    I = MatrixXcd::Identity(6,6);

    MatrixXcd TM(6,6);
    MatrixXcd TM_full(6,6);
    TM_full = (-d*TM_).exp();
    TM = I-d*TM_;

    double d2 = 0.0006;
    MatrixXcd TM2(6,6);
    TM2 = I-d2*TM_2;

    MatrixXcd TMf(6,6);
    // TMf = TM2*TM*TM2;
    TMf = TM;
    // cout<<"TM(0,1)"<<TM(0,1)<<endl;
    // cout<<"TM(3,1)"<<TM(3,1)<<endl;
    // cout<<"TM(4,1)"<<TM(4,1)<<endl;
    Matrix4cd B;
    B << -TMf(0,0), -TMf(0,3), -TMf(0,2), 0.,
        -TMf(3,0), -TMf(3,3), -TMf(3,2), 0.,
        -TMf(4,0), -TMf(4,3), -TMf(4,2), 0.,
        -TMf(2,0), -TMf(2,3), -TMf(2,2), 1.;
    Matrix4cd A;
    A << TMf(0,1), TMf(0,5), TMf(0,4), 0.,
        TMf(3,1), TMf(3,5), TMf(3,4), 0.,
        TMf(4,1), TMf(4,5), TMf(4,4), -1.,
        TMf(2,1), TMf(2,5), TMf(2,4), 0.;
    Matrix4cd SM;
    SM = B.completeOrthogonalDecomposition().pseudoInverse()*A;
    // cout<<"Impedance matrix (full)"<<TM_full<<endl;
    // cout<<"Impedance matrix"<<TM<<endl;
    // cout<<"SM(0,1): "<<SM(0,1)<<endl;

    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<0>>,
            xValOperatorNeg<xExtractCompVectorT<0>>, valtype > vx_ux; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<0>>,
            xValOperatorNeg<xExtractCompVectorT<1>>, valtype > vx_uy; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<0>>,
            xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > vx_pp; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<0>>,
            xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > vx_pn; 

    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<1>>,
            xValOperatorNeg<xExtractCompVectorT<0>>, valtype > vy_ux; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<1>>,
            xValOperatorNeg<xExtractCompVectorT<1>>, valtype > vy_uy; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<1>>,
            xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > vy_pp; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xExtractCompVectorT<1>>,
            xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > vy_pn; 

    xFormBilinearWithoutLaw<xValOperatorNeg<xtool::xIdentity < double >>,
            xValOperatorNeg<xExtractCompVectorT<0>>, valtype > qn_ux; 
    xFormBilinearWithoutLaw<xValOperatorNeg<xtool::xIdentity < double >>,
            xValOperatorNeg<xExtractCompVectorT<1>>, valtype > qn_uy; 

    xFormBilinearWithoutLaw<xValOperatorPos<xtool::xIdentity < double >>,
            xValOperatorNeg<xExtractCompVectorT<0>>, valtype > qp_ux; 
    xFormBilinearWithoutLaw<xValOperatorPos<xtool::xIdentity < double >>,
            xValOperatorNeg<xExtractCompVectorT<1>>, valtype > qp_uy; 



    assembler.setCoeff(SM(1,0));
    Assemble(vx_ux, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(1,1));
    Assemble(vx_uy, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(1,3));
    Assemble(vx_pp, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(1,2));
    Assemble(vx_pn, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(0,0));
    Assemble(vy_ux, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(0,1));
    Assemble(vy_uy, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(0,3));
    Assemble(vy_pp, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(0,2));
    Assemble(vy_pn, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(2,0));
    Assemble(qn_ux, assembler, integration_rule, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(2,1));
    Assemble(qn_uy, assembler, integration_rule, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(2,3));
    Assemble(qnpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(SM(2,2));
    Assemble(qnpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(-SM(3,0));
    Assemble(qp_ux, assembler, integration_rule, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(-SM(3,1));
    Assemble(qp_uy, assembler, integration_rule, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(-SM(3,3));
    Assemble(qppp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  

    assembler.setCoeff(-SM(3,2));
    Assemble(qppn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  
#endif

#endif
    cout << "End bilinear direct terms" << endl;
}
else {
    assembler.setCoeff(lambda);
    // Assemble(penalty_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(penalty_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(-1.*lambda*alpha/(omega*omega));
    // Assemble(jump_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(jump_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1),xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(-1.*lambda*alpha/(omega*omega));
    // Assemble(average_jump_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(average_jump_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
       
    assembler.setCoeff(lambda*alpha*alpha/(omega*omega*omega*omega));
    // Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    
    // [q]*average
    assembler.setCoeff(-1./(omega*omega));
    // Assemble(jump_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(jump_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(-1./(omega*omega));
    // Assemble(average_jump_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(average_jump_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
// // //  - delta*double_average
    assembler.setCoeff(1.*alpha/(omega*omega*omega*omega));
    // Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(1.*alpha/(omega*omega*omega*omega));
    // Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

// // +delta * double_average
    assembler.setCoeff(-1.*alpha/(omega*omega*omega*omega));
    // Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(double_average_press, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    // coupling terms with u^s
    cout<<"coupling terms for Nitsche"<<endl;
#ifdef BIOTBIOT
    assembler.setCoeff(1.);
    Assemble(disp_jumpp, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(1.);
    // Assemble(jumpp_disp, assembler, integration_rule_neg, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha/(omega*omega));
    // Assemble(averagep_disp, assembler, integration_rule_neg, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha/(omega*omega));
    Assemble(disp_averagep, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(alpha);
    Assemble(disp_disp, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(epsilon*alpha);
    Assemble(disp_jumpp, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(epsilon*alpha);
    // Assemble(jumpp_disp, assembler, integration_rule_neg, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*epsilon*alpha*alpha/(omega*omega));
    // Assemble(averagep_disp, assembler, integration_rule_neg, press_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*epsilon*alpha*alpha/(omega*omega));
    Assemble(disp_averagep, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(epsilon*alpha*alpha);
    Assemble(disp_disp, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);  
#endif
#if defined(AIRBIOT)||defined(CYL)

    xFormBilinearWithoutLaw<xValOperatorWithParamSolidPart<xComputeNormalDisplacement>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > disp_jumpp_1(opdispnormal);      
    xFormBilinearWithoutLaw<xValOperatorWithParamSolidPart<xComputeNormalDisplacement>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > disp_averagep_1(opdispnormal, opMeanNormalVelocity); 
    xFormBilinearWithoutLaw<xValOperatorWithParamSolidPart<xComputeNormalDisplacement>,
            xValOperatorWithParamSolidPart<xComputeNormalDisplacement>, valtype > disp_disp_1(opdispnormal); 
    // assembler.setCoeff(1.);
    TagSideNegative(interface);
    if (gamma !=1.){
    cout<<"coupling terms for Nitsche in neg part"<<endl;
    double gamma_2 = 1-gamma;
    assembler.setCoeff(gamma_2);
    // Assemble(disp_jumpp_1, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_jumpp_1, assembler, integration_rule_neg, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    
    assembler.setCoeff(1.);

    assembler.setCoeff(-gamma_2*alpha/(omega*omega));
    // Assemble(disp_averagep_1, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_averagep_1, assembler, integration_rule_neg, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(alpha*gamma_2*gamma_2);
    // Assemble(disp_disp_1, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_disp_1, assembler, integration_rule_neg, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(gamma_2*lambda*alpha);
    // Assemble(disp_jumpp_1, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_jumpp_1, assembler, integration_rule_neg, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(-gamma_2*lambda*alpha*alpha/(omega*omega));
    // Assemble(disp_averagep_1, assembler, integration_rule, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_averagep_1, assembler, integration_rule_neg, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(lambda*alpha*alpha*gamma_2*gamma_2);
    // Assemble(disp_disp_1, assembler, integration_rule, disp_l, disp_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(disp_disp_1, assembler, integration_rule_neg, disp_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);  
    }
#endif

    cout << "End bilinear Nitsche terms" << endl;

}
    
#endif
    // Treatment of Robin BCs
    // TreatmentOfMixedBCs(data, assembler, integration_rule_env, press_l, omega);
    cout << "End Building the Linear System" << endl;
    if(parsedinfos.getInt("debug") == 1){
        cout << "export matrix" << endl;
       std::ofstream mata_("matA_.mm");  
       A.printMatrixMarket(mata_);
    //    mmPrettyPrint(&mata_);

       std::ofstream matb_("matB_.mm");
       A.printMatrixMarket(matb_);

       //Export the matrix and vectors for checking the testcase
       string outNameb = "vecb_"+ to_string(proc_id)+".txt";
       std::ofstream ofsvecb(outNameb.c_str());
       b.Printdata(ofsvecb);
       string outNamea = "matA_"+ to_string(proc_id)+".mm";
       std::ofstream ofsmata(outNamea);
       A.printMatrixMarket(ofsmata);
    }
cout << "End export matrix" << endl;
// #ifdef CHECK_WITH_EIGEN
//     typedef Eigen::SparseMatrix<valtype, Eigen::ColMajor > SpMat;
//     SpMat Acomplex(A.getM(),A.getN());
//     typedef Eigen::Triplet<valtype > SparseTriplet;
//     std::vector<SparseTriplet> coefficients(graphsym.getNNZ());
//     fillTripletSym(A, coefficients);
//     Acomplex.setFromTriplets(coefficients.begin(), coefficients.end());
//     Acomplex.makeCompressed();
//     //    cout<<"ACOMPLEX\n"<<Acomplex<<endl;

//     b.switchInsertModeOff();
//     b.switchInsertModeOn();
//     Eigen::VectorXcd bcomplex(dist_index->getGlobalIndexSize());
//     for(unsigned int i = 0; i < bcomplex.rows(); ++i) bcomplex(i) = b[i];
// #endif

#ifdef USE_ZMUMPS
    struct timespec begin, end; 
    clock_gettime(CLOCK_REALTIME, &begin);
   solver_t_ solver;

    //#ifdef WITH_DMUMPS
    solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,14,400);
    //#endif
    // if check the solved system like conditioning number 
    //decomment them, verify the Ringo(6)

    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::INTER,0,1);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,1,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,2,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,3,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,4,2);

    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,11, 1);
    solver.connectMatrix(A);
    b.switchInsertModeOff();
    solver.solve(b, sol);
    sol.switchInsertModeOff();
    sol.switchToGlobalValue();

    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    printf("Time measured: %.3f seconds.\n", elapsed);
#endif

// #ifdef CHECK_WITH_EIGEN
//     Eigen::SparseLU<SpMat> esolver(Acomplex);
//     auto solcomplex = esolver.solve(bcomplex);
//     sol.switchInsertModeOn();
//     for(unsigned int i = 0; i < bcomplex.rows(); ++i) sol[i] = solcomplex(i);
//     //cout<<"EigenSol:\n";
//     //sol.Printdata(cout);
// #endif

    Visit(xWriteSolutionVisitor < xlinalg::xDistVector<valtype> >(sol),
          complex_manager.begin("dofs"),
          complex_manager.end("dofs"));
    cout << "End Solving the Linear System" << endl;

//#####################################post processing############################################//
    cout << "Start Post Treatment" << endl;
    complex_manager.PrintForDebug("res"+pid+".dbg");
    xexport::xExportGmshMPIIO pexport;
// displacement evaluation
    xGetRealPartVComponent realV0(0);
    xGetNormVComponent normV0(0);
    xGetRealPartVComponent realV1(1);
    xEvalField<xtool::xIdentity<xtensor::xVector<valtype> >, xField<valtype> > eval_disp(disp_l);
    xEvalField < xGetRealPartVComponent, xField<valtype> > eval_disp_realx(disp_l, realV0);
    xEvalField < xGetNormVComponent, xField<valtype> > eval_disp_normx(disp_l, normV0);
    xEvalField < xGetRealPartVComponent, xField<valtype> > eval_disp_realy(disp_l, realV1);
    // pexport.setNbSplit(2.*degree);
#if defined(AIRBIOT)||defined(CYL)
    pexport.setNbSplit(3.*degree);
    // Export(eval_disp_realx, pexport, "DISPLACEMENT_REAL_X", integration_rule_neg, porous.begin(), porous.end());
    // Export(eval_disp_realy, pexport, "DISPLACEMENT_REAL_Y", integration_rule_neg, porous.begin(), porous.end());
#endif

#ifdef BIOTBIOT
    Export(eval_disp_realx, pexport, "DISPLACEMENT_REAL_X", integration_rule, all.begin(), all.end());
    Export(eval_disp_realy, pexport, "DISPLACEMENT_REAL_Y", integration_rule, all.begin(), all.end());
#endif
// pressure evaluation
    xEvalField <xGetComplex, xField<valtype>> eval_press(press_l);
    xEvalField < xGetRealPart, xField<valtype> > eval_press_real(press_l);
    xEvalField < xGetImagPart, xField<valtype> > eval_press_imag(press_l);
    xEvalField < xGetModu<valtype>, xField<valtype> > eval_press_norm(press_l);
    xEvalField <xGetPressureIndB, xField<valtype>> eval_press_indb(press_l);
    xEvalUnary<xGetModu<double> > pressureindb(eval_press_indb);
    // Export(eval_press_norm, pexport, "PRESSURE_NORM", integration_rule, all.begin(), all.end());
    // pexport.setNbSplit(3.*degree);
    Export(eval_press_real, pexport, "XFEM_PRESSURE_REAL", integration_rule_env, all.begin(), all.end());
#ifndef CYL
    xEvalBinary<xMult<valtype, xtensor::xVector<valtype>, xtensor::xVector<valtype>> > exact_velocity_flux(fluid_inv_density, Exact_velocity);
    xEvalUnary<xGetRealPartV> exact_velocity_real(exact_velocity_flux);
    // Export(exact_velocity_real, pexport, "EXACT_SOL_Grad_press", integration_rule, all.begin(), all.end()); 
#endif
#ifdef CAR
// pexport.setNbSplit(2*degree);
// Export(eval_press_indb, pexport, "PRESSURE_REAL", integration_rule_env, all.begin(), all.end());
    // double intsolution{0.};
    // xIntegrateEvalCommand< xEvalUnary<xGetModu<double> > > integrateSolution(pressureindb, intsolution);  // int one 
    // ApplyCommandOnIntegrationRule(integrateSolution, integration_rule, all.begin(), all.end());
    // cout<<"solution global: "<<intsolution<<endl;
    // double intsolutionP{0.};
    // xIntegrateEvalCommand< xEvalUnary<xGetModu<double> > > integrateSolutionP(pressureindb, intsolutionP);  // int one 
    // ApplyCommandOnIntegrationRule(integrateSolutionP, integration_rule_neg, porous.begin(), porous.end());
    // cout<<"solution porous: "<<intsolutionP<<endl;
#endif

// cylinder case
#ifdef CYL
    xEvalUnary<xGetRealPartV> exact_disp_real(Cylinder_displacement);
    // Export(exact_disp_real, pexport, "EXACT_SOL_CYL_DISP", integration_rule_neg, porous.begin(), porous.end());
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > biot_exact_disp_norm(Cylinder_displacement); 
    xEvalUnary<xtool::xIdentity<xtensor::xVector<valtype> > > exact_disp(Cylinder_displacement);   

    xEvalUnary<xGetRealPart> cylinder_exact_real(Cylinder_pressure);
    pexport.setNbSplit(1.*degree);
    Export(cylinder_exact_real, pexport, "EXACT_SOL_PRESS_REAL", integration_rule, all.begin(), all.end()); 
    xEvalUnary<xGetComplex> exact_press(Cylinder_pressure);
#endif

#ifndef THREE_DIM
#ifdef AIRBIOT
//   air-biot square
    // xEvalGradField<xtensor::xSymmetrize<valtype>, xField<valtype> > eval_strain(disp_l);
    // xEvalBinary<xtool::xMult<xtensor::xTensor4<valtype>, xtensor::xTensor2<valtype>, xtensor::xTensor2<valtype>>> stress(hooke, eval_strain);
    // xEvalUnary<xGetNormTensor> stress_norm(stress);
 
    // Export(stress_norm, pexport, "STRESS", integration_rule_neg, porous.begin(), porous.end());

    // xEvalUnary<xGetNormTensor> cylinder_exact_stress_norm(Cylinder_stress); 
    // Export(cylinder_exact_stress_norm, pexport, "EXACT_SOL_STRESS", integration_rule_neg, porous.begin(), porous.end());

    xEvalUnary<xtool::xIdentity<xtensor::xVector<valtype> > > exact_disp(Biot_exact_disp);
    xEvalUnary<xGetRealPartV> biot_exact_disp_real(exact_disp);

    // Export(biot_exact_disp_real, pexport, "EXACT_SOL_DISP", integration_rule_neg, porous.begin(), porous.end());
    
    xEvalUnary<xGetRealPart> biot_exact_real(Biot_exact_press);
#ifndef FEM_STANDARD
    pexport.setNbSplit(1*degree);
#endif
    Export(biot_exact_real, pexport, "EXACT_SOL_PRESS_REAL", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetComplex> exact_press(Biot_exact_press);
    xEvalUnary<xGetRealPartV> disp_total_real(Biot_disp_total);
    pexport.setNbSplit(2*degree);
    // Export(disp_total_real, pexport, "EXACT_DISP_TOTAL", integration_rule, porous.begin(), porous.end());

#endif

#ifdef BIOTBIOT
    xEvalUnary<xtool::xIdentity<xtensor::xVector<valtype> > > exact_disp(Coup_Biot_exact_disp);
    xEvalUnary<xGetRealPartV> coup_biot_exact_disp_real(Coup_Biot_exact_disp);

    // Export(coup_biot_exact_disp_real, pexport, "EXACT_SOL_DISP", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetRealPart> Coup_biot_exact_real(Coup_Biot_exact_press);
    Export(Coup_biot_exact_real, pexport, "EXACT_SOL_PRESS_REAL", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetComplex> exact_press(Coup_Biot_exact_press);
#endif

// common displacement post processing
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > exact_disp_norm(exact_disp);
    xEvalBinary<std::minus<xtensor::xVector<valtype>>> diff_disp(eval_disp, exact_disp);
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > diff_disp_norm(diff_disp);


#ifdef AIRBIOT
    // Export(diff_disp_norm, pexport, "ERROR_ABSOLUTE_DISP", integration_rule_neg, porous.begin(), porous.end());
#endif
#ifdef BIOTBIOT
    // Export(diff_disp_norm, pexport, "ERROR_ABSOLUTE_DISP", integration_rule, all.begin(), all.end());
#endif

//  common  pressure post processing
    xEvalUnary<xGetModu<valtype> > exact_norm(exact_press);
    xEvalBinary<xMinusComplex <valtype, valtype, valtype>> diff_press(eval_press, exact_press);
    xEvalUnary < xGetModu <valtype>> diff_press_norm(diff_press);
    // Export(diff_press_norm, pexport, "ERROR_ABSOLUTE_PRESS", integration_rule_env, all.begin(), all.end());
    // xEvalUnary<xGetRealPart> biot_exact_real(Coup_Biot_exact_press)

// Linear combination error indicator
    // xEvalBinary<xPlusComplexVS <valtype, xtensor::xVector<valtype>, valtype>> diff_sum(diff_press, diff_disp);
    // xEvalBinary<xPlusComplexVS <valtype, xtensor::xVector<valtype>, valtype>> sol_sum(exact_press, exact_disp);
    // xEvalUnary < xGetModu <valtype>> diff_sum_norm(diff_sum);
    // xEvalUnary < xGetModu <valtype>> sol_sum_norm(sol_sum);
    // xEvalUnary <xGetSquare<double>> squa_diff(diff_sum_norm);
    // xEvalUnary <xGetSquare<double>> squa_sol(sol_sum_norm);
 #if defined(XFEM_CONTINUE) || defined(XFEM_DISCONTINUE)
    // xEvalApproxFunction<double> evalappro(enrichment_function);
    // Export(evalappro, pexport, "shape", integration_rule, all.begin(), all.end()); 
#endif


    cout<<"end the post pro"<<endl;
// Acoustic indicator,i.e. reflexion, absorption
    double fluid_density = 1.213;
    xEvalGradField <xtool::xIdentity<xtensor::xVector<valtype> > , xField<valtype> > eval_grad_cplx(press_l);
    xEvalVelocity<xGetComplexV> eval_velocity(eval_grad_cplx, fluid_density, omega);
    // xGetNormalvelocity opNormalvelocity(eval_normal_, fluid_density, classNeg, classPos);
    xEvalBinary <xDivisComplex <valtype, valtype, valtype>> impedence(eval_press, eval_velocity);
    xEvalUnary<xReflection> ref(impedence);
    xEvalUnary<xGetRealPart> imp_real(impedence);
    xEvalUnary<xGetImagPart> imp_imag(impedence);
    // Export(imp_real, pexport, "REF_REAL", integration_rule, all.begin(), all.end());
    // Export(imp_imag, pexport, "REF_IMAG", integration_rule, all.begin(), all.end());
    // xEvalUnary<xAbsorption<valtype>> absorpt(ref);
    // Export(absorpt, pexport, "ABSORPTION", integration_rule, all.begin(), all.end());

#if 0
   cout << "-----test eval jump------" << endl;
    xEvalJumpB< valtype > eval_press_jump(eval_press, &press_l);
    xEvalJumpB< valtype > eval_exact_jump(exact_press);
    // xEvalJump< double > eval_exact_jump(exact_solution_real);

    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_jump(eval_exact_jump, eval_press_jump);
    xEvalUnary < xGetModu <valtype>> mod_diff_jump(diff_jump);
    xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    xEvalUnary < xGetModu <valtype>> mod_press_jump(eval_press_jump);
    xEvalUnary <xGetSquare<double>> sqart_err_jump(mod_diff_jump);
    xEvalUnary <xGetSquare<double>> sqart_exa_jump(mod_exact_jump);
    xIntegrationRulePartition integration_rule_exa(xMesh::getPartition, 5*(degree+1));
    // Export(mod_exact_jump, pexport, "EVAL_EXACT_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(mod_press_jump, pexport, "EVAL_PRESS_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    double intJump1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntJump1(sqart_err_jump, intJump1);  // int one 
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateIntJump1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateIntJump1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    double intJump2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntJump2(sqart_exa_jump, intJump2);
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateIntJump2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateIntJump2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    xEvalBinary< xDivisComplex<double, double, double>> err_rela_jump(sqart_err_jump, sqart_exa_jump);
    // Export(err_rela_jump, pexport, "ERROR_RELA_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    double errorJump = sqrt(intJump1/intJump2)*100.;
    cout << "absolute jump error : " <<intJump1<< endl;
    cout << "exact jump : " <<intJump2<< endl;
    cout << "jump error : " <<errorJump<< endl;
    
    // Export(mod_diff_jump, pexport, "MODULE_ERROR_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    cout << "-----End eval jump-------" << endl;

#if 0
    double intDiffSum{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntDiffSum(squa_diff, intDiffSum);  // int one 
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateIntJump1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateIntDiffSum, integration_rule_exa,  all.begin(), all.end());
    double intSum{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntSum(squa_sol, intSum);
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateIntJump2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateIntSum, integration_rule_exa,  all.begin(), all.end());
    // xEvalBinary< xDivisComplex<double, double, double>> err_rela_jump(sqart_err_jump, sqart_exa_jump);
    // Export(err_rela_jump, pexport, "ERROR_RELA_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    double errorSum = sqrt(intDiffSum/intSum)*100.;
    cout << "relative sum error : " <<errorSum<< endl;
#endif

#endif


    /* ****************************error computation*******************************
    created several operator class in xComplexUtils.h for complex and for err computation, facilite to manager*****
    **********************************************************/
    //start
    cout << "Start Error computaiton" << endl;
#ifndef FEM_STANDARD
// if defined(XFEM_CONTNUE)|defined(XFEM_DISCONTINUE)
//     xIntegrationRulePartition integration_rule_err(xMesh::getPartition, 2*(degree+1));
// else
//     xIntegrationRuleStoredPartition integration_rule_exa(xMesh::getPartition, 5*(degree+1))
    // xValueManagerDist<double> double_manager;
    xField<> err_info_l(&double_manager, xSpaceConstant("ERROR_INFO_PRESS"));
    xValueCreator<xValueError>  creator_err;
    DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());

    xEvalUnary < xGetSquare <double>> square_mod_err(diff_press_norm);
    xEvalUnary < xGetSquare <double>> square_exact_sol(exact_norm);
    xEvalUnary < xGetSquare <double>> square_mod_err_disp(diff_disp_norm);
    xEvalUnary < xGetSquare <double>> square_exact_sol_disp(exact_disp_norm);

#if 0
    double int1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1(square_mod_err, int1);  // int one 
    ApplyCommandOnIntegrationRule(integrateInt1, integration_rule, all.begin(), all.end());
    // ApplyCommandOnIntegrationRule(integrateInt1, integration_rule_exa, all.begin(), all.end());
    // cout<<"absolute error: "<<sqrt(int1)<<endl;
    double int2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2(square_exact_sol, int2);  // int one 
    ApplyCommandOnIntegrationRule(integrateInt2, integration_rule, all.begin(), all.end());
    // ApplyCommandOnIntegrationRule(integrateInt2, integration_rule_exa, all.begin(), all.end());
    double error = 100.*sqrt(int1/int2);
    cout<<"pressure relative error: "<<std::setprecision(10)<<error<<endl;
    xEvalBinary< xDivisComplex<double, double, double>> err_rela(square_mod_err, square_exact_sol);
    // Export(err_rela, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());

#if defined(AIRBIOT)||defined(CYL)
    double intd1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntd1(square_mod_err_disp, intd1);  // int one 
    ApplyCommandOnIntegrationRule(integrateIntd1, integration_rule_neg, porous.begin(), porous.end());
    // cout<<"disp absolute error: "<<sqrt(intd1)<<endl;
    double intd2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntd2(square_exact_sol_disp, intd2);  // int one 
    ApplyCommandOnIntegrationRule(integrateIntd2, integration_rule_neg, porous.begin(), porous.end());
    // double errord = 100.*sqrt(intd1/intd2);
    double errorsum = 100.*sqrt(intd1)/sqrt(intd2);
    cout<<"displacement relative error: "<<std::setprecision(10)<<errorsum<<endl;
    // std::ofstream out("error.txt",ios::app);
    // out << degree<<" "<<complex_manager.size("dofs")<<" "<<"error pressure "<<error<<" error jump "<<errorJump<<endl;
    // out.close();
#endif

    
#ifdef BIOTBIOT
    double intd1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntd1(square_mod_err_disp, intd1);  // int one 
    ApplyCommandOnIntegrationRule(integrateIntd1, integration_rule, all.begin(), all.end());
    // cout<<"disp absolute error: "<<sqrt(intd1)<<endl;
    double intd2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntd2(square_exact_sol_disp, intd2);  // int one 
    ApplyCommandOnIntegrationRule(integrateIntd2, integration_rule, all.begin(), all.end());
    double errord = sqrt(intd1/intd2)*100.;
#endif

#endif

    // Export(diff_press_norm, pexport, "ERROR_ABSOLUTE_PRESS", integration_rule, all.begin(), all.end()); 
    // // integration
    // xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
    //         xEvalUnary < xGetModu <valtype>> > abs_exa_form(exact_norm);
    // xValueError::choice("ENG_EXA");
    // xFillFieldFromZeroForm<> fill_abs_exa(err_info_l, abs_exa_form);
    // ApplyCommandOnIntegrationRule(fill_abs_exa, integration_rule_exa, all.begin(), all.end());

    // xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
    //         xEvalUnary < xGetModu <valtype>> > abs_err_form(diff_press_norm);
    // xValueError::choice("ABS2_EXA");
    // xFillFieldFromZeroForm<> fill_err_exa(err_info_l, abs_err_form);
    // ApplyCommandOnIntegrationRule(fill_err_exa, integration_rule_exa, all.begin(), all.end());

    // xValueError::finalize(mesh->getPartitionManager().getComm());
    // xValueError::choice("REL_EXA");

    // xEvalField<xtool::xIdentity<double> > val_err_info(err_info_l);
    // Export(val_err_info, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());

    // cout<<"Wrinting error file"<<endl;
    // std::ofstream out;
    // if (!proc_id)
    // {
    //   out.open("error.txt",std::ofstream::out);
    //   out << "l2 error of pressure  " << endl;
    //   out << "solution exact is   " << xValueError::total("ENG_EXA") << endl;
    //   out << "err abs exact is  " << xValueError::total("ABS_EXA") << endl;
    //   out << "err rel exact is  " << xValueError::total("REL_EXA") << endl;
    // }
    // xValueError::clear();
    
    // xField<> err_info_inter(&double_manager, xSpaceConstant("ERROR_INFO_INTER"));
    
    // DeclareInterpolation(err_info_inter, creator_err, interface->begin(1), interface->end(1));
    
    // xEvalField<xtool::xIdentity<double> > val_err_interf_info(err_info_inter);

    // // integration
    // xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
    //         xEvalUnary < xGetModu <valtype>> > abs_exa_jump_form(mod_exact_jump);
    // xValueError::choice("ENG_EXA");
    // xFillFieldFromZeroForm<> fill_abs_exa_jump(err_info_inter, abs_exa_jump_form);
    // TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(fill_abs_exa_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    // xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
    //         xEvalUnary < xGetModu <valtype>> > abs_err_jump_form(mod_diff_jump);
    // xValueError::choice("ABS2_EXA");
    // xFillFieldFromZeroForm<> fill_err_jump(err_info_inter, abs_err_jump_form);
    // TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(fill_err_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    // xValueError::finalize(mesh->getPartitionManager().getComm());
    // xValueError::choice("REL_EXA");

    // // xEvalField<xtool::xIdentity<double> > val_err_jump_info(err_info_inter);
    // // Export(val_err_info, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());
    // cout<<"err rel jump: "<<xValueError::total("REL_EXA")*100. << endl;

#endif
#if 0
    xValueError::clear();
    // xValueManagerDist<double> double_manager;
    // xValueCreator<xValueError>  creator_err;
    xField<> err_info_disp(&double_manager, xSpaceConstant("ERROR_INFO_DISP"));
    // xIntegrationRulePartition integration_rule_exa(xMesh::getPartition, 2*(degree+1));
    DeclareInterpolation(err_info_disp, creator_err, porous.begin(), porous.end());
    // DeclareInterpolation(err_info_disp, creator_err, all.begin(), all.end());
    xEvalField<xtool::xIdentity<double> > val_err_disp_info(err_info_disp);

    // integration
    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModuV <xtensor::xVector<valtype> > >,
            xEvalUnary < xGetModuV <xtensor::xVector<valtype> > > > abs_exa_disp_form(biot_exact_disp_norm);
    xValueError::choice("ENG_EXA");
    xFillFieldFromZeroForm<> fill_abs_exa_disp(err_info_disp, abs_exa_disp_form);
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_abs_exa_disp, integration_rule_neg, porous.begin(), porous.end());
    // ApplyCommandOnIntegrationRule(fill_abs_exa_disp, integration_rule_exa, all.begin(), all.end());

    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModuV <xtensor::xVector<valtype>>>,
            xEvalUnary < xGetModuV <xtensor::xVector<valtype> >> > abs_err_disp_form(diff_disp_norm);
    xValueError::choice("ABS2_EXA");
    xFillFieldFromZeroForm<> fill_err_disp(err_info_disp, abs_err_disp_form);
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_err_disp, integration_rule_neg, porous.begin(), porous.end());
    // ApplyCommandOnIntegrationRule(fill_err_disp, integration_rule_exa, all.begin(), all.end());
    xValueError::finalize(mesh->getPartitionManager().getComm());

    xValueError::choice("REL_EXA");

    // xEvalField<xtool::xIdentity<double> > val_err_disp_info(err_info_disp);
    // Export(val_err_disp_info, pexport, "ERROR_RELATIVE_DISP", integration_rule, porous.begin(), porous.end());
    
    // std::ofstream out;
    if (!proc_id)
    {
        // out.open("error.txt",std::ofstream::out);
      out << "l2 error  of dispalcement  " << endl;
      out << "solution exact is   " << xValueError::total("ENG_EXA") << endl;
      out << "err abs exact is  " << xValueError::total("ABS_EXA") << endl;
      out << "err rel exact is  " << xValueError::total("REL_EXA") << endl;
    }

    if (!proc_id)
      out.close();

#endif

#endif
// Division operator for complex, compute the relative_error
    // xEvalBinary< xDivisComplex<double, double, double>> err_rela(square_mod_err, square_exact_sol);
    // xEvalUnary <xGetSqrt<double>> sqart_err_rela(err_rela);

  
    // TagSidePositive(interface);
    // Export(eval_press_jump, pexport, "ERROR_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // UnTagSide(interface);
    
     // Binaryopeator for substraction, compute the difference between exact and FEM
    // xEvalBinary <xMinusComplex <valtype, valtype, valtype>> errJ(eval_press_jump, eval_press_jump_exact);
    //  // xGetModu class get the absollute value
    // xEvalUnary < xGetModu <valtype>> mod_errJ(errJ);
    //  // xGetSquare struc : get Square
    // xEvalUnary < xGetSquare <double>> square_mod_errJ(mod_errJ);

    // xEvalUnary < xGetModu <valtype>> mod_exact_solJ(eval_press_jump_exact);
    // xEvalUnary < xGetSquare <double>> square_exact_solJ(mod_exact_solJ);

    //  // Division operator for complex, compute the relative_error
    // xEvalBinary< xDivisComplex<double, double, double>> err_relaJ(square_mod_errJ, square_exact_solJ);
    // xEvalUnary <xGetSqrt<double>> sqart_err_relaJ(err_relaJ);

    // double int1{0.};
    // xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1(square_mod_errJ, int1);  // int one 
    // ApplyCommandOnIntegrationRule(integrateInt1, integration_rule, interface->begin(1), interface->end(1), xUpperCreator());

    

    // export the line result in text tPRESSURE_REAL.posmax);
    // xtensor::xPoint BBmin, BBmax;
    // mesh->compute_bounding_box(BBmin, BBmax);
    // Trellis_Util::mPoint Plbegin(0.5*(BBmin(0)+BBmax(0))-0.002, BBmin(1), 0.);
    // // Trellis_Util::mPoint Plend(0.5*(BBmin(0)+BBmax(0))-0.002, BBmax(1), 0.);
    // // Trellis_Util::mPoint Prbegin(0.5*(BBmin(0)+BBmax(0))+0.002, BBmin(1), 0.);
    // // Trellis_Util::mPoint Prend(0.5*(BBmin(0)+BBmax(0))+0.002, BBmax(1), 0.);
    // xtensor::xPoint Pbegin(BBmin(0), 0.5*(BBmin(1)+BBmax(1)), 0.);
    // xtensor::xPoint Pend(BBmax(0), 0.5*(BBmin(1)+BBmax(1)), 0.);


    // // lineExporterbt(*mesh, Pbegin, Pmend, biot_exact_stress_norm, 2000, "STRESS_NORM_LINE.txt");
    // // lineExporterbt(*mesh, Pbegin, Pend, eval_press_norm, 2000, "PRESS_NORM_LINE.txt");
    // lineExporterbt(*mesh, Pbegin, Pend, eval_press_real, 2000, "PRESS_XFEMR_LINE.txt");
    // // lineExporterbt(*mesh, Pbegin, Pend, biot_exact_real, 2000, "PRESS_EXACT_LINE.txt");
    // lineExporterbt(*mesh, Pbegin, Pend, eval_press_imag, 2000, "PRESS_XFEMI_LINE.txt");
    // lineExporterbt(*mesh, Pbegin, Pend, eval_press_real, 2000, "PRESS_LINE.txt");
    // // lineExporterbt(*mesh, Pbegin, Pend, biot_exact_real, 2000, "PRESS_LINE.txt");
    // lineExporterbt(*mesh, Pbegin, Pend, eval_disp_realx, 2000, "DISP_REALX_LINE.txt");
        // xtensor::xPoint BBmin, BBmax;
    // auto BB = data.getMesh()->compute_bounding_box();
#ifdef CAR
    cout << "start export solution in line" << endl;
    auto BB = mesh->compute_bounding_box();
    cout<<"Frequency is :"<<freq<<endl;
    xtensor::xPoint Peval(0.897, 0.5*(BB.min(1)+BB.max(1)), 0.);
    std::ofstream out_press("press_point.txt",ios::app);
    out_press <<"freq:"<<" "<<freq<<" "<<"pressure_value"<<" "<<PointExporterbt(*mesh, *mesh, Peval, eval_press_indb)<<endl;

    xtensor::xPoint Pbegin(BB.min(0), 0.5*(BB.min(1)+BB.max(1)), 0.);
    xtensor::xPoint Pend(BB.max(0), 0.5*(BB.min(1)+BB.max(1)), 0.);
    // xtensor::xPoint Pbegin(0.2, -0.6, 0.);
    // xtensor::xPoint Pend(2.6, 0.05, 0.);
    cout << "exporting ..." << endl;
    // lineExporterbt(*mesh, *mesh, Pbegin, Pend, eval_press_indb, 2000, "PRESS_LINE.txt");
#endif
    cout << "End Post Treatment" << endl;
    cout<<"--------------------------------------------------------\n";
//  }
    return;
   
}


template <typename FIELD >
void helmholtzNdbiot :: TreatmentOfEssEnv (const FIELD& listFunctionSpace, xData &data){

  for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it) {
    const xEnv& env = *it;
    string phys = env.Phys; 
    int type = env.Type; 
    int entity = env.Entity;
    //

    //
    if (phys == "DISPLACEMENT_X" || phys == "DISPLACEMENT_Y" || phys == "DISPLACEMENT_Z" || phys == "PRESSURE" ) {       
      if (type == FIX) {
        cout<<"Treat essential BCs for physical entity "<<env.Entity<< " of dim " << env.getDimension()<<endl;
        xClassRegion bc(data.getMesh(), entity, env.getDimension());
        std::complex<double> val {env.getValue(), 0.};
	    DirichletBoundaryCondition (listFunctionSpace, phys, bc.begin(), bc.end(), val);
      }
      else assert(1 == 0);      
    }
} 
return;
}

template <class ASSEMBLER, class FIELD >
void  helmholtzNdbiot :: TreatmentOfMixedBCs(xData &data, ASSEMBLER &assembler, const xIntegrationRule& integration_rule,
                                       const FIELD& press_l, double omega){


    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it) {
        const xEnv& env = *it;
        if (env.Phys == "ELASTIC_SUPPORT") {
            assert(env.Type == FIX);
            cout<<"Treat Robin BCs for physical entity "<<env.Entity<< " of dim " << env.getDimension()<<endl;

            xUniformMaterialSensitivity<valtype> inv_celerity("inv_celerity");
            xFormBilinearWithLaw<xValOperator<xtensor::xCastToComplex<double> >,
                    xEval<valtype>,
                    xValOperator<xtensor::xCastToComplex<double> > ,valtype > mass(inv_celerity);
            assembler.setCoeff(env.getValue() * omega * j);

            xClassRegion gammaR(data.getMesh(), env.Entity, env.getDimension());
            Assemble(mass, assembler, integration_rule, press_l, press_l, gammaR.begin(), gammaR.end(), xUpperAdjacency());
            assembler.setCoeff(1.);
        }
    }

}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNdbiot :: TreatmentOfNatEnvAcoustic   (xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule,
                                          xfem::xEval<xtensor::xVector<std::complex<double> > > &exact_velocity,
                                          double omega_sq, double beta)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;
    

        if(env.Phys == "DISPLACEMENT"){
            // total displacement (one of the surface terms in weak form for pressure field) 
            assert(env.Type == FIX);
            double val = env.getValue();
            xtensor::xVector<valtype>  vec_val({val*cos(beta), 0.}, {val*sin(beta), 0.}, {0.,0.});
            xEvalConstant<xtensor::xVector<valtype>> flux(vec_val);
            // std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<" val: "<<val<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > scaled_flux(flux, eval_normal);
            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > >, std::complex<double> > lin(scaled_flux);
            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), xUpperAdjacency());
        }

        else if(env.Phys == "NORMAL_VELOCITY"){
            assert(env.Type == FIX);
            double val = env.getValue();
            xtensor::xVector<valtype>  vec_val{val};
            xEvalConstant<xtensor::xVector<valtype>> flux(vec_val);
//            xEvalUnary<xtensor::xCastToVectorComplex<double> > flux(vec_val);
            std::cout<<"In TreatmentOfNatEnv, treating Bc for Normal velocity "<<env.Entity<<" of dim "<<env.getDimension()<<" val: "<<val<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);

            xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > scaled_flux(flux, eval_normal);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > >, std::complex<double> > lin(scaled_flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            assembler.setCoeff(-1./(1.213*omega_sq));
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), xUpperAdjacency());
            assembler.setCoeff(1.);
            std::cout<<"End TreatmentOfNatEnv, treating Bc for Normal velocity "<<endl;
        }

        else if(env.Phys == "VELOCITY_X"){
            // exact gradiant of pressure
            assert(env.Type == FIX);

            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            // std::cout<<"In TreatmentOfNatEnv, treating velocity x "<<env.Entity<<" of dim "<<env.getDimension()<<std::endl;
       
            xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);
            
            // env.Entity.print();
            xUniformMaterialSensitivity<valtype> inv_density("inv_density");

            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, flux);
            // xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, exact_velocity);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);
            
            // xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());
            xClassRegion bc(mesh, env.Entity, env.getDimension());
            assembler.setCoeff(1./(omega_sq));
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.);
         
        }
    }
    return;
}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNdbiot :: TreatmentOfNatEnvPEMSolid   (const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, 
                                         xfem::xEval<xtensor::xTensor2<std::complex<double> > > &stress)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;

        // stress BCs one of the surface term for displacement field
        if(env.Phys == "TRACTION_X" || env.Phys == "TRACTION_Y" ){

            assert(env.Type == FIX);
            xtensor::xVector<valtype> val;
            if (env.Phys == "TRACTION_X") val(0) = env.getValue();
            if (env.Phys == "TRACTION_Y") val(1) =  env.getValue();
            std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<" val: "<<val<<std::endl;
            xEvalConstant<xtensor::xVector<valtype> >  flux(val);
            xFormLinearWithLoad<xValOperator<xtensor::xCastToVectorComplex<double> >, xEvalConstant<xtensor::xVector<valtype> >, valtype > lin(flux); 
            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());
            assembler.setCoeff(-1.);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), 
	        xUpperAdjacency()); 
            assembler.setCoeff(1.);
        }

        else if (env.Phys == "STRESS")
        // exact solution of stress for BCs
        {
            assert(env.Type == FIX);
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xTensor2<valtype >, xtensor::xVector<valtype >, xtensor::xVector<valtype > > > flux(stress, eval_normal);
            xFormLinearWithLoad < xValOperator < xtensor::xCastToVectorComplex<double> >,
                                  xEvalBinary < xtool::xMult < xtensor::xTensor2<valtype>, xtensor::xVector<valtype >, xtensor::xVector<valtype> > >, valtype > lin(flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());
            // assembler.setCoeff(-1.);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            // assembler.setCoeff(1.);
        }
    }
    return;
}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNdbiot :: TreatmentOfNatEnvPEMFluid   (const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, 
                                         xfem::xEval<xtensor::xVector<std::complex<double> > > &disp_total)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;

            // exact solution BCs of total displacement
        if (env.Phys == "ELECTRIC_DISPLACEMENT")
        {
            assert(env.Type == FIX);
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(disp_total, eval_normal);

            xFormLinearWithLoad < xValOperator < xtensor::xCastToComplex<double> >,
                                  xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > >, valtype > lin(flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());
            // assembler.setCoeff(-1.);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            // assembler.setCoeff(1.);     
        }
    }
    return;
}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNdbiot :: TreatmentOfNatEnvAcoustic   (xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                           ASSEMBLER & assembler, xfem::xIntegrationRule & integration_rule,
                            xfem::xEval<std::complex<double> > &exact_velocity, 
                            double omega_sq)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;

            // exact solution BCs of total displacement
        if (env.Phys == "VELOCITY_X")
        {
            assert(env.Type == FIX);
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            
  
            // xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);
            
            // env.Entity.print();
            xUniformMaterialSensitivity<valtype> inv_density("inv_density");

            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, exact_velocity);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);
            

            // xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
            //                         xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<valtype >, valtype > >, valtype> lin(flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            // cout<<"Treat VELOCITY_X for "<<env.Entity<<" of dim "<<env.getDimension()<<endl;

            // assembler.setCoeff(1.);
            assembler.setCoeff(1./omega_sq);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.); 
        }
    }
    return;
}







