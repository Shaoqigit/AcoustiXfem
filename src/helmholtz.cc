#include "helmholtz.h"

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
#include "xComplexUtils.h"
#include "exactSolutions.h"
#include "BiotexactSolutions.h"
#include "BiotCoupexactSolutions.h"
#include "CylinderExactSolutions.h"
#include "exactSolutionsGenericInterface.h"
#include "AnalyticalMaterial.h"
#include "acousticMaterial.h"
#include "heParser.h"
#include "nitscheOperators.h"
#include "EnrichmentUtils.h"
#include "xNitscheFile.h"

// PSG

// from NEXFEM
#include "nurbsUtilities.h"
#include "sisl.h"
#include "xSpacePolynomialOLD.h"
#include "orthoEnrichFunction2.h"

#include "xSimpleGeometry.h"

#include <Eigen/Dense>

using namespace std;
using namespace xfem;
using namespace xtool;
using namespace xtensor;
using namespace AOMD;
using namespace xlinalg;
using namespace xcut;
using namespace Eigen;

// #define GENERIC_INTERFACE
#define XFEM_CONTINUE
// #define XFEM_DISCONTINUE
#define SQUARE
// #define NURBS
// #define CYL
// #define THREE_DIM
// #define OCTREE

constexpr double cx{0.};
constexpr double cy{0.};


#ifdef CHECK_WITH_EIGEN
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "xGenericSparseMatrixTraitPolicy.h"
#endif

void lineExporterFE(xMesh &m, xtensor::xPoint p1, xtensor::xPoint p2, xEval<double> &eval, int nbPoints, string filename){

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
      auto elt_uvw_list = m.locateElement(pt);



      if(!elt_uvw_list.empty()){
          auto elt_uvw = elt_uvw_list.front();
        AOMD::mEntity* pe = const_cast<AOMD::mEntity*>(elt_uvw.first);

          xGeomElem geo(pe);
          geo.setUVWForXYZ(pt);
          double value;
          eval(&geo,&geo,value);

          file << curvAbs << " " << pt(0)<< " "<<  pt(1)<< " "<<  pt(2)<< " "<< value<<endl;

        }

      curvAbs += step / length;
    }

  file.close();
}

helmholtzNd::helmholtzNd (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos) : data(data_) {

    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);
    // read materials
    //Give the pulsation to the porous material
    xfem::xAcousticEquivalentFluid::setOmega(parsedinfos.getDouble("Frequency_Hz"));
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_fluid_porous", xfem::xMaterialCreator < xfem::xAcousticEquivalentFluid >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_air", xfem::xMaterialCreator < xfem::xAcousticAir >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_liquid", xfem::xMaterialCreator < xfem::xAcousticLiquid >() );
//    xexport::xExportGmshAsciiSort pexport;
    xexport::xExportGmshAsciiSort pexport;

    //Redirection of all output to individual files
#if 0
    pid = std::to_string(proc_id);
    string fo = "proc_"+pid+"_output.txt";
    FILE * ok = freopen (fo.c_str(),"w",stdout);
    if (!ok) {std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<< std::endl; throw; }
#endif

    xfem::xMesh* mesh=data.getMesh();
    AOMD::mMesh* mmesh=&mesh->getMesh();
    // set  partition with parmetis
    xinterface::aomd::xAttachedDataManagerAOMD < int > target_proc = xmeshtool::setParmetisPartition(*mmesh, mesh->getPartitionManager(),  nb_proc  );
    xmeshtool::migrateEntities(*mmesh, mesh->getPartitionManager(), target_proc);

    xmeshtool::exportGMSH_V2_Dist(*mmesh, mesh->getPartitionManager(), "partmesh");

    //Read zones
    data.ReadZones();
}

void helmholtzNd::TreatmentOfFormulation (xParseData &parsedinfos) {



    cout << "Starting Formulation for Complex2d" << endl;
    cout << "Start creating the function space" << endl;
    xfem::xMesh* mesh=data.getMesh();
    AOMD::mMesh* mmesh=&mesh->getMesh();
    xRegion all(mesh);

    const int degree = parsedinfos.getInt("order_approx_min");
    // xSpacePolynomialBernstein approxSpace("PRESSURE", xSpace::SCALAR, degree);
    xSpacePolynomialLagrange approxSpace("PRESSURE", xSpace::SCALAR, degree);
    // xSpaceLagrange::lag_degree_t order = xSpaceLagrange::DEGREE_TWO;
    // xSpaceLagrange approxSpace("PRESSURE", xSpace::SCALAR, order);


    double freq = parsedinfos.getDouble("Frequency_Hz");
    double omega = 2.*M_PI*freq;
    double beta = parsedinfos.getDouble("H_wave_angle");  // degree
    double wave_angle = beta*M_PI/180.;  // radians
    double radius = parsedinfos.getDouble("Cylinder_radius");
    double sigma_film  = parsedinfos.getDouble("Resistivity");
    double d = parsedinfos.getDouble("Mat_Thickness");
    double alpha_ = parsedinfos.getDouble("Nitsche_fluid");
    double gamma = parsedinfos.getDouble("Average_gamma");
    double h = parsedinfos.getDouble("Element_size");
    double air_density = 1.213;
    valtype v0 = exp(j);
    

// Generic interface parameters
    // double M11 = parsedinfos.getDouble("A");
    // // double M12 = 1. * sigma_film * d;
    // double M12 = parsedinfos.getDouble("B");
    // double M21 = parsedinfos.getDouble("C");
    // double M22 = parsedinfos.getDouble("D");
    double phi_F = parsedinfos.getDouble("Porosity");
    double alpha_F = parsedinfos.getDouble("Tortuosity");
    double lambda_prime_F = parsedinfos.getDouble("Thermal_length");
    double lambda_F = parsedinfos.getDouble("Viscous_length");


    Analytical_MatEF mat_F(omega, beta, phi_F, sigma_film, alpha_F, lambda_prime_F, lambda_F);
    valtype Z = mat_F.rho_eq_til * mat_F.c_eq_til;
    valtype k_a = mat_F.k_a;
    valtype k_ay = k_a*sin(wave_angle);
    valtype k_eqx = sqrt(pow(mat_F.k_eq_til,2)-pow(k_a,2)*pow(sin(wave_angle),2));

    Matrix2cd M1;
    M1(0,0) = cos(k_eqx * d); 
    M1(0,1) = j * omega * mat_F.rho_eq_til/k_eqx * sin(k_eqx * d); 
    M1(1,0) = j * k_eqx/(omega * mat_F.rho_eq_til) * sin(k_eqx * d);
    M1(1,1) = cos(k_eqx * d); 

    // Matrix2cd M1;
    // M1(0,0) = 1.; 
    // M1(0,1) = j * omega * mat_F.rho_eq_til/k_eqx * sin(k_eqx * d);
    // // M1(0,1) = sigma_film * d;
    // M1(1,0) = 1e-6+j*1e-6;
    // M1(1,1) = 1.; 

    Analytical_MatEF mat_F2(omega, beta, 0.72, 87e3, 1.02, 480e-6, 480e-6);
    valtype Z2 = mat_F2.rho_eq_til * mat_F2.c_eq_til;
    valtype k_a2 = mat_F2.k_a;
    valtype k_ay2 = k_a2*sin(beta);
    valtype k_eqx2 = sqrt(pow(mat_F2.k_eq_til,2)-pow(k_a2,2)*pow(sin(beta),2));
    double d2 = 0.0003;

    Matrix2cd M2;
    M2(0,0) = cos(k_eqx2 * d2); 
    M2(0,1) = j * omega * mat_F2.rho_eq_til/k_eqx2 * sin(k_eqx2 * d2); 
    M2(1,0) = j * k_eqx2/(omega * mat_F2.rho_eq_til) * sin(k_eqx2 * d2);
    M2(1,1) = cos(k_eqx2 * d2); 

    Matrix2cd M;
    // M = M2*M1*M2;
    M = M1;
    // cout<<"Film wave number: "<<k_eqx<<endl;
    valtype M11 = M(0,0); 
    valtype M12 = M(0,1); 
    // valtype M12 = 0.;
    
    valtype M21 = M(1,0); 
    valtype M22 = M(1,1); 
    // cout<<"Coeff M_11: "<<M11<<endl;
    // cout<<"Coeff M_12: "<<M12<<endl;
    // cout<<"Coeff M_21: "<<M21<<endl;
    // cout<<"Coeff M_22: "<<M22<<endl;
    cout<<"Pressure jump: "<<-sigma_film *d<<endl;
    // valtype M11 = 1; cout<<"Coeff M_11: "<<M11<<endl;
    // valtype M12 = mat_F.rho_eq_til*omega*d/j; cout<<"Coeff M_12: "<<M12<<endl;
    // valtype M21 = -j*omega*d/mat_F.K_eq_til-j*k_eqx*k_eqx*d/(mat_F.rho_eq_til*omega); cout<<"Coeff M_21: "<<M21<<endl;
    // valtype M22 = 1.; cout<<"Coeff M_22: "<<M22<<endl;



#ifdef NURBS
           SISLCurve *pc=nullptr;

            int kdim=2, kn, kk;
            //  double aepsge=0.000001; /* geometric tolerance */
            double aepsge=1.e-12;
            double aepsco=1.e-10; /* computational tolerance */


            const double radiusCircle=radius;
            xtensor::xPoint centerCircle(cx,cy,0.);

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
            // throw;

            AOMD_Util::Instance()->ex_port("surface.msh",&surfaceAOMD);
            cout<<"create surface"<<endl;

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
            cout<<"create slice"<<endl;

// throw;
            xIntegrationRulePartition integration_ruleC(xMesh::getPartition, 4*degree);
            xIntegrationRuleBasic integration_ruleB(2*degree);
            xexport::xExportGmshAsciiSort pexportA;
            xEvalConstant<double > one(1.);
            // Export(one, pexportA, "ONE", integration_ruleC, all.begin(), all.end());
            // Export(one, pexportA, "ONEB", integration_ruleB, all.begin(), all.end());

            double result=0.;
            xIntegrateEvalCommand<xEvalConstant<double> > integOne(one, result);

            xAcceptAll acc;
            xRegion matrix(data.getMesh());
            // Export(one, pexportA, "ONEM", integration_ruleB, matrix.begin(), matrix.end());
            xIntegrationRulePartition integration_ruleA(xMesh::getPartition, 4*degree);
            // Export(one, pexportA, "AONE", integration_ruleA, matrix.begin(), matrix.end());

            //    xIntegrationRuleStoredBasic integration_ruleSB(2*degp);
            xIntegrationRuleStoredPartition integration_ruleSB(xMesh::getPartition, 6*(degree+1));
            //  xIntegrationRulePartition integration_ruleSB(2*degp);
            //  xIntegrationRulePartition integration_ruleSB(2*degp);

            // ApplyCommandOnIntegrationRule(integOne,integration_ruleSB,all.begin(),all.end());

            // cout<<"Integration result is "<<result<<endl;

            // integOne.reset_result();
            // ApplyCommandOnIntegrationRule(integOne,integration_ruleSB,matrix.begin(),matrix.end());
            // cout<<"Integration result is "<<result<<endl;


            xIntegrationRulePartition integration_rule_envT(xMesh::getPartition,  50);
            xIntegrationRulePartition integration_ruleTest(xMesh::getPartition,  30);

            // attachSubMeshToBnd(data.getMesh());

            // cout<<"construct the spline"<<endl;
            // throw;

            
#endif

    cout<<"Approximation order is "<<degree<<endl;  
    // exact solution for pressure and velocity
#ifdef GENERIC_INTERFACE
    // xPressureTubeNormal exact_pressure(omega, 1., M11, M12, M21, M22);
    xSolObliquePlaneGenericInterface exact_pressure(omega, wave_angle, M11, M12, M21, M22);
    xGradSolObliquePlaneGenericInterface exact_velocity(omega, wave_angle, M11, M12, M21, M22);
    
#endif
#ifdef SQUARE
    xEvalObliqueExactSolutionBouillard2DCplx exact_pressure(omega, wave_angle, sigma_film, d);
    xEvalObliqueSpeedBouillard2DCplx exact_velocity(omega, wave_angle, sigma_film, d);
#endif


#ifdef CYL
    xEvalCylinderExactSolution exact_pressure(omega, radius, sigma_film, d);
    xEvalCylinderExactspeed exact_velocity(omega, radius, sigma_film, d);
#endif

    xValueManagerDist<double> double_manager;
    xfem::xValueManagerDist<valtype> complex_manager;
    xField<valtype> press_l(&complex_manager, approxSpace);

#ifdef NURBS

    OLD::xSpacePolynomialOLDLagrange pressE("PRESSURE", xSpace::SCALAR, degree);
    xSpaceFiltered::filter_t filter(acc);
    xSpaceFiltered::filter_t filter_E = xAcceptSupportHasSubmesh();
    xValKeyExtend key_modifier_E("_Enrichment");
    xScalarFunctionDerivDiscAnalytical enrichment_function_E(xPoint(cx,cy,0), radius);  

    xSpacePolynomialLagrange lagScalar("FCTQ", xSpace::SCALAR, 1);
    xField<> alphaField(&double_manager, lagScalar);
    createDomainRamp cramp(all.begin(),all.end(), data.getMesh(), filter_E, &double_manager, alphaField);

    xEvalField<xIdentity<double> > eval_ramp(alphaField);
    xEvalGradField<xIdentity<xVector<>> > eval_Gramp(alphaField);


    correctedEnrichmentFunction enrichment_function_C(enrichment_function_E, alphaField);//Corrected function
    xSpaceFiltered::filter_t filterC=bind1st(mem_fun(&createDomainRamp::isInDomain), &cramp);

//            xSpaceFiltered enriched_tip_C(xSpaceXFEM(lagrangeE, enrichment_function_E,key_modifier_E), filterC);

//            xSpaceFiltered enriched_tip_C(xSpaceXFEMShifted(lagrangeE, enrichment_function_E,key_modifier_E), filterC);

    //OK
    xSpaceFiltered enriched_tip_C(xSpaceXFEMShiftedCorrectedOLD(pressE, enrichment_function_E,key_modifier_E, &eval_ramp, &eval_Gramp), filterC);
    press_l.insert(enriched_tip_C);

    // xFilteredRegion<xIter, xSpaceFiltered::filter_t>  enrichzone0(all.begin(0), all.end(0), filterC);
    // xFilteredRegion<xIter, xSpaceFiltered::filter_t>  enrichzone1(all.begin(1), all.end(1), filterC);
    // xFilteredRegion<xIter, xSpaceFiltered::filter_t>  enrichzone2(all.begin(2), all.end(2), filterC);
    // Export(one, pexportA, "ENRICHZONEV", integration_ruleC, enrichzone0.begin(), enrichzone0.end());
    // Export(one, pexportA, "ENRICHZONEE", integration_ruleC, enrichzone1.begin(), enrichzone1.end());
    // Export(one, pexportA, "ENRICHZONES", integration_ruleC, enrichzone2.begin(), enrichzone2.end());
    // xEvalApproxFunction<double> evenr(enrichment_function_C);
    // xEvalGradApproxFunction<xVector<>> evGenr(enrichment_function_C);
    // pexportA.setNbSplit(5);
    // Export(evenr, pexportA, "EVENR_", integration_ruleC, all.begin(), all.end());
    // Export(evGenr, pexportA, "EVGENR_", integration_ruleC, all.begin(), all.end());
cout<<"End enrichement"<<endl;
#endif

    //If a levelset was provided
    xClassifyOn classNeg("material_neg");
    xClassifyOn classPos("material_pos");

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
                    else classNeg(elem_sub);

                    // cout<<"CLASS\n";
                }

                // cout<<"----\n";

            }

            // xIntegrationRuleStoredPartition integration_ruleSBM(xMesh::getPartition, xAccept("material_pos"), 10*degree);
            // Export(one, pexportA, "MATRIX", integration_ruleSBM, all.begin(), all.end());

            // xIntegrationRuleStoredPartition integration_ruleSBI(xMesh::getPartition, xAccept("material_neg"), 10*degree);
            // Export(one, pexportA, "INCLUSION", integration_ruleSBI, all.begin(), all.end());
            cout<<"classify\n";
#endif
#ifdef XFEM_CONTINUE
    
        if(lset){
            //--- Create physical surface object
            xcut::xPhysSurfParameter param_surf(classNeg, classPos);
            inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));

            // //--- Setup enriched space
            // //a/ keyword to create the enriched keys
            xValKeyExtend key_modifier("_MAT_ENRICH");
            //b*/ enrichment function : Ridge
            xScalarFunctionDerivDiscXFEM enrichment_function(inclusion->getLevelSet());
            //c*/ filter for ridge function "supportCutStrictlyEltWise", do not cut at bord of element
            xSpaceFiltered::filter_t filter(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictlyEltWise, inclusion, std::placeholders::_1));
            //d/ create enriched space (but ALL the dofs are enriched...)
            xSpaceXFEM space_full(approxSpace, enrichment_function, key_modifier);
            //e/ filter the space in order to keep only the relevant nodes
            xSpaceFiltered enriched(space_full, filter);
            //f/ add the enriched space to the approximation
            press_l.insert(enriched);
        }
#endif
#ifdef XFEM_DISCONTINUE
        if(lset){
            //--- Create physical surface object
            xcut::xPhysSurfParameter param_surf(classNeg, classPos);
            inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));

            //--- Setup enriched space
            //a/ keyword to create the enriched keys
            xValKeyExtend key_modifier("_MAT_ENRICH");
            //b/ enrichment function : Heaviside
            xScalarFunctionDiscXFEM enrichment_function(inclusion->getLevelSet());
            //c/ filter to enrich the necessary dofs only (supportCutStrictlyEltWise is defined in xPhysSurfByTagging)
            xSpaceFiltered::filter_t filter(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictly, inclusion, std::placeholders::_1));
            {
            xEvalConstant<double> one(1.);
            xIntegrationRulePartition integr(xMesh::getPartition,xAccept("material_pos"),1);
            xexport::xExportGmshMPIIO pexport;
            pexport.setNbSplit(2*degree);
            Export(one, pexport, "POROUS2", integr, all.begin(), all.end());
            Export(*lset, pexport,"LSET");
            }


            //d/ create enriched space (but ALL the dofs are enriched...)
            xSpaceXFEM space_full(approxSpace, enrichment_function, key_modifier);
            //e/ filter the space in order to keep only the relevant nodes
            xSpaceFiltered enriched(space_full, filter);
            //f/ add the enriched space to the approximation
            press_l.insert(enriched);

        }
#endif



    xValueCreator < xSingleValue<valtype> >  creator;
    DeclareInterpolation<xField<valtype> >(press_l, creator, all.begin(), all.end());
    // DeclareInterpolation<xField<valtype> >(press_l, creator, matrix.begin(), matrix.end());


    std::vector < string > vars_ess {"PRESSURE"};
    TreatmentOfEssEnv(press_l, data.PhysicalEnv->begin(), data.PhysicalEnv->end(),
                      vars_ess, mesh);

    complex_manager.PrintForDebug("cplx_dcl_"+pid+".dbg");

    // In case some one use this test with Dirichlet BC we need to uniformize fixed values
    finalizeFixedUniform(mesh->getPartitionManager(),complex_manager);
    complex_manager.PrintForDebug("cplx_ess_"+pid+".dbg");

    xDistStateDofCreatorUniform<xMesh::partman_t, xTrue, valtype> snh(mesh->getPartitionManager(),complex_manager, "dofs");
    DeclareState<xField<valtype> >(press_l, snh, all.begin(), all.end());

    std::shared_ptr < xlinalg::xDistIndex > dist_index = finalizeDofUniform(snh);
    cout << "End creating the function space" << endl;

    cout << "Start Building the Linear System" << endl;
    complex_manager.PrintForDebug("cplx_sym_"+pid+".dbg");
    xIntegrationRulePartition integration_rule_env(xMesh::getPartition, 50);
    xIntegrationRulePartition integration_rule(xMesh::getPartition, 2*(degree+1));    //2*degree is the degree to integrate
    
    xlinalg::xDistVector<valtype> b(*dist_index);
    xlinalg::xDistVector<valtype> sol(*dist_index);
    xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),1);
    // xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),false);
    AssembleGraph(graphsym,press_l, all.begin(), all.end() );
    // AssembleGraph(graphsym,press_l, matrix.begin(), matrix.end() );

#ifdef USE_ZMUMPS
    typedef xlinalg::xLinearSystemSolverMumpsDistributed < valtype > solver_t_;
    typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
    // typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
#endif
    Matrix_SDP_t_ A(graphsym);
    xAssemblerBasic < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler(A, b);
    xAssemblerBasicAndTranspose < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler_transpose(A, b);

    TreatmentOfNatEnv(press_l, assembler, integration_rule_env, exact_velocity, 1, -j*exp(j), wave_angle);
    // Treatment of Robin BCs
    // TreatmentOfMixedBCs(data, assembler, integration_rule_env, press_l, omega);


// =====================================Bilinear form===================================================
    // first term 1/rho
    // xEvalConstant<valtype> one(1./1.213*j);
    xUniformMaterialSensitivity<valtype > fluid_inv_density("inv_density");
    xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
            xEval<valtype>,
            xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > diffusive(fluid_inv_density);  // result type is valtype
     // second term -omega^2/c^c * rho
    xUniformMaterialSensitivity<valtype> inv_celerity_density("inv_celerity_square_density");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xValOperator<xtensor::xCastToComplex<double> > ,valtype > mass(inv_celerity_density);
            
#ifdef NURBS
    Assemble(diffusive, assembler, integration_ruleSB, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(-omega*omega);
    Assemble(mass, assembler, integration_ruleSB, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);
#else
    Assemble(diffusive, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(-omega*omega); 
    Assemble(mass, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);
#endif
// source volume point
    // xEvalPointSource pointsource(omega);
    // xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(fluid_inv_density, pointsource);

    // xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
    //                               xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

    // assembler.setCoeff(1.);
    // Assemble(lin, assembler, integration_rule, press_l, all.begin(), all.end());
    // assembler.setCoeff(1.);

  int Dim = mesh->dim();
  cout<<"dimension of problem: "<<Dim<<endl;
    
    xMesh* interface = inclusion->getMesh_bnd();
    if(!lset) throw;
    xEvalNormalFronLevelSet eval_normal_(*lset);
// interface forms
#ifdef XFEM_DISCONTINUE
    xComputeNormalVelocity opNormalvelocity(eval_normal_, fluid_inv_density, classNeg, classPos, 1.);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(gamma, opNormalvelocity);
    xValJumpOperator<valtype > jump();

    valtype alpha = sigma_film*d/(j*omega);
        cout<<"stabilization"<<alpha_<<endl;
    
// interface coefficient

// lambda * [q][p]-delta p
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > penalty;
    // Penalty terms
// [[q]]<p/n>
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > jump_average(opMeanNormalVelocity);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > average_jump(opMeanNormalVelocity);
// <1/rho*q/n><1/rho*p/n>
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average(opMeanNormalVelocity, opMeanNormalVelocity); 

// ========================Genric interface terms===================================
    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
        xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > qnpn;
    
    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
        xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > qnpp;

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
        xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > qppn;

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
        xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > qppp;
if (alpha_ == 0.){
    if (parsedinfos.getInt("Generic_Interface") == 1){

// ========================================Generic interface terms================================
    cout<<"Direct form for Generic interface"<<endl;
    valtype F11 = M22/M12; cout<<"F11: "<<F11*j*omega<<endl;
    valtype F12 = (M12*M21-M11*M22)/M12; cout<<"F12: "<<F12*j*omega<<endl;
    valtype F21 = 1./M12; cout<<"F21: "<<F21*j*omega<<endl;
    valtype F22 = -M11/M12; cout<<"F22: "<<F22*j*omega<<endl;

//  inteface dimension should be adapted to the problem dimension
    assembler.setCoeff(j*omega*F12);
    Assemble(qnpp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(j*omega*F11);
    Assemble(qnpn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(-j*omega*F22);
    Assemble(qppp, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(-j*omega*F21);
    Assemble(qppn, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    }

    else{
    cout<<"Direct form"<<endl;

    assembler.setCoeff(1./alpha);
    Assemble(penalty, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    }


}
else {

    if (parsedinfos.getInt("Generic_Interface") == 1){

    cout<<"Nitsche generic interface terms"<<endl;

    valtype I_11 = M11/M21; cout<<"I_11: "<<I_11/(j*omega)<<endl;
    // valtype I_11 = 1./M21;
    valtype I_12 = -(M11*M22-M21*M12)/M21; cout<<"I_12: "<<I_12/(j*omega)<<endl;
    // valtype I_12 = -M22/M21;
    valtype I_21 = 1./M21; cout<<"I_21: "<<I_21/(j*omega)<<endl;
    valtype I_22 = -M22/M21; cout<<"I_22: "<<I_22/(j*omega)<<endl;

    valtype gamma_ = I_11/(j*omega);
    // valtype lambda_ = 1./(gamma_+h/alpha_);
    valtype lambda_ = alpha_;

    // xGradOperatorWithParamNeg<xComputeNormalVelocity> opNormalVelocityNeg(opNormalvelocity);
    // xGradOperatorWithParamPos<xComputeNormalVelocity> opNormalVelocityPos(opNormalvelocity);

    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opNormalVelocityPos(1., opNormalvelocity);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opNormalVelocityNeg(0., opNormalvelocity);

    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qn_dpn(opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > dqn_pn(opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValOperatorNeg<xtensor::xCastToComplex<double> >, valtype > dqp_pn(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > dqn_dpp(opNormalVelocityNeg, opNormalVelocityPos);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > dqp_dpp(opNormalVelocityPos, opNormalVelocityPos);
    
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > dqp_dpn(opNormalVelocityPos, opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > dqn_dpn(opNormalVelocityNeg, opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xValOperatorNeg<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qn_dpp(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > dqn_pp(opNormalVelocityNeg);
    
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValOperatorPos<xtensor::xCastToComplex<double> >, valtype > dqp_pp(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qp_dpp(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xValOperatorPos<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qp_dpn(opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qj_dpp(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > qj_dpn(opNormalVelocityNeg);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > dqp_pj(opNormalVelocityPos);

    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > dqn_pj(opNormalVelocityNeg);
// ===================================


// qn_dpn
    assembler.setCoeff(-1.);
    Assemble(qn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
// dqn_pn
    assembler.setCoeff(-1.);
    Assemble(dqn_pn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    // assembler.setCoeff(-1.);
    // Assemble(dqn_pj, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    assembler.setCoeff((-I_12)/(j*omega));
    Assemble(dqn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    // dqn_dpn
    assembler.setCoeff((-I_11)/(j*omega));
    Assemble(dqn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    cout<<"coeff physical nitshce: "<<I_12-I_22+I_11-I_21<<endl;

// dqn_dpp
    // assembler.setCoeff((I_12-I_22)/(j*omega));
    // Assemble(dqn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);
    // // dqn_dpn
    // assembler.setCoeff((I_11-I_21)/(j*omega));
    // Assemble(dqn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);


    assembler.setCoeff(1.);
    Assemble(qp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
// dqp_pp
    assembler.setCoeff(1.);
    Assemble(dqp_pp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

// dqn_dpp
    assembler.setCoeff((I_22)/(j*omega));
    Assemble(dqp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    // dqp_dpn
    assembler.setCoeff((I_21)/(j*omega));
    Assemble(dqp_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);



// stabilization terms:
// lambda* qn_pn
//     assembler.setCoeff(1.*lambda_);
//     Assemble(qnpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda qn_dpp
//     assembler.setCoeff(-1.*lambda_*I_12/(j*omega));
//     Assemble(qn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda qn_dpn
//     assembler.setCoeff(-1.*lambda_*I_11/(j*omega));
//     Assemble(qn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_pn
//     assembler.setCoeff(-1.*lambda_*I_12/(j*omega));
//     Assemble(dqp_pn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_dpp
//     assembler.setCoeff(1.*lambda_*I_12*I_12/(j*omega*j*omega));
//     Assemble(dqp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_dpn
//     assembler.setCoeff(1.*lambda_*I_12*I_11/(j*omega*j*omega));
//     Assemble(dqp_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_pn
//     assembler.setCoeff(-1.*lambda_*I_11/(j*omega));
//     Assemble(dqn_pn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_dpp
//     assembler.setCoeff(1.*lambda_*I_11*I_12/(j*omega*j*omega));
//     Assemble(dqn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_dpn
//     assembler.setCoeff(1.*lambda_*I_11*I_11/(j*omega*j*omega));
//     Assemble(dqn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);

// // // ===================================================== Positive part ============================================================


// // // stabilization terms:
// // // lambda* qp_pp
//     assembler.setCoeff(1.*lambda_);
//     Assemble(qppp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda qp_dpp
//     assembler.setCoeff(-1.*lambda_*I_22/(j*omega));
//     Assemble(qp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda qp_dpn
//     assembler.setCoeff(-1.*lambda_*I_21/(j*omega));
//     Assemble(qp_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_pp
//     assembler.setCoeff(-1.*lambda_*I_22/(j*omega));
//     Assemble(dqp_pp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_dpp
//     assembler.setCoeff(1.*lambda_*I_22*I_22/(j*omega*j*omega));
//     Assemble(dqp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqp_dpn
//     assembler.setCoeff(1.*lambda_*I_22*I_21/(j*omega*j*omega));
//     Assemble(dqp_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_pp
//     assembler.setCoeff(-1.*lambda_*I_21/(j*omega));
//     Assemble(dqn_pp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_dpp
//     assembler.setCoeff(1.*lambda_*I_21*I_22/(j*omega*j*omega));
//     Assemble(dqn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);
// // lambda dqn_dpn
//     assembler.setCoeff(1.*lambda_*I_21*I_21/(j*omega*j*omega));
//     Assemble(dqn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//     assembler.setCoeff(1.);



//  common stabilization terms
    // assembler.setCoeff(alpha_);
    // Assemble(penalty, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_12-I_22)/(j*omega));
    // Assemble(qj_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_11-I_21)/(j*omega));
    // Assemble(qj_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_12-I_22)/(j*omega));
    // Assemble(dqp_pj, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_12-I_22)*(I_12-I_22)/(j*omega*j*omega));
    // Assemble(dqp_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_12-I_22)*(I_11-I_21)/(j*omega*j*omega));
    // Assemble(dqp_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_11-I_21)/(j*omega));
    // Assemble(dqn_pj, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_11-I_21)*(I_12-I_22)/(j*omega*j*omega));
    // Assemble(dqn_dpp, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);

    // assembler.setCoeff(alpha_*(I_11-I_21)*(I_11-I_21)/(j*omega*j*omega));
    // Assemble(dqn_dpn, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    // assembler.setCoeff(1.);
}
else{
    // ==================== Nitsche pressure jump==========================
    cout<<"Nitsche pressure jump"<<endl;
    cout<<"pressure jump: "<<alpha/(j*omega)<<endl;
    cout<<"stabilization lambda is used"<<endl;
        
    // alpha_ = gamma_max * 16.;
    valtype lambda = 1./(alpha+1./alpha_);
    assembler.setCoeff(lambda);
    Assemble(penalty, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(1.*lambda*alpha);
    Assemble(jump_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(1.*lambda*alpha);
    Assemble(average_jump, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
       
    assembler.setCoeff(lambda*alpha*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    
    // // [q]*average
    assembler.setCoeff(-1.);
    Assemble(jump_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    Assemble(average_jump, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
// // //  - delta*double_average
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

// +delta * double_average
    assembler.setCoeff(1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(Dim-1), interface->end(Dim-1), xUpperCreator());
    assembler.setCoeff(1.);

    // ============================================= element wise stabilization========================================================================//

}


    cout << "End bilinear Nitsche terms" << endl;
}
    

#endif


    cout << "End Building the Linear System" << endl;

    if(parsedinfos.getInt("debug") == 1){
       std::ofstream mata_("matA_.mm");  
       A.printMatrixMarket(mata_);

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

// throw;
#ifdef USE_ZMUMPS
   solver_t_ solver;

    //#ifdef WITH_DMUMPS
    //    solver_.setParam(xlinalg::xLinearSystemSolverMumpsBase::INTER,0,1);
    solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,14,200);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::INTER,0,1);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,1,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,2,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,3,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,4,2);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,9,1);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,11, 1);
    //#endif
    solver.connectMatrix(A);
    b.switchInsertModeOff();
    solver.solve(b, sol);
    sol.switchInsertModeOff();
    sol.switchToGlobalValue();
#endif
    //    cout<<"------sol------\n";
    //    sol.Printdata(cout);

    Visit(xWriteSolutionVisitor < xlinalg::xDistVector<valtype> >(sol),
          complex_manager.begin("dofs"),
          complex_manager.end("dofs"));
    cout << "End Solving the Linear System" << endl;

// ==========================================================Post treatmeant=======================================
    cout << "Start Post Treatment" << endl;
    complex_manager.PrintForDebug("res"+pid+".dbg");
    xexport::xExportGmshMPIIO pexport;
    xEvalField < xGetRealPart, xField<valtype> > eval_press_real(press_l);
    xEvalField <xGetComplex, xField<valtype>> eval_press(press_l);
#ifdef XFEM_DISCONTINUE
    pexport.setNbSplit(2.*degree);
#endif
    Export(eval_press_real, pexport, "PRESSURE_REAL", integration_rule, all.begin(), all.end());
//     xEvalField < xGetImagPart, xField<valtype> > eval_press_imag(press_l);
//    Export(eval_press_imag, pexport, "PRESSURE_IMAG", integration_rule, all.begin(), all.end());
    xEvalField <xGetModu<valtype>, xField<valtype>> eval_press_norm(press_l);
    Export(eval_press_norm, pexport, "PRESSURE_MODULE", integration_rule, all.begin(), all.end());
//     xEvalField <xGetPhase<valtype>, xField<valtype>> eval_press_phase(press_l);
//     Export(eval_press_phase, pexport, "PRESSURE_ARG", integration_rule, all.begin(), all.end());

#ifndef THREE_DIM

    xEvalGradField <xtool::xIdentity<xtensor::xVector<valtype> > , xField<valtype> > eval_grad_cplx(press_l);
    xEvalVelocity<xGetComplexV> eval_velocity(eval_grad_cplx, air_density, omega);
//    xEvalUnary<xGetRealPartV> velocity_real(eval_grad_cplx);

    xEvalBinary<xMult<valtype, xtensor::xVector<valtype>, xtensor::xVector<valtype>> > velocity_flux(fluid_inv_density, eval_grad_cplx);
    xEvalBinary<xMult<valtype, xtensor::xVector<valtype>, xtensor::xVector<valtype>> > exact_velocity_flux(fluid_inv_density, exact_velocity);
    xEvalUnary<xGetRealPartV> velocity_real(velocity_flux);
    xEvalUnary<xGetRealPartV> exact_velocity_real(exact_velocity_flux);
    xEvalUnary<xGetModuV<xtensor::xVector<valtype>>> velocity_module(velocity_flux);
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > exact_velocity_module(exact_velocity_flux);
    // xGetNormalvelocity 
    // xEvalUnary<xGetNormalvelocity> velocity_real(velocity_flux);
    // pexport.setNbSplit(degree);
    // TagSideNegative(interface);
    // Export(velocity_module, pexport, "NORMAL_VELOCITY", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // TagSideNegative(interface);
    // Export(exact_velocity_module, pexport, "EXACT_NORMAL_VELOCITY", integration_rule, all.begin(), all.end());
    // xEvalUnary<xGetRealPartV> velocity_exact(exact_velocity);

    
//    xEvalVelocity<xGetRealPartV> eval_velocity_real(eval_grad_cplx, air_density, omega);
   Export(exact_velocity_real, pexport, "EXACT_VELOCITY_REAL", integration_rule, all.begin(), all.end());
   Export(velocity_real, pexport, "VELOCITY_REAL", integration_rule, all.begin(), all.end());

    // pexport.setNbSplit(2*degree);
    xEvalUnary<xGetRealPart> exact_solution_real(exact_pressure);
    Export(exact_solution_real, pexport, "EXACT_SOL_REAL", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetModu<valtype>> exact_solution_norm(exact_pressure);
    Export(exact_solution_norm, pexport, "EXACT_SOL_MODULE", integration_rule, all.begin(), all.end());

    xEvalUnary<xGetComplex> exact_solution(exact_pressure);

//     xEvalUnary<xGetImagPart> exact_solution_imag(exact_pressure);
// //    Export(exact_solution_imag, pexport, "EXACT_SOL_IMAG", integration_rule, all.begin(), all.end());
    // xEvalUnary<xGetModu<valtype>> exact_solution_norm(exact_pressure);
    // Export(exact_solution_norm, pexport, "EXACT_SOL_MODULE", integration_rule, all.begin(), all.end());
//     xEvalUnary<xGetPhase<valtype>> exact_solution_phase(exact_pressure);
//     Export(exact_solution_phase, pexport, "EXACT_SOL_ARG", integration_rule, all.begin(), all.end());

    
    // Export(exact_solution_real, pexport, "EXACT_SOL_BND_POS_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(eval_press_real, pexport, "PRESSURE_BND_POS_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());

    // TagSideNegative(interface);
    // Export(exact_solution_real, pexport, "EXACT_SOL_BND_NEG_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(eval_press_real, pexport, "PRESSURE_BND_NEG_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // UnTagSide(interface);

    // herited xtool::identity to get complex value xGetComplex
    // xEvalField <xGetComplex, xField<valtype>> eval_press(press_l);
    // xEvalUnary<xGetComplex> exact_solution(exact_pressure);

    xEvalBinary <xDivisComplex <valtype, valtype, valtype>> impedence(eval_press, eval_velocity);
    xEvalUnary<xGetRealPart> impedence_real(impedence);
    // Export(impedence_real, pexport, "IMPEDENCE", integration_rule, all.begin(), all.end());
    // just in case where has a normal incidence
    xEvalUnary<xReflection> ref(impedence);
    xEvalUnary<xGetModu<valtype>> ref_abs(ref);
    // xEvalUnary<xGetPhase<valtype>> ref_angle(ref);
    // Export(ref_abs, pexport, "REFLECTION_REAL", integration_rule, all.begin(), all.end());
    // Export(ref_angle, pexport, "REFLECTION_ANGLE", integration_rule, all.begin(), all.end());
    xEvalUnary<xAbsorption<valtype>> absorpt(ref);
    pexport.setNbSplit(degree);
    // Export(absorpt, pexport, "ABSORPTION", integration_rule, all.begin(), all.end());

    // cout << "-----test eval jump------" << endl;

    xEvalJump< valtype > eval_press_jump(eval_press, &press_l);
    // // xEvalJump< double > eval_press_jump(eval_press_real);
    xEvalJump< valtype > eval_exact_jump(exact_pressure);
    // // xEvalJump< double > eval_exact_jump(exact_solution_real);

    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_jump(eval_press_jump, eval_exact_jump);
    xEvalUnary < xGetModu <valtype>> mod_diff_jump(diff_jump);
    xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    xEvalUnary <xGetSquare<double>> sqart_err_jump(mod_diff_jump);
    xEvalUnary <xGetSquare<double>> sqart_exa_jump(mod_exact_jump);
    // xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    xEvalUnary < xGetModu <valtype>> mod_press_jump(eval_press_jump);
    // Export(mod_diff_jump, pexport, "EVAL_MOD_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // // Export(mod_diff_jump, pexport, "EVAL_PRESS_JUMP_TEST", integration_rule, all.begin(), all.end());
    // Export(mod_exact_jump, pexport, "EVAL_EXACT_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
// 
    // Export(mod_press_jump, pexport, "EVAL_PRESS_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());



    /* ****************************error computation*******************************
    created several operator class in xComplexUtils.h for complex and for errr computation, facilite to manager*****
    **********************************************************/
    //start
#if 1
    cout << "Start Error computation" << endl;

    
    xField<> err_info_l(&double_manager, xSpaceConstant("ERROR_INFO"));
    xValueCreator<xValueError>  creator_err;
    DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());
    xIntegrationRulePartition integration_rule_exa(xMesh::getPartition, 4*(degree+1));

    // Binaryopeator for substraction, compute the difference between exact and FEM
    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff(eval_press, exact_solution);
    // xGetModu class get the absollute value
    xEvalUnary < xGetModu <valtype>> mod_diff(diff);
    // xGetSquare struc : get Square
    xEvalUnary < xGetSquare <double>> square_mod_err(mod_diff);
    xEvalUnary < xGetModu <valtype>> mod_exact_sol(exact_solution);
    xEvalUnary < xGetSquare <double>> square_exact_sol(mod_exact_sol);



    double int1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1(square_mod_err, int1);  // int one 
#ifdef NURBS
    ApplyCommandOnIntegrationRule(integrateInt1, integration_ruleSB, all.begin(), all.end());
#else
    ApplyCommandOnIntegrationRule(integrateInt1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
#endif

    cout<<"absolute error: "<<sqrt(int1)<<endl;
    double int2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2(square_exact_sol, int2);  // int one 
#ifdef NURBS
    ApplyCommandOnIntegrationRule(integrateInt2, integration_ruleSB, all.begin(), all.end());
#else
     ApplyCommandOnIntegrationRule(integrateInt2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
#endif
    double error = 100.*sqrt(int1/int2);
    // cout<<"exact solution: "<<sqrt(int2)<<endl;
    cout<<"relative error: "<<error<<endl;

// *************************************computation of normal flux on two sides

    xEvalNormal eval_normal_real;
    xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal_inter(eval_normal_);
    xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > normal_exact_flux(exact_velocity, eval_normal_inter);
    xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > normal_flux(eval_grad_cplx, eval_normal_inter);
    
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > flux_norm(eval_grad_cplx);
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > exact_flux_norm(exact_velocity);
    xEvalBinary<std::minus<xtensor::xVector<valtype>>> diff_flux(eval_grad_cplx, exact_velocity);
    xEvalUnary<xGetModuV<xtensor::xVector<valtype> > > diff_flux_norm(diff_flux);
    xEvalUnary < xGetSquare <double>> square_diff_flux_norm(diff_flux_norm);
    xEvalUnary < xGetSquare <double>> square_exact_flux_norm(exact_flux_norm);
    
    Export(diff_flux_norm, pexport, "diff_flux", integration_rule, all.begin(), all.end());
    // xEvalUnary < xGetModu <valtype>> mod_exa_flux(normal_exact_flux);
    // xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_flux(normal_flux, normal_exact_flux);
    // xEvalUnary < xGetModu <valtype>> mod_diff_flux(diff_flux);
    // xEvalUnary < xGetSquare <double>> square_exact_flux(mod_exa_flux);
    // xEvalUnary < xGetSquare <double>> square_diff_flux(mod_diff_flux);

    // integration_rule_partition_in
    double int_fluxp1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxp1(square_diff_flux_norm, int_fluxp1);
    TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(integrateInt_fluxp1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    double int_fluxp2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxp2(square_exact_flux_norm, int_fluxp2);
    TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(integrateInt_fluxp2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    // 
    double int_fluxn1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxn1(square_diff_flux_norm, int_fluxn1);
    TagSideNegative(interface);
    ApplyCommandOnIntegrationRule(integrateInt_fluxn1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    double int_fluxn2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxn2(square_exact_flux_norm, int_fluxn2);
    TagSideNegative(interface);
    ApplyCommandOnIntegrationRule(integrateInt_fluxn2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    double error_flux_pos = 100.*sqrt(int_fluxp1/int_fluxp2);
    double error_flux_neg = 100.*sqrt(int_fluxn1/int_fluxn2);
    cout<<"relative error flux pos side: "<<error_flux_pos<<endl;
    cout<<"relative error  flux neg side: "<<error_flux_neg<<endl;
    cout<<"relative average flux error: "<<(error_flux_pos+error_flux_neg)/2.<<endl;
    std::ofstream out_press("error_press.txt",ios::app);
    out_press << degree<<" "<<complex_manager.size("dofs")<<" "<<"error pressure "<<error<<endl;

#ifdef XFEM_DISCONTINUE
    double intJump1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntJump1(sqart_err_jump, intJump1);  // int one 
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(integrateIntJump1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    
    double intJump2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateIntJump2(sqart_exa_jump, intJump2);
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(integrateIntJump2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    // xEvalBinary< xDivisComplex<double, double, double>> err_rela_jump(sqart_err_jump, sqart_exa_jump);
    // Export(err_rela_jump, pexport, "ERROR_RELA_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    double errorJump = sqrt(intJump1/intJump2)*100.;
    cout << "jump error : " <<errorJump<< endl;
    // cout << "absolute jump : " <<intJump1<< endl;
    
    // Export(mod_errJ, pexport, "MODULE_ERROR_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    cout << "-----End eval jump-------" << endl;
#endif

#endif

#if 1
    xValueError::clear();
    
    xField<> err_info_inter(&double_manager, xSpaceConstant("ERROR_INFO_INTER"));
    
    DeclareInterpolation(err_info_inter, creator_err, interface->begin(1), interface->end(1));
    
    xEvalField<xtool::xIdentity<double> > val_err_interf_info(err_info_inter);

    // integration
    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModuV <xtensor::xVector<valtype>>>,
            xEvalUnary < xGetModuV <xtensor::xVector<valtype>>> > abs_exa_jump_form(exact_flux_norm);
    xValueError::choice("ENG_EXA");
    xFillFieldFromZeroForm<> fill_abs_exa_jump(err_info_inter, abs_exa_jump_form);
    TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_abs_exa_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModuV <xtensor::xVector<valtype>>>,
            xEvalUnary < xGetModuV <xtensor::xVector<valtype>>> > abs_err_jump_form(diff_flux_norm);
    xValueError::choice("ABS2_EXA");
    xFillFieldFromZeroForm<> fill_err_jump(err_info_inter, abs_err_jump_form);
    TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_err_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    xValueError::finalize(mesh->getPartitionManager().getComm());
    xValueError::choice("REL_EXA");

    // xEvalField<xtool::xIdentity<double> > val_err_jump_info(err_info_inter);
    // Export(val_err_info, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());
    cout<<"err rel jump: "<<xValueError::total("REL_EXA")*100. << endl;
    std::ofstream out;
    if (!proc_id)
    {
      out.open("error.txt",std::ofstream::out);
      out << "energy error  of interface   " << endl;
      out << "solution exact is   " << xValueError::total("ENG_EXA") << endl;
      out << "err abs exact is  " << xValueError::total("ABS_EXA") << endl;
      out << "err rel exact is  " << xValueError::total("REL_EXA") << endl;
    }


    if (!proc_id){
      out.close();
    //   out2.close();
    }
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

    //export the line result in text to compare with 1D results
    xtensor::xPoint BBmin, BBmax;
    auto BB = data.getMesh()->compute_bounding_box();


    xtensor::xPoint Pbegin(BB.min(0), 0.5*(BB.min(1)+BB.max(1)), 0.);
    xtensor::xPoint Pend(BB.max(0), 0.5*(BB.min(1)+BB.max(1)), 0.);

    lineExporterFE(*mesh, Pbegin, Pend, eval_press_real, 2000, "PRESS_LINE.txt");
#endif
    cout << "End Post Treatment" << endl;


    cout<<"--------------------------------------------------------\n";



    return;
}



template <class FIELD , class ITER>
void helmholtzNd :: TreatmentOfEssEnv(const FIELD& fct, ITER itenv, ITER endenv,
                                    const std::vector<string>& vars, xMesh* mesh)
{
    //	const bool debug = xdebug_flag;
    for (; itenv != endenv; ++itenv)
    {
        const xEnv& env = *itenv;
        string phys = env.Phys;
        int type = env.Type;
        int entity = env.Entity;
        if (find(vars.begin(), vars.end(), phys) != vars.end())
        {
            if (type == FIX) {
                xClassRegion bc(mesh, entity, env.getDimension());
                std::complex<double> val = exp(-j*100.*env.getValue());
                DirichletBoundaryCondition (fct, phys, bc.begin(), bc.end(), val);
            }
        }
    }
    return;
}



template <class ASSEMBLER, class FIELD >
void  helmholtzNd :: TreatmentOfMixedBCs(xData &data, ASSEMBLER &assembler, const xIntegrationRule& integration_rule,
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
            assembler.setCoeff(env.getValue() * omega * j/340.);

            xClassRegion gammaR(data.getMesh(), env.Entity, env.getDimension());
            Assemble(mass, assembler, integration_rule, press_l, press_l, gammaR.begin(), gammaR.end(), xUpperAdjacency());
            assembler.setCoeff(1.);
        }
    }

}



template < typename ASSEMBLER, typename FIELD >
void helmholtzNd :: TreatmentOfNatEnv   (const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, xEval<xtensor::xVector<std::complex<double> > > &exact_velocity,
                                         std::complex<double> minusjrhock, std::complex<double> minusjomegaexp, double beta)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;
        if (env.Phys == "ELECTRIC_DISPLACEMENT")
        {
            assert(env.Type == FIX);
            // std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);
            // cout<<"assemble exact"<<endl;
            xUniformMaterialSensitivity<valtype> inv_density("inv_density");
            // xEvalInvDensity<valtype> ev_inv_density(inv_density);
            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, flux);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            assembler.setCoeff(minusjrhock);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.);

        }else if(env.Phys == "NORMAL_VELOCITY"){

            assert(env.Type == FIX);
            double val = env.getValue();
            xtensor::xVector<valtype>  vec_val({val*cos(beta)}, {val*sin(beta)}, {0.,0.});
            xEvalConstant<xtensor::xVector<valtype>> flux(vec_val);
//            xEvalUnary<xtensor::xCastToVectorComplex<double> > flux(vec_val);
            // std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<" val: "<<val<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);

            xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > scaled_flux(flux, eval_normal);
            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > >, std::complex<double> > lin(scaled_flux);
//            xEvalUnary<xtool::xIdentity<double> > scaled_flux(flux);
//            xFormLinearWithLoad<xValOperator<xtool::xIdentity<double> >,
//                      xEvalUnary<xtool::xIdentity<double> >, std::complex<double>> lin(scaled_flux); //Shape functions are assumed double
            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            assembler.setCoeff(1./1.213);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), xUpperAdjacency());
            assembler.setCoeff(1.);

        }
    }

    return;
}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNd :: TreatmentOfNatEnv   (const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, xEval<std::complex<double> > &exact_velocity,
                                         std::complex<double> minusjrhock, std::complex<double> minusjomegaexp, double beta)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;
        if (env.Phys == "ELECTRIC_DISPLACEMENT")
        {
            assert(env.Type == FIX);
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            // xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);

            xUniformMaterialSensitivity<valtype> inv_density("inv_density");
            xEvalInvDensity<valtype> ev_inv_density(inv_density);
            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, exact_velocity);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

            xClassRegion bc(data.getMesh(), env.Entity, env.getDimension());

            assembler.setCoeff(minusjrhock);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.);

        }
    }

    return;
}


xEvalPointSource::xEvalPointSource(double omega_): omega(omega_), o(0.,0.,0.){}
void xEvalPointSource::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const {
         xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    std::complex<double> j{0.,1.};
    res = 1e3*exp(-100.*100.*((x-0.)*(x-0.)+(y-0.1)*(y-0.1)));
    return;
 }
