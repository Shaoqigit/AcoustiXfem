#include "helmholtzPSG.h"
#include <cmath>
#include <iomanip>

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
#include "acousticMaterial.h"
#include "AnalyticalMaterial.h"
#include "heParser.h"
#include "nitscheOperators.h"
#include "EnrichmentUtils.h"
#include "xNitscheFile.h"

// PSG
#include "hoGeometricalModel_hoOctree.h"
// #include "hoGeometricalModel.h"
#include "xPhysSurfByTagging.h"

// from NEXFEM
#include "nurbsUtilities.h"
#include "sisl.h"
#include "xSpacePolynomialOLD.h"
#include "orthoEnrichFunction2.h"

#include "xSimpleGeometry.h"

using namespace std;
using namespace xfem;
using namespace xtool;
using namespace xtensor;
using namespace AOMD;
using namespace xlinalg;
using namespace xcut;
using namespace xgeom;


#define XFEM_DISCONTINUE
#define SQUARE
// #define OCTREE
// #define CYL
#define HEL
#define GNITSCHE
// #define CAR
// #define FEM_STANDARD
// #define CHECK_WITH_EIGEN


constexpr double cx{0.};
constexpr double cy{0.};


// #ifdef CHECK_WITH_EIGEN
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "xGenericSparseMatrixTraitPolicy.h"
#include <Eigen/IterativeLinearSolvers>

// #endif


double PointExporter(xMesh &m, xMesh &mG, xtensor::xPoint p, xEval<double> &eval){


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


void lineExporter(xMesh &m, xMesh &mG, xtensor::xPoint p1, xtensor::xPoint p2, xEval<double> &eval, int nbPoints, string filename){

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


helmholtzNdPSG::helmholtzNdPSG (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos) : data(data_){


    // std::cout<<"initial the construct"<<std::endl;
    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);
    // read materials
    //Give the pulsation to the porous material
    xfem::xAcousticEquivalentFluid::setOmega(parsedinfos.getDouble("Frequency_Hz"));
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_fluid_porous", xfem::xMaterialCreator < xfem::xAcousticEquivalentFluid >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_air", xfem::xMaterialCreator < xfem::xAcousticAir >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_liquid", xfem::xMaterialCreator < xfem::xAcousticLiquid >() );
    xfem::xMaterialManagerSingleton::instance().registerMaterial("acoustic_biot", xfem::xMaterialCreator< xfem::xAcousticPEMBiot>() );
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

void helmholtzNdPSG::TreatmentOfFormulation (xParseData &parsedinfos, xMinus m2) {

#ifdef OCTREE

    cout<<"Octree mesh construct" <<endl;
    xfem::xMesh* meshO=data.getMesh(); // read the original mesh from data
    cout<<"Read original mesh" <<endl;
    // AOMD::mMesh* mmeshO=&meshO->getMesh();

    const int Dim = 2;
    hoGeometricalModel_hoOctree geoModel(Dim);

    xIntegrationRuleBasic integ_basic(2);
    xIntegrationRulePartition integ_part(xMesh::getPartition, 2);
    xIntegrationRuleBasic integ_basicPos(xAccept("material_pos"), 2);
    xIntegrationRulePartition integ_partPos(xMesh::getPartition, xAccept("material_pos"), 2);
    xIntegrationRulePartition integ_partNeg(xMesh::getPartition, xAccept("material_neg"), 2);

    const int nLevels = parsedinfos.getInt("levelG");;
    // auto lset_ = [](std::array<double,3> xyz) -> double { return sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) - 0.12;};
    // auto lset_ = [](std::array<double,3> xyz) -> double { return sqrt(xyz[0]*xyz[0] + (xyz[1]-0.25)*(xyz[1]-0.25)) - 0.12;};
    // auto lset_ = [](std::array<double,3> xyz) -> double { return -xyz[0]+0.;};
    // auto lset_ = [](std::array<double,3> xyz) -> double { return -(xyz[0]*xyz[0] + xyz[1]*xyz[1] - 0.1*0.1);};
    auto lset_ = [](std::array<double,3> xyz) -> double { return xyz[0]-0. ;};
        // car
    // auto lset_ = [&m2] (std::array<double,3> xyz) -> double { xPoint p(xyz[0], xyz[1], xyz[2]); return m2(p);};
    geoModel.createMeshGFromLs(*meshO, lset_, nLevels);
    geoModel.createLsG(lset_);
    // geoModel.createPhysicalSurface("material_pos", "material_neg");
    geoModel.createPhysicalSurface("material_neg", "material_pos");
    geoModel.classifyElementsC();
    geoModel.preprocessFilters();

{
    xEvalConstant<double> one(1.);
    xexport::xExportGmshAsciiSort pexport;

    Export(one, pexport, "PARTPOS", integ_partPos, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));
    Export(one, pexport, "PARTNEG", integ_partNeg, geoModel.getMeshC()->begin(Dim), geoModel.getMeshC()->end(Dim));



    // Export(one, pexport, "EDGE", integ_part, geoModel.getMeshG()->begin(1), geoModel.getMeshG()->end(1));
    cout<<"export boundary"<<endl;

}

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

    cout<<"initial the problem parameter"<<endl;
    double freq = parsedinfos.getDouble("Frequency_Hz");
    double omega = 2.*M_PI*freq;
// for (double freq=100.; freq<=750.; freq=freq+5.){
//     double omega = 2.*M_PI*freq;

    
    double beta = parsedinfos.getDouble("H_wave_angle");  // degree
    double wave_angle = beta*M_PI/180.;  // radians
    double radius = parsedinfos.getDouble("Cylinder_radius");
    double sigma_film  = parsedinfos.getDouble("Resistivity");
    double d = parsedinfos.getDouble("Mat_Thickness");
    double alpha_ = parsedinfos.getDouble("Nitsche_fluid");
    double gamma = parsedinfos.getDouble("Average_gamma");
    double h = parsedinfos.getDouble("Element_size");
    double air_density = 1.213;
    // valtype v0 = exp(j);


#ifdef CYL
// exact solution initialization
    xEvalCylinderExactSolution exact_pressure(omega, radius, sigma_film, d);
    xEvalCylinderExactspeed exact_velocity(omega, radius, sigma_film, d);
#endif

#ifdef SQUARE
    xEvalObliqueExactSolutionBouillard2DCplx exact_pressure(omega, wave_angle, sigma_film, d);
    xEvalObliqueSpeedBouillard2DCplx exact_velocity(omega, wave_angle, sigma_film, d);
#endif

    cout << "Start creating the function space" << endl;
    
    const int degree = parsedinfos.getInt("order_approx_min");
    xSpacePolynomialBernstein press("PRESSURE", xSpace::SCALAR, degree);
    // xSpacePolynomialLagrange press("PRESSURE", xSpace::SCALAR, degree);

    cout<<"Approximation order is "<<degree<<endl;  

    xValueManagerDist<double> double_manager;
    xfem::xValueManagerDist<valtype> complex_manager;
    xField<valtype> press_l(&complex_manager, press);

#ifndef FEM_STANDARD 
    //a/ keyword to create the enriched keys
    xValKeyExtend key_modifier("_MAT_ENRICH");


#ifndef OCTREE
    xClassifyOn classNeg("material_neg");
    xClassifyOn classPos("material_pos");
    //     for(xIter ita = all.begin(1); ita != all.end(1); ++ita){
    //      classNeg(*ita);
    //  }
    xcut::xPhysSurfParameter param_surf(classNeg, classPos); // in and out
    inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));
    xScalarFunctionDiscXFEM enrichment_function(inclusion->getLevelSet());
    xSpaceFiltered::filter_t filter(std::bind(&xcut::xPhysSurfByTagging::supportCutStrictly, inclusion, std::placeholders::_1));
#else
    xScalarFunctionDiscXFEM_ enrichment_function(*geoModel.getLsG());
    //c/ filter to enrich the necessary dofs only (supportCutStrictlyEltWise is defined in xPhysSurfByTagging)
    auto filter = bind1st(mem_fun(&hoGeometricalModel_hoOctree::filterSupportCutStrictly), &geoModel);

#endif

    //b/ enrichment function : Heaviside
    // a chercher element calcul, partition *side of->geo_integ*
    // export the function. 
    // 
    {
    // xEvalConstant<double> one(1.);
    // xexport::xExportGmshAsciiSort pexport;
    xexport::xExportGmshMPIIO pexport;
    Export(*lset, pexport,"LSET");
    // // xEvalApproxFunction<double> evalappro(enrichment_function);
    // // Export(evalappro, pexport, "Heaviside", integ_part, geoModel.getMeshG()->begin(2), geoModel.getMeshG()->end(2)); 

    }
    //d/ create enriched space (but ALL the dofs are enriched...)
    xSpaceXFEM space_full(press, enrichment_function, key_modifier);
    //e/ filter the space in order to keep only the relevant nodes
    xSpaceFiltered enriched(space_full, filter);
    //f/ add the enriched space to the approximation
    press_l.insert(enriched);
       
        cout << "End enrichment" << endl;
#endif
// throw;

    xValueCreator < xSingleValue<valtype> >  creator;
    DeclareInterpolation<xField<valtype> >(press_l, creator, all.begin(), all.end());

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
#ifdef OCTREE
    xClassifyOn classNeg(geoModel.getInString());
    xClassifyOn classPos(geoModel.getOutString());
    // xClassifyOn classNeg("material_neg");
    // xClassifyOn classPos("material_pos");
    //     for(xIter ita = all.begin(1); ita != all.end(1); ++ita){
    //      classPos(*ita);
    //  }
//     cout<<"re classifiy the elements"<<endl;
//     for(auto it=all.begin() ; it!= all.end(); ++it){
//                 mEntity *elem = *it;
//                 xGeomElem geo(elem);
//                 xPartition partition;
//                 mesh->getPartition(elem,partition);


//                 for(xPartition::iterator itp = partition.begin() ; itp!= partition.end(); ++itp){
//                     mEntity *elem_sub = *itp;
//                     xGeomElem geo_sub(elem_sub);
//                     xPoint cdg_sub = geo_sub.getCDGxyz();
//                     geo.setUVWForXYZ(cdg_sub);
//                     xPoint cdg = geo.getUVW();
//                     double val=0.;
// //                    lsHole.getVal(elem,cdg, val);
//                     val = sqrt((cdg_sub(0)-0.)*(cdg_sub(0)-0.) + (cdg_sub(1)-0.25)*(cdg_sub(1)-0.25))-0.12;

//                     if(val>0.) classPos(elem_sub);
//                     else{
//                         classNeg(elem_sub);
//                      classNeg(elem);
//                     }

//                     // cout<<"CLASS\n";
//                 }

//                 // cout<<"----\n";

//             }


#endif

#ifdef HEL
    xIntegrationRulePartition integration_rule_env(xMesh::getPartition, 6*(degree+1));
    xIntegrationRulePartition integration_rule(xMesh::getPartition, 2*(degree+1));    //2*degree is the degree to integrate
    xIntegrationRulePartitionBoundary integration_bc_pos(xMesh::getPartition, xAccept("material_pos"), 4*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc_neg(xMesh::getPartition, xAccept("material_neg"), 4*(degree+1));
    xIntegrationRulePartitionBoundary integration_bc(xMesh::getPartition, 2*(degree+1));
#endif

// {
//     xEvalConstant<double> one(1.);
//     xexport::xExportGmshMPIIO pexport;
//     Export(one, pexport, "POROUS", integration_bc_pos, all.begin(), all.end());
// }


    
    xlinalg::xDistVector<valtype> b(*dist_index);
    xlinalg::xDistVector<valtype> sol(*dist_index);
    xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),1);
    

#ifdef HEL
    AssembleGraph(graphsym,press_l, all.begin(), all.end() );
#endif

#ifdef USE_ZMUMPS
    typedef xlinalg::xLinearSystemSolverMumpsDistributed < valtype > solver_t_;
    typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
    // typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
#endif
    Matrix_SDP_t_ A(graphsym);
    xAssemblerBasic < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler(A, b);
    xAssemblerBasicAndTranspose < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler_transpose(A, b);
    cout<<"Prescribe Neumann BCs"<<endl;
#ifdef OCTREE
    TreatmentOfNatEnv(mesh, press_l, assembler, integration_bc_pos, exact_velocity);
    xMesh* interface = geoModel.getBnd();
    cout<<"read interface boundary"<<endl;
    xEvalNormalFronLevelSet eval_normal_(*geoModel.getLsG());
    cout<<"calculation of interface normal"<<endl;
#else
    TreatmentOfNatEnv(mesh, press_l, assembler, integration_rule_env, exact_velocity);
    xMesh* interface = inclusion->getMesh_bnd();
    if(!lset) throw;
    xEvalNormalFronLevelSet eval_normal_(*lset);
#endif

    
    cout<<"Bilinear form construct "<<endl;
#ifdef HEL
    // first term 1/rho
    xUniformMaterialSensitivity<valtype > fluid_inv_density("inv_density");
    xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
            xEval<valtype>,
            xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > diffusive(fluid_inv_density);  // result type is valtype
     // second term -omega^2/c^c * rho
    xUniformMaterialSensitivity<valtype> inv_celerity_density("inv_celerity_square_density");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xValOperator<xtensor::xCastToComplex<double> > ,valtype > mass(inv_celerity_density);
    

    Assemble(diffusive, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(-omega*omega); 
    Assemble(mass, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);
#endif

    cout<<"end volumic Bilinear form  "<<endl;

#ifndef FEM_STANDARD
    // =============================================interface forms==========================================//
    // Analytical_MatEF mat_P(omega, beta, 0.97, 57e3, 1.54, 73.8e-6, 24.6e-6);
    // Analytical_MatEF mat_P(omega, beta, 0.98, 45e3, 1., 250e-6, 110e-6);
    Analytical_MatEF mat_P(omega, beta, 0.98, 13.5e3, 1.7, 160e-6, 80e-6);  // XFM
    valtype inv_rho_eq = 1./mat_P.rho_eq_til;
    // valtype inv_rho_eq = 1./1.213;
    valtype alpha = sigma_film*d/(j*omega);
    valtype lambda = 1./(alpha+1/alpha_);

    cout<<"pressure jump: "<<std::setprecision(10)<<1./alpha<<endl;
    // lambda * [q][p]-delta p
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > penalty;

if (alpha_ == 0.){
    cout<<"Direct form"<<endl;
    assembler.setCoeff(1./alpha);
    // Assemble(penalty, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    Assemble(penalty, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
}
else {

    xComputeNormalVelocity opNormalvelocity(eval_normal_, fluid_inv_density, classNeg, classPos, 1.);
    


#ifndef GNITSCHE
    cout<<"Classical Nitsche stabilization"<<endl;
    // Penalty terms
    // Nitsche_param_XFEM nitsche(mesh, press_l, *lset, classNeg, classPos);
    // nitsche.ComputeParameter(*interface);
    // xEvalStabParam<valtype> evallambdaC(nitsche, *interface, alpha, *lset);
    
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(gamma, opNormalvelocity);
    
    xEvalConstant<valtype> evallambdaC(lambda);
// [[q]]<p/n>
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > jump_average(opMeanNormalVelocity);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > average_jump(opMeanNormalVelocity);
// <1/rho*q/n><1/rho*p/n>
    xFormBilinearWithoutLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average(opMeanNormalVelocity, opMeanNormalVelocity); 

    xFormBilinearWithLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
        xEval<valtype>,
        xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > lambda_penalty(evallambdaC);

// [[q]]<p/n>
    xFormBilinearWithLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > lambda_jump_average(evallambdaC, opMeanNormalVelocity);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xEval<valtype>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > lambda_average_jump(opMeanNormalVelocity, evallambdaC);

    xFormBilinearWithLaw<xGradMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xEval<valtype>,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > lambda_double_average(opMeanNormalVelocity, evallambdaC, opMeanNormalVelocity); 

#else
    cout<<"gamma Nitsche stabilization"<<endl;
// ================================= Weighting mean =====================================================//
#ifdef OCTREE
    Nitsche_param_XFEM nitsche(mesh, press_l, *geoModel.getLsG(), classNeg, classPos);
    Nitsche_param nitscheJ(mesh, press_l, *geoModel.getLsG(), classNeg, classPos);
    cout<<"calculation of nitsche stability"<<endl;
    nitscheJ.ComputeParameter(*interface, inv_rho_eq);
    cout<<"eval nitsche stability"<<endl;
    xEvalStabParamJiang<valtype> evallambda(nitscheJ, *interface, alpha, *geoModel.getLsG(), inv_rho_eq, alpha_);
#else
    // Nitsche_param_XFEM nitsche(mesh, press_l, *lset, classNeg, classPos);
    // nitsche.ComputeParameter(*interface);
    Nitsche_param nitscheJ(mesh, press_l, *lset, classNeg, classPos);
    nitscheJ.ComputeParameter(*interface, inv_rho_eq);
    xEvalStabParamJiang<valtype> evallambda(nitscheJ, *interface, alpha, *lset, inv_rho_eq, alpha_);
#endif


#if 1
{
    cout<<"export stabilization parameter"<<endl;
    Nitsche_param_XFEM nitsche(mesh, press_l, *lset, classNeg, classPos);
    nitsche.ComputeParameter(*interface);
#ifdef OCTREE
    xEvalGamma<double> evalgamma(nitsche, *interface, alpha, *geoModel.getLsG());
    xEvalBetaJiang<double> evalbeta(nitscheJ, *interface, alpha, *geoModel.getLsG(), inv_rho_eq);
#else
    xEvalGamma<double> evalgamma(nitsche, *interface, alpha, *lset);
    xEvalBetaJiang<double> evalbeta(nitscheJ, *interface, alpha, *lset, inv_rho_eq);
#endif
    xexport::xExportGmshMPIIO pexport;
    xIntegrationRulePartition integ_part_e(xMesh::getPartition, 20);
    Export(evalgamma, pexport, "BETA", integration_rule, all.begin(), all.end());
    Export(evalbeta, pexport, "BETAJ", integration_rule, all.begin(), all.end());
    // Export(evalgamma_real_test, pexport, "STABLE_TEST", integration_rule, all.begin(), all.end());
}   
#endif

    xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocityJ(*interface, nitscheJ, inv_rho_eq, opNormalvelocity);
    cout<<"end stabilization calculation"<<endl;

    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > jump_average(opMeanNormalVelocityJ);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithoutLaw<xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > average_jump(opMeanNormalVelocityJ);
// <1/rho*q/n><1/rho*p/n>
    xFormBilinearWithoutLaw<xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average(opMeanNormalVelocityJ, opMeanNormalVelocityJ); 

    xFormBilinearWithLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
        xEval<valtype>,
        xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > lambda_penalty(evallambda);

    xFormBilinearWithLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > lambda_jump_average(evallambda, opMeanNormalVelocityJ);
//  <1/rho*q/n>[[p]]
    xFormBilinearWithLaw<xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xEval<valtype>,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > lambda_average_jump(opMeanNormalVelocityJ, evallambda);

    xFormBilinearWithLaw<xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>,
            xEval<valtype>,
            xGradWeightMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > lambda_double_average(opMeanNormalVelocityJ, evallambda, opMeanNormalVelocityJ); 
#endif


#ifndef OCTREE
    cout<<"Nitsche form"<<endl;
    assembler.setCoeff(1.);
    Assemble(lambda_penalty, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(1.*alpha);
    Assemble(lambda_jump_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(1.*alpha);
    Assemble(lambda_average_jump, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
       
    assembler.setCoeff(alpha*alpha);
    Assemble(lambda_double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    
    // // [q]*average
    assembler.setCoeff(-1.);
    Assemble(jump_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    Assemble(average_jump, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
// // //  - delta*double_average
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

// +delta * double_average
    assembler.setCoeff(1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
#else
    cout<<"Octree Nitsche form"<<endl;
    assembler.setCoeff(1.);
    Assemble(lambda_penalty, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(1.*alpha);
    Assemble(lambda_jump_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

    assembler.setCoeff(1.*alpha);
    Assemble(lambda_average_jump, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
       
    assembler.setCoeff(alpha*alpha);
    Assemble(lambda_double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    
    // // [q]*average
    assembler.setCoeff(-1.);
    Assemble(jump_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.);
    Assemble(average_jump, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
// // //  - delta*double_average
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

// +delta * double_average
    assembler.setCoeff(1.*alpha);
    Assemble(double_average, assembler, integration_rule_env, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    assembler.setCoeff(1.);

#endif
    cout << "End bilinear Nitsche terms" << endl;

}
#endif
    cout << "End Building the Linear System" << endl;

    if(parsedinfos.getInt("debug") == 1){
       std::ofstream mata_("matA_.mm");  
       A.printMatrixMarket(mata_);

       std::ofstream matb_("matB_.mm");
       A.printMatrixMarket(matb_);

// throw;
       //Export the matrix and vectors for checking the testcase
       string outNameb = "vecb_"+ to_string(proc_id)+".txt";
       std::ofstream ofsvecb(outNameb.c_str());
       b.Printdata(ofsvecb);
       string outNamea = "matA_"+ to_string(proc_id)+".mm";
       std::ofstream ofsmata(outNamea);
       A.printMatrixMarket(ofsmata);
    }


#ifdef USE_ZMUMPS
   solver_t_ solver;

    //#ifdef WITH_DMUMPS
    solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,14,400);

    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::INTER,0,1);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,1,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,2,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,3,6);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,4,2);
    // // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,9,1);
    // solver.setParam(xlinalg::xLinearSystemSolverMumpsBase::ICNTL,11, 1);
    //#endif
    solver.connectMatrix(A);
    b.switchInsertModeOff();
    solver.solve(b, sol);
    sol.switchInsertModeOff();
    sol.switchToGlobalValue();
#endif

#if 0
// interative solver
   typedef Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> SpMat;
   SpMat Acomplex(A.getM(), A.getN());
   typedef Eigen::Triplet<std::complex<double>> SparseTriplet;
   std::vector<SparseTriplet> coefficients(graphsym.getNNZ());
   fillTripletSym(A, coefficients);
   Acomplex.setFromTriplets(coefficients.begin(), coefficients.end());
   Acomplex.makeCompressed();
   //    cout<<"ACOMPLEX\n"<<Acomplex<<endl;

   Eigen::VectorXcd bcomplex(dist_index->getGlobalIndexSize());
   for (unsigned int i = 0; i < bcomplex.rows(); ++i) bcomplex(i) = b[i];

//    Eigen::SparseLU<SpMat> esolver(Acomplex);
   ConjugateGradient<SpMat> cg;
   cg.compute(Acomplex);
   Eigen::VectorXcd solcomplex(dist_index->getGlobalIndexSize());
   solcomplex = cg.solve(bcomplex);
   std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    // update b, and solve again
    solcomplex = cg.solve(bcomplex);
   for (unsigned int i = 0; i < bcomplex.rows(); ++i) sol[i] = solcomplex(i);
   Visit(xWriteSolutionVisitor<xlinalg::xDistVector<valtype>>(sol), complex_manager.begin("dofs"), complex_manager.end("dofs"));
#endif




    Visit(xWriteSolutionVisitor < xlinalg::xDistVector<valtype> >(sol),
          complex_manager.begin("dofs"),
          complex_manager.end("dofs"));
    cout << "End Solving the Linear System" << endl;

cout<<"-----------------------new---------------------------\n";
    cout << "Start Post Treatment" << endl;
    complex_manager.PrintForDebug("res"+pid+".dbg");
    xexport::xExportGmshMPIIO pexport;
    xEvalField < xGetRealPart, xField<valtype> > eval_press_real(press_l);
    xEvalField < xGetImagPart, xField<valtype> > eval_press_imag(press_l);
    xEvalField <xGetComplex, xField<valtype>> eval_press(press_l);
    xEvalField <xGetModu<valtype>, xField<valtype>> eval_press_norm(press_l);

    xEvalField <xGetPressureIndB, xField<valtype>> eval_press_indb(press_l);
    xEvalUnary<xGetModu<double> > pressureindb(eval_press_indb);

    // pexport.setNbSplit(2.*degree);
    Export(eval_press_real, pexport, "PRESSURE_REAL", integration_rule, all.begin(), all.end());
    


#ifdef CAR
    // pexport.setNbSplit(2.*degree);
    // Export(eval_press_indb, pexport, "PRESSURE_REAL", integration_rule, all.begin(), all.end());
    // Export(eval_press_imag, pexport, "PRESSURE_IMAG", integration_rule, all.begin(), all.end());

    auto BB = mesh->compute_bounding_box();
    xtensor::xPoint Pbegin(BB.min(0), 0.6*(BB.min(1)+BB.max(1)), 0.);
    xtensor::xPoint Pend(BB.max(0), 0.6*(BB.min(1)+BB.max(1)), 0.);

    xtensor::xPoint Peval(1.21, 0.6*(BB.min(1)+BB.max(1)), 0.);
    double press_pt = PointExporter(*mesh, *mesh, Peval, pressureindb);
    // cout<<"pressure at point (1.17, 0.024) octree:" <<press_pt<<endl;
    cout<<"pressure at point (1.3, -0.045) octree:" <<press_pt<<endl;
    std::ofstream out_press("press_point.txt",ios::app);
    out_press <<"freq:"<<" "<<freq<<" "<<"pressure_value:"<<" "<<press_pt<<endl;
    lineExporter(*mesh, *mesh, Pbegin, Pend, eval_press_indb, 2000, "PRESS_LINE.txt");
#endif

#if 1
   xEvalGradField <xtool::xIdentity<xtensor::xVector<valtype> > , xField<valtype> > eval_grad_cplx(press_l);

   xEvalVelocity<xGetComplexV> eval_velocity(eval_grad_cplx, air_density, omega);
   xEvalUnary<xGetRealPartV> velocity_real(eval_grad_cplx);
    // Export(velocity_real, pexport, "NORMAL_VELOCITY", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetRealPartV> velocity_exact(exact_velocity);
    // Export(velocity_exact, pexport, "EXACT_VELOCITY", integration_rule, all.begin(), all.end());
    // pexport.setNbSplit(2*degree);
    xEvalUnary<xGetRealPart> exact_solution_real(exact_pressure);
    Export(exact_solution_real, pexport, "EXACT_SOL_REAL", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetModu<valtype>> exact_solution_norm(exact_pressure);
    // Export(exact_solution_norm, pexport, "EXACT_SOL_ABS", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetComplex> exact_solution(exact_pressure);

#endif
#if 1
    cout << "Start Error computation" << endl;
    xIntegrationRulePartition integration_rule_exa(xMesh::getPartition, 4*(degree+1));

    // Binaryopeator for substraction, compute the difference between exact and FEM
    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff(eval_press, exact_pressure);
    // xGetModu class get the absollute value
    xEvalUnary < xGetModu <valtype>> mod_diff(diff);

    xEvalUnary < xGetSquare <double>> square_mod_err(mod_diff);
    xEvalUnary < xGetModu <valtype>> mod_exact_sol(exact_pressure);
    xEvalUnary < xGetSquare <double>> square_exact_sol(mod_exact_sol);
    xEvalBinary<xDivisComplex <double, double, double>> rela_error(square_mod_err, square_exact_sol);

    // Export(mod_diff, pexport, "ABS_ERROR", integration_rule, all.begin(), all.end());

    double int1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1(square_mod_err, int1);  // int one 
    ApplyCommandOnIntegrationRule(integrateInt1, integration_rule_exa, all.begin(), all.end());

    cout<<"absolute error: "<<std::setprecision(10)<<sqrt(int1)<<endl;
    double int2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2(square_exact_sol, int2);  // int one 
     ApplyCommandOnIntegrationRule(integrateInt2, integration_rule_exa, all.begin(), all.end());

    double error = 100.*sqrt(int1/int2);
    // cout<<"exact solution: "<<sqrt(int2)<<endl;
    cout<<"relative error: "<<std::setprecision(10)<<error<<endl;

#if 0
    double int1_l{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1_left(square_mod_err, int1_l);  // int one 
    TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(integrateInt1_left, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt1_left, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));

    double int1_r{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1_right(square_mod_err, int1_r);  // int one 
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateInt1_right, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt1_right, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));

    double int2_l{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2_left(square_exact_sol, int2_l);  // int one 
    TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(integrateInt2_left, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt2_left, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));

    double int2_r{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2_right(square_exact_sol, int2_r);  // int one 
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateInt2_right, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt2_right, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));

    cout<<"pressure error at left interface "<<100.*sqrt(int1_l/int2_l)<<endl;
    cout<<"pressure error at right interface "<<100.*sqrt(int1_r/int2_r)<<endl;
    cout<<"pressure error at right interface "<<0.5*(100.*sqrt(int1_r/int2_r)+100.*sqrt(int1_l/int2_l))<<endl;

#endif
//  ================================ Jump error and interface flux error==========================================//
    xEvalJump< valtype > eval_press_jump(eval_press, &press_l);

    xEvalJump< valtype > eval_exact_jump(exact_pressure);

    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_jump(eval_press_jump, eval_exact_jump);
    xEvalUnary < xGetModu <valtype>> mod_diff_jump(diff_jump);
    xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    xEvalUnary <xGetSquare<double>> sqart_err_jump(mod_diff_jump);
    xEvalUnary <xGetSquare<double>> sqart_exa_jump(mod_exact_jump);
    // xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    xEvalUnary < xGetModu <valtype>> mod_press_jump(eval_press_jump);

    xEvalNormal eval_normal_real;
    xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal_inter(eval_normal_);
    xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > normal_exact_flux(exact_velocity, eval_normal_inter);
    xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > normal_flux(eval_grad_cplx, eval_normal_inter);
    xEvalUnary < xGetModu <valtype>> mod_exa_flux(normal_exact_flux);
    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_flux(normal_flux, normal_exact_flux);
    xEvalUnary < xGetModu <valtype>> mod_diff_flux(diff_flux);
    xEvalUnary < xGetSquare <double>> square_exact_flux(mod_exa_flux);
    xEvalUnary < xGetSquare <double>> square_diff_flux(mod_diff_flux);

    double int_fluxp1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxp1(square_diff_flux, int_fluxp1);
    TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(integrateInt_fluxp1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt_fluxp1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    double int_fluxp2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxp2(square_exact_flux, int_fluxp2);
    TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(integrateInt_fluxp2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt_fluxp2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    // 
    double int_fluxn1{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxn1(square_diff_flux, int_fluxn1);
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateInt_fluxn1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt_fluxn1, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    double int_fluxn2{0.};
    xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt_fluxn2(square_exact_flux, int_fluxn2);
    TagSideNegative(interface);
    // ApplyCommandOnIntegrationRule(integrateInt_fluxn2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    ApplyCommandOnIntegrationRule(integrateInt_fluxn2, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreatorOnMeshC(mesh));
    double error_flux_pos = 100.*sqrt(int_fluxp1/int_fluxp2);
    double error_flux_neg = 100.*sqrt(int_fluxn1/int_fluxn2);
    cout<<"relative error flux pos side: "<<error_flux_pos<<endl;
    cout<<"relative error  flux neg side: "<<error_flux_neg<<endl;
    // cout<<"flux neg side: "<<int_fluxn2<<endl;
    cout<<"relative average flux error: "<<(error_flux_pos+error_flux_neg)/2.<<endl;
    // std::ofstream out_press("error_press.txt",ios::app);
    // out_press << degree<<" "<<complex_manager.size("dofs")<<" "<<"error pressure "<<error<<endl;

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
    // xEvalBinary< xDivisComplex<double, double, double>> err_rela_jump(sqart_err_jump, sqart_exa_jump);
    // Export(sqart_exa_jump, pexport, "ERROR_RELA_JUMP", integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());
    double errorJump = sqrt(intJump1/intJump2)*100.;
    double abserrorJump = sqrt(intJump1);
    double Jump = sqrt(intJump2)*100.;
    cout << "jump error : " <<errorJump<< endl;
    cout << "abs jump error : " <<abserrorJump<< endl;
#endif

    cout << "Overval" << endl;
    cout << "End Post Treatment" << endl;


    cout<<"--------------------------------------------------------\n";
// #ifndef OCTREE
// }
// #endif

    return;
}



template <class FIELD , class ITER>
void helmholtzNdPSG :: TreatmentOfEssEnv(const FIELD& fct, ITER itenv, ITER endenv,
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
                std::complex<double> val {env.getValue(), 0.};
                DirichletBoundaryCondition (fct, phys, bc.begin(), bc.end(), val);
            }
        }
    }
    return;
}



template <class ASSEMBLER, class FIELD >
void  helmholtzNdPSG :: TreatmentOfMixedBCs(xData &data, ASSEMBLER &assembler, const xIntegrationRule& integration_rule,
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
void helmholtzNdPSG :: TreatmentOfNatEnv   (xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, xEval<xtensor::xVector<std::complex<double> > > &exact_velocity)
{

    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it)
    {
        const xEnv & env = *it;
        if (env.Phys == "ELECTRIC_DISPLACEMENT")
        {
            assert(env.Type == FIX);
            std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);
            // cout<<"assemble exact"<<endl;
            xUniformMaterialSensitivity<valtype> inv_density("inv_density");
            xEvalInvDensity<valtype> eval_inv_density(inv_density);
            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, flux);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

            // xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
            //                       xEvalBinary<xMult<xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double>> >, valtype > lin(flux);
// //          
            xClassRegion bc(mesh, env.Entity, env.getDimension());

// {   xexport::xExportGmshAsciiSort pexport;
//     xEvalConstant<double> one(1.);
//     Export(one, pexport, "EDGE"+env.Entity, integration_rule, bc.begin(), bc.end());
// }
            // assembler.setCoeff(1./1.213);
            assembler.setCoeff(1.);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.);
        }
        else if(env.Phys == "NORMAL_VELOCITY"){

            assert(env.Type == FIX);
            double val = env.getValue();
            xtensor::xVector<valtype>  vec_val({val*cos(0.)}, {val*sin(0.)}, {0.,0.});
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
            xClassRegion bc(mesh, env.Entity, env.getDimension());

            // assembler.setCoeff(minusjomegaexp);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), xUpperAdjacency());
            // assembler.setCoeff(1.);

        }
    }

    return;
}

template < typename ASSEMBLER, typename FIELD >
void helmholtzNdPSG :: TreatmentOfNatEnv   (xfem::xMesh* mesh, const FIELD& listFunctionSpace,
                                         ASSEMBLER & assembler, xIntegrationRule & integration_rule, xEval<xtensor::xVector<std::complex<double> > > &exact_velocity, double omega_sq)
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
            // xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(inv_density, flux);

            // xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
            //                       xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double>> >, valtype > lin(flux);

            xClassRegion bc(mesh, env.Entity, env.getDimension());

            assembler.setCoeff(1./(1.213*omega_sq));
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(),
                     xUpperAdjacency());
            assembler.setCoeff(1.);

        }
    }

    return;
}
