#include "helmholtz.h"

// Export
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"


//xfem::distmesh
#include "xLoadBalanceTools.h"

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

// From helmholtzNd
#include "xComplexUtils.h"
#include "exactSolutions.h"
#include "acousticMaterial.h"
#include "heParser.h"
#include "nitscheOperators.h"

using namespace std;
using namespace xfem;
using namespace xtool;
using namespace xtensor;
using namespace AOMD;
using namespace xlinalg;
using namespace xcut;

// #define CHECK_WITH_EIGEN 1
// #define XFEM_CONTINUE
// #define XFEM_DISCONTINUE
// #define USE_NITSCHE_DISCONTINUE
// #define USE_PENALTY_DISCONTIUE
// #define USE_NITSCHE_CONTINUE
// #define FEM_PENALTY


#ifdef CHECK_WITH_EIGEN
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "xGenericSparseMatrixTraitPolicy.h"
#endif


void lineExporter(xMesh &m, Trellis_Util::mPoint p1, Trellis_Util::mPoint p2, xEval<double> &eval, int nbPoints, string filename){

  xtensor::xVector<> v12(p1,p2);
  double length = v12.normValue();
  double step = length / (nbPoints-1);


  //Evaluation points
  std::vector<Trellis_Util::mPoint> points;
  points.reserve(nbPoints);
  //   std::vector<double> curvAbs(nbPoints);
  points.push_back(p1);
  //   curvAbs.push_back(0);

  for(int iPt=1; iPt < nbPoints; ++iPt){
      points.push_back(Trellis_Util::mPoint(p1(0) + iPt * step * v12(0) , p1(1) + iPt * step * v12(1), p1(2) + iPt * step * v12(2)));
      //     curvAbs.push_back(iPt * step / length);
    }

  ofstream file(filename, ofstream::out);

  double curvAbs{0.};
  for(Trellis_Util::mPoint pt: points){
      std::set<mEntity *> elts;
      m.locateElement(pt,elts);


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


helmholtzNd::helmholtzNd (xfem::xData &data_, MPI_Comm world, xParseData &parsedinfos) : data(data_) {

    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);
    // read materials
    //Give the pulsation to the porous material
    xfem::xAcousticEquivalentFluid::setOmega(parsedinfos.getDouble("H_angular_frequency"));
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

    // set  partition with parmetis
    xinterface::aomd::xAttachedDataManagerAOMD < int > target_proc = xmeshtool::setParmetisPartition(*( data.mesh ), (data.mesh)->getPartitionManager(),  nb_proc  );
    xmeshtool::migrateEntities(*( data.mesh ),(data.mesh)->getPartitionManager(), target_proc);

    xmeshtool::exportGMSH_V2_Dist(*( data.mesh ),( data.mesh )->getPartitionManager(), "partmesh");


    //Read zones
    data.ReadZones();


}


void helmholtzNd::TreatmentOfFormulation (xParseData &parsedinfos) {


    cout << "Starting Formulation for Complex2d" << endl;
    cout << "Start creating the function space" << endl;
    xRegion all(data.mesh);
    //xSpaceLagrange lagrange("PRESSURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_TWO);

//    const int degree = 1;
    const int degree = parsedinfos.getInt("order_approx_min");
    xSpacePolynomialBernstein approxSpace("PRESSURE", xSpace::SCALAR, degree);

    cout<<"Approximation order is "<<degree<<endl;


    double omega = parsedinfos.getDouble("H_angular_frequency");
    double beta = parsedinfos.getDouble("H_wave_angle");  // degree
    double wave_angle = beta*M_PI/180.;  // radians
    double sigma_film  = parsedinfos.getDouble("Resistivity");
    double d = parsedinfos.getDouble("Mat_Thickness");
    double alpha_ = parsedinfos.getDouble("Nitsche");
    double h = parsedinfos.getDouble("Element_size");
    

    double air_density = 1.213;
    double sound_celerity =340;
    //double k = omega / sound_celerity;
    valtype v0 = exp(j);
//    double wave_angle = 0.;//radians

    // xEvalObliqueExactSolutionBouillard2DCplx exact_pressure(omega, wave_angle, d, v0);
    xEvalObliqueExactSolutionBouillard2DCplx exact_pressure(omega, wave_angle, sigma_film, d);
    xEvalObliqueSpeedBouillard2DCplx exact_velocity(omega, wave_angle, sigma_film, d);


    xfem::xValueManagerDist<valtype> complex_manager;
    xField<valtype> press_l(&complex_manager, approxSpace);

    //If a levelset was provided

        // heaviside enrichment
#ifdef XFEM_CONTINUE
    xClassifyOn classNeg("material_neg");
    xClassifyOn classPos("material_pos");
        if(lset){
            //--- Create physical surface object
            xcut::xPhysSurfParameter param_surf(classNeg, classPos);
            inclusion.reset(new xcut::xPhysSurfByTagging(*lset, param_surf));

            //--- Setup enriched space
            //a/ keyword to create the enriched keys
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
    xClassifyOn classNeg("material_neg");
    xClassifyOn classPos("material_pos");

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
    std::vector < string > vars_ess {"PRESSURE"};
    TreatmentOfEssEnv(press_l, data.PhysicalEnv->begin(), data.PhysicalEnv->end(),
                      vars_ess, data.mesh);


    complex_manager.PrintForDebug("cplx_dcl_"+pid+".dbg");

    // In case some one use this test with Dirichlet BC we need to uniformize fixed values
    finalizeFixedUniform(data.mesh->getPartitionManager(),complex_manager);
    complex_manager.PrintForDebug("cplx_ess_"+pid+".dbg");

    xDistStateDofCreatorUniform<xMesh::partman_t, xTrue, valtype> snh(data.mesh->getPartitionManager(),complex_manager, "dofs");
    DeclareState<xField<valtype> >(press_l, snh, all.begin(), all.end());


    std::shared_ptr < xlinalg::xDistIndex > dist_index = finalizeDofUniform(snh);
    cout << "End creating the function space" << endl;



    cout << "Start Building the Linear System" << endl;
    complex_manager.PrintForDebug("cplx_sym_"+pid+".dbg");
    xIntegrationRulePartition integration_rule_env(2*(degree+1));
    xIntegrationRulePartition integration_rule(2*(degree+1));    //2*degree is the degree to integrate
    xlinalg::xDistVector<valtype> b(*dist_index);
    xlinalg::xDistVector<valtype> sol(*dist_index);
    xlinalg::xGraphMatrix graphsym(dist_index->getGlobalIndexSize(),dist_index->getGlobalIndexSize(),1);
    AssembleGraph(graphsym,press_l, all.begin(), all.end() );

#ifdef USE_ZMUMPS
    typedef xlinalg::xLinearSystemSolverMumpsDistributed < valtype > solver_t_;
    typedef xlinalg::xGenericSparseMatrix < valtype,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,solver_t_::matrix_storage,solver_t_::matrix_indexing > Matrix_SDP_t_;
#endif
    Matrix_SDP_t_ A(graphsym);
    xAssemblerBasic < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler(A, b);
    xAssemblerBasicAndTranspose < Matrix_SDP_t_,xlinalg::xDistVector<valtype>, valtype > assembler_transpose(A, b);


    TreatmentOfNatEnv(press_l, assembler, integration_rule_env, exact_velocity, 1, -j*omega*exp(j), wave_angle);
    cout << "End treating natural BCs" << endl;


    //    std::ofstream vecb_("vecb_.txt");//Export per-process vectors...
    //    b.Printdata(vecb_);

    // first term 1/rho
    xUniformMaterialSensitivity<valtype > fluid_inv_density("inv_density");
    xFormBilinearWithLaw<xGradOperator< xtensor::xCastToVectorComplex<double> >, // xtensor convert to vector complex double type
            xEval<valtype>,
            xGradOperator<xtensor::xCastToVectorComplex<double> >, valtype > diffusive(fluid_inv_density);  // result type is valtype


    Assemble(diffusive, assembler, integration_rule, press_l, press_l, all.begin(), all.end());

    // second term -omega^2/c^c * rho
    xUniformMaterialSensitivity<valtype> inv_celerity_density("inv_celerity_square_density");
    xFormBilinearWithLaw<xValOperator<xtensor::xCastToComplex<double> >,
            xEval<valtype>,
            xValOperator<xtensor::xCastToComplex<double> > ,valtype > mass(inv_celerity_density);


    assembler.setCoeff(-omega*omega);
    Assemble(mass, assembler, integration_rule, press_l, press_l, all.begin(), all.end());
    assembler.setCoeff(1.);

#ifdef USE_NITSCHE_DISCONTINUE
// interface forms

    xMesh* interface = inclusion->getMesh_bnd();
    xEvalNormalFronLevelSet eval_normal_(*(lset.get()));
    xComputeNormalVelocity opNormalvelocity(eval_normal_, fluid_inv_density, classNeg, classPos);
    xGradMeanOperatorWithParamInter<xComputeNormalVelocity> opMeanNormalVelocity(opNormalvelocity);

// interface coefficient
    valtype alpha = sigma*d/(j*omega);
    valtype lambda = 1./(alpha+h/alpha_);

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
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > double_average(opMeanNormalVelocity); 

    assembler.setCoeff(lambda);
    Assemble(penalty, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(-1.*lambda*alpha);
    Assemble(jump_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    assembler.setCoeff(-1.*lambda*alpha);
    Assemble(average_jump, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
       
    assembler.setCoeff(lambda*alpha*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    
    // [q]*average
    Assemble(jump_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    Assemble(average_jump, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
//  - delta*double_average
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);
    assembler.setCoeff(-1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

// +delta * double_average
    assembler.setCoeff(1.*alpha);
    Assemble(double_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    cout << "End bilinear Nitsche terms" << endl;
#endif

#ifdef FEM_PENALTY
    valtype alpha = j*omega / (sigma * d);
    xFormBilinearWithoutLaw<xValJumpOperator<xtensor::xCastToComplex<double> >,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > penalty;

    xClassRegion interface(data.mesh, 5, 2);
    assembler.setCoeff(alpha);
    Assemble(penalty, assembler, integration_rule, press_l, press_l, interface.begin(), interface.end(),  xUpperCreator());
    assembler.setCoeff(1.);
    cout << "End bilinear FEM_Penalty terms" << endl;
#endif

#ifdef CHARAC_WAVE_METHOD
    xFormBilinearWithoutLaw<xValJumpOperatorwithInpedance<xtensor::xCastToComplex<double> >,
            xValJumpOperator<xtensor::xCastToComplex<double> >, valtype > impedance_jump;
    assembler.setCoeff(-1.*j*omega);
    Assemble(impedance_jump, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    xFormBilinearWithoutLaw<xValJumpOperatorwithInpedance<xtensor::xCastToComplex<double> >,
            xGradMeanOperatorWithParamInter<xComputeNormalVelocity>, valtype > impedance_average(opMeanNormalVelocity);

    
    assembler.setCoeff(-1.*sigma*d);
    Assemble(impedance_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());
    assembler.setCoeff(1.);

    Assemble(jump_average, assembler, integration_rule, press_l, press_l, interface->begin(1), interface->end(1), xUpperCreator());


#endif

    // Treatment of Robin BCs
    TreatmentOfMixedBCs(data, assembler, integration_rule_env, press_l, omega);
    cout << "End Building the Linear System" << endl;
    //    std::ofstream mata_("matA_.mm");  
    //    A.printMatrixMarket(mata_);

    //    std::ofstream matb_("matB_.mm");
    //    A.printMatrixMarket(matb_);


    //    //Export the matrix and vectors for checking the testcase
    //    string outNameb = "vecb_"+ to_string(proc_id)+".txt";
    //    std::ofstream ofsvecb(outNameb.c_str());
    //    b.Printdata(ofsvecb);
    //    string outNamea = "matA_"+ to_string(proc_id)+".mm";
    //    std::ofstream ofsmata(outNamea);
    //    A.printMatrixMarket(ofsmata);

#ifdef CHECK_WITH_EIGEN
    typedef Eigen::SparseMatrix<valtype, Eigen::ColMajor > SpMat;
    SpMat Acomplex(A.getM(),A.getN());
    typedef Eigen::Triplet<valtype > SparseTriplet;
    std::vector<SparseTriplet> coefficients(graphsym.getNNZ());
    fillTripletSym(A, coefficients);
    Acomplex.setFromTriplets(coefficients.begin(), coefficients.end());
    Acomplex.makeCompressed();
    //    cout<<"ACOMPLEX\n"<<Acomplex<<endl;

    b.switchInsertModeOff();
    b.switchInsertModeOn();
    Eigen::VectorXcd bcomplex(dist_index->getGlobalIndexSize());
    for(unsigned int i = 0; i < bcomplex.rows(); ++i) bcomplex(i) = b[i];
#endif

#ifdef USE_ZMUMPS
   solver_t_ solver;

    //#ifdef WITH_DMUMPS
    //    solver_.setParam(xlinalg::xLinearSystemSolverMumpsBase::INTER,0,1);
    //#endif
    solver.connectMatrix(A);
    b.switchInsertModeOff();
    solver.solve(b, sol);
    sol.switchInsertModeOff();
    sol.switchToGlobalValue();
#endif
    //    cout<<"------sol------\n";
    //    sol.Printdata(cout);

#ifdef CHECK_WITH_EIGEN
    Eigen::SparseLU<SpMat> esolver(Acomplex);
    auto solcomplex = esolver.solve(bcomplex);
    sol.switchInsertModeOn();
    for(unsigned int i = 0; i < bcomplex.rows(); ++i) sol[i] = solcomplex(i);
    //cout<<"EigenSol:\n";
    //sol.Printdata(cout);
#endif


    Visit(xWriteSolutionVisitor < xlinalg::xDistVector<valtype> >(sol),
          complex_manager.begin("dofs"),
          complex_manager.end("dofs"));
    cout << "End Solving the Linear System" << endl;


    cout << "Start Post Treatment" << endl;
    complex_manager.PrintForDebug("res"+pid+".dbg");
    xexport::xExportGmshMPIIO pexport;
    xEvalField < xGetRealPart, xField<valtype> > eval_press_real(press_l);
    Export(eval_press_real, pexport, "PRESSURE_REAL", integration_rule, all.begin(), all.end());
    xEvalField < xGetImagPart, xField<valtype> > eval_press_imag(press_l);
//    Export(eval_press_imag, pexport, "PRESSURE_IMAG", integration_rule, all.begin(), all.end());
    xEvalField <xGetModu<valtype>, xField<valtype>> eval_press_norm(press_l);
    Export(eval_press_norm, pexport, "PRESSURE_MODULE", integration_rule, all.begin(), all.end());
    xEvalField <xGetPhase<valtype>, xField<valtype>> eval_press_phase(press_l);
    Export(eval_press_phase, pexport, "PRESSURE_ARG", integration_rule, all.begin(), all.end());

    
    
   xEvalGradField <xtool::xIdentity<xtensor::xVector<valtype> > , xField<valtype> > eval_grad_cplx(press_l);
   xEvalVelocity<xGetComplexV> eval_velocity(eval_grad_cplx, air_density, omega);
   xEvalUnary<xGetRealPart> velocity_real(eval_velocity);
    // Export(velocity_real, pexport, "NORMAL_VELOCITY", integration_rule, all.begin(), all.end());
//    xEvalVelocity<xGetImagPartV> eval_velocity_imag(eval_grad_cplx, air_density, omega);
//    Export(eval_velocity_real, pexport, "VELOCITY_REAL", integration_rule, all.begin(), all.end());
//    Export(eval_velocity_imag, pexport, "VELOCITY_IMAG", integration_rule, all.begin(), all.end());

    // xEvalUnary<xGetRealPart> exact_solution_real(exact_pressure);
//     Export(exact_solution_real, pexport, "EXACT_SOL_REAL", integration_rule, all.begin(), all.end());
//     xEvalUnary<xGetImagPart> exact_solution_imag(exact_pressure);
// //    Export(exact_solution_imag, pexport, "EXACT_SOL_IMAG", integration_rule, all.begin(), all.end());
    xEvalUnary<xGetModu<valtype>> exact_solution_norm(exact_pressure);
    Export(exact_solution_norm, pexport, "EXACT_SOL_MODULE", integration_rule, all.begin(), all.end());
//     xEvalUnary<xGetPhase<valtype>> exact_solution_phase(exact_pressure);
//     Export(exact_solution_phase, pexport, "EXACT_SOL_ARG", integration_rule, all.begin(), all.end());

    
    // Export(exact_solution_real, pexport, "EXACT_SOL_BND_POS_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(eval_press_real, pexport, "PRESSURE_BND_POS_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());

    // TagSideNegative(interface);
    // Export(exact_solution_real, pexport, "EXACT_SOL_BND_NEG_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(eval_press_real, pexport, "PRESSURE_BND_NEG_REAL", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // UnTagSide(interface);

    // herited xtool::identity to get complex value xGetComplex
    xEvalField <xGetComplex, xField<valtype>> eval_press(press_l);
    xEvalUnary<xGetComplex> exact_solution(exact_pressure);

    xEvalBinary <xDivisComplex <valtype, valtype, valtype>> impedence(eval_press, eval_velocity);
    // xEvalUnary<xGetRealPart> impedence_real(impedence);
    // Export(impedence_real, pexport, "IMPEDENCE", integration_rule, all.begin(), all.end());
    xEvalUnary<xReflection> ref(impedence);
    xEvalUnary<xGetModu<valtype>> ref_abs(ref);
    xEvalUnary<xGetPhase<valtype>> ref_angle(ref);
    Export(ref_abs, pexport, "REFLECTION_REAL", integration_rule, all.begin(), all.end());
    Export(ref_angle, pexport, "REFLECTION_ANGLE", integration_rule, all.begin(), all.end());
    // xEvalUnary<xAbsorption<valtype>> absorpt(ref);
    // Export(absorpt, pexport, "ABSORPTION", integration_rule, all.begin(), all.end());

    cout << "-----test eval jump------" << endl;

<<<<<<< HEAD
    // xEvalJump< valtype > eval_press_jump(eval_press, &press_l);
=======
    xEvalJump< valtype > eval_press_jump(eval_press, &press_l);
>>>>>>> ab8f4cb1a1d015e2f35ee7fa6ce640397484d33c
    // xEvalJump< double > eval_press_jump(eval_press_real);
    // xEvalJump< valtype > eval_exact_jump(exact_solution);
    // xEvalJump< double > eval_exact_jump(exact_solution_real);

<<<<<<< HEAD
    // xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff_jump(eval_press_jump, eval_exact_jump);
    // xEvalUnary < xGetModu <valtype>> mod_diff_jump(diff_jump);
    // xEvalUnary < xGetModu <valtype>> mod_exact_jump(eval_exact_jump);
    // xEvalUnary <xGetSqrt<double>> sqart_err_jump(mod_diff_jump);
    // xEvalUnary <xGetSqrt<double>> sqart_exa_jump(mod_exact_jump);
    
    // Export(mod_diff_jump, pexport, "EVAL_PRESS_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // Export(mod_diff_jump, pexport, "EVAL_PRESS_JUMP_TEST", integration_rule, all.begin(), all.end());

    // double int1{0.};
    // xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt1(sqart_err_jump, int1);  // int one 
    // TagSidePositive(interface);
    // ApplyCommandOnIntegrationRule(integrateInt1, integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // double int2{0.};
    // xIntegrateEvalCommand<xEvalUnary < xGetSquare <double>> > integrateInt2(sqart_exa_jump, int2);
    // ApplyCommandOnIntegrationRule(integrateInt2, integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // xEvalBinary< xDivisComplex<double, double, double>> err_rela_jump(sqart_err_jump, sqart_exa_jump);
    // Export(err_rela_jump, pexport, "ERROR_RELA_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    // double error = int1/int2;
    // cout << "jump error : " <<error<< endl;
    
    
    
    // Export(mod_errJ, pexport, "MODULE_ERROR_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
=======
    xEvalUnary < xGetModu <valtype>> eval_press_jump_(eval_press_jump);
    xEvalUnary < xGetModu <valtype>> eval_exact_jump_(eval_exact_jump);
     Export(eval_press_jump_, pexport, "EVAL_PRESS_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
     Export(eval_exact_jump_, pexport, "EVAL_EXACT_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
//    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> errJ(eval_press_jump, eval_exact_jump);
//    xEvalUnary < xGetModu <valtype>> mod_errJ(errJ);
    
//    Export(mod_errJ, pexport, "MODULE_ERROR_JUMP", integration_rule, interface->begin(1), interface->end(1), xUpperCreator());
    throw;
>>>>>>> ab8f4cb1a1d015e2f35ee7fa6ce640397484d33c
    cout << "-----End eval jump-------" << endl;


    /* ****************************error computation*******************************
    created several operator class in xComplexUtils.h for complex and for errr computation, facilite to manager*****
    **********************************************************/
    //start
#if 0
    cout << "Start Error computaiton" << endl;

    xValueManagerDist<double> double_manager;
    xField<> err_info_l(&double_manager, xSpaceConstant("ERROR_INFO"));
    xValueCreator<xValueError>  creator_err;
    DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());
    xIntegrationRulePartition integration_rule_exa(2*(degree+1));

    // Binaryopeator for substraction, compute the difference between exact and FEM
    xEvalBinary <xMinusComplex <valtype, valtype, valtype>> diff(eval_press, exact_solution);
    // xGetModu class get the absollute value
    xEvalUnary < xGetModu <valtype>> mod_diff(diff);
    // xGetSquare struc : get Square
    xEvalUnary < xGetSquare <double>> square_mod_err(mod_diff);
    xEvalUnary < xGetModu <valtype>> mod_exact_sol(exact_solution);
    xEvalUnary < xGetSquare <double>> square_exact_sol(mod_exact_sol);

    // integration
    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
            xEvalUnary < xGetModu <valtype>> > abs_exa_form(mod_exact_sol);
    xValueError::choice("ENG_EXA");
    xFillFieldFromZeroForm<> fill_abs_exa(err_info_l, abs_exa_form);
    ApplyCommandOnIntegrationRule(fill_abs_exa, integration_rule_exa, all.begin(), all.end());

    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
            xEvalUnary < xGetModu <valtype>> > abs_err_form(mod_diff);
    xValueError::choice("ABS2_EXA");
    xFillFieldFromZeroForm<> fill_err_exa(err_info_l, abs_err_form);
    ApplyCommandOnIntegrationRule(fill_err_exa, integration_rule_exa, all.begin(), all.end());

    xValueError::finalize(data.mesh->getPartitionManager().getComm());
    xValueError::choice("REL_EXA");

    xEvalField<xtool::xIdentity<double> > val_err_info(err_info_l);
    Export(val_err_info, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());


    std::ofstream out;
    if (!proc_id)
    {
      out.open("error.txt",std::ofstream::out);
      out << "energy error      " << endl;
      out << "solution exact is   " << xValueError::total("ENG_EXA") << endl;
      out << "err abs exact is  " << xValueError::total("ABS_EXA") << endl;
      out << "err rel exact is  " << xValueError::total("REL_EXA") << endl;
    }
#endif

#if 0
    xValueError::clear();
    
    xField<> err_info_inter(&double_manager, xSpaceConstant("ERROR_INFO_INTER"));
    
    DeclareInterpolation(err_info_inter, creator_err, interface->begin(1), interface->end(1));
    
    xEvalField<xtool::xIdentity<double> > val_err_interf_info(err_info_inter);

    // integration
    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
            xEvalUnary < xGetModu <valtype>> > abs_exa_jump_form(mod_exact_jump);
    xValueError::choice("ENG_EXA");
    xFillFieldFromZeroForm<> fill_abs_exa_jump(err_info_inter, abs_exa_jump_form);
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_abs_exa_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    xFormZeroEvalBilinearFormWithoutLaw<xEvalUnary < xGetModu <valtype>>,
            xEvalUnary < xGetModu <valtype>> > abs_err_jump_form(mod_diff_jump);
    xValueError::choice("ABS2_EXA");
    xFillFieldFromZeroForm<> fill_err_jump(err_info_inter, abs_err_jump_form);
    // TagSidePositive(interface);
    ApplyCommandOnIntegrationRule(fill_err_jump, integration_rule_exa, interface->begin(1), interface->end(1), xUpperCreator());

    xValueError::finalize(data.mesh->getPartitionManager().getComm());
    xValueError::choice("REL_EXA");

    // xEvalField<xtool::xIdentity<double> > val_err_jump_info(err_info_inter);
    // Export(val_err_info, pexport, "ERROR_RELATIVE", integration_rule, all.begin(), all.end());

    if (!proc_id)
    {
      
      out << "energy error  of interface   " << endl;
      out << "solution exact is   " << xValueError::total("ENG_EXA") << endl;
      out << "err abs exact is  " << xValueError::total("ABS_EXA") << endl;
      out << "err rel exact is  " << xValueError::total("REL_EXA") << endl;
    }


    if (!proc_id)
      out.close();
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
    Trellis_Util::mPoint BBmin, BBmax;
    data.mesh->compute_bounding_box(BBmin, BBmax);
    Trellis_Util::mPoint Plbegin(0.5*(BBmin(0)+BBmax(0))-0.002, BBmin(1), 0.);
    Trellis_Util::mPoint Plend(0.5*(BBmin(0)+BBmax(0))-0.002, BBmax(1), 0.);
    Trellis_Util::mPoint Prbegin(0.5*(BBmin(0)+BBmax(0))+0.002, BBmin(1), 0.);
    Trellis_Util::mPoint Prend(0.5*(BBmin(0)+BBmax(0))+0.002, BBmax(1), 0.);
    Trellis_Util::mPoint Pmbegin(0.5*(BBmin(0)+BBmax(0)), BBmin(1), 0.);
    Trellis_Util::mPoint Pmend(0.5*(BBmin(0)+BBmax(0)), BBmax(1), 0.);
    Trellis_Util::mPoint Pbegin(BBmin(0), 0.5*(BBmin(1)+BBmax(1)), 0.);
    Trellis_Util::mPoint Pend(BBmax(0), 0.5*(BBmin(1)+BBmax(1)), 0.);

    lineExporter(*data.mesh, Pbegin, Pend, eval_press_real, 5000, "PRESSURE_REAL_LINE.txt");
    lineExporter(*data.mesh, Pbegin, Pend, eval_press_norm, 5000, "PRESSURE_NORM_LINE.txt");
    // lineExporter(*data.mesh, Pbegin, Pend, exact_solution_norm, 5000, "PRESSURE_EXACT_LINE.txt");
    // lineExporter(*data.mesh, Prbegin, Prend, eval_press_real, 1000, "PRESSURE_RIGHT_LINE.txt");
    // lineExporter(*data.mesh, Prbegin, Prend, exact_solution_real, 1000, "EXAC_RIGHT_LINE.txt");
    // lineExporter(*data.mesh, Pmbegin, Pmend, sqart_err_rela, 1000, "ERROR_LINE.txt");
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
                std::complex<double> val {env.getValue(), 0.};
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
            assembler.setCoeff(env.getValue() * omega * j);

            xClassRegion gammaR(data.mesh, env.Entity, env.getDimension());
            Assemble(mass, assembler, integration_rule, press_l, press_l, gammaR.begin(), gammaR.end(), xUpperAdjacency());
            assembler.setCoeff(1.);
        }
    }

}



//template <class ASSEMBLER, class FIELD, class BILINEARFORM >
//void  helmholtzNd :: TreatmentOfMixedBCs(xData &data, ASSEMBLER &assembler, const xIntegrationRule& integration_rule,
//                                       const FIELD& press_l, BILINEARFORM &mass, double k){


//    for (xPhysicalEnv::const_iterator it = data.PhysicalEnv->begin(); it != data.PhysicalEnv->end(); ++it) {
//        const xEnv& env = *it;
//        if (env.Phys == "ELASTIC_SUPPORT") {
//            assert(env.Type == FIX);
//            cout<<"Treat Robin BCs for physical entity "<<env.Entity<< " of dim " << env.getDimension()<<endl;

//
//            assembler.setCoeff(env.getValue() * omega * j);

//            xClassRegion gammaR(data.mesh, env.Entity, env.getDimension());
//            Assemble(mass, assembler, integration_rule, press_l, press_l, gammaR.begin(), gammaR.end(), xUpperAdjacency());
//        }
//    }

//    assembler.setCoeff(1.);
//}



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
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);
            xEvalBinary < xtool::xMult < xtensor::xVector<std::complex<double> >, xtensor::xVector<std::complex<double> >, std::complex<double> > > flux(exact_velocity, eval_normal);

            xUniformMaterialSensitivity<valtype> inv_density("inv_density");
            xEvalInvDensity<valtype> ev_inv_density(inv_density);
            xEvalBinary<xMult<valtype, valtype, valtype> > scaled_flux(ev_inv_density, flux);

            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary<xMult<valtype, valtype, valtype> >, valtype > lin(scaled_flux);

            xClassRegion bc(data.mesh, env.Entity, env.getDimension());

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
            std::cout<<"In TreatmentOfNatEnv, treating Bc "<<env.Entity<<" of dim "<<env.getDimension()<<" val: "<<val<<std::endl;
            xEvalNormal eval_normal_real;
            xEvalUnary<xtensor::xCastToVectorComplex<double> > eval_normal(eval_normal_real);

            xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > > scaled_flux(flux, eval_normal);
            xFormLinearWithLoad < xValOperator < xtool::xIdentity < double > >,
                                  xEvalBinary < xtool::xMult < xtensor::xVector<valtype>, xtensor::xVector<std::complex<double> >, std::complex<double> > >, std::complex<double> > lin(scaled_flux);
//            xEvalUnary<xtool::xIdentity<double> > scaled_flux(flux);
//            xFormLinearWithLoad<xValOperator<xtool::xIdentity<double> >,
//                      xEvalUnary<xtool::xIdentity<double> >, std::complex<double>> lin(scaled_flux); //Shape functions are assumed double
            xClassRegion bc(data.mesh, env.Entity, env.getDimension());

            assembler.setCoeff(minusjomegaexp);
            Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), xUpperAdjacency());
            assembler.setCoeff(1.);




        }
    }


    return;
}
