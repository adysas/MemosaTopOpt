// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWMODEL_H_
#define _FLOWMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FlowFields.h"

#include "Mesh.h"
#include "LinearSolver.h"

#include "FlowBC.h"

template<class T>
class FlowModel : public Model
{
 public:

  typedef std::map<int,FlowBC<T>*> FlowBCMap;
  typedef std::map<int,FlowVC<T>*> FlowVCMap;
  
  class Impl;
  
  
  FlowModel(const GeomFields& geomFields,
	    FlowFields& thermalFields, const MeshList& meshes);
  
  virtual ~FlowModel();

  virtual void init();
  
  virtual map<string,shared_ptr<ArrayBase> >&
    getPersistenceData();

  virtual void restart();
  
  FlowBCMap& getBCMap();
  FlowVCMap& getVCMap();

  FlowModelOptions<T>& getOptions();

  void printBCs();

  // do the specified number of iterations, return true if converged 
  bool advance(const int niter);
  void advanceMomentum();
  void advanceContinuity();

#ifdef PV_COUPLED
  bool advanceCoupled(const int niter);
#endif
  
  void updateTime();

  Vector<T,3> getPressureIntegral(const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getMomentumDerivativeIntegral(const Mesh& mesh);


  Vector<T,3> getPressureIntegralonIBFaces(const Mesh& mesh);
  Vector<T,3> getMomentumFluxIntegralonIBFaces(const Mesh& mesh);
  void getTraction(const Mesh& mesh);

  boost::shared_ptr<ArrayBase> getStressTensor(const Mesh& mesh, const ArrayBase& cellIds);

  void computeIBFaceVelocity(const StorageSite& particles);
  void computeSolidSurfaceForce(const StorageSite& particles);
  void computeSolidSurfaceForcePerUnitArea(const StorageSite& particles);
  void ComputeStressTensorES(const StorageSite& particles);
  void checkReverseFlow(const Mesh& mesh, const int faceGroupId);

  //////////////////////////////////////// Flow topology optimization //////////////////////////////////////////

  /*
#include "FlowTopologyOptimizationFunctionPrototypes.h";
#include "ThermalFlowLaminarTopologyOptimizationFunctionPrototypes.h";
  */

  //////////////////////////////////////// Flow topology optimization //////////////////////////////////////////
  void interpolateMassFluxApproximate();
  void interpolateMassFluxInteriorFacesSetFaceCellPressureOnPressureBoundarySetMassFluxOtherBoundaries();
  void updateBetaFieldGhostCells();
  void initMomApField();
  void initpreviousVelocityField();
  void interpolateSourceFieldSIMP();
  void interpolateSourceFieldSIMPRamp();
  void interpolateViscosityFieldSIMP();
  void initResidualFieldsFlow();
  void copyVelocityPressureField(const Field& velocityDoubleField, const Field& pressureDoubleField, const MeshList& meshesDouble);
  void copyVelocityBoundaryFace(const Field& velocityComponentDoubleField, Field& velocityComponentField, const StorageSite& facesDouble, const StorageSite& faces);
  void copyMomApField(const Field& momApDoubleField, const MeshList& meshesDouble);
  void copyMomentumResidual(const Field& momentumRes, const MeshList& meshesFrom);
  void copyContinuityResidual(const Field& continuityRes, const MeshList& meshesFrom);
  void copypreviousVelocityField(const Field& previousVelocityDoubleField, const MeshList& meshesDouble);
  void copypressureGradientField(const Field& pressureGradientDoubleField, const MeshList& meshesDouble);
  void copyfacePressureField(const Field& facePressureDoubleField, const MeshList& meshesDouble);
  void copyMassFluxField(const Field& massFluxFaceDoubleField, const MeshList& meshesDouble);
  void copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble);
  void copySourceField(const Field& sourceDoubleField, const MeshList& meshesDouble);
  void ComputeMomApCoefficients();
  void updateFacePressure();
  void updatePressureGradient();
  void ComputeContinuityResidual();
  void ComputeMomentumResidual();
  void ComputeMomentumResidualSaveMomAp();
  void setObjectiveFunctionFlow(const T objFunction);
  void setObjectiveFunctionFlow(const Vector<T,3> objFunction, const int component);

#if ( defined(USING_ATYPE_RAPID) )
  void setIndependentStateVariablesFlow();
  void setIndependentDesignVariablesFlow();
  void setViscosityAsIndependentVariable();
  void setSourceAsIndependentVariable();
  void setSourceAsIndependentVariable(const int cellId);
  void setSourceFieldAsIndependentVariables();
  void setDensityAsIndependentVariable();
  void WriteResidualGradientsToPetscFlow();
  void WriteResidualGradientsToMatlabFlow();
  void WriteObjectiveFunctionGradientToPetscFlow();
  void WriteObjectiveFunctionGradientToMatlabFlow();
  void WritePressureDropGradientToPetsc(const int faceId, const int direction);
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID) )
  void dumpBetaSensistivityFields();
  void dumpVelocityVector(const string fileBase);
  void dumpContinuityMatrix(const string fileBase);
  void ComputeVelocityGradient();
#endif

  // Functions for developing the code
  void printPressureIntegrals();
  void printMomentumFluxIntegrals();
  void printMassFluxIntegrals();
  void printCellPressure(string fileName);
  void printCellVelocity(string fileName);
  void printMomApField();
  void printMomApField(string fileName);
  void printpreviousVelocityField(string fileName);
  void printpressureGradientField();
  void printpressureGradientField(string fileName);
  void printfacePressureField();
  void printfacePressureField(string fileName);
  void printMassFlux();
  void printMassFlux(string fileName);
  void printMomentumResidual();
  void printMomentumResidual(string fileName);
  void printContinuityResidual();
  void printContinuityResidual(string fileName);
  void printSourceField();

#if ( defined(USING_ATYPE_RAPID) )
  void printIndependentStateVariables();
  void printIndependentDesignVariables();
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID) )
  void buildNeighborsList(const T rMin);
  void applyFilter();
#endif

#if ( defined(USING_ATYPE_RAPID) )
  void setIndependentStateVariablesThermoFluidLaminar();
  void setIndependentDesignVariablesThermoFluidLaminar();
#endif


  //Turbulence
  void ComputeVorticityMagnitude();
  void ComputeTotalViscosity();
  void ComputeShearStress(const Mesh& mesh, const int faceGroupId);
  void ComputeWallStress(const Mesh& mesh, const int faceGroupId);

  void copyEddyViscosityField(const Field& eddyViscosityField);

  Vector<T,3> getPressureGradientIntegral();
  Vector<T,3> getTotalPressureLossesBetweenInletOutlet(const Mesh& mesh, const int inlet, const int outlet);
  Vector<T,3> getDynamicPressureIntegral(const Mesh& mesh, const int faceId);

  //Vector<T,3> getDynamicPressureIntegralOutlet(const Mesh& mesh, const int outlet);
  //Vector<T,3> getDynamicPressureIntegralInlet(const Mesh& mesh, const int inlet);

  //Turbulent topology optimization
#if ( defined(USING_ATYPE_RAPID) )
  void setIndependentStateVariablesTurbulentFlow();
  void setIndependentDesignVariablesTurbulentFlow();
#endif

 private:
  shared_ptr<Impl> _impl;
};

#endif

