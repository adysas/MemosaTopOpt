// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ThermalFields.h"
#include "Mesh.h"
#include "ThermalBC.h"

template<class T>
class ThermalModel : public Model
{
public:

  typedef std::map<int,ThermalBC<T>*> ThermalBCMap;
  typedef std::map<int,ThermalVC<T>*> ThermalVCMap;
  typedef Vector<T,3> VectorT3;
  class Impl;
  
  
  ThermalModel(const GeomFields& geomFields,
               ThermalFields& thermalFields, const MeshList& meshes);
  
  virtual ~ThermalModel();

  virtual void init();
  
  ThermalBCMap& getBCMap();
  ThermalVCMap& getVCMap();
  
  ThermalBC<T>& getBC(const int id);

  ThermalModelOptions<T>& getOptions();
  

  void advance(const int niter);
  void updateTime();
  T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
  void printBCs();

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void dumpMatrix(const string fileBase);
#endif

  //#include "ThermalTopologyOptimizationFunctionPrototypes.h"

  //////////////////////////Thermal topology optimization///////////////////////////////////////////////////////////

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void buildNeighborsList(const T rMin); 
  void applyFilter(); 
  void applyFilterConstraintFunction();
#endif


  void initResidualFieldThermal();
  void updateBetaFieldGhostCells();
  void copyTemperatureField(const Field& tempDoubleField, const MeshList& meshesDouble);
  void copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble);
  void updateSourceFieldType1();
  void updateSourceFieldType2();
  void updateSourceFieldType3();
  void interpolateThermalConductivityFieldSIMP();
  void interpolateMaterialThermalPropertiesFieldSIMP();
  void interpolateThermalConductivityFieldSIMPReverse();
  void interpolateThermalConductivityFieldSIMPRamp();
  void interpolateHeatTransferCoefficient();
  void interpolateRadiallyVaryingHeatTransferCoefficient();
  void ComputeTemperatureResidual();
  T getCellTemperature(const Mesh& mesh, const int cellId);
  T getAverageTemperature(const Mesh& mesh);
  void setObjectiveFunctionThermal(const T objFunction);

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void dumpBetaSensistivityFields();
  void dumpSensistivityConstraintFields();
  void dumpVolumeField();

#endif

#if (defined(USING_ATYPE_RAPID))
  void setIndependentStateVariablesThermal();
  void setIndependentDesignVariablesThermal();
  void WriteResidualGradientsToPetscThermal();
  void WriteObjectiveFunctionGradientToPetscThermal();
#endif

  void initializeBetaField();
  void printResidualField();
  void printTemperatureField();
  void printBetaField();
  void printConductivityField();
  void printBoundaryCellNumbers(const Mesh& mesh, const int faceGroupId);


  /////////////////////////////////////////////////////////////////////////////////

  T getDeltaEnthalpyIntegralThermoFluid(const Mesh& mesh, const int faceGroupId1, const int faceGroupId2);
  T getBulkTemperature(const Mesh& mesh, const int faceGroupId);
  T getRootMeanSquareTemperatureOutletTemperatureIdeal(const Mesh& mesh, const int faceGroupIdIdeal, const int faceGroupIdOut);

  void copyMomentumResidual(const Field& momentumRes);
  void copyContinuityResidual(const Field& continuityRes);
  void copyMassFlux(const Field& massFlux);
  void setObjectiveFunctionFlow(const VectorT3 objFunction, const int direction);
  void setPressureDropComponent(const VectorT3 pDropComp, const int direction);

#if ( defined(USING_ATYPE_RAPID) )
  void setIndependentStateVariablesThermoFluidLaminar();
  void setIndependentDesignVariablesThermoFluidLaminar();
  void WriteResidualGradientsToPetscThermoFluidLaminar();
  void WriteResidualGradientsToMatlabThermoFluidLaminar();
  void WriteObjectiveFunctionGradientToPetscThermoFluidLaminar();
  void WritePressureDropGradientToPetscThermoFluidLaminar();
  void WritePressureDropGradientToMatlabThermoFluidLaminar();
#endif

  void copyBetaField(const Field& betaField);
  void copyTemperatureBoundaryFace(const Field& temperatureDoubleField, Field& temperatureField, const StorageSite& facesDouble, const StorageSite& faces);
  ////////////////////////////////////////////////////////////////////////////////


private:
  shared_ptr<Impl> _impl;
};

#endif

