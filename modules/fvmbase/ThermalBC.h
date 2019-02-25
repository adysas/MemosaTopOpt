// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct ThermalBC : public FloatVarDict<T>
{
  ThermalBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedHeatFlux",T(0.0));
      this->defineVar("convectiveCoefficient", T(0.0));
      this->defineVar("farFieldTemperature", T(300.0));
      this->defineVar("surfaceEmissivity",T(1.0));
  }
  string bcType;
};

template<class T>
struct ThermalVC : public FloatVarDict<T>
{
  ThermalVC()
    {
      this->defineVar("thermalConductivity",T(10.0));
      this->defineVar("thermalConductivityHigh",T(1.0));
      this->defineVar("thermalConductivityLow",T(0.001));

      this->defineVar("density", T(1.0));
      this->defineVar("densityHigh",T(1.0));
      this->defineVar("densityLow",T(0.001));

      this->defineVar("specificHeat", T(1.0));
      this->defineVar("specificHeatHigh", T(1.0));
      this->defineVar("specificHeatLow", T(0.001));
      	
    }
  string vcType;
};

template<class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
  ThermalModelOptions()
  {
    this->defineVar("initialTemperature",T(300.0));
    this->defineVar("timeStep", T(1e-7));

    /////////////////////Top Opt///////////
    this->defineVar("alpha",T(0.));
    this->defineVar("volumeFraction",T(0.5));
    this->defineVar("p", T(3.0)); 
    this->defineVar("maxHeatFluxFaceId", T(1)); 
    this->defineVar("niter", T(1));
    this->defineVar("scaleThermalModelObjFunction", T(1.));
    this->defineVar("scaleFlowModelObjFunction", T(1.));
    this->defineVar("heatTransferCoefficient", T(0.));
    //The following has to be in volume condition; keeping here for simplicity
    this->defineVar("TemperatureAmbient", T(0.));
    ///////////////////////////////////////

    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
    this->linearSolver = 0;
    this->useCentralDifference=false;
    this->transient=false;
    this->timeDiscretizationOrder = 1;
  }
  double relativeTolerance;
  double absoluteTolerance;
  bool useCentralDifference;
  LinearSolver *linearSolver;
  bool transient;
  int timeDiscretizationOrder;
#ifndef SWIG
  LinearSolver& getLinearSolver()
  {
    if (this->linearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolver = ls;
    }
    return *this->linearSolver ;
  }
#endif
};

