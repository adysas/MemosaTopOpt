// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct SpalartAllmarasBC : public FloatVarDict<T>
{
  SpalartAllmarasBC()
    {
      this->defineVar("specifiedNuTilde",T(0.0001));
      this->defineVar("specifiedNuTildeFlux",T(0.0));
    }
  string bcType;
};

template<class T>
struct SpalartAllmarasVC : public FloatVarDict<T>
{
  SpalartAllmarasVC()
    {
      this->defineVar("viscosity",T(1.0));
      this->defineVar("density", T(1.0));
      this->defineVar("specificHeat", T(1.0));	
    }
  string vcType;
};

template<class T>
struct SpalartAllmarasModelOptions : public FloatVarDict<T>
{
  SpalartAllmarasModelOptions()
    {
      this->defineVar("initialNuTilde", T(0.0001));
      this->defineVar("timeStep", T(1e-7));

      ///////////////////// SA model///////////
      this->defineVar("sigma",T(0.6666666666666666));
      this->defineVar("Cb1",T(0.1355));
      this->defineVar("Cb2",T(0.622));
      this->defineVar("kappa",T(0.41));
      this->defineVar("Cw2",T(0.3));
      this->defineVar("Cw3", T(2.0));
      this->defineVar("Cv1",T(7.1));
      this->defineVar("Ct4", T(2.0));  
  
      this->defineVar("nuTildeURF", T(0.5));
      
      this->defineVar("alphaSolid", T(1.0e5));
      this->defineVar("alphaFluid", T(0.));
      this->defineVar("p", T(3.0));
      this->defineVar("volumeFraction", T(0.4));

      this->defineVar("scaleSpalartAllmarasModelObjFunction", T(0.0));
      this->defineVar("scaleFlowModelObjFunction",T(1.0));
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

