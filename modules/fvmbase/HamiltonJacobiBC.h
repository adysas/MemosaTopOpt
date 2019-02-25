// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct HamiltonJacobiBC : public FloatVarDict<T>
{
  HamiltonJacobiBC()
    {
      this->defineVar("specifiedPhi",T(0.0));
      this->defineVar("specifiedPhiFlux",T(0.0));
    }
  string bcType;
};

template<class T>
struct HamiltonJacobiVC : public FloatVarDict<T>
{
  HamiltonJacobiVC()
    {
      this->defineVar("gamma",T(-1.0));
      this->defineVar("density", T(1.0));
      this->defineVar("specificHeat", T(1.0));      
    }
  string vcType;
};

template<class T>
struct HamiltonJacobiModelOptions : public FloatVarDict<T>
{
  HamiltonJacobiModelOptions()
    {
      this->defineVar("initialPhi",T(0.0));
      this->defineVar("volumeFraction",T(0.5));
      this->defineVar("epsilon",T(0.1));
      this->defineVar("phiURF",T(0.7));
      this->relativeTolerance=1e-8;
      this->absoluteTolerance=1e-16;
      this->linearSolver = 0;
      this->useCentralDifference=false;
    }
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
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

