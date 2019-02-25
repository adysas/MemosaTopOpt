// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct WallDistanceBC : public FloatVarDict<T>
{
  WallDistanceBC()
    {
      this->defineVar("specifiedPhi",T(0.0));
      this->defineVar("specifiedPhiFlux",T(0.0));
    }
  string bcType;
};


template<class T>
struct WallDistanceModelOptions : public FloatVarDict<T>
{
  WallDistanceModelOptions()
    {
      this->defineVar("initialPhi",T(0.0));
      this->relativeTolerance=1e-8;
      this->absoluteTolerance=1e-16;
      this->linearSolver = 0;

      this->defineVar("alphaSolid", T(1.0e5));
      this->defineVar("alphaFluid", T(0.));
      this->defineVar("p", T(3.0));
      this->defineVar("volumeFraction", T(0.4));
    }
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
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

