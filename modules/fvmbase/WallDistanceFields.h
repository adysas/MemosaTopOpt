// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PHIFIELDS_H_
#define _PHIFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct WallDistanceFields
{
  WallDistanceFields(const string baseName);

  Field phi;
  Field phiFlux;
  Field phiGradient;
  Field source;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
  Field wallDistance;

  //Topology optimization
  Field sourceSolid;
  Field beta;
  Field wallDistanceResidual;
  
};

#endif
