// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _HAMILTONJACOBIFIELDS_H_
#define _HAMILTONJACOBIFIELDS_H_


#include "Field.h"

struct HamiltonJacobiFields
{
  HamiltonJacobiFields(const string baseName);

  Field phi;
  Field phiGradient;
  Field convectionFlux;
  Field gamma;
  Field density;
  Field zero;
  Field one;
  Field continuityResidual;

  Field phiFlux;
  
  Field wallDistance;
  // Source field to include Ergun model or Brinkman's model one field
  Field source;
  Field sourceSolid;

  //For RAPID
  Field beta;

};

#endif
