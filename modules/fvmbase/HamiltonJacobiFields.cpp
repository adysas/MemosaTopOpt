// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "HamiltonJacobiFields.h"

HamiltonJacobiFields::HamiltonJacobiFields(const string baseName) :
  phi(baseName + ".phi"),
  phiGradient(baseName + ".phiGradient"),
  convectionFlux(baseName + ".convectionFlux"),
  gamma(baseName + ".gamma"),
  density(baseName + ".density"),
  zero(baseName + ".zero"),
  one(baseName + ".one"),
  continuityResidual(baseName + ".continuityResidual"),
  phiFlux(baseName + ".phiFlux"),
  wallDistance(baseName + ".wallDistance"),
  source(baseName + ".source"),
  sourceSolid(baseName + ".sourceSolid"),
  beta(baseName + ".beta")
{}

