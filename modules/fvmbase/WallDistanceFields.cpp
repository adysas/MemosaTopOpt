// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "WallDistanceFields.h"

WallDistanceFields::WallDistanceFields(const string baseName) :
  phi(baseName + ".phi"),
  phiGradient(baseName + ".phiGradient"),
  phiFlux(baseName + ".phiFlux"),
  source(baseName + ".source"),
  zero(baseName + ".zero"),
  one(baseName + ".one"),
  wallDistance(baseName + ".wallDistance"),
  sourceSolid(baseName + ".sourceSolid"),
  beta(baseName + ".beta"),
  wallDistanceResidual(baseName + ".wallDistanceResidual")
{}


