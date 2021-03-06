// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "SpalartAllmarasFields.h"

SpalartAllmarasFields::SpalartAllmarasFields(const string baseName) :
  nuTilde(baseName + ".nuTilde"),
  nuTildeN1(baseName + ".nuTildeN1"),
  nuTildeN2(baseName + ".nuTildeN2"),
  nuTildeFlux(baseName + ".nuTildeFlux"),
  density(baseName + ".density"),
  nu(baseName + ".nu"),
  nuTildeGradient(baseName + ".nuTildeGradient"),
  sourceTerm1(baseName + ".sourceFieldTerm1"),
  sourceTerm2(baseName + ".sourceFieldTerm2"),
  convectionFlux(baseName + ".convectionFlux"),
  zero(baseName + ".zero"),
  one(baseName + ".one"),
  diffusivity(baseName + ".diffusivity"),
  vorticityMagnitude(baseName + ".vorticityMagnitude"),
  chi(baseName + ".chi"),
  fv1(baseName + ".fv1"),
  fv2(baseName + ".fv2"),
  STilde(baseName + ".STilde"),
  r(baseName + ".r"),
  g(baseName + ".g"),
  fw(baseName + ".fw"),
  ft2(baseName + ".ft2"),
  eddyViscosity(baseName + ".eddyViscosity"),
  wallDistance(baseName + ".wallDistance"),
  source(baseName + ".source"),
  spalartAllmarasResidual(baseName + ".spalartAllmarasResidual"),
  momentumResidual(baseName + ".momentumResidual"),
  continuityResidual(baseName + ".continuityResidual"),
  wallDistanceResidual(baseName + ".wallDistanceResidual"),  
  beta(baseName + ".beta"),
  sensitivity(baseName + ".sensitivity")
{}


