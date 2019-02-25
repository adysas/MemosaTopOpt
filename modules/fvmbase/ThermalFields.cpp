// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "ThermalFields.h"

ThermalFields::ThermalFields(const string baseName) :
  temperature(baseName + ".temperature"),
  temperatureN1(baseName + ".temperatureN1"),
  temperatureN2(baseName + ".temperatureN2"),
  heatFlux(baseName + ".heatFlux"),
  temperatureGradient(baseName + ".temperatureGradient"),
  conductivity(baseName + ".conductivity"),
  source(baseName + ".source"),
  convectionFlux(baseName + ".convectionFlux"),
  zero(baseName + "zero"),
  one(baseName + "one"),
  specificHeat(baseName + ".specificHeat"),
  density(baseName + ".density"),
  volumetricHeatCapacity(baseName + ".volumetricHeatCapacity"),
  beta(baseName + ".beta"),
  sensitivity(baseName + ".sensitivity"),
  momentumResidual(baseName + ".momResidual"),
  continuityResidual(baseName + ".continuityResidual"),
  temperatureResidual(baseName + ".temperatureResidual"),
  sensitivityConstraint(baseName + ".sensitivityConstraint"),
  heatTransferCoefficient(baseName + ".heatTransferCoefficient")
{}


