// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FlowFields.h"

FlowFields::FlowFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  massFlux(baseName + ".massFlux"),
  velocityGradient(baseName + ".velocityGradient"),
  pressureGradient(baseName + ".pressureGradient"),
  momentumFlux(baseName + ".momentumFlux"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  continuityResidual(baseName + ".continuityResidual"),
  velocityN1(baseName + ".velocityN1"),
  velocityN2(baseName + ".velocityN2"),
  tractionX(baseName + ".tractionX"),
  tractionY(baseName + ".tractionY"),
  tractionZ(baseName + ".tractionZ"),
  stress(baseName + ".stress"),
  force(baseName + ".force"),
  eddyviscosity(baseName + ".eddyviscosity"),
  totalviscosity(baseName + ".totalviscosity"),
  utau(baseName + ".utau"),
  uparallel(baseName + ".uparallel"),
  tau(baseName + ".tau"),
  tauwall(baseName + ".tauwall"),

  //ESInterface
  InterfaceVelocity(baseName+".InterfaceVelocity"),
  InterfacePressure(baseName+".InterfacePressure"), 
  InterfaceStress(baseName+".InterfaceStress"),
  InterfaceDensity(baseName+".InterfaceDensity"),

//Source field for Ergun model / Brinkman model
  source(baseName+".source"),
  residualMomentum(baseName + ".residualMomentum"),
  residualContinuity(baseName + ".residualContinuity"),
  beta(baseName + ".beta"),
  sensitivity(baseName + ".sensitivity"),
  momentumAp(baseName + ".momentumAp"),
  previousVelocity(baseName + ".previousVelocity"),

// for SA turbulence
  vorticityMagnitude(baseName + ".vorticityMagnitude"),
  tauWall(baseName + ".tauWall"),
  wallDistance(baseName + ".wallDistance")
{}

