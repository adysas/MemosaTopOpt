// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _THERMALFIELDS_H_
#define _THERMALFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ThermalFields
{
  ThermalFields(const string baseName);

  Field temperature;
  Field temperatureN1;
  Field temperatureN2;
  Field specificHeat;           //(J/(kg.K))
  Field density;                //(kg/m^3)
  Field volumetricHeatCapacity; //(J/(m^3.K))
  Field heatFlux;
  Field temperatureGradient;
  Field conductivity;
  Field source;
  Field convectionFlux;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
  //Added for RAPID
  Field beta;
  Field sensitivity;

  //Added for thermo-fluid applications
  Field momentumResidual;
  Field continuityResidual;
  Field temperatureResidual;

  //this is for pressure drop constraint
  Field sensitivityConstraint;
  Field heatTransferCoefficient;
};

#endif
