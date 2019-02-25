// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <atype.h>

#include "MeshMetricsCalculator.h"
#include "MeshMetricsCalculator_impl.h"

template class MeshMetricsCalculator<ATYPE>;


#include "ThermalModel.h"
#include "ThermalModel_impl.h"
template class ThermalModel<ATYPE>;



#ifndef USING_ATYPE_PC

#endif



#include "FlowModel.h"
#include "FlowModel_impl.h"
template class FlowModel<ATYPE>;


#include "SpalartAllmarasModel.h"
#include "SpalartAllmarasModel_impl.h"
template class SpalartAllmarasModel<ATYPE>;

#include "WallDistanceModel.h"
#include "WallDistanceModel_impl.h"
template class WallDistanceModel<ATYPE>;

#include "HamiltonJacobiModel.h"
#include "HamiltonJacobiModel_impl.h"
template class HamiltonJacobiModel<ATYPE>;


#ifndef USING_ATYPE_PC

#endif



