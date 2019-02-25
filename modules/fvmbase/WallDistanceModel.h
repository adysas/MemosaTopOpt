// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _WALLDISTANCEMODEL_H_
#define _WALLDISTANCEMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "WallDistanceFields.h"

#include "Mesh.h"
#include "WallDistanceBC.h"

template<class T>
class WallDistanceModel : public Model
{
public:

  typedef std::map<int,WallDistanceBC<T>*> WallDistanceBCMap;

  class Impl;
  
  
  WallDistanceModel(const GeomFields& geomFields,
               WallDistanceFields& wallDistanceFields, const MeshList& meshes);
  
  virtual ~WallDistanceModel();

  virtual void init();
  
  WallDistanceBCMap& getBCMap();
  
  WallDistanceBC<T>& getBC(const int id);

  WallDistanceModelOptions<T>& getOptions();
  

  T getPhiFluxIntegral(const Mesh& mesh, const int faceGroupId);
  
  void printBCs();

  void advance(const int niter);
  void ComputeWallDistance();

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void dumpMatrix(const string fileBase);
#endif

  void interpolateSourceFieldSIMP();
  void copyPhiField(const Field& phiDoubleField, const MeshList& meshesDouble);
  void copyBetaField(const Field& betaField);
  void ComputeWallDistanceModelResidual();
  void initResidualFieldWallDistanceModel();
  
#if defined(USING_ATYPE_RAPID)
  void setIndependentStateVariablesFlowTurbulent();
  //void setIndependentDesignVariablesFlowTurbulent();
#endif
  //void updateBetaFieldGhostCells();
  
private:
  shared_ptr<Impl> _impl;
};

#endif

