// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _HAMILTONJACOBIMODEL_H_
#define _HAMILTONJACOBIMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "HamiltonJacobiFields.h"

#include "Mesh.h"
#include "HamiltonJacobiBC.h"

template<class T>
class HamiltonJacobiModel : public Model
{
public:

  typedef std::map<int,HamiltonJacobiBC<T>*> HamiltonJacobiBCMap;

  class Impl;
  
  
  HamiltonJacobiModel(const GeomFields& geomFields,
               HamiltonJacobiFields& hamiltonJacobiFields, const MeshList& meshes);
  
  virtual ~HamiltonJacobiModel();

  virtual void init();
  
  HamiltonJacobiBCMap& getBCMap();
  
  HamiltonJacobiBC<T>& getBC(const int id);

  HamiltonJacobiModelOptions<T>& getOptions();
  

  T getPhiFluxIntegral(const Mesh& mesh, const int faceGroupId);
  
  void printBCs();

  void advance(const int niter);
  void advanceDiffusionOnly(const int niter);
  void ComputeWallDistance();

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void dumpMatrix(const string fileBase);
#endif

private:
  shared_ptr<Impl> _impl;
};

#endif

