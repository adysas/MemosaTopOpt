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
  typedef std::map<int,HamiltonJacobiVC<T>*> HamiltonJacobiVCMap;
  typedef Vector<T,3> VectorT3;
  class Impl;
  
  
  HamiltonJacobiModel(const GeomFields& geomFields,
               HamiltonJacobiFields& hamiltonJacobiFields, const MeshList& meshes);
  
  virtual ~HamiltonJacobiModel();

  virtual void init();
  
  HamiltonJacobiBCMap& getBCMap();
  HamiltonJacobiVCMap& getVCMap();
  
  HamiltonJacobiBC<T>& getBC(const int id);

  HamiltonJacobiModelOptions<T>& getOptions();
  
  void advance(const int niter);
  T getPhiFluxIntegral(const Mesh& mesh, const int faceGroupId);
  void ComputeWallDistance();

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void dumpMatrix(const string fileBase);
#endif


private:
  shared_ptr<Impl> _impl;
};

#endif

