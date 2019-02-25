// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SPALARTALLMARASGRADIENTDISCRETIZATION_H_
#define _SPALARTALLMARASGRADIENTDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "Gradient.h"
#include "DiagonalTensor.h"


template<class X>
  class SpalartAllmarasGradientDiscretization : public Discretization
{
 public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef Gradient<X> XGrad;
  typedef Array<X> XArray;
  typedef Array<XGrad> GradArray;
  typedef Array<T_Scalar> TArray;

 SpalartAllmarasGradientDiscretization(const MeshList& meshes,
				       const GeomFields& geomFields,
				       Field& varField,
				       const Field& varGradientField,
				       const Field& densityField, 
				       const T_Scalar Cb2,
				       const T_Scalar sigma) :
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varGradientField(varGradientField),
    _densityField(densityField),
    _Cb2(Cb2),
    _sigma(sigma)
    {}

  void discretize(const Mesh& mesh, MultiFieldMatrix&,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
   
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
  
    const TArray& density =
      dynamic_cast<const TArray&>(_densityField[cells]);

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    const GradArray& xGradCell =
      dynamic_cast<const GradArray&>(_varGradientField[cells]);

    const int nCells = cells.getSelfCount();
    
    for(int c=0; c<nCells; c++)
      {
	T_Scalar gradNuSquare = (pow(xGradCell[c][0],2) + pow(xGradCell[c][1],2) + pow(xGradCell[c][2],2));
        rCell[c] += _Cb2*_sigma*density[c]*gradNuSquare*cellVolume[c];
      }

  }
 private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _varGradientField;
  const Field& _densityField;
  const T_Scalar _Cb2;
  const T_Scalar _sigma;
};

#endif
