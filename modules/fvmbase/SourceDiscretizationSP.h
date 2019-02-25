// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SOURCEDISCRETIZATIONSP_H_
#define _SOURCEDISCRETIZATIONSP_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "CRMatrixRect.h"

template<class X, class Diag, class OffDiag>
  class SourceDiscretizationSP : public Discretization
{
 public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;  
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;

 SourceDiscretizationSP(const MeshList& meshes,
			const GeomFields& geomFields,
			const Field& varField,
			const Field& sourceField)  :
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _sourceField(sourceField)
    {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    const TArray& source = dynamic_cast<const TArray&>(_sourceField[cells]);
    const MultiField::ArrayIndex cVarIndex(&_varField, &cells);
    /*
#if !(defined(USING_ATYPE_RAPID))
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex, cVarIndex));   
    DiagArray& diag = matrix.getDiag();
    OffDiagArray& offDiag = matrix.getOffDiag();
    const CRConnectivity& conn=matrix.getConnectivity();
    const Array<int>& row=conn.getRow();
    const Array<int>& col=conn.getCol();
#endif
    */

#if !(defined(USING_ATYPE_RAPID))
    DiagArray& diag = matrix.getDiag();
#endif

    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);   
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    const int nCells = cells.getSelfCount();
    /*
      for(int c=0; c<nCells; c++)
      {
      rCell[c] += cellVolume[c]*source[c]*xCell[c];
      #if !(defined(USING_ATYPE_RAPID))
      diag[c] += cellVolume[c]*source[c];
      #endif

      }
    */
    for(int c=0; c<nCells; c++)
      {
	if (source[c] == 1.0e5)
	  {
	    rCell[c] = 0;
#if !(defined(USING_ATYPE_RAPID))
	    diag[c] = 1;
#endif
	  }
      }
    /*
#if !(defined(USING_ATYPE_RAPID))
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;
	const StorageSite& faces = fg.site;
	const int nFaces = faces.getCount();
	const CRConnectivity& faceCells = mesh.getFaceCells(faces);	
	CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
	for(int f=0; f<nFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    if (source[c0] == 1.0e5)
	      {
		assembler.getCoeff01(f) = 0;
	      }
	    else if (source[c1] == 1.0e5)
	      {
		assembler.getCoeff10(f) = 0;
	      }
	  }
      }
  }
#endif
    */
 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _sourceField;
};

#endif
