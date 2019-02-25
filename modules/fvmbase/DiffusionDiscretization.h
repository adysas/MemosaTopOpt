// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _DIFFUSIONDISCRETIZATION_H_
#define _DIFFUSIONDISCRETIZATION_H_

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

template<class T>
inline T harmonicAverage(const T& x0, const T& x1)
{
  const T sum = x0+x1;
  if (x0+x1 != NumTypeTraits<T>::getZero())
    return 2.0*x0*x1/sum;
  else
    return sum;
}

  
template<class X, class Diag, class OffDiag>
  class DiffusionDiscretization : public Discretization
{
 public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Gradient<X> XGrad;
  typedef Array<int> IntArray;
  
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;
  
 DiffusionDiscretization(const MeshList& meshes,
			 const GeomFields& geomFields,
			 Field& varField,
			 const Field& diffusivityField,
			 const Field& varGradientField,
			 const T_Scalar thickness=T_Scalar(0.0)) :
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _diffusivityField(diffusivityField),
    _varGradientField(varGradientField),
    _thickness(thickness)
    {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix, MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

#if !(defined(USING_ATYPE_RAPID))
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex, cVarIndex));   
    DiagArray& diag = matrix.getDiag();
    OffDiagArray& offDiag = matrix.getOffDiag();
    const CRConnectivity& conn=matrix.getConnectivity();
    const Array<int>& row=conn.getRow();
    const Array<int>& col=conn.getCol();
#endif
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
   
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    const GradArray& xGradCell =
      dynamic_cast<const GradArray&>(_varGradientField[cells]);

    const TArray& diffCell =
      dynamic_cast<const TArray&>(_diffusivityField[cells]);

   
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;

	const StorageSite& faces = fg.site;
	const int nFaces = faces.getCount();
	const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	const VectorT3Array& faceArea =
	  dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);    
	const TArray& faceAreaMag =
	  dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	const VectorT3Array& faceCentroid =
	  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);

#if !(defined(USING_ATYPE_RAPID))
	//DiagArray& diag = matrix.getDiag();	
	CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
#endif

	for(int f=0; f<nFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);

	    T_Scalar vol0 = cellVolume[c0];
	    T_Scalar vol1 = cellVolume[c1];
        
	    VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
        
	    T_Scalar faceDiffusivity(1.0);
	    
	    if (vol0 == 0.)
	      faceDiffusivity = diffCell[c1];
	    else if (vol1 == 0.)
	      faceDiffusivity = diffCell[c0];
	    else
	      faceDiffusivity = harmonicAverage(diffCell[c0],diffCell[c1]);
	
	    const T_Scalar diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);
	    const T_Scalar diffCoeff = faceDiffusivity*diffMetric;
	    const VectorT3 secondaryCoeff = faceDiffusivity*(faceArea[f]-ds*diffMetric);
        
	    const XGrad gradF = (xGradCell[c0]*vol0+xGradCell[c1]*vol1)/(vol0+vol1);

	    const X dFluxSecondary = gradF*secondaryCoeff;

#if (defined(USING_ATYPE_RAPID))
	    const X dFlux = diffCoeff*(xCell[c1]-xCell[c0]) + dFluxSecondary;
	    rCell[c0] += dFlux;
	    rCell[c1] -= dFlux;
#endif

#if !(defined(USING_ATYPE_RAPID))
	    
	    const X dFlux = diffCoeff*(xCell[c1]-xCell[c0]) + dFluxSecondary;
	    /*
	    cout << "----------------- Check ----------------------------------"<<endl;
	    dumpMatrixandRfield(mesh, mfmatrix, rField, "Kali");
	    cout << "----------------- Check ----------------------------------"<<endl;
	    */
	    rCell[c0] += dFlux;
	    rCell[c1] -= dFlux;


	    assembler.getCoeff01(f) +=diffCoeff;

	    /*
	    cout << "----------------- Check ----------------------------------"<<endl;
	    dumpMatrixandRfield(mesh, mfmatrix, rField, "DD");
	    cout << "----------------- Check ----------------------------------"<<endl;
	    */
	    assembler.getCoeff10(f) +=diffCoeff;

	    /*
	    cout << "----------------- Check ----------------------------------"<<endl;
	    dumpMatrixandRfield(mesh, mfmatrix, rField, "DD");
	    cout << "----------------- Check ----------------------------------"<<endl;
	    */

	    diag[c0] -= diffCoeff;
	    /*
	    cout << "----------------- Check ----------------------------------"<<endl;
	    dumpMatrixandRfield(mesh, mfmatrix, rField, "DD");
	    cout << "----------------- Check ----------------------------------"<<endl;
	    */
	    
	    diag[c1] -= diffCoeff;
	    /*
	    cout << "----------------- Check ----------------------------------"<<endl;
	    dumpMatrixandRfield(mesh, mfmatrix, rField, "DD");
	    cout << "----------------- Check ----------------------------------"<<endl;
	    */
#endif
	    //cout <<"face " << f << ": dflux" << dFlux << endl;
	    
	    //cout << " c0 and c1 " << c0 << " " << c1 << endl;
	    //cout << " temperatures" << xCell[c0] << " " << xCell[c1] << endl;
	    //cout << " diffcell " << diffCell[c0] << " " << diffCell[c1] << " " << faceDiffusivity << endl;
	    //cout << " diffmetric  " << diffMetric << endl;
	    //cout << " diffCoeff  " << diffCoeff << endl;
	    //cout << " ds " << ds << endl;
	    //cout << " Cell centroid " << cellCentroid[c0] << "           " << cellCentroid[c1] << endl;
	    
	  }
      }
#if !(defined(USING_ATYPE_RAPID))
    /*
    cout << "-----------------POST DD----------------------------------"<<endl;
    const int nCells = cells.getCountLevel1();
    for(int i=0; i<nCells; i++)
      {
	cout << i+1 << "," << diag[i] << "; " ;
	//	<< row[i] << ", " << row[i+1] << ", " << col[i] << endl;
	for (int j=row[i]; j<row[i+1]; j++)
	  cout << i+1 << "," << col[j]+1 << "=" << offDiag[j] << ";";
	cout << " R=" << rCell[i] << endl;
      }
      cout << "-----------------POST DD----------------------------------"<<endl;*/
#endif
  }

  void dumpMatrixandRfield(const Mesh& mesh, MultiFieldMatrix& mfmatrix, MultiField& b, string fileBase)
  {
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))

    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getCount();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

    //////////////////////////////////////////////////////////////////
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex, cVarIndex));   
    DiagArray& diag = matrix.getDiag();
    OffDiagArray& offDiag = matrix.getOffDiag();


    const CRConnectivity& cr = matrix.getConnectivity();

    const Array<int>& row = cr.getRow();
    const Array<int>& col = cr.getCol();


    string coeffFileName = fileBase + ".coeff";
	
    FILE *matFile = fopen(coeffFileName.c_str(),"wb");

    cout << "********************************************************************"<<endl;
    for(int i=0; i<nCells; i++)
      {

	//cout << vDiag[i] << vCoeff[i] << endl;
	//fprintf(matFile, "%d %d %le %le\n", i+1, i+1, vDiag[i][0], vDiag[i][1]); 
	cout << i+1 << "," << diag[i] << "; " ;
	fprintf(matFile, "%d %d %22.15le \n", i+1, i+1, diag[i][0]); 
	//	<< row[i] << ", " << row[i+1] << ", " << col[i] << endl;
	for (int j=row[i]; j<row[i+1]; j++)
	  {
	    cout << i+1 << "," << col[j]+1 << "=" << offDiag[j] << "; ";
	    fprintf(matFile, "%d %d %22.15le \n", i+1, col[j]+1, offDiag[j]); 
	  }
	cout << endl;

      }
    cout << "********************************************************************"<<endl;
    fclose(matFile);

#endif    
  }


 private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _diffusivityField; 
  const Field& _varGradientField;
  const T_Scalar _thickness;
};

#endif
