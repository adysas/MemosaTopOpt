// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GENERICBCS_H_
#define _GENERICBCS_H_

#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"

template<class X, class Diag, class OffDiag>
class BaseGenericBCS
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<T_Scalar> TArray;
  typedef Array<int> IntArray;
  
  typedef Vector<T_Scalar,3> VectorT3;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef FluxJacobianMatrix<Diag,X> FMatrix;
  typedef DiagonalMatrix<Diag,X> BBMatrix;

  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;
  
  typedef Array<X> XArray;
  typedef Array<VectorT3> VectorT3Array;

  
  BaseGenericBCS(const StorageSite& faces,
             const Mesh& mesh,
             const GeomFields& geomFields,
             Field& varField,
             Field& fluxField,
             MultiFieldMatrix& matrix,
             MultiField& xField, MultiField& rField) :
    _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _fluxField(fluxField),
    _xIndex(&_varField,&_cells),
    _fluxIndex(&_fluxField,&_faces),
    _dRdX(dynamic_cast<CCMatrix&>(matrix.getMatrix(_xIndex,_xIndex))),
    _dFluxdX(dynamic_cast<FMatrix&>(matrix.getMatrix(_fluxIndex,_xIndex))),
    _dFluxdFlux(dynamic_cast<BBMatrix&>(matrix.getMatrix(_fluxIndex,_fluxIndex))),
    _assembler(_dRdX.getPairWiseAssembler(_faceCells)),
    _dRdXDiag(_dRdX.getDiag()),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _flux(dynamic_cast<XArray&>(xField[_fluxIndex])),
    _rFlux(dynamic_cast<XArray&>(rField[_fluxIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _is2D(mesh.getDimension()==2)    
  {}
    
  void applyDirichletBC(int f, const X& bValue) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    // the current value of flux and its Jacobians
    const X fluxB = -_r[c1];
    const OffDiag dFluxdXC0 = -_assembler.getCoeff10(f);
    const Diag dFluxdXC1 = -_dRdXDiag[c1];
    const OffDiag dRC0dXC1 = _assembler.getCoeff01(f);
    
    // since we know the boundary value, compute the boundary
    // x correction and it's contribution to the residual for c0; we
    // can then eliminate the coefficient to the boundary cell
    
    const X dXC1 = bValue - _x[c1];
    const X dFlux = dFluxdXC1*dXC1;
    const X dRC0 = dRC0dXC1*dXC1;
    _r[c0] += dRC0;  
    
    _assembler.getCoeff01(f) = NumTypeTraits<OffDiag>::getZero();

    // set the boundary value and make its correction equation an
    // identity
    _x[c1] = bValue;
    _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getZero();
    _r[c1] = NumTypeTraits<X>::getZero();
    _dRdXDiag[c1] = NumTypeTraits<Diag>::getNegativeUnity();

    //setup the equation for the boundary flux correction
    _dFluxdX.setCoeffL(f,dFluxdXC0);
    _dFluxdX.setCoeffR(f,NumTypeTraits<OffDiag>::getZero());
    _flux[f] = fluxB;
    _rFlux[f] = dFlux;
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void applyDirichletBC(const X& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue);
  }
  
  void applyDirichletBC(const FloatValEvaluator<X>& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue[i]);
  }
  
  void applyNeumannBC(const int f,
                      const X& specifiedFlux) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    // the current value of flux and its Jacobians
    const X fluxB = -_r[c1];
        
        
    // since we know the boundary flux, compute the boundary flux
    // correction and add it to the c0 residual; also eliminate
    // coeff to the boundary cell and remove it from the ap coeff
    
    const X dFlux = specifiedFlux*_faceAreaMag[f] - fluxB;
    // setup the equation for the boundary value; the coefficients
    // are already computed so just need to set the rhs
    _r[c1] = dFlux;

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    // fix the boundary flux to the specified value
    _flux[f] = specifiedFlux*_faceAreaMag[f];
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }


  void applyNeumannBC(const X& bFlux) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNeumannBC(i,bFlux);
  }
  
  void applyNeumannBC(const FloatValEvaluator<X>& bFlux) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNeumannBC(i,bFlux[i]);
  }

  void applyExtrapolationBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
	applyExtrapolationBC(i);
  }

  // boundary value = cell value, flux as defined by interior discretization
  
  void applyExtrapolationBC(const int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
        
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // the current value of flux and its Jacobians
    const X fluxB = -_r[c1];
    const OffDiag dFluxdXC0 = -_assembler.getCoeff10(f);
    const Diag dFluxdXC1 = -_dRdXDiag[c1];

    const X xc0mxc1 = _x[c0]-_x[c1];
        
    // eliminate boundary dependency from cell equation
    _dRdXDiag[c0] += dFluxdXC1;
    _r[c0] += dFluxdXC1*xc0mxc1;
    _assembler.getCoeff01(f) = 0;

    // boundary value equation
    _dRdXDiag[c1] = NumTypeTraits<Diag>::getNegativeUnity();
    _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getUnity();
    _r[c1] = xc0mxc1;
    _dRdX.setBoundary(c1);

    //setup the equation for the boundary flux correction
    _dFluxdX.setCoeffL(f,dFluxdXC0);
    _dFluxdX.setCoeffR(f,dFluxdXC0); // should really be dFluxdXC1 
    _flux[f] = fluxB;
    _rFlux[f] = T_Scalar(0);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }
  
  void applyConvectionBC(const int f,
                         const X& hCoeff, const X& Xinf) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // the current value of flux and its Jacobians
    const X fluxInterior = -_r[c1];

    // flux based on current boundary value
    const X fluxBoundary = -hCoeff*(_x[c1]-Xinf)*_faceAreaMag[f];

    const X dFlux = fluxBoundary-fluxInterior;

    _r[c1] = dFlux;

    // add this to complete the Jacobian wrt boundar value
    _dRdXDiag[c1] -= hCoeff*_faceAreaMag[f];

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    _flux[f] = fluxBoundary;
    _rFlux[f] = 0;
    _dFluxdX.setCoeffL(f,NumTypeTraits<X>::getZero());
    _dFluxdX.setCoeffR(f,-hCoeff*_faceAreaMag[f]);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void applyConvectionBC(const X& hCoeff, const X& Xinf) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyConvectionBC(i,hCoeff,Xinf);
  }
  
  void applyRadiationBC(const int f,
                         const X& emissivity, const X& Xinf) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    //The value of the Stefan-Boltzman constant
    double s_b_const = 5.670373E-8;

    // the current value of flux and its Jacobians
    const X fluxInterior = -_r[c1];

    // flux based on current boundary value
    const X fluxBoundary = -emissivity*s_b_const*\
      (_x[c1]*_x[c1]*_x[c1]*_x[c1]-Xinf*Xinf*Xinf*Xinf)*_faceAreaMag[f];

    const X dFlux = fluxBoundary-fluxInterior;

    _r[c1] = dFlux;

    // add this to complete the Jacobian wrt boundary value
    _dRdXDiag[c1] -= \
      4*emissivity*s_b_const*_x[c1]*_x[c1]*_x[c1]*_faceAreaMag[f];

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    _flux[f] = fluxBoundary;
    _rFlux[f] = 0;
    _dFluxdX.setCoeffL(f,NumTypeTraits<X>::getZero());
    _dFluxdX.setCoeffR(f,-4*emissivity*s_b_const*_x[c1]*_x[c1]*_x[c1]*_faceAreaMag[f]);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void applyMixedBC(const int f, const X& hCoeff,
                         const X& emissivity, const X& Xinf) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    //The value of the Stefan-Boltzman constant
    double s_b_const = 5.670373E-8;

    // the current value of flux and its Jacobians
    const X fluxInterior = -_r[c1];

    // flux based on current boundary value
    const X fluxBoundary = (-emissivity*s_b_const*(_x[c1]*_x[c1]*_x[c1]*_x[c1]-Xinf*Xinf*Xinf*Xinf)-hCoeff*(_x[c1]-Xinf))*_faceAreaMag[f];

    const X dFlux = fluxBoundary-fluxInterior;

    _r[c1] = dFlux;

    // add this to complete the Jacobian wrt boundary value
    _dRdXDiag[c1] -= (4*emissivity*s_b_const*_x[c1]*_x[c1]*_x[c1]+hCoeff)*_faceAreaMag[f];

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    _flux[f] = fluxBoundary;
    _rFlux[f] = 0;
    _dFluxdX.setCoeffL(f,NumTypeTraits<X>::getZero());
    _dFluxdX.setCoeffR(f,-4*emissivity*s_b_const*_x[c1]*_x[c1]*_x[c1]*_faceAreaMag[f]);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void applyInterfaceBC(const int f) const
  {
    // the boundary cell could be either c0 or c1 at an interface
    int cb = _faceCells(f,1);
    T_Scalar sign(NumTypeTraits<T_Scalar>::getUnity());
    if (cb < _cells.getSelfCount())
    {
        cb = _faceCells(f,0);
        sign *= -1.0;
    }
    

    // the current value of flux and its Jacobians
    const X fluxInterior = -_r[cb];
    const OffDiag dFluxdXC0 = -sign*_assembler.getCoeff10(f);
    const OffDiag dFluxdXC1 = sign*_assembler.getCoeff01(f);

   
    _r[cb] = T_Scalar(0);

    if (sign>0)
      _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getZero();
    else
      _assembler.getCoeff01(f) = NumTypeTraits<OffDiag>::getZero();
    
    //setup the equation for the boundary flux correction
    _dFluxdX.setCoeffL(f,dFluxdXC0);
    _dFluxdX.setCoeffR(f,dFluxdXC1);
    _flux[f] = fluxInterior;
    _rFlux[f] = T_Scalar(0);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void applyInterfaceBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyInterfaceBC(i);
  }

  //*********************************************************************
  //special interface boundary condition for dielectric layer
  
  void applyDielectricInterfaceBC(const int f, const X& hCoeff, 
  				const X& Xinf, const X& source) const
  {
  // here the hCoeff = dielectric_constant / dielectric_thickness
  // source = totalcharge * dielectric_thickness / 4. for 2D
  // source = totalcharge * dielectric_thickness / 6. for 3D
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // the current value of flux and its Jacobians
    const X fluxInterior = -_r[c1];

    // flux based on current boundary value
   
    X fluxSource = source * _faceAreaMag[f];   
    if (_is2D)
    	fluxSource /= 2.0;
    else 
    	fluxSource /= 2.0; 
    const X fluxBoundary = -hCoeff*(_x[c1]-Xinf)*_faceAreaMag[f] + fluxSource;
    const X dFlux = fluxBoundary-fluxInterior;

    _r[c1] = dFlux;

    // add this to complete the Jacobian wrt boundar value
    _dRdXDiag[c1] -= hCoeff*_faceAreaMag[f];

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    _flux[f] = fluxBoundary;
    _rFlux[f] = 0;
    _dFluxdX.setCoeffL(f,NumTypeTraits<X>::getZero());
    _dFluxdX.setCoeffR(f,-hCoeff*_faceAreaMag[f]);
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
    
  }
    
  void applyDielectricInterfaceBC(const X& hCoeff, 
  				const X& Xinf, const X& source) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDielectricInterfaceBC(i,hCoeff, Xinf, source );
  }
 


  void applyFlowBC(const TArray& convFlux, const X& bValue) const
  {
    for(int f=0; f<_faces.getCount(); f++)
      if (convFlux[f] < 0)
        applyDirichletBC(f,bValue);
      else
        applyExtrapolationBC(f);
  }

  void applyNonzeroDiagBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNonzeroDiagBC(i);
  }

  void applyNonzeroDiagBC(int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;    
   
    _dRdXDiag[c1][0] = T_Scalar(-1.0);   
  }

  
protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const IntArray& _ibType;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const Field& _fluxField;
  const MultiField::ArrayIndex _xIndex;
  const MultiField::ArrayIndex _fluxIndex;
  CCMatrix& _dRdX;
  FMatrix& _dFluxdX;
  BBMatrix& _dFluxdFlux;
  CCAssembler& _assembler;
  DiagArray& _dRdXDiag;
  XArray& _x;
  XArray& _r;
  XArray& _flux;
  XArray& _rFlux;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const bool _is2D;
  
};

#endif

