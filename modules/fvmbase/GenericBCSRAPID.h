// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GENERICBCSRAPID_H_
#define _GENERICBCSRAPID_H_

#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"

template<class X>
  class BaseGenericBCSRapid
{
 public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<T_Scalar> TArray;
  typedef Array<int> IntArray;
  
  typedef Vector<T_Scalar,3> VectorT3;
  
  
  typedef Array<X> XArray;
  typedef Array<VectorT3> VectorT3Array;

  
 BaseGenericBCSRapid(const StorageSite& faces,
		const Mesh& mesh,
		const GeomFields& geomFields,
		Field& varField,
		Field& fluxField,
		MultiField& xField, 
		MultiField& rField) :
  _faces(faces),
    _cells(mesh.getCells()),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _fluxField(fluxField),
    _xIndex(&_varField,&_cells),
    _fluxIndex(&_fluxField,&_faces),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _flux(dynamic_cast<XArray&>(xField[_fluxIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _is2D(mesh.getDimension()==2)    
      {}
    
  void applyDirichletBC(int f, const X& bValue) const
  {
    const int c1 = _faceCells(f,1);
    //cout << "Before Flux - " << c1 << "-->" << _flux[f] << endl;
    const X fluxB = -_r[c1];
    _r[c1] = bValue - _x[c1];
    _flux[f] = fluxB;
    
    //_flux[f] = -_r[c1];
    //cout << "Flux - " << c1 << "-->" << _flux[f] << endl;
    //cout << "rc1 - " << c1 << "-->" << _r[c1] << endl; 
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
    const int c1 = _faceCells(f,1);
    const X fluxB = -_r[c1];
    _r[c1] = specifiedFlux*_faceAreaMag[f] - fluxB;
    _flux[f] = fluxB;
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
    const int c1 = _faceCells(f,1);
    const int c0 = _faceCells(f,0);
    //cout << "Before Flux - " << c1 << "-->" << _flux[f] << endl;
    const X fluxB = -_r[c1];
    _r[c1] = _x[c0] - _x[c1];
    _flux[f] = fluxB;
    
  }
  
  void applyConvectionBC(const int f,
                         const X& hCoeff, const X& Xinf) const
  {
    /*
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
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();*/
    throw;
  }

  void applyConvectionBC(const X& hCoeff, const X& Xinf) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyConvectionBC(i,hCoeff,Xinf);
  }
  
  void applyRadiationBC(const int f,
			const X& emissivity, const X& Xinf) const
  {
    /*
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
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();*/
    throw;
  }

  void applyMixedBC(const int f, const X& hCoeff,
		    const X& emissivity, const X& Xinf) const
  {
    /*
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
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();*/
    throw;
  }

  void applyInterfaceBC(const int f) const
  {
    /*
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
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();*/
    throw;
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
    /*
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
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();*/
    throw;
    
  }
    
  void applyDielectricInterfaceBC(const X& hCoeff, 
				  const X& Xinf, const X& source) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDielectricInterfaceBC(i,hCoeff, Xinf, source );
  }
 


  void applyFlowBC(const TArray& convFlux, const X& bValue) const
  {
    /*
    for(int f=0; f<_faces.getCount(); f++)
      if (convFlux[f] < 0)
        applyDirichletBC(f,bValue);
      else
      applyExtrapolationBC(f);*/
    throw;
  }

  void applyNonzeroDiagBC() const
  {
    /*
    for(int i=0; i<_faces.getCount(); i++)
      applyNonzeroDiagBC(i);
    */
    throw;
  }

  void applyNonzeroDiagBC(int f) const
  {
    /*
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;    
   
      _dRdXDiag[c1][0] = T_Scalar(-1.0);   */
    throw;
  }

  
 protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const Field& _fluxField;
  const MultiField::ArrayIndex _xIndex;
  const MultiField::ArrayIndex _fluxIndex;
  XArray& _x;
  XArray& _r;
  XArray& _flux;
  //XArray& _rFlux;   // _rFlux(dynamic_cast<XArray&>(rField[_fluxIndex])),
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const bool _is2D;
  
};


template<class X>
class GenericBCSRapid : public BaseGenericBCSRapid<X>
{
 public:

  typedef BaseGenericBCSRapid<X> T_Parent;
  
 GenericBCSRapid(const StorageSite& faces,
	    const Mesh& mesh,
	    const GeomFields& geomFields,
	    Field& varField,
	    Field& fluxField,
	    MultiField& xField, MultiField& rField) :
  T_Parent(faces,mesh,geomFields,varField,fluxField,xField,rField)
    {}
  

  // see the specialization for Vectors below
  void applySymmetryBC() const
  {
    
    for(int f=0; f<this->_faces.getCount(); f++)
      {
    
        const int c0 = this->_faceCells(f,0);
        const int c1 = this->_faceCells(f,1);
        
	this->_r[c1] = this->_x[c0] - this->_x[c1];
        this->_flux[f] = X(0);
	//Anything for _rFlux?
      }
  }
  
};

#endif
