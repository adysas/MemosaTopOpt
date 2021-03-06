// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWMODELPRESSUREBCS_H_
#define _FLOWMODELPRESSUREBCS_H_

// this file is meant to be included inside the FlowModel::Impl class
// the code is here just because FlowModel_impl.h has grown too big

void fixedPressureMomentumBC(const StorageSite& faces,
                             const Mesh& mesh,
                             MultiFieldMatrix& matrix,
                             const MultiField& xField, MultiField& rField)
{
  const StorageSite& cells = mesh.getCells();
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
  const TArray& faceAreaMag =
    dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

  VVMatrix& vvMatrix =
    dynamic_cast<VVMatrix&>(matrix.getMatrix(vIndex,vIndex));

  VVDiagArray& vvDiag = vvMatrix.getDiag();

  const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);

  const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);


  TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

  const T momURF(_options["momentumURF"]);
    
  const int nFaces = faces.getCount();

  for(int f=0; f<nFaces; f++)
    {
      if (massFlux[f] < 0.)
	{
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          const VectorT3& Af = faceArea[f];         
	  const T dpdV = -rho[c0]*mag2(V[c1])/momURF;
          vvDiag[c0][0] += dpdV*Af[0]*Af[0]/faceAreaMag[f];
          vvDiag[c0][1] += dpdV*Af[1]*Af[1]/faceAreaMag[f];
          vvDiag[c0][2] += dpdV*Af[2]*Af[2]/faceAreaMag[f];
	}
    }
}


T fixedPressureContinuityBC(const StorageSite& faces,
                            const Mesh& mesh,
                            MultiFieldMatrix& matrix,
                            const MultiField& xField, 
			    MultiField& rField, 
			    const FlowBC<T>& bc)
{
  const StorageSite& cells = mesh.getCells();
  MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
  MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
  const TArray& cellVolume =
    dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
  const VectorT3Array& cellCentroid =
    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

#if !(defined(USING_ATYPE_RAPID))            
  PPMatrix& ppMatrix =
    dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);

#endif

#if (defined(USING_ATYPE_RAPID)) 
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
  //const VectorT3Array& momAp = dynamic_cast<const VectorT3Array&>((*_momApField)[cells]);
#endif  
  
  const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
#if !(defined(USING_ATYPE_RAPID)) 
  const VectorT3Array& Vprev = dynamic_cast<const VectorT3Array&>((*_previousVelocity)[cells]);
#endif
  const TArray& p = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
  const PGradArray& pGrad = dynamic_cast<const PGradArray&>(_flowFields.pressureGradient[cells]);

  const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);

#if !(defined(USING_ATYPE_RAPID))            
  PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
  PPDiagArray& ppDiag = ppMatrix.getDiag();
#endif

  TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
  TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

#if !(defined(USING_ATYPE_RAPID))            
  FMatrix& dFluxdP = dynamic_cast<FMatrix&>(matrix.getMatrix(mfIndex,pIndex));
#endif



  FloatValEvaluator<T>  bp(bc.getVal("specifiedPressure"),  faces);
  
  const T momURF(_options["momentumURF"]);
#if !(defined(USING_ATYPE_RAPID))
  const T OneMinusmomURF(T(1.0)-momURF);
#endif

#if (defined(USING_ATYPE_RAPID))
  const T OneMinusmomURF(T(0.0));
  //const T OneMinusmomURF(T(1.0)-momURF);
#endif
    
  const int nFaces = faces.getCount();

  T netFlux(0.);
    
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);
      const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
      const VectorT3& Af = faceArea[f];
        
      const T dpf = pGrad[c0]*ds - p[c1] + p[c0] ;
      const T rhoF = rho[c0];
      
		     
      // Q < 0
      const T Q = rhoF*(Af[0]*Af[0] / momAp[c0][0]  +
                        Af[1]*Af[1] / momAp[c0][1]  +
                        Af[2]*Af[2] / momAp[c0][2])*cellVolume[c0]/(dot(Af,ds));

#if !(defined(USING_ATYPE_RAPID))      
      const T massFluxI = rhoF*(dot(V[c0],Af) - OneMinusmomURF*dot(Vprev[c0],Af)) - Q*dpf +
        OneMinusmomURF*massFlux[f];
#endif

#if (defined(USING_ATYPE_RAPID))      
      const T massFluxI = rhoF*dot(V[c0],Af) - Q*dpf;
#endif

      const VectorT3& Vb = V[c1];
      const T massFluxB = rhoF*dot(Vb,Af);
        
      netFlux += massFluxI;
        
      massFlux[f] = massFluxI;

      T Vb_dpdVb(0);
      if (massFluxB < 0)
        Vb_dpdVb = -mag2(Vb)*rhoF; 

      //dynamic pressure: Vb_dpdVb is the dynamic pressure
      //The cell pressure p[c0] or p[c1] is the static pressure 

#if !(defined(USING_ATYPE_RAPID))            
      const T denom = massFluxI-Q*Vb_dpdVb;

      if (denom != 0)
	{
	  const T dMassFluxdp0 = -Q*massFluxI/denom;
	  dFluxdP.setCoeffL(f,dMassFluxdp0);
	  dFluxdP.setCoeffR(f,T(0.));
          const T dpbdp0 = -Q*Vb_dpdVb/denom;
	  rCell[c0] -= massFlux[f];        

          // contribution to cell equation
	  ppDiag[c0] -= dMassFluxdp0;
          ppDiag[c1] = -1;
          ppAssembler.getCoeff01(f) = 0.;
          ppAssembler.getCoeff10(f) = dpbdp0;

	}
      else
	{
          // treat as fixed pressure
	  dFluxdP.setCoeffL(f,-Q);
	  dFluxdP.setCoeffR(f,T(0.));
	  ppDiag[c0] += Q;
	  ppDiag[c1] = -1;
	  rCell[c0] -= massFlux[f];
	  rCell[c1] = 0;
     	  ppAssembler.getCoeff10(f) = 0.;
	  ppAssembler.getCoeff01(f) = 0.;       
	}

#endif


#if (defined(USING_ATYPE_RAPID)) 
      if (massFluxB < 0) //flow coming in
	{
	  //rCell[c1] = boundary pressure - total cell pressure
	  //Here boundary pressure is the sum of specified pressure and the dynamic pressure
	  rCell[c0] -= massFlux[f];
	  rCell[c1] = (bp[f] + Vb_dpdVb) - p[c1];
	}
      else //flow leaving
	{
  	  //rCell[c1] = boundary pressure - total cell pressure
	  //here the boundary pressure is the specified pressure 
	  rCell[c0] -= massFlux[f];
	  rCell[c1] = bp[f]  - p[c1];

	  /*
	  cout <<"\n p[" << c1 << "] = " << p[c1] << endl;
	  cout <<"\n rCell[" << c1 << "] = " << rCell[c1] << endl;
	  cout <<"\n rCell[" << c1 << "] = " << rCell[c1]._v << endl;*/
	}

#endif

#if !(defined(USING_ATYPE_RAPID))            
      ppMatrix.setBoundary(c1);
#endif          
    }
  return netFlux;
}

void pressureBoundaryPostContinuitySolve(const StorageSite& faces,
                                         const Mesh& mesh,
                                         const T bp)
{
  const StorageSite& cells = mesh.getCells();

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
  const TArray& faceAreaMag =
    dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
  VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
  TArray& p = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
  const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);
  TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

  const int nFaces = faces.getCount();
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);
      const T rhoF = rho[c0];

      if (massFlux[f] > 0)
	{
          //const T Vn = massFlux[f]/(rhoF*faceAreaMag[f]);
          //const VectorT3 Vt = V[c0]
          //  - dot(V[c0],faceArea[f])/(faceAreaMag[f]*faceAreaMag[f])*faceArea[f];
          // V[c1] = Vn/faceAreaMag[f]*faceArea[f]+Vt ;
          V[c1] = V[c0];
          p[c1]=bp;
	}
      else
	{
          const T Vn = -massFlux[f]/(rhoF*faceAreaMag[f]);

          V[c1] = -Vn*faceArea[f]/faceAreaMag[f];
          p[c1] = bp - 0.5*rhoF*mag2(V[c1]);
	}
    }
}
// the real face pressure is calculated by bc's and stored in p[c1]
// for boundary faces but if we are also storing the pressure at all
// faces we need to copy it here
  
void updateFacePressureBoundary(const Mesh& mesh,
                                const StorageSite& faces)
{
  const StorageSite& cells = mesh.getCells();

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  const TArray& pCell = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
  TArray& pFace = dynamic_cast<TArray&>(_flowFields.pressure[faces]);

  const int nFaces = faces.getCount();
  for(int f=0; f<nFaces; f++)
    {
      //        const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);
        
      pFace[f] = pCell[c1];
    }
}

#endif
/*
	  //rCell[c1] = boundary pressure - total cell pressure
	  //rCell[c1] = rhoF*(dot(V[c0],Af)) - dpf;
	  //rCell[c1] = rhoF*(dot(V[c0],Af)) - pGrad[c0]*ds - p[c1] + p[c0];
	  //rCell[c1] = dpbdp0*p[c0] - p[c1];
	  rCell[c1] = -Q*Vb_dpdVb/denom*p[c0] - p[c1];
	  massFluxI = rhoF*(dot(V[c0],Af) - Q*dpf; 
*/
