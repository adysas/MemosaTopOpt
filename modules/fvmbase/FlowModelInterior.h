// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWMODELINTERIOR_H_
#define _FLOWMODELINTERIOR_H_

T discretizeMassFluxInterior(const Mesh& mesh,
                             const StorageSite& faces,
                             MultiFieldMatrix& mfmatrix,
                             const MultiField& xField, MultiField& rField,
                             const bool isSymmetry=false)
{
  const StorageSite& cells = mesh.getCells();
  const int nCellsInterior = cells.getSelfCount();
  
  MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
  const TArray& cellVolume =
    dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
  const TArray& faceAreaMag =
    dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    
  const VectorT3Array& cellCentroid =
    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  const Array<int>& row = faceCells.getRow();
  const Array<int>& col = faceCells.getCol();

#if !(defined(USING_ATYPE_RAPID))
  PPMatrix& ppMatrix =
    dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
#endif

#if (defined(USING_ATYPE_RAPID))
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
  //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
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
  TArray& ppCoeff = ppMatrix.getOffDiag();
#endif

  TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
  TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

  const T momURF(_options["momentumURF"]);
#if !(defined(USING_ATYPE_RAPID))
  const T OneMinusmomURF(T(1.0)-momURF);
#endif

#if (defined(USING_ATYPE_RAPID))
  const T OneMinusmomURF(T(0.0));
  //const T OneMinusmomURF(T(1.0)-momURF);
#endif
    
  const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
  const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);

  const StorageSite& ibFaces = mesh.getIBFaces();
    
  const VectorT3Array* ibVelocity = (ibFaces.getCount() > 0) ?
    &(dynamic_cast<const VectorT3Array&>(_flowFields.velocity[mesh.getIBFaces()])) : 0;

  // the net flux from ib faces
  T boundaryFlux(0);

  //temp variable; delete this later
  T temp(0);
      
  const int nFaces = faces.getCount();
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
      const VectorT3& Af = faceArea[f];
        
      const T diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(Af,ds);

      const T momApBar0 = (momAp[c0][0]+momAp[c0][1]+momAp[c0][2])/3.0;
      const T momApBar1 = isSymmetry ? momApBar0 : (momAp[c1][0]+momAp[c1][1]+momAp[c1][2])/3.0 ;
      const T momApBarFace = momApBar0 + momApBar1;

#if !(defined(USING_ATYPE_RAPID))
      const T VdotA0 = dot(V[c0],Af) - OneMinusmomURF*dot(Vprev[c0],Af);
      const T VdotA1 = dot(V[c1],Af) - OneMinusmomURF*dot(Vprev[c1],Af);
#endif

#if (defined(USING_ATYPE_RAPID))
      const T VdotA0 = dot(V[c0],Af);
      const T VdotA1 = dot(V[c1],Af);
#endif
        
      const T dpf = cellVolume[c0]*(pGrad[c0]*ds) + cellVolume[c1]*(pGrad[c1]*ds);
      // const T dpf = (cellVolume[c0]*(pGrad[c0]*ds) + cellVolume[c1]*(pGrad[c1]*ds))/
      //(cellVolume[c0]+cellVolume[c1]) - p[c1] + p[c0];

      const T Vn = (VdotA0*momApBar0 + VdotA1*momApBar1 -dpf*diffMetric) / momApBarFace;
      //const T Vn = (VdotA0*momApBar0 + VdotA1*momApBar1 -
      //      (cellVolume[c0]+cellVolume[c1])*dpf*diffMetric) / momApBarFace;

      const T rhoF = 0.5*(rho[c0]+rho[c1]);
      const T aByMomAp = Af[0]*Af[0] / (momAp[c0][0] + momAp[c1][0]) +
        Af[1]*Af[1] / (momAp[c0][1] + momAp[c1][1]) +
        Af[2]*Af[2] / (momAp[c0][2] + momAp[c1][2]);


      const T pCoeff = rhoF*aByMomAp*(cellVolume[c0]+cellVolume[c1])/(dot(Af,ds));

      if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
          (ibType[c1] == Mesh::IBTYPE_FLUID))
	{
#if !(defined(USING_ATYPE_RAPID)) 
          massFlux[f] = rhoF*Vn - pCoeff*(p[c0]-p[c1]) + (1-momURF)*massFlux[f];
#endif

#if (defined(USING_ATYPE_RAPID)) 
          massFlux[f] = rhoF*Vn - pCoeff*(p[c0]-p[c1]);
	  //massFlux[f] = rhoF*Vn - pCoeff*(p[c0]-p[c1]) + (1-momURF)*massFlux[f];
#endif
	  /*
	  if ((c0 == 4) || (c1 == 4))
	    {
	      cout << "\nCheck\n";
	      cout << "\nMassFlux [4] = " << massFlux[4] << endl; 
	      cout << "\nrCell [4] = " << rCell[4] << endl; 
	      cout << "\nc0,c1 = " << c0 << c1 << endl;
	      }*/
	  /*
#if !(defined(USING_ATYPE_RAPID)) 
	  if ((c0 == 6) || (c1 == 6))
	    {	      
	      T T1 = rhoF*Vn;
	      cout << "\nT1 = " << T1 << "\t";
	      T T2 = pCoeff*(p[c0]-p[c1]);
	      cout << "T2 = " << T2 << "\t";
	      T T3 = (1-momURF)*massFlux[f];
	      cout << "T3 = " << T3 << "\t";
	      cout << "Sum = " << massFlux[f] << endl;
	      if (c0 == 6)
		temp = temp - massFlux[f];
	      else
		temp = temp + massFlux[f];
	      cout << "\nMassflux " << temp << endl;
	    }
	  
	    
#endif

#if (defined(USING_ATYPE_RAPID)) 
	  if ((c0 == 6) || (c1 == 6))
	    {	      
	      T T1 = rhoF*Vn;
	      cout << "\nT1 = " << T1._v << "\t";
	      T T2 = pCoeff*(p[c0]-p[c1]);
	      cout << "T2 = " << T2._v << "\t";
	      T T3 = (1-momURF)*massFlux[f];
	      cout << "T3 = " << T3._v << "\t";
	      cout << "Sum = " << massFlux[f]._v << endl;
	      if (c0 == 6)
		temp = temp - massFlux[f];
	      else
		temp = temp + massFlux[f];
	      cout << "\nMassflux " << temp._v << endl;
	    }
	  
	    
#endif
	  */
          if (isSymmetry) massFlux[f] = 0.;
          
          //massFlux[f] = rhoF*Vn  + OneMinusmomURF*massFlux[f];
            
          rCell[c0] -= massFlux[f];
          rCell[c1] += massFlux[f];
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
	  //preDumpPressureLinearSystem(mfmatrix, rField, "FMI");
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
#endif
#if !(defined(USING_ATYPE_RAPID))            
          ppAssembler.getCoeff01(f) -=pCoeff;
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
	  //preDumpPressureLinearSystem(mfmatrix, rField, "FMI");
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;

          ppAssembler.getCoeff10(f) -=pCoeff;
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
	  //preDumpPressureLinearSystem(mfmatrix, rField, "FMI");
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
            
          ppDiag[c0] += pCoeff;
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
	  //preDumpPressureLinearSystem(mfmatrix, rField, "FMI");
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
          ppDiag[c1] += pCoeff;
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
	  //preDumpPressureLinearSystem(mfmatrix, rField, "FMI");
	  //cout << "----------------- Post Pressure ----------------------------------"<<endl;
#endif
          if (isSymmetry)
	    {
#if !(defined(USING_ATYPE_RAPID))            
              ppDiag[c0] -= ppAssembler.getCoeff01(f);
              ppAssembler.getCoeff01(f) =0;
#endif
              rCell[c1] = 0;
#if !(defined(USING_ATYPE_RAPID))
              ppMatrix.setBoundary(c1);
#endif
	    }
	}
      else if (((ibType[c0] == Mesh::IBTYPE_FLUID)
                && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
               ((ibType[c1] == Mesh::IBTYPE_FLUID)
                && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	{
          const int ibFace = ibFaceIndex[f];
          if (ibFace < 0)
            throw CException("invalid ib face index");
          
          const VectorT3& faceVelocity = (*ibVelocity)[ibFace];

          // this is an iBFace, determine which cell is interior and
          // which boundary. Treat as a fixed flux boundary. 
          if (ibType[c0] == Mesh::IBTYPE_FLUID)
	    {
              massFlux[f]= rho[c0]*dot(Af,faceVelocity);
              rCell[c0] -= massFlux[f];
              rCell[c1] = 0;
#if !(defined(USING_ATYPE_RAPID))
              ppMatrix.setDirichlet(c1);
#endif
              boundaryFlux += massFlux[f];
	    }
          else
	    {
              rCell[c0] = 0;
#if !(defined(USING_ATYPE_RAPID))
              ppMatrix.setDirichlet(c0);
#endif

              // skip c1 if this is a ghost cell since the other mesh will
              // compute the right face flux
              if (c1 < nCellsInterior)
		{
                  massFlux[f]= rho[c1]*dot(Af,faceVelocity);
                  rCell[c1] += massFlux[f];
                  boundaryFlux -= massFlux[f];
		}
              else
		{
                  rCell[c1] = 0;
#if !(defined(USING_ATYPE_RAPID))
                  ppAssembler.getCoeff10(f) =0;
                  ppDiag[c1] = -1;
#endif
		}
	    }
                
	}
      else 
	{
          if ((ibType[c0] == Mesh::IBTYPE_FLUID) ||
              (ibType[c1] == Mesh::IBTYPE_FLUID))
            throw CException("invalid face to skip");
            
          // setup to get zero corrections
          massFlux[f]=0;
#if !(defined(USING_ATYPE_RAPID))
          ppDiag[c0] = -1;
          ppDiag[c1] = -1;
          ppMatrix.setDirichlet(c0);
          ppMatrix.setDirichlet(c1);
#endif
	}
    }

#ifdef PV_COUPLED
  if (mfmatrix.hasMatrix(pIndex,vIndex))
    {
      PVMatrix& pvMatrix =
        dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
      PVAssembler& pvAssembler = pvMatrix.getPairWiseAssembler(faceCells);
      PVDiagArray& pvDiag = pvMatrix.getDiag();

      for(int f=0; f<nFaces; f++)
	{
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          const VectorT3& Af = faceArea[f];
            
          const T momApBar0 = (momAp[c0][0]+momAp[c0][1]+momAp[c0][2])/3.0;
          const T momApBar1 = (momAp[c1][0]+momAp[c1][1]+momAp[c1][2])/3.0;
          const T momApBarFace = momApBar0 + momApBar1;
          const T rhoF = 0.5*(rho[c0]+rho[c1]);

          VectorT3T coeff0(rhoF*momApBar0/momApBarFace*Af);
          VectorT3T coeff1(rhoF*momApBar1/momApBarFace*Af);
            
        
          pvAssembler.getCoeff01(f) -=coeff1;
          pvAssembler.getCoeff10(f) +=coeff0;
            
          pvDiag[c0] -= coeff0;
          pvDiag[c1] += coeff1;
	}
    }
#endif
  /*
  cout << "\nMassFlux [4] = " << massFlux[4] << endl; 
  cout << "\nrCell [4] = " << rCell[4] << endl; */
  return boundaryFlux;
}


void correctVelocityInterior(const Mesh& mesh,
                             const StorageSite& faces,
                             const MultiField& ppField,
                             const bool isSymmetry = false)                               
{
  const StorageSite& cells = mesh.getCells();
  const int nCellsInterior = cells.getSelfCount();

  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
  const VectorT3Array& cellCentroid =  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
  VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
  const TArray& pp = dynamic_cast<const TArray&>(ppField[pIndex]);
  const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
  const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

  const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
  const int nFaces = faces.getCount();
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
          (ibType[c1] == Mesh::IBTYPE_FLUID))
	{
          const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
          const VectorT3& Af = faceArea[f];
            
          const T aByMomAp0 = Af[0]*Af[0] / momAp[c0][0] +
            Af[1]*Af[1] / momAp[c0][1] +
            Af[2]*Af[2] / momAp[c0][2];
            
          const T aByMomAp1 = Af[0]*Af[0] / momAp[c1][0] +
            Af[1]*Af[1] / momAp[c1][1] +
            Af[2]*Af[2] / momAp[c1][2];
            
          const T Adotes = dot(Af,ds)/mag(ds);
          const T coeff0  = cellVolume[c0]*rho[c0]*aByMomAp0/Adotes;
          const T coeff1  = isSymmetry ? coeff0 :
            cellVolume[c1]*rho[c1]*aByMomAp1/Adotes;
            
          const T ppFace = (coeff0*pp[c0]+coeff1*pp[c1])/(coeff0+coeff1);
          const VectorT3 ppA = ppFace*faceArea[f];
            
          V[c0] += ppA/momAp[c0];
          if (!isSymmetry)
            V[c1] -= ppA/momAp[c1];
	}
      else if (((ibType[c0] == Mesh::IBTYPE_FLUID)
                && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
               ((ibType[c1] == Mesh::IBTYPE_FLUID)
                && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	{
          // this is an iBFace, determine which cell is interior and
          // which boundary. Correct the interior cell's velocity
          // using the cell pressure correction as the face pressure
          // correction
          if (ibType[c0] == Mesh::IBTYPE_FLUID)
	    {
              const T ppFace = pp[c0];
              const VectorT3 ppA = ppFace*faceArea[f];
                
              V[c0] += ppA/momAp[c0];
	    }
          else if (c1 < nCellsInterior)
	    {
              const T ppFace = pp[c1];
              const VectorT3 ppA = ppFace*faceArea[f];
                
              V[c1] -= ppA/momAp[c1];
	    }
	}
      // nothing needs to be done for the solid/solid or solid/ib faces
    }
}

void updateFacePressureInterior(const Mesh& mesh,
                                const StorageSite& faces,
                                const bool isSymmetry=false)
{
  const StorageSite& cells = mesh.getCells();

  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
  const VectorT3Array& cellCentroid =  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

#if !(defined(USING_ATYPE_RAPID))
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
#endif
#if (defined(USING_ATYPE_RAPID))
  //const VectorT3Array& momAp = dynamic_cast<const VectorT3Array&>((*_momApField)[cells]);
  const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
#endif    
  const TArray& pCell = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
  TArray& pFace = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
  const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
  const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

  const int nFaces = faces.getCount();
  const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
          (ibType[c1] == Mesh::IBTYPE_FLUID))
	{
          const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
          const VectorT3& Af = faceArea[f];
            
          const T aByMomAp0 = Af[0]*Af[0] / momAp[c0][0] +
            Af[1]*Af[1] / momAp[c0][1] +
            Af[2]*Af[2] / momAp[c0][2];
            
          const T aByMomAp1 = Af[0]*Af[0] / momAp[c1][0] + 
            Af[1]*Af[1] / momAp[c1][1] +
            Af[2]*Af[2] / momAp[c1][2];
            
          const T Adotes = dot(Af,ds)/mag(ds);
          const T coeff0  = cellVolume[c0]*rho[c0]*aByMomAp0/Adotes;
          const T coeff1  = isSymmetry ? coeff0 :
            cellVolume[c1]*rho[c1]*aByMomAp1/Adotes;
            
          pFace[f] = (coeff0*pCell[c0]+coeff1*pCell[c1])/(coeff0+coeff1);
	}
      else if (((ibType[c0] == Mesh::IBTYPE_FLUID)
                && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
               ((ibType[c1] == Mesh::IBTYPE_FLUID)
                && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	{
          // this is an iBFace, determine which cell is interior and
          // which boundary. copy pressure from the fluid cell
          if (ibType[c0] == Mesh::IBTYPE_FLUID)
	    {
              pFace[f] = pCell[c0];
	    }
          else
	    {
              pFace[f] = pCell[c1];
	    }
	}
      else
        // for solid/solid and solid/ib faces pressure is never used
        pFace[f]=0;
    }
}

void correctMassFluxInterior(const Mesh& mesh,
                             const StorageSite& faces,
                             MultiFieldMatrix& mfmatrix,
                             const MultiField& xField)
{
  const StorageSite& cells = mesh.getCells();
  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
  const TArray& pp = dynamic_cast<const TArray&>(xField[pIndex]);

  PPMatrix& ppMatrix =
    dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));
    
  PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);

  TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    
  const int nFaces = faces.getCount();
  for(int f=0; f<nFaces; f++)
    {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      // should work for ib and ib/solid etc faces as well since the
      // coefficients at such faces are all zero
      massFlux[f] -= ppAssembler.getCoeff01(f)*pp[c1] -
        ppAssembler.getCoeff10(f)*pp[c0];
    }

#ifdef PV_COUPLED
  MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
  if (mfmatrix.hasMatrix(pIndex,vIndex))
    {
      PVMatrix& pvMatrix =
        dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
      PVAssembler& pvAssembler = pvMatrix.getPairWiseAssembler(faceCells);
      const VectorT3Array& Vp = dynamic_cast<const VectorT3Array&>(xField[vIndex]);

      for(int f=0; f<nFaces; f++)
	{
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
            
          massFlux[f] += pvAssembler.getCoeff01(f)*Vp[c1] +
            pvAssembler.getCoeff10(f)*Vp[c0];
	}
    }
#endif
}

#endif
