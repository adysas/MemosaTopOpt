// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWMODELVELOCITYBC_H_
#define _FLOWMODELVELOCITYBC_H_

// this file is meant to be included inside the FlowModel::Impl class
// the code is here just because FlowModel_impl.h has grown too big

  T fixedFluxContinuityBC(const StorageSite& faces,
                          const Mesh& mesh,
                          MultiFieldMatrix& matrix,
                          MultiField& xField,
                          MultiField& rField,
                          const FlowBC<T>& bc)
  {
    const StorageSite& cells = mesh.getCells();

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

#if !(defined(USING_ATYPE_RAPID))            
    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

    FMatrix& dFluxdP = dynamic_cast<FMatrix&>(matrix.getMatrix(mfIndex,pIndex));

    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
    PPDiagArray& ppDiag = ppMatrix.getDiag();
#endif

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
    TArray& xCell = dynamic_cast<TArray&>(xField[pIndex]);

    const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);

    FloatValEvaluator<VectorT3>
      bVelocity(bc.getVal("specifiedXVelocity"),
                bc.getVal("specifiedYVelocity"),
                bc.getVal("specifiedZVelocity"),
                faces);
    
    FloatValEvaluator<T>  bp(bc.getVal("specifiedPressure"),  faces);

    const bool fixPressure = (bc.bcType == "FixedBoundary");
    
    const int nFaces = faces.getCount();

    T netFlux(0.);
    
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

	//cout << "\nBefore" << endl;
	//cout << "\nrCell[" << c0 << "] = " << rCell[c0] << endl;
	//cout << "\nxCell[" << c0 << "] = " << xCell[c0] << endl;
	//cout << "\nrCell[" << c1 << "] = " << rCell[c1] << endl;
	//cout << "\nxCell[" << c1 << "] = " << xCell[c1] << endl;

        massFlux[f] = density[c0]*dot(bVelocity[f],faceArea[f]);

        rCell[c0] -= massFlux[f];
	//implies through RAPID
	//d(rCell_c0)/d(pCell_c0) -= 0; 
	//and
	//d(rCell_c0)/d(pCell_c1) -= 0; Ref. Eq. 1 recovered
        netFlux += massFlux[f];

#if !(defined(USING_ATYPE_RAPID))            
        ppAssembler.getCoeff01(f) = 0;
	//implies
	//d(rCell_c0)/d(pCell_c1) = 0; Define Ref. Eq. 1
        ppDiag[c1] = -1; 
	//implies
	//d(rCell_c1)/d(pCell_c1) = -1; Define Ref. Eq. 2
#endif

        if (fixPressure)
        {
	  //Don't care of this condition
	  //Fix needed
            
#if !(defined(USING_ATYPE_RAPID))     
	  rCell[c1] = bp[f] - pCell[c1];       
	  ppAssembler.getCoeff10(f) = 0;
#endif

#if (defined(USING_ATYPE_RAPID))     
	  //Fix needed
	  rCell[c1] = bp[f] - xCell[c1];
#endif

        }
        else
        {          

#if !(defined(USING_ATYPE_RAPID))   
	  rCell[c1] = 0.;         
          ppAssembler.getCoeff10(f) = 1;
	  //implies
	  //d(rCell_c1)/d(pCell_c0) = 1; Define Ref. Eq. 3
#endif


#if (defined(USING_ATYPE_RAPID))            
	  rCell[c1] = xCell[c0] - xCell[c1];
	  //implies through RAPID
	  //d(rCell_c1)/d(pCell_c0) = 1; Ref. Eq. 3 recovered
	  //d(rCell_c1)/d(pCell_c1) = -1; Ref. Eq. 2 recovered
	  //Question: Anything needed for flux? Something like below
	  //What is the original value for rCell[c1]?
	  //const T fluxBoundary = rCell[c1];
	  //_flux[f] = fluxBoundary;
#endif
        }

#if !(defined(USING_ATYPE_RAPID))                    
        ppMatrix.setBoundary(c1);
        dFluxdP.setCoeffL(f,T(0.));
        dFluxdP.setCoeffR(f,T(0.));
#endif   
	//cout << "\nAfter" << endl;
	//cout << "\nrCell[" << c0 << "] = " << rCell[c0] << endl;
	//cout << "\nxCell[" << c0 << "] = " << xCell[c0] << endl;
	//cout << "\nrCell[" << c1 << "] = " << rCell[c1] << endl;
	//cout << "\nxCell[" << c1 << "] = " << xCell[c1] << endl;
    }
#ifdef PV_COUPLED
    if (matrix.hasMatrix(vIndex,pIndex))
    {
        VPMatrix& vpMatrix =
          dynamic_cast<VPMatrix&>(matrix.getMatrix(vIndex,pIndex));
        
        VPAssembler& vpAssembler = vpMatrix.getPairWiseAssembler(faceCells);
        VPDiagArray& vpDiag = vpMatrix.getDiag();

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            vpDiag[c0] += vpAssembler.getCoeff01(f);
            vpAssembler.getCoeff01(f) = NumTypeTraits<T>::getZero();
        }
    }
#endif
    return netFlux;
  }

#endif
