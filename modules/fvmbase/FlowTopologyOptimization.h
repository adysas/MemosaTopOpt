/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void interpolateMassFluxApproximate()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const FlowVC<T>& vc = *_vcMap[mesh.getID()];
            
      const StorageSite& cells = mesh.getCells();
      const StorageSite& faces = mesh.getFaces();
      const int nFaces = faces.getCount();
      const CRConnectivity& faceCells = mesh.getAllFaceCells();

      const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

      const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);
      const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

      for(int f=0; f<nFaces; f++)
	{
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);
	  massFlux[f] = 0.5*(density[c0]*dot(V[c0],faceArea[f]) + density[c1]*dot(V[c1],faceArea[f]));
	}
    }


  for (int n=0; n<numMeshes; n++)
    {	
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;

	  const FlowBC<T>& bc = *_bcMap[fg.id];

	  FloatValEvaluator<VectorT3>
	    bVelocity(bc.getVal("specifiedXVelocity"),
		      bc.getVal("specifiedYVelocity"),
		      bc.getVal("specifiedZVelocity"),
		      faces);

	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	  if ((bc.bcType == "NoSlipWall") ||
	      (bc.bcType == "SlipJump") ||
	      (bc.bcType == "VelocityBoundary"))
	    {
	      // arrays for this face group
	      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
                
	      const VectorT3Array& faceArea =
		dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
                
	      for(int f=0; f<nFaces; f++)
		{
		  const int c0 = faceCells(f,0);
		  massFlux[f] = density[c0]*dot(bVelocity[f],faceArea[f]);
		}
	    }
	  //else if (bc.bcType == "SpecifiedPressure")
	  else if (bc.bcType == "PressureBoundary")
	    {
	      T bp = bc["specifiedPressure"];
	      TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
	      TArray& cellPressure = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
	      for(int f=0; f<nFaces; f++)
		{
		  const int c1 = faceCells(f,1);
		  //facePressure[f] = cellPressure[c1] = bp;
		  facePressure[f] = bp;
		}
	    }
	  else if ((bc.bcType == "Symmetry"))
	    {
	      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
	      massFlux.zero();
	    }
	}
    }
  computeContinuityResidual();
}

void interpolateMassFluxInteriorFacesSetFaceCellPressureOnPressureBoundarySetMassFluxOtherBoundaries()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const FlowVC<T>& vc = *_vcMap[mesh.getID()];
            
      const StorageSite& cells = mesh.getCells();
      const StorageSite& faces = mesh.getFaces();
      const int nFaces = faces.getCount();
      const CRConnectivity& faceCells = mesh.getAllFaceCells();

      const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

      const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);
      const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

      for(int f=0; f<nFaces; f++)
	{
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);
	  massFlux[f] = 0.5*(density[c0]*dot(V[c0],faceArea[f]) + density[c1]*dot(V[c1],faceArea[f]));
	}

	
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;

	  const FlowBC<T>& bc = *_bcMap[fg.id];

	  FloatValEvaluator<VectorT3>
	    bVelocity(bc.getVal("specifiedXVelocity"),
		      bc.getVal("specifiedYVelocity"),
		      bc.getVal("specifiedZVelocity"),
		      faces);

	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	  if ((bc.bcType == "NoSlipWall") ||
	      (bc.bcType == "SlipJump") ||
	      (bc.bcType == "VelocityBoundary"))
	    {
	      // arrays for this face group
	      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
                
	      const VectorT3Array& faceArea =
		dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
                
	      for(int f=0; f<nFaces; f++)
		{
		  const int c0 = faceCells(f,0);
		  massFlux[f] = density[c0]*dot(bVelocity[f],faceArea[f]);
		}
	    }
	  //else if (bc.bcType == "SpecifiedPressure")
	  else if (bc.bcType == "PressureBoundary")
	    {
	      T bp = bc["specifiedPressure"];
	      TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
	      TArray& cellPressure = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
	      for(int f=0; f<nFaces; f++)
		{
		  const int c1 = faceCells(f,1);
		  facePressure[f] = cellPressure[c1] = bp;
		  //facePressure[f] = bp;
		}
	    }
	  else if ((bc.bcType == "Symmetry"))
	    {
	      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
	      massFlux.zero();
	    }
	}
    }
  computeContinuityResidual();
}

void updateBetaFieldGhostCells()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      TArray& beta = dynamic_cast<TArray&>(_flowFields.beta[cells]);	

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);
	      beta[c1] = beta[c0];
	    }
	}

    }    
}

  
void initMomApField()
{
  _momApField = shared_ptr<Field>(new Field("momAp"));
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();	
      //shared_ptr<VectorT3Array> momAp(new VectorT3Array(cells.getCountLevel1()));
      shared_ptr<VVDiagArray> momAp(new VVDiagArray(cells.getCountLevel1()));
      momAp->zero();
      _momApField->addArray(cells,momAp);
    }
}

void initpreviousVelocityField()
{
  _previousVelocity = shared_ptr<Field>(new Field("previousVelocity"));
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();	
      shared_ptr<VectorT3Array> prevVel(new VectorT3Array(cells.getCountLevel1()));
      prevVel->zero();
      _previousVelocity->addArray(cells,prevVel);
    }
}

void interpolateSourceFieldSIMP()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      TArray& source = dynamic_cast<TArray&>(_flowFields.source[cells]);
      const TArray& beta = dynamic_cast<const TArray&>(_flowFields.beta[cells]);	
	
      const int nCells = cells.getCountLevel1();
		
      /*
	for ( int c=0; c<nCells; c++)
	{	    
	source[c] = 1.*(_options["alphaSolid"] + (_options["alphaFluid"] - _options["alphaSolid"])*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	}
      */
	  
      const T alphaSolid = _options["alphaSolid"];
      const T alphaFluid = _options["alphaFluid"];
      const T p = _options["p"];
      for ( int c=0; c<nCells; c++)
	{

	  //source[c] = 1.*(_options["alphaSolid"] + (_options["alphaFluid"] - _options["alphaSolid"])*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	    
	  source[c] = -1.0*((alphaSolid - alphaFluid)*pow(beta[c],p) + alphaFluid);
	    
	  //cout << c << ":" << alphaSolid << "," << alphaFluid << "," << beta[c] << "," << source[c] << endl;
	}

    }    
}

void interpolateSourceFieldSIMPRamp()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      TArray& source = dynamic_cast<TArray&>(_flowFields.source[cells]);
      const TArray& beta = dynamic_cast<const TArray&>(_flowFields.beta[cells]);	
	
      const int nCells = cells.getCountLevel1();
		
      /*
	for ( int c=0; c<nCells; c++)
	{	    
	source[c] = 1.*(_options["alphaSolid"] + (_options["alphaFluid"] - _options["alphaSolid"])*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	}
      */
	  
      const T alphaSolid = _options["alphaSolid"];
      const T alphaFluid = _options["alphaFluid"];
      const T p = _options["p"];
      for ( int c=0; c<nCells; c++)
	{

	  //source[c] = -1.*(_options["alphaSolid"] + (_options["alphaFluid"] - _options["alphaSolid"])*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	    	    
	  //cout << c << ":" << alphaSolid << "," << alphaFluid << "," << beta[c] << "," << source[c] << endl;
	  source[c] = -1.*(_options["alphaFluid"] + (_options["alphaSolid"] - _options["alphaFluid"])*(1-beta[c])/(1 + beta[c]*_options["p"]));
	}

    }    
}

void interpolateViscosityFieldSIMP()
{
  const int numMeshes = _meshes.size();
  T muSolid(0.);
  T muFluid(0.);
  T p(0.);

  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const FlowVC<T>& vc = *_vcMap[mesh.getID()];
      TArray& mu = 
	dynamic_cast<TArray&>(_flowFields.viscosity[cells]);

      const TArray& beta = 
	dynamic_cast<const TArray&>(_flowFields.beta[cells]);
	
      muSolid = vc["viscositySolid"];
      muFluid = vc["viscosityFluid"];
	
      p = _options["p"];
	
      for ( int c=0; c<nCells; c++)
	{
	  mu[c] = (muSolid - muFluid)*pow(beta[c],p) + muFluid; 
	  //cout << c << ":" << cond[c] << endl;
	}
    }
}

void initResidualFieldsFlow()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& resMomCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      TArray& resContCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);

      //shared_ptr<VectorT3Array> resMomCell(new VectorT3Array(cells.getCountLevel1()));
      resMomCell.zero();
      //_flowFields.residualMomentum.addArray(cells,resMomCell);
	
      //shared_ptr<TArray> resContCell(new TArray(cells.getCountLevel1()));
      resContCell.zero();
      //_flowFields.residualContinuity.addArray(cells,resContCell);
    }  
}

void copyVelocityPressureField(const Field& velocityDoubleField, const Field& pressureDoubleField, const MeshList& meshesDouble)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
      MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

      VectorT3Array& vCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
      const Array<VectorDouble3>& velocityCellDouble = dynamic_cast<const Array<VectorDouble3>&>(velocityDoubleField[cellsDouble]);

      TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
      const Array<double>& pressureCellDouble = dynamic_cast<const Array<double>&>(pressureDoubleField[cellsDouble]);
    
      for(int i=0; i<nCells; i++)
	{
	  vCell[i][0] = velocityCellDouble[i][0];
	  vCell[i][1] = velocityCellDouble[i][1];
	  vCell[i][2] = velocityCellDouble[i][2];
	  pCell[i] = pressureCellDouble[i];
	}
    }    
}

void copyVelocityBoundaryFace(const Field& velocityComponentDoubleField, Field& velocityComponentField, const StorageSite& facesDouble, const StorageSite& faces)
{
  //We are copying only one component of velocity
  const int nFaces = faces.getCount();
  const Array<double>& velocityComponentFaceDouble = dynamic_cast<const Array<double>&>(velocityComponentDoubleField[facesDouble]);
  TArray& velocityComponentFace = dynamic_cast<TArray&>(velocityComponentField[faces]);
  for(int f=0; f<nFaces; f++)
    {
      velocityComponentFace[f] = velocityComponentFaceDouble[f];
    }
}

void copyMomApField(const Field& momApDoubleField, const MeshList& meshesDouble)
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);
      const VVDoubleDiagArray& momApDouble = dynamic_cast<const VVDoubleDiagArray&>(momApDoubleField[cellsDouble]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  momAp[i][0] = momApDouble[i][0];
	  momAp[i][1] = momApDouble[i][1];
	  momAp[i][2] = momApDouble[i][2];
	}
    }    

  //_previousVelocity = dynamic_pointer_cast<Field>(_flowFields.velocity.newCopy());
}

void copyMomentumResidual(const Field& momentumRes, const MeshList& meshesFrom)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      const Mesh& meshFrom = *meshesFrom[n];
      const StorageSite& cellsFrom = meshFrom.getCells();	

      const int nCells = cells.getCountLevel1();

      VectorT3Array& momResTo = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      const VectorT3Array& momResFrom = dynamic_cast<const VectorT3Array&>(momentumRes[cellsFrom]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  momResTo[i][0] = momResFrom[i][0];
	  momResTo[i][1] = momResFrom[i][1];
	  momResTo[i][2] = momResFrom[i][2];
	}
    }
}

void copyContinuityResidual(const Field& continuityRes, const MeshList& meshesFrom)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      const Mesh& meshFrom = *meshesFrom[n];
      const StorageSite& cellsFrom = meshFrom.getCells();	
		
      const int nCells = cells.getCountLevel1();

      TArray& contResTo = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
      const TArray& contResFrom = dynamic_cast<const TArray&>(continuityRes[cellsFrom]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  contResTo[i] = contResFrom[i];
	}
    }
}

void copypreviousVelocityField(const Field& previousVelocityDoubleField, const MeshList& meshesDouble)
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      VectorT3Array& prevVel = dynamic_cast<VectorT3Array&>((*_previousVelocity)[cells]);
      const VectorDouble3Array& prevVelDouble = dynamic_cast<const VectorDouble3Array&>(previousVelocityDoubleField[cellsDouble]);
	    
      for(int i=0; i<nCells; i++)
	{
	  prevVel[i][0] = prevVelDouble[i][0];
	  prevVel[i][1] = prevVelDouble[i][1];
	  prevVel[i][2] = prevVelDouble[i][2];
	}
    }    

}

void copypressureGradientField(const Field& pressureGradientDoubleField, const MeshList& meshesDouble)
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      PGradArray& pGradCell = dynamic_cast<PGradArray&>(_flowFields.pressureGradient[cells]);
	
      const PDoubleGradArray& pGradCellDouble = dynamic_cast<const PDoubleGradArray&>(pressureGradientDoubleField[cellsDouble]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  pGradCell[i][0] = pGradCellDouble[i][0];
	  pGradCell[i][1] = pGradCellDouble[i][1];
	  pGradCell[i][2] = pGradCellDouble[i][2];
	  //pGradCell[i][0] = pressureGradientDoubleField[i][0];
	  //pGradCell[i][1] = pressureGradientDoubleField[i][1];
	  //pGradCell[i][2] = pressureGradientDoubleField[i][2];
	}
    }    
}

void copyfacePressureField(const Field& facePressureDoubleField, const MeshList& meshesDouble)
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& facesDouble = meshDouble.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
      const Array<double>& facePressureDouble = dynamic_cast<const Array<double>&>(facePressureDoubleField[facesDouble]);
	
      for(int f=0; f<nFaces; f++)
	{
	  facePressure[f] = facePressureDouble[f];
	}

    } 
}

void copyMassFluxField(const Field& massFluxFaceDoubleField, const MeshList& meshesDouble)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& facesDouble = meshDouble.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& massFluxFace = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
      const Array<double>& massFluxFaceDouble = dynamic_cast<const Array<double>&>(massFluxFaceDoubleField[facesDouble]);
	
      for(int f=0; f<nFaces; f++)
	{
	  massFluxFace[f] = massFluxFaceDouble[f];
	}

    } 
}
 
void copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
      const Array<double>& betaCellDouble = dynamic_cast<const Array<double>&>(betaDoubleField[cellsDouble]);
	
      for(int i=0; i<nCells; i++)
	{
	  betaCell[i] = betaCellDouble[i];
	}

    }    
}

void copySourceField(const Field& sourceDoubleField, const MeshList& meshesDouble)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      TArray& sourceCell = dynamic_cast<TArray&>(_flowFields.source[cells]);
      const Array<double>& sourceCellDouble = dynamic_cast<const Array<double>&>(sourceDoubleField[cellsDouble]);
	
      for(int i=0; i<nCells; i++)
	{
	  sourceCell[i] = sourceCellDouble[i];
	}

    }    
}

void ComputeMomApCoefficients()
{
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const StorageSite& cells = mesh.getCells();
      const int nCellsSelfCount = cells.getSelfCount();
      const StorageSite& faces = mesh.getFaces();
      const int nFaces = faces.getCount();
      const CRConnectivity& faceCells = mesh.getFaceCells(faces);
      const VectorT3Array& cellCentroid =
	dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
      const TArray& cellVolume =
	dynamic_cast<const TArray&>(_geomFields.volume[cells]);
      const VectorT3Array& faceArea =
	dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);  
      const TArray& faceAreaMag =
	dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

      const TArray& diffCell =
	dynamic_cast<const TArray&>(getViscosityField()[cells]);

      const TArray& massFluxFace =
	dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);

      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);

      const TArray& continuityResidual =
	dynamic_cast<const TArray&>(_flowFields.continuityResidual[cells]);

      const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
      _previousVelocity = dynamic_pointer_cast<Field>(_flowFields.velocity.newCopy());

      const TArray& source = dynamic_cast<const TArray&>(_flowFields.source[cells]);
      const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);
      /////////////////////////Diffusion and convection Discretization /////////////
      for(int f=0; f<nFaces; f++)
	{
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);

	  //Diffusion discretization
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
	  momAp[c0][0] -= diffCoeff;
	  momAp[c0][1] -= diffCoeff;
	  momAp[c0][2] -= diffCoeff;
	    
	  momAp[c1][0] -= diffCoeff;
	  momAp[c1][1] -= diffCoeff;
	  momAp[c1][2] -= diffCoeff;
	    
	  //Convection discretization
	  const T_Scalar faceCFlux = massFluxFace[f];
	    
	  if (faceCFlux > T_Scalar(0))
	    {
	      momAp[c0][0] -= faceCFlux;
	      momAp[c0][1] -= faceCFlux;
	      momAp[c0][2] -= faceCFlux;
	    }
	  else
	    {
	      momAp[c1][0] += faceCFlux;
	      momAp[c1][1] += faceCFlux;
	      momAp[c1][2] += faceCFlux;
	    }
	}

      //Handle source linearization of Brinkman
	  
      for(int c=0;c<nCellsSelfCount;c++)
	{
	  const T_Scalar cImb = continuityResidual[c];
	  momAp[c][0] += cImb;
	  momAp[c][1] += cImb;
	  momAp[c][2] += cImb;
	}


      for(int c=0; c<nCellsSelfCount; c++)
	{
	  momAp[c][0] += cellVolume[c]*source[c];
	  momAp[c][1] += cellVolume[c]*source[c];
	  momAp[c][2] += cellVolume[c]*source[c];
	}
    }


  // pressure boundary

  //I haven't handled any boundary conditions except pressure boundary


  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);
      const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);

      const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);
      //Handle pressure boundary
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const T momURF(_options["momentumURF"]);
	  const FlowBC<T>& bc = *_bcMap[fg.id];

	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	  const VectorT3Array& faceArea =
	    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
	  const TArray& faceAreaMag =
	    dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

	  const TArray& massFluxFace =
	    dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);

	  if (bc.bcType == "PressureBoundary")
	    {
	      for(int f=0; f<nFaces; f++)
		{
		  if (massFluxFace[f] < 0.)
		    {
		      const int c0 = faceCells(f,0);
		      const int c1 = faceCells(f,1);
		      const VectorT3& Af = faceArea[f];
		      const T dpdV = -rho[c0]*mag2(V[c1])/momURF;
		      momAp[c0][0] += dpdV*Af[0]*Af[0]/faceAreaMag[f];
		      momAp[c0][1] += dpdV*Af[1]*Af[1]/faceAreaMag[f];
		      momAp[c0][2] += dpdV*Af[2]*Af[2]/faceAreaMag[f];
		    }
		}
	    }	    
	}
    }

  // Under-relaxation: before or after?
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);	
      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);
      //Under-relaxing
      const T_Scalar momURF(_options["momentumURF"]);
      const int nCells = cells.getSelfCount();   
      for(int c=0; c<nCells; c++)
	{
	  momAp[c][0] /= momURF;
	  momAp[c][1] /= momURF;
	  momAp[c][2] /= momURF;
	}
    }

  //Copying the mom ap coefficients of c0 to c1 cells
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);	
      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);
	
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);
	      momAp[c1] = momAp[c0];
	    }
	}
    }

  _momApField->syncLocal();
}

void updateFacePressure()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const StorageSite& cells = mesh.getCells();
      const StorageSite& iFaces = mesh.getInteriorFaceGroup().site;
        
      updateFacePressureInterior(mesh,iFaces);

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const FlowBC<T>& bc = *_bcMap[fg.id];
	    
	  updateFacePressureBoundary(mesh,faces);
	}

    }
}

void updatePressureGradient()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getSelfCount();

      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      VVDiagArray& momAp = dynamic_cast<VVDiagArray&>((*_momApField)[cells]);

      PGradArray& pGradCell = dynamic_cast<PGradArray&>(_flowFields.pressureGradient[cells]);
	
      const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
	
      const StorageSite& faces = mesh.getFaces();
      const CRConnectivity& faceCells = mesh.getAllFaceCells();
      const VectorT3Array& faceArea =
	dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
      const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
      pGradCell.zero();

      const int nFaces = faces.getCount();
      for(int f=0; f<nFaces; f++)
	{
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);
	  pGradCell[c0].accumulate(faceArea[f],facePressure[f]);
	  pGradCell[c1].accumulate(faceArea[f],-facePressure[f]);
	}

      for(int c=0; c<nCells; c++)
	pGradCell[c] /= cellVolume[c];

      // copy boundary values from adjacent cells
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const int faceCount = faces.getCount();

	  if (fg.groupType == "symmetry")
	    {
	      const VectorT3Array& faceArea =
		dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	      const TArray& faceAreaMag =
		dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	      for(int f=0; f<faceCount; f++)
		{
		  const int c0 = faceCells(f,0);
		  const int c1 = faceCells(f,1);
		  const VectorT3 en = faceArea[f]/faceAreaMag[f];
		  reflectGradient(pGradCell[c1], pGradCell[c0], en);
		}
	    }
	  else
	    {
	      for(int f=0; f<faceCount; f++)
		{
		  const int c0 = faceCells(f,0);
		  const int c1 = faceCells(f,1);
		    
		  pGradCell[c1] = pGradCell[c0];
		}
	    }	      	       
	}    
    }
}

void ComputeContinuityResidual()
//MultiField& ComputeContinuityResidual() 
{
  //Continuity residual
  LinearSystem lsContinuity;
  initContinuityLinearization(lsContinuity);
  lsContinuity.initAssembly();
  linearizeContinuity(lsContinuity);
  _previousVelocity = shared_ptr<Field>();
  _momApField = shared_ptr<Field>();
  //postContinuitySolve(lsContinuity);

  //Copying into the field residualContinuity
  //Ask Dr. Mathur whether there is an easier way to copy from a multifield to a field
  MultiField& resContinuitymf = lsContinuity.getB();

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
      TArray& resContinuitymfCell = dynamic_cast<TArray&>(resContinuitymf[pIndex]);
      TArray& residualContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
      for(int i=0; i<nCells; i++)
	{
	  residualContinuityCell[i] = resContinuitymfCell[i];
#if !( defined(USING_ATYPE_RAPID) )
	  //cout << "\nCont. Res["<<i<<"]= " << residualContinuityCell[i];
#endif
#if ( defined(USING_ATYPE_RAPID) )
	  //cout << "\nCont. Res["<<i<<"]= " << residualContinuityCell[i];
#endif
	}
    }
  //return resContinuitymf;
}

void ComputeMomentumResidual()
{
  //Momentum residual
  LinearSystem lsMomentum;
  initMomentumLinearization(lsMomentum);
  lsMomentum.initAssembly();
  linearizeMomentum(lsMomentum);
  //_momApField = shared_ptr<Field>();
    
  MultiField& resMomentummf = lsMomentum.getB();

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
      VectorT3Array& resMomentummfCell = dynamic_cast<VectorT3Array&>(resMomentummf[vIndex]);
      VectorT3Array& residualMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
	
      for(int i=0; i<nCells; i++)
	{
	  residualMomentumCell[i] = resMomentummfCell[i];
	  //cout << "\nU Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][0] << endl;
	  //cout << "\nV Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][1] << endl;
	  //cout << "\nW Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][2] << endl;
	}
    }
}

void ComputeMomentumResidualSaveMomAp()
{
  //Momentum residual
  LinearSystem lsMomentum;
  initMomentumLinearization(lsMomentum);
  lsMomentum.initAssembly();
  linearizeMomentum(lsMomentum);
  //_momApField = shared_ptr<Field>();
    
  MultiField& resMomentummf = lsMomentum.getB();

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
      VectorT3Array& resMomentummfCell = dynamic_cast<VectorT3Array&>(resMomentummf[vIndex]);
      VectorT3Array& residualMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
	
      for(int i=0; i<nCells; i++)
	{
	  residualMomentumCell[i] = resMomentummfCell[i];
	  //cout << "\nU Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][0] << endl;
	  //cout << "\nV Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][1] << endl;
	  //cout << "\nW Mom. Res. ["<<i<<"]= " << residualMomentumCell[i][2] << endl;
	}
    }

    
  // save current velocity for use in continuity discretization
  _previousVelocity = dynamic_pointer_cast<Field>(_flowFields.velocity.newCopy());

  // save the momentum ap coeffficients for use in continuity discretization
  _momApField = shared_ptr<Field>(new Field("momAp"));
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
      const VVMatrix& vvMatrix =
	dynamic_cast<const VVMatrix&>(lsMomentum.getMatrix().getMatrix(vIndex,vIndex));
      const VVDiagArray& momAp = vvMatrix.getDiag();
      _momApField->addArray(cells,dynamic_pointer_cast<ArrayBase>(momAp.newCopy()));
    }
  _momApField->syncLocal();
    
}

void setObjectiveFunctionFlow(const VectorT3 objFunction, const int component)
{
  _objFunctionFlowModel = objFunction[component];
}

void setObjectiveFunctionFlow(const T objFunction)
{
  _objFunctionFlowModel = objFunction;
}

#if ( defined(USING_ATYPE_RAPID) )
 
void setIndependentStateVariablesFlow()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      VectorT3Array& vCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
      TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
      const int nCells = cells.getCountLevel1(); 
      int index = 0;
      for(int c=0; c<nCells; c++)
	{
	  index = 4*c;
	  vCell[c][0].setIndex(index);
	  vCell[c][1].setIndex(index+1);
	  vCell[c][2].setIndex(index+2);
	  pCell[c].setIndex(index+3);
	}
    }        
}

void setIndependentDesignVariablesFlow()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  betaCell[c].setIndex(4*nCells+c);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}

void setViscosityAsIndependentVariable()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();  
      const int index = 4*nCells;
      TArray& viscosityCell = dynamic_cast<TArray&>(_flowFields.viscosity[cells]);
      for (int c=0; c<nCells; c++)
	{
	  viscosityCell[c].setIndex(index);
	}
      _nDesignVariables = 1;
    }
}


   
void setSourceAsIndependentVariable()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();  
      const int index = 4*nCells;
      TArray& sourceCell = dynamic_cast<TArray&>(_flowFields.source[cells]);
      for (int c=0; c<nCells; c++)
	{
	  sourceCell[c].setIndex(index);
	}
      _nDesignVariables = 1;
    }
}

void setSourceAsIndependentVariable(const int cellId)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();  
      const int index = 4*nCells;
      TArray& sourceCell = dynamic_cast<TArray&>(_flowFields.source[cells]);
      sourceCell[cellId].setIndex(index);
      _nDesignVariables = 1;
    }
}

void setSourceFieldAsIndependentVariables()
{
  const int numMeshes = _meshes.size();
	
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& sourceCell = dynamic_cast<TArray&>(_flowFields.source[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  sourceCell[c].setIndex(4*nCells+c);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}

void setDensityAsIndependentVariable()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();  
      const int index = 4*nCells;
      TArray& densityCell = dynamic_cast<TArray&>(_flowFields.density[cells]);
      for (int c=0; c<nCells; c++)
	{
	  densityCell[c].setIndex(index);
	}
      _nDesignVariables = 1;
    }
}


void WriteResidualGradientsToPetscFlow()
{
    
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& rMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      TArray& rContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
	
      const int nCells = cells.getCountLevel1();
      //const int nCellsWithoutGhosts = cells.getSelfCount();

      int nCoeffUVWP = 0;
      int nCoeffBeta = 0;

      //Let the total number of cells be N
      //Let the total number of interior cells be M
      //N = M + Number of boundary faces

      //There are 4 state variables U, V, W, P; Therefore 4N state variables
      //There are M design variables

      //There are 4N residuals
      //dRdUVWP is a matrix of size 4N X 4N

      //nnzElementsRowdRdUVWP is a double array that stores the number of non-zero values in each row. Or the number of independent variables on which the residual of the cell depends upon. 
      double nnzElementsRowdRdUVWP[4*nCells];
      memset(nnzElementsRowdRdUVWP, 0, (4*nCells)*sizeof(double));
	
      //dRdBeta is a matrix of size 4N X M
      double nnzElementsRowdRdBeta[4*nCells];
      memset(nnzElementsRowdRdBeta, 0, (4*nCells)*sizeof(double));

      cout << endl;
      for(int i=0; i<nCells; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  //First check whether the derivatives are with respect to X momentum (or U)
		  if (ii.first < 4*nCells)
		    {
		      nnzElementsRowdRdUVWP[4*i+j]++;
		      nCoeffUVWP++;
		    }
		  else
		    {			  
		      nnzElementsRowdRdBeta[4*i+j]++;
		      nCoeffBeta++;
		    }
		}
	    }
	  //cout << "P-Res = " << rContinuityCell[i] << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 4*nCells)
		{
		  nnzElementsRowdRdUVWP[4*i+3]++;
		  nCoeffUVWP++;
		}
	      else
		{
		  nnzElementsRowdRdBeta[4*i+3]++;
		  nCoeffBeta++;
		}
	    }
	}
      
      //Write in binary format so that Petsc can read directly
      string dRdUVWPFileNameBinary = "dRdUVWP.binary";
      FILE *dRdUVWPFileBinary = fopen(dRdUVWPFileNameBinary.c_str(),"wb"); 
      string dRdBetaFileNameBinary = "dRdBeta.binary";
      FILE *dRdBetaFileBinary = fopen(dRdBetaFileNameBinary.c_str(),"wb"); 

      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;

      //some id number for petsc
      bigEndianInteger = 1211216;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      //Number of rows for dRdUVWP = 4N
      bigEndianInteger = 4*nCells;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	  

      //Number of columns for dRdUVWP = 4N
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);

	  
      //Number of rows for dRdBeta = 4N
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      //Number of columns for dRdBeta = M
      //bigEndianInteger = nCellsWithoutGhosts;
      bigEndianInteger = _nDesignVariables;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      // Number of non-zero elements in dRdUVWP
      bigEndianInteger = nCoeffUVWP;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);

      // Number of non-zero elements in dRdBeta
      bigEndianInteger = nCoeffBeta;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);



      // Number of non-zeros in each row in both dRdUVWP and dRdBeta

      for(int i=0; i<4*nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdUVWP[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	}

      for(int i=0; i<4*nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdBeta[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
	}



      // Column indices of all non-zeros elements in both dRdUVWP and dRdBeta (starting index is zero)
	
      for(int i=0; i<nCells; i++)
	{
	  for (int j=0; j<3; j++)
	    {
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  if (ii.first < 4*nCells)
		    {
		      bigEndianInteger = ii.first;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		    }
		  else
		    {
		      bigEndianInteger = ii.first-4*nCells;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		    }
		}

	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 4*nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-4*nCells;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		}
	    }
	}

	  
      // Values of all non-zeros for both dRdUVWP and dRdBeta
      for(int i=0; i<nCells; i++)
	{
	  for (int j=0; j<3; j++)
	    {
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  bigEndianDouble = ii.second;
		  LowEndian2BigEndian(bigEndianDouble);
		  if (ii.first < 4*nCells)
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
		  else
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 4*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }
	    
	}

      fclose(dRdUVWPFileBinary);
      fclose(dRdBetaFileBinary);
    }
}

void WriteResidualGradientsToMatlabFlow()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& rMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      TArray& rContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
	
      const int nCells = cells.getCountLevel1();
      //const int nCellsWithoutGhosts = cells.getSelfCount();

      string SV = "DRDSV.mat";
      string DV = "DRDDV.mat";
      FILE *DRDSVFile = fopen(SV.c_str(),"wb");
      FILE *DRDDVFile = fopen(DV.c_str(),"wb");
	
      cout << endl;
      for(int i=0; i<nCells; i++)
	{
	  for (int j=0; j<3; j++)
	    {
	      //cout << rMomentumCell[i][j]._v << endl;
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  if (ii.first < 4*nCells)	       
		    fprintf(DRDSVFile, "%d %d %22.15le \n", 4*i+j+1, ii.first+1, ii.second); 
		  else
		    fprintf(DRDDVFile, "%d %d %22.15le \n", 4*i+j+1, ii.first-4*nCells+1, ii.second); 
		}
	    }
	  //cout << rContinuityCell[i]._v << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 4*nCells)
		{
		  fprintf(DRDSVFile, "%d %d %22.15le \n", 4*i+3+1, ii.first+1, ii.second); 
		}
	      else
		{
		  fprintf(DRDDVFile, "%d %d %22.15le \n", 4*i+3+1, ii.first-4*nCells+1, ii.second);
		}
	    }	    
	}
      fclose(DRDSVFile);
      fclose(DRDDVFile);
    }
}

void WriteObjectiveFunctionGradientToPetscFlow()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();


      double temporaryArray[4*nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (4*nCells+_nDesignVariables)*sizeof(double));
    
      double scaleFlowModel = _options["scaleFlowModelObjFunction"]._v;            
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionFlowModel._dv)
	{
	  temporaryArray[ii.first] = scaleFlowModel*ii.second;
	}

      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing dcdPhi
      string dcdUVWPFileNameBinary = "dcdUVWP.binary";
      FILE *dcdUVWPFileBinary = fopen(dcdUVWPFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdUVWPFileBinary);
                
      bigEndianInteger = 4*nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdUVWPFileBinary);
      
      for(int i=0; i<4*nCells; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdUVWPFileBinary);
	}
      fclose(dcdUVWPFileBinary);

      //Writing dcdBeta
      string dcdBetaFileNameBinary = "dcdBeta.binary";
      FILE *dcdBetaFileBinary = fopen(dcdBetaFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdBetaFileBinary);
      
      bigEndianInteger = _nDesignVariables;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdBetaFileBinary);

      for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdBetaFileBinary);
	}
      fclose(dcdBetaFileBinary);
    }
}

void WritePressureDropGradientToPetsc(const int faceId, const int direction)
{

  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();

      VectorT3 pressureDrop = getPressureIntegral(mesh, faceId);
      T pressureDropComponent = pressureDrop[direction];


      double temporaryArray[4*nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (4*nCells+_nDesignVariables)*sizeof(double));
    
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, pressureDropComponent._dv)
	{
	  temporaryArray[ii.first] = ii.second;
	}

      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing dcdPhi
      string dPDdUVWPFileNameBinary = "dPDdUVWP.binary";
      FILE *dPDdUVWPFileBinary = fopen(dPDdUVWPFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdUVWPFileBinary);
                
      bigEndianInteger = 4*nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdUVWPFileBinary);
      
      for(int i=0; i<4*nCells; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dPDdUVWPFileBinary);
	}
      fclose(dPDdUVWPFileBinary);

      //Writing dcdBeta
      string dPDdBetaFileNameBinary = "dPDdBeta.binary";
      FILE *dPDdBetaFileBinary = fopen(dPDdBetaFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdBetaFileBinary);
      
      bigEndianInteger = _nDesignVariables;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdBetaFileBinary);

      for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dPDdBetaFileBinary);
	}
      fclose(dPDdBetaFileBinary);
    }
}


void WriteObjectiveFunctionGradientToMatlabFlow()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      
      double temporaryArray[4*nCells+_nDesignVariables];
      memset(temporaryArray, 0, (4*nCells+_nDesignVariables)*sizeof(double));    
       
      double scaleFlowModel = _options["scaleFlowModelObjFunction"]._v;   
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionFlowModel._dv)
	{
	  temporaryArray[ii.first] = scaleFlowModel*ii.second;
	}

      string SV = "DcDSV.mat";
      string DV = "DcDDV.mat";
      FILE *DcDSVFile = fopen(SV.c_str(),"wb");
      FILE *DcDDVFile = fopen(DV.c_str(),"wb");

      for(int i=0; i<4*nCells; i++)
	{
	  fprintf(DcDSVFile, "%22.15le \n", temporaryArray[i]); 
	}
      fclose(DcDSVFile);

      //for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
      for(int i=4*nCells; i<4*nCells+_nDesignVariables; i++)
	{
	  fprintf(DcDDVFile, "%22.15le \n", temporaryArray[i]); 
	}
      fclose(DcDDVFile);
    }
}
#endif


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID) )

void dumpBetaSensistivityFields()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
	
      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
      TArray& sensCell = dynamic_cast<TArray&>(_flowFields.sensitivity[cells]);
      const int nCellsWithoutGhosts = cells.getSelfCount();

	
      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing Sensitivity
      string sensitivityFileNameBinary = "Sensitivity.bin";
      FILE *sensitivityFileBinary = fopen(sensitivityFileNameBinary.c_str(),"wb"); 

      string betaFileNameBinary = "Beta.bin";
      FILE *betaFileBinary = fopen(betaFileNameBinary.c_str(),"wb"); 

      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, sensitivityFileBinary);
      fwrite(&bigEndianInteger, sizeof(int), 1, betaFileBinary);

      bigEndianInteger = nCellsWithoutGhosts;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, sensitivityFileBinary);
      fwrite(&bigEndianInteger, sizeof(int), 1, betaFileBinary);

      for(int i=0; i<nCellsWithoutGhosts; i++)
	{
	  bigEndianDouble = sensCell[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, sensitivityFileBinary);

	  bigEndianDouble = betaCell[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, betaFileBinary);
	  //cout << "\nKui: " << betaCell[i] << "; " << sensCell[i] << endl;
	}
      fclose(sensitivityFileBinary);
      fclose(betaFileBinary);
	
    }
}


void dumpVelocityVector(const string fileBase)
{

  const Mesh& mesh = *_meshes[0];
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  VectorT3Array& velocityCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);

  string velFileName = fileBase + ".vel";
  FILE *velFile = fopen(velFileName.c_str(),"wb");
    
  for(int i=0; i<nCells; i++)
    fprintf(velFile, "%22.15le \t %22.15le \n", velocityCell[i][0], velocityCell[i][1]);

  fclose(velFile);
}


void ComputeVelocityGradient()
{
  _velocityGradientModel.compute();    
}

#endif



// Functions for developing the code

void printPressureIntegrals()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
            
	  VectorT3 r(VectorT3::getZero());
            
	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea =
	    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const int nFaces = faces.getCount();
	  const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
	  for(int f=0; f<nFaces; f++)
	    r += faceArea[f]*facePressure[f];

	  cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
	}
    }
}
  
void printMomentumFluxIntegrals()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
            
	  VectorT3 r(VectorT3::getZero());
            
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const VectorT3Array& momFlux =
	    dynamic_cast<const VectorT3Array&>(_flowFields.momentumFlux[faces]);
	  for(int f=0; f<nFaces; f++)
	    r += momFlux[f];

	  cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
	}
    }
}
  
void printMassFluxIntegrals()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
            
	  T r(0.);
            
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
	  for(int f=0; f<nFaces; f++)
	    r += massFlux[f];

	  cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
	}
    }
}



void printCellPressure(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const TArray& pressure = dynamic_cast<const TArray&>((_flowFields.pressure)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nCell pressure \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \n", n+1, pressure[n]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \n", n+1, pressure[n]._v);
#endif	
	}
    }   
}


void printCellVelocity(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>((_flowFields.velocity)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nCell velocity \n";
      for(int n=0; n<nCells; n++)
	{	
#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, velocity[n][0], velocity[n][1], velocity[n][2]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, velocity[n][0]._v, velocity[n][1]._v, velocity[n][2]._v);
#endif	
	}
    }   
}


void printMomApField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nMom Ap Coefficient \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << n << " Mom Ap Coefficient (Double)" << momAp[n][0] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  //cout << n << " Mom Ap Coefficient (Rapid) " << momAp[n][0]._v << endl;
	  cout << n << " Mom Ap Coefficient (Rapid) " << momAp[n][0] << endl;
#endif	
	}
    } 
    
}


void printMomApField(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nMom Ap Coefficient \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, momAp[n][0], momAp[n][1], momAp[n][2]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, momAp[n][0]._v, momAp[n][1]._v, momAp[n][2]._v);
#endif	
	}
    } 
    
}


void printpreviousVelocityField(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      const VectorT3Array& previousVelocity = dynamic_cast<const VectorT3Array&>((*_previousVelocity)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nPrevious Velocity \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, previousVelocity[n][0], previousVelocity[n][1], previousVelocity[n][2]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, previousVelocity[n][0]._v, previousVelocity[n][1]._v, previousVelocity[n][2]._v);
#endif	
	}
    } 
    
}

void printpressureGradientField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      //VectorT3Array& momAp = dynamic_cast<VectorT3Array&>((*_momApField)[cells]);
      const PGradArray& pgradCell = dynamic_cast<const PGradArray&>((_flowFields.pressureGradient)[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nPressure Gradient  \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << n << " Pressure Gradient (Double) " << pgradCell[n][0] << ", " << pgradCell[n][1] << ", " << pgradCell[n][2] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  cout << n << " Pressure Gradient (Rapid) " << pgradCell[n][0]._v << ", " << pgradCell[n][1]._v << ", " << pgradCell[n][2]._v << endl;
#endif	
	}
    } 
    
}

void printpressureGradientField(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const PGradArray& pgradCell = dynamic_cast<const PGradArray&>((_flowFields.pressureGradient)[cells]);
      const int nCells = cells.getCount();  

      cout << "\nPressure Gradient  \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, pgradCell[n][0], pgradCell[n][1], pgradCell[n][2]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, pgradCell[n][0]._v, pgradCell[n][1]._v, pgradCell[n][2]._v);
#endif	
	}

    } 
}



void printfacePressureField()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
	
      for(int f=0; f<nFaces; f++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << f << " face Pressure (Double) " << facePressure[f] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  cout << f << " face Pressure (Rapid) " << facePressure[f]._v << endl;
#endif		    
	}
    } 
}


void printfacePressureField(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
	
      for(int f=0; f<nFaces; f++)
	{

#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \n", f+1, facePressure[f]); 
#endif	

#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \n", f+1, facePressure[f]._v); 
#endif		    
	}
    } 
}


void printMassFlux()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
	
      for(int f=0; f<nFaces; f++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << f << " Mass flux (Double) " << massFlux[f] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  cout << f << " Mass flux (Rapid) " << massFlux[f]._v << endl;
#endif	
	    
	}

    } 
}

void printMassFlux(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
	
      for(int f=0; f<nFaces; f++)
	{

#if !( defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_TANGENT) )
	  fprintf(file, "%d \t %22.15le \n", f+1, massFlux[f]); 
#endif	

#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \n", f+1, massFlux[f]._v); 
#endif		    
	}
    } 
}

void printMomentumResidual()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex vIndex(&_flowFields.velocity, &cells);

      //TArray& resContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
      VectorT3Array& resMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\Momentum residual \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << n << "th cell Mom Residual = " << resMomentumCell[n][0] <<", " << resMomentumCell[n][1] << ", " << resMomentumCell[n][2] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  //cout << n << "th cell Residual = " << resMomentumCell[n][0]._v << endl;
	  cout << n << "th cell Mom Residual = " << resMomentumCell[n][0]._v <<", " << resMomentumCell[n][1]._v << ", " << resMomentumCell[n][2]._v << endl;
	  /*
	    cout << n << "th cell Mom Residual = " << resMomentumCell[n][0] << endl;
	    cout << n << "th cell Mom Residual = " << resMomentumCell[n][1] << endl;
	    cout << n << "th cell Mom Residual = " << resMomentumCell[n][2] << endl;
	  */
	  /*
	    foreach(const typename Rapid::PartialDerivatives::value_type& ii, resMomentumCell[n][0]._dv)
	    {
	    if (ii.first == n)
	    cout << n << "th Ap Coeff = " << ii.second << endl;
	    }*/
#endif
	    

	}	
    }    
}


void printMomentumResidual(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex vIndex(&_flowFields.velocity, &cells);
      VectorT3Array& resMomentumCell = dynamic_cast<VectorT3Array&>(_flowFields.residualMomentum[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nMomentum residual \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || (defined(USING_ATYPE_TANGENT)) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, resMomentumCell[n][0], resMomentumCell[n][1], resMomentumCell[n][2]);
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \t %22.15le \t %22.15le \n", n+1, resMomentumCell[n][0]._v, resMomentumCell[n][1]._v, resMomentumCell[n][2]._v);
#endif	    
	}	
    }    
}

void printContinuityResidual()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex vIndex(&_flowFields.velocity, &cells);

      TArray& resContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nContinuity residual \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << n << "th cell Con Residual(double) = " << resContinuityCell[n] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  //cout << n << "th cell Residual = " << resMomentumCell[n][0]._v << endl;
	  //cout << n << "th cell Con Residual(Rapid) = " << resContinuityCell[n]._v << endl;
	  cout << n << "th cell Con Residual(Rapid) = " << resContinuityCell[n] << endl;
	  /*
	    foreach(const typename Rapid::PartialDerivatives::value_type& ii, resMomentumCell[n][0]._dv)
	    {
	    if (ii.first == n)
	    cout << n << "th Ap Coeff = " << ii.second << endl;
	    }*/
#endif
	    

	}	
    }    
}

	
void printContinuityResidual(string fileName)
{
  const int numMeshes = _meshes.size();
  FILE *file = fopen(fileName.c_str(),"wb"); 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex vIndex(&_flowFields.velocity, &cells);

      TArray& resContinuityCell = dynamic_cast<TArray&>(_flowFields.residualContinuity[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\nContinuity residual \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) || (defined(USING_ATYPE_TANGENT)) )
	  fprintf(file, "%d \t %22.15le \n", n+1, resContinuityCell[n]); 
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  fprintf(file, "%d \t %22.15le \n", n+1, resContinuityCell[n]._v); 
#endif	    
	}	
    }    
}


void printSourceField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      TArray& sourceCell = dynamic_cast<TArray&>(_flowFields.source[cells]);
      const int nCells = cells.getCount();  
     
      cout << "\n Source field \n";
      for(int n=0; n<nCells; n++)
	{
#if !( defined(USING_ATYPE_RAPID) )
	  cout << n << "th cell source(double) = " << sourceCell[n] << endl;
#endif	
#if ( defined(USING_ATYPE_RAPID) )
	  //cout << n << "th cell Residual = " << resMomentumCell[n][0]._v << endl;
	  //cout << n << "th cell Con Residual(Rapid) = " << resContinuityCell[n]._v << endl;
	  cout << n << "th cell source(Rapid) = " << sourceCell[n] << endl;
#endif	  
	}	
    }
}

#if ( defined(USING_ATYPE_RAPID) )

void printIndependentStateVariables()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      VectorT3Array& vCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
      const int nCells = cells.getCount();   
      cout << endl;
      for(int c=0; c<nCells; c++)
	{
	  cout << "U[" << c << "]" << vCell[c][0] << endl;
	  cout << "V[" << c << "]" << vCell[c][1] << endl;
	  cout << "W[" << c << "]" << vCell[c][2] << endl;
	}
    }
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
      const int nCells = cells.getCount();  
      for(int c=0; c<nCells; c++)
	{
	  cout << "P[" << c << "]" << pCell[c] << endl;
	}
    }
}
  
void printIndependentDesignVariables()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
      const int nCellsWithoutGhosts = cells.getSelfCount();
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  cout << "beta[" << c << "]" << betaCell[c] << endl;
	}
    }
}
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
void buildNeighborsList(const T rMin)
{    
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      _filter.buildNeighborListWeight(mesh, rMin);
    }    
}

void applyFilter()
{
    
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      _filter.applyFilter(mesh, _flowFields.sensitivity, _flowFields.beta);
    }
}
#endif

////////////////////////////////////////////////////////////////////////////////////
#if ( defined(USING_ATYPE_RAPID) )
void setIndependentStateVariablesThermoFluidLaminar()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      VectorT3Array& vCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
      TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
      const int nCells = cells.getCountLevel1(); 
      int index = 0;
      for(int c=0; c<nCells; c++)
	{
	  index = 5*c;
	  vCell[c][0].setIndex(index);
	  vCell[c][1].setIndex(index+1);
	  vCell[c][2].setIndex(index+2);
	  pCell[c].setIndex(index+3);
	}
    }        
}

void setIndependentDesignVariablesThermoFluidLaminar()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  betaCell[c].setIndex(5*nCells+c);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}

void setIndependentStateVariablesTurbulentFlow()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      //u,v,w=0,1,2
      //p=3
      //nuTilde=4
      //phi=5 
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      VectorT3Array& vCell = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
      TArray& pCell = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
      const int nCells = cells.getCountLevel1(); 
      int index = 0;
      for(int c=0; c<nCells; c++)
	{
	  index = 6*c;
	  vCell[c][0].setIndex(index);
	  vCell[c][1].setIndex(index+1);
	  vCell[c][2].setIndex(index+2);
	  pCell[c].setIndex(index+3);
	}
    }        
}

void setIndependentDesignVariablesTurbulentFlow()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& betaCell = dynamic_cast<TArray&>(_flowFields.beta[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  betaCell[c].setIndex(6*nCells+c);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}


#endif


void copyEddyViscosityField(const Field& eddyViscosityField)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();

      TArray& eddyViscTo = dynamic_cast<TArray&>(_flowFields.eddyviscosity[cells]);
      const TArray& eddyViscFrom = dynamic_cast<const TArray&> (eddyViscosityField[cells]);

      for(int i=0; i<nCells; i++)
	{
	  eddyViscTo[i] = eddyViscFrom[i];
	}

    }    
}


VectorT3 getPressureGradientIntegral()
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const PGradArray& pGrad = dynamic_cast<const PGradArray&>(_flowFields.pressureGradient[cells]);
      const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
      const TArray& beta = dynamic_cast<const TArray&>(_flowFields.beta[cells]);

      //Here I am weighting the pressure gradient with the volume fraction of the liquid; the intend is to capture only the pressure drop in the fluid. 
      /*
	for(int i=0; i<nCells; i++)
	{
	r[0] += (T(1.0) - beta[i])*pGrad[i][0]*cellVolume[i];
	r[1] += (T(1.0) - beta[i])*pGrad[i][1]*cellVolume[i];
	r[2] += (T(1.0) - beta[i])*pGrad[i][2]*cellVolume[i];
	} 
      */

      
#if ( defined(USING_ATYPE_RAPID) )
      for(int i=0; i<nCells; i++)
	{
	  r[0] += (T(1.0) - beta[i]._v)*pGrad[i][0]*cellVolume[i];
	  r[1] += (T(1.0) - beta[i]._v)*pGrad[i][1]*cellVolume[i];
	  r[2] += (T(1.0) - beta[i]._v)*pGrad[i][2]*cellVolume[i];
	} 
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID) )
      for(int i=0; i<nCells; i++)
	{
	  r[0] += (1.0 - beta[i])*pGrad[i][0]*cellVolume[i];
	  r[1] += (1.0 - beta[i])*pGrad[i][1]*cellVolume[i];
	  r[2] += (1.0 - beta[i])*pGrad[i][2]*cellVolume[i];
	} 
#endif
      
    }
      
  return r;
}

/*
VectorT3 getDynamicPressureIntegralInlet(const Mesh& mesh, const int inlet)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const StorageSite& cells = mesh.getCells();

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == inlet)
	{

	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const int nFaces = faces.getCount();
	  const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
	  const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
	  const FlowBC<T>& bc = *_bcMap[fg.id];
	  FloatValEvaluator<VectorT3> bVelocity(bc.getVal("specifiedXVelocity"),
						bc.getVal("specifiedYVelocity"),
						bc.getVal("specifiedZVelocity"),
						faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const T rhoF = rho[c0];
	      const T dp = -0.5*mag2(bVelocity[f])*rhoF; 
	      r += dp*faceArea[f];
	      //r += 0.5*massFlux[f]*bVelocity[f];
	    }

	  found=true;
	}
    }
  if (!found)
    throw CException("getPressureIntegral: invalid faceGroupID");
  return r;
}


VectorT3 getDynamicPressureIntegralOutlet(const Mesh& mesh, const int outlet)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const StorageSite& cells = mesh.getCells();

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == outlet)
	{
	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
	  const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
	  const int nFaces = faces.getCount();
	  const FlowBC<T>& bc = *_bcMap[fg.id];
	  FloatValEvaluator<T> bpressure(bc.getVal("specifiedPressure"),faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);
	      const T rhoF = rho[c0];
	      const T dp = -0.5*mag2(velocity[c1])*rhoF; 
	      r += dp*faceArea[f];
	      //r += 0.5*massFlux[f]*velocity[c0];
	    }

	  found=true;
	}
    }
  if (!found)
    throw CException("getPressureIntegral: invalid faceGroupID");
  return r;
}

*/

VectorT3 getDynamicPressureIntegral(const Mesh& mesh, const int faceId)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const StorageSite& cells = mesh.getCells();

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceId)
	{
	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
	  const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
	  const int nFaces = faces.getCount();
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);
	      const T rhoF = rho[c0];
	      const T dp = -0.5*mag2(velocity[c1])*rhoF; 
	      r += dp*faceArea[f];
	    }

	  found=true;
	}
    }
  if (!found)
    throw CException("getPressureIntegral: invalid faceGroupID");
  return r;
}

VectorT3 getTotalPressureLossesBetweenInletOutlet(const Mesh& mesh, const int inlet, const int outlet)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const StorageSite& cells = mesh.getCells();

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == inlet)
	{
	  r += (getPressureIntegral(mesh, inlet) + getDynamicPressureIntegral(mesh, inlet));
	  found=true;
	}

      if (fg.id == outlet)
	{
	  r -= (getPressureIntegral(mesh, outlet) + getDynamicPressureIntegral(mesh, outlet));
	  found=true;
	}
    }
  if (!found)
    throw CException("getPressureIntegral: invalid faceGroupID");
  return r;
}

/*
VectorT3 getTotalPressureLossesBetweenInletOutlet(const Mesh& mesh, const int inlet, const int outlet)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;

  const StorageSite& cells = mesh.getCells();

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == inlet)
	{
	  r += (getPressureIntegral(mesh, inlet) + getDynamicPressureIntegral(mesh, inlet));
	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const int nFaces = faces.getCount();
	  const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
	  const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
	  const FlowBC<T>& bc = *_bcMap[fg.id];
	  FloatValEvaluator<VectorT3> bVelocity(bc.getVal("specifiedXVelocity"),
						bc.getVal("specifiedYVelocity"),
						bc.getVal("specifiedZVelocity"),
						faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      r += (faceArea[f]*facePressure[f] + 0.5*massFlux[f]*bVelocity[f]);
	    }

	  found=true;
	}

      if (fg.id == outlet)
	{
	  const StorageSite& faces = fg.site;
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
	  const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
	  const int nFaces = faces.getCount();
	  const FlowBC<T>& bc = *_bcMap[fg.id];
	  FloatValEvaluator<T> bpressure(bc.getVal("specifiedPressure"),faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      r -= (faceArea[f]*bpressure[f] + 0.5*massFlux[f]*velocity[c0]);
	    }

	  found=true;
	}
    }
  if (!found)
    throw CException("getPressureIntegral: invalid faceGroupID");
  return r;
}
*/
