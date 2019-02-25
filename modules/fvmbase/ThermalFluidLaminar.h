#ifndef _THERMALFLUIDLAMINAR_H_
#define _THERMALFLUIDLAMINAR_H_

void initResidualField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();    	
      TArray& resTempCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);	
      resTempCell.zero();
    }  
}

void copyTemperatureField(const Field& tempDoubleField, const MeshList& meshesDouble)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const Mesh& meshDouble = *meshesDouble[n];
      const StorageSite& cellsDouble = meshDouble.getCells();
	
      const int nCells = cells.getCountLevel1();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& tCell = dynamic_cast<TArray&>(_thermalFields.temperature[cells]);
      const Array<double>& tempCellDouble = dynamic_cast<const Array<double>&>(tempDoubleField[cellsDouble]);

      for(int i=0; i<nCells; i++)
	{
	  tCell[i] = tempCellDouble[i];
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

      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const Array<double>& betaCellDouble = dynamic_cast<const Array<double>&>(betaDoubleField[cellsDouble]);
	
      for(int i=0; i<nCells; i++)
	{
	  betaCell[i] = betaCellDouble[i];
	}

    }    
}


void interpolateThermalConductivityFieldSIMP()
{
  const int numMeshes = _meshes.size();
  T kHigh(0.);
  T kLow(0.);
  T p(0.);

  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const ThermalVC<T>& vc = *_vcMap[mesh.getID()];
      TArray& cond = 
	dynamic_cast<TArray&>(_thermalFields.conductivity[cells]);

      const TArray& beta = 
	dynamic_cast<const TArray&>(_thermalFields.beta[cells]);
	
      kHigh = vc["thermalConductivityHigh"];
      kLow = vc["thermalConductivityLow"];
	
      p = _options["p"];
	
      //cout << endl << "Conductivity" << endl; 
      for ( int c=0; c<nCells; c++)
	{
	  cond[c] = (kHigh - kLow)*pow(beta[c],p) + kLow; 
	  //cout << c << ":" << cond[c] << endl;
	}
    }
}

void updateBetaFieldGhostCells()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
	

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

void updateSourceFieldType1()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& source = dynamic_cast<TArray&>(_thermalFields.source[cells]);

      const int nCells = cells.getCount();
		
      for ( int c=0; c<nCells; c++)
	{
	  source[c] = _options["alpha"];
	}

    }    
}


void updateSourceFieldType2()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& source = dynamic_cast<TArray&>(_thermalFields.source[cells]);
      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const int nCells = cells.getCount();
		
      for ( int c=0; c<nCells; c++)
	{
	  source[c] = _options["alpha"]*beta[c];
	}

    }    
}

void updateSourceFieldType3()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& source = dynamic_cast<TArray&>(_thermalFields.source[cells]);
      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const int nCells = cells.getCount();
		
      for ( int c=0; c<nCells; c++)
	{
	  source[c] = _options["alpha"]*(1. - beta[c]);
	}

    }    
}

void ComputeTemperatureResidual()
{
  //Temperature residual
  LinearSystem ls;
  initLinearization(ls);        
  ls.initAssembly();
  linearize(ls);
  MultiFieldMatrix& matrix = ls.getMatrix();
  MultiField& resTemperaturemf = ls.getB();

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& resTemperatureField = dynamic_cast<TArray&>(resTemperaturemf[tIndex]);
      TArray& residualTemperatureCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);
      for(int i=0; i<nCells; i++)
	{
	  residualTemperatureCell[i] = resTemperatureField[i];
	}
    }
}

void printResidualField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& rCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);
      const int nCells = cells.getCount();    
      cout << "\n";
      for(int n=0; n<nCells; n++)
	{
	  cout << n << "th cell Residual" << rCell[n] << endl;
	}
    }    
}


void copyMomentumResidual(const Field& momentumRes)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
	
      const int nCells = cells.getCountLevel1();

      VectorT3Array& momentumResidualThermalModel = dynamic_cast<VectorT3Array&>(_thermalFields.momentumResidual[cells]);
      const VectorT3Array& momentumResidualFlowModel = dynamic_cast<const VectorT3Array&>(momentumRes[cells]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  momentumResidualThermalModel[i][0] = momentumResidualFlowModel[i][0];
	  momentumResidualThermalModel[i][1] = momentumResidualFlowModel[i][1];
	  momentumResidualThermalModel[i][2] = momentumResidualFlowModel[i][2];
	}
    }
}

void copyContinuityResidual(const Field& continuityRes)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
			
      const int nCells = cells.getCountLevel1();

      TArray& continuityResidualThermalModel = dynamic_cast<TArray&>(_thermalFields.continuityResidual[cells]);
      const TArray& continuityResidualFlowModel = dynamic_cast<const TArray&>(continuityRes[cells]);
	
    
      for(int i=0; i<nCells; i++)
	{
	  continuityResidualThermalModel[i] = continuityResidualFlowModel[i];
	}
    }
}

void copyMassFlux(const Field& massFluxFromFlowModel)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& faces = mesh.getFaces();
	
      const int nFaces = faces.getCount();

      TArray& convectionFlux = dynamic_cast<TArray&>(_thermalFields.convectionFlux[faces]);
      const TArray& massFlux = dynamic_cast<const TArray&>(massFluxFromFlowModel[faces]);
	
    
      for(int f=0; f<nFaces; f++)
	{
	  convectionFlux[f] = massFlux[f];
	}
    }
}


void copyPressureIntegral(const VectorT3 pressureIntegral, int direction)
{
  _objFunctionFlowModel = pressureIntegral[direction];
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  cout << "\nPressureIntegral double " << _objFunctionFlowModel << endl;
#endif
#if ( defined(USING_ATYPE_RAPID) )
  cout << "\nPressureIntegral rapid " << _objFunctionFlowModel._v << endl;
#endif
}

double getDeltaEnthalpyIntegralRAPIDThermoFluid(const Mesh& mesh, const int faceGroupId1, const int faceGroupId2)
{
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  const int nCellsWithoutGhosts = cells.getSelfCount();
  T r1(0.0);
  T r2(0.0);
  const TArray& temperature = dynamic_cast<const TArray&>(_thermalFields.temperature[cells]);
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId1)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& convectionFlux = dynamic_cast<const TArray&>(_thermalFields.convectionFlux[faces]);	  
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c1 = faceCells(f,1);
	      r1 += temperature[c1]*convectionFlux[f];
	    }
	}
      if (fg.id == faceGroupId2)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& convectionFlux = dynamic_cast<const TArray&>(_thermalFields.convectionFlux[faces]);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c1 = faceCells(f,1);
	      r2 += temperature[c1]*convectionFlux[f];
	    }
	}      
    }
  _objFunctionThermalModel = r2 - r1;

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  return _objFunctionThermalModel;
#endif
#if defined(USING_ATYPE_RAPID)
  return _objFunctionThermalModel._v;
#endif
}


#if ( defined(USING_ATYPE_RAPID) )

void setIndependentStateVariablesThermal()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      TArray& tCell = dynamic_cast<TArray&>(_thermalFields.temperature[cells]);
      const int nCells = cells.getCountLevel1(); 
      int index = 0;
      for(int c=0; c<nCells; c++)
	{
	  tCell[c].setIndex(c);
	}
    }        
}

void setIndependentDesignVariablesThermal()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  betaCell[c].setIndex(c + nCells);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}

void WriteResidualGradientsToPetscThermal()
{
  for ( int id = 0; i < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& rMomentumCell = dynamic_cast<VectorT3Array&>(_thermalFields.momentumResidual[cells]);
      TArray& rContinuityCell = dynamic_cast<TArray&>(_thermalFields.continuityResidual[cells]);
      TArray& rTemperatureCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);
	
      const int nCells = cells.getCountLevel1();
      //const int nCellsWithoutGhosts = cells.getSelfCount();

      int nCoeffUVWP = 0;
      int nCoeffBeta = 0;

      //Let the total number of cells be N
      //Let the total number of interior cells be M
      //N = M + Number of boundary faces

      //There are 5 state variables U, V, W, P, T; Therefore 5N state variables
      //There are M design variables

      //There are 5N residuals
      //dRdUVWP is a matrix of size 5N X 5N

      //nnzElementsRowdRdUVWP is a double array that stores the number of non-zero values in each row. Or the number of independent variables on which the residual of the cell depends upon. 
      double nnzElementsRowdRdUVWP[5*nCells];
      memset(nnzElementsRowdRdUVWP, 0, (5*nCells)*sizeof(double));
	
      //dRdBeta is a matrix of size 5N X M
      double nnzElementsRowdRdBeta[5*nCells];
      memset(nnzElementsRowdRdBeta, 0, (5*nCells)*sizeof(double));

      cout << endl;
      for(int i=0; i<nCells; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
		  
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  //First check whether the derivatives are with respect to X momentum (or U)
		  if (ii.first < 5*nCells)
		    {
		      nnzElementsRowdRdUVWP[5*i+j]++;
		      nCoeffUVWP++;
		    }
		  else
		    {			  
		      nnzElementsRowdRdBeta[5*i+j]++;
		      nCoeffBeta++;
		    }
		}
	    }
	  //cout << "P-Res = " << rContinuityCell[i] << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  nnzElementsRowdRdUVWP[5*i+3]++;
		  nCoeffUVWP++;
		}
	      else
		{
		  nnzElementsRowdRdBeta[5*i+3]++;
		  nCoeffBeta++;
		}
	    }
	  //Adding temperature residual
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  nnzElementsRowdRdUVWP[5*i+4]++;
		  nCoeffUVWP++;
		}
	      else
		{
		  nnzElementsRowdRdBeta[5*i+4]++;
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

      //Number of rows for dRdUVWP = 5N
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);	  

      //Number of columns for dRdUVWP = 5N
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	  
      //Number of rows for dRdBeta = 5N
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
      for(int i=0; i<5*nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdUVWP[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	}
      for(int i=0; i<5*nCells; i++)
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
		  if (ii.first < 5*nCells)
		    {
		      bigEndianInteger = ii.first;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		    }
		  else
		    {
		      bigEndianInteger = ii.first-5*nCells;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		    }
		}

	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-5*nCells;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-5*nCells;
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
		  if (ii.first < 5*nCells)
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
		  else
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }	    
	}

      fclose(dRdUVWPFileBinary);
      fclose(dRdBetaFileBinary);
    }
}

void setIndependentStateVariablesThermoFluid()
{
  const int numMeshes = _meshes.size();
    
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      TArray& tCell = dynamic_cast<TArray&>(_thermalFields.temperature[cells]);
      const int nCells = cells.getCountLevel1(); 
      int index = 0;
      for(int c=0; c<nCells; c++)
	{
	  index = 5*c;
	  tCell[c].setIndex(index+4);
	}
    }        
}

void setIndependentDesignVariablesThermoFluid()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount(); 
      const int nCellsWithoutGhosts = cells.getSelfCount();
      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
	
      for(int c=0; c<nCellsWithoutGhosts; c++)
	{
	  betaCell[c].setIndex(5*nCells+c);
	}
      _nDesignVariables = nCellsWithoutGhosts;
    }
}


void getHeatFluxIntegralRAPIDThermoFluid(const Mesh& mesh, const int faceGroupId)
{
  bool found = false;
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  const int nCellsWithoutGhosts = cells.getSelfCount();
  T r(0.0);

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const TArray& heatFlux =
	    dynamic_cast<const TArray&>(_thermalFields.heatFlux[faces]);
	  for(int f=0; f<nFaces; f++)
	    {
	      r += heatFlux[f];
	    }
	  found=true;
	}
    }
  _objFunctionThermalModel = r;
}


double WriteObjectiveFunctionGradientToPetsc()
{
  double objFunction(0.);
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();


      double temporaryArray[5*nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (5*nCells+_nDesignVariables)*sizeof(double));
    
      double scaleFlowModel = _options["scaleObjectiveFunctionFlowModel"]._v;            
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionFlowModel._dv)
	{
	  temporaryArray[ii.first] = scaleFlowModel*ii.second;
	}

      double scaleThermalModel = _options["scaleObjectiveFunctionThermalModel"]._v;
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionThermalModel._dv)
	{
	  temporaryArray[ii.first] += scaleThermalModel*ii.second;
	}

      objFunction = (scaleFlowModel*_objFunctionFlowModel._v + scaleThermalModel*_objFunctionThermalModel._v);
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
                
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdUVWPFileBinary);
      
      for(int i=0; i<5*nCells; i++)
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

      //for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
      for(int i=5*nCells; i<5*nCells+_nDesignVariables; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdBetaFileBinary);
	}
      fclose(dcdBetaFileBinary);

    }
  return objFunction;
}

void WriteThermoFluidLaminarResidualGradientsToPetsc()
{
  for ( int id = 0; i < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& rMomentumCell = dynamic_cast<VectorT3Array&>(_thermalFields.momentumResidual[cells]);
      TArray& rContinuityCell = dynamic_cast<TArray&>(_thermalFields.continuityResidual[cells]);
      TArray& rTemperatureCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);
	
      const int nCells = cells.getCountLevel1();
      //const int nCellsWithoutGhosts = cells.getSelfCount();

      int nCoeffUVWP = 0;
      int nCoeffBeta = 0;

      //Let the total number of cells be N
      //Let the total number of interior cells be M
      //N = M + Number of boundary faces

      //There are 5 state variables U, V, W, P, T; Therefore 5N state variables
      //There are M design variables

      //There are 5N residuals
      //dRdUVWP is a matrix of size 5N X 5N

      //nnzElementsRowdRdUVWP is a double array that stores the number of non-zero values in each row. Or the number of independent variables on which the residual of the cell depends upon. 
      double nnzElementsRowdRdUVWP[5*nCells];
      memset(nnzElementsRowdRdUVWP, 0, (5*nCells)*sizeof(double));
	
      //dRdBeta is a matrix of size 5N X M
      double nnzElementsRowdRdBeta[5*nCells];
      memset(nnzElementsRowdRdBeta, 0, (5*nCells)*sizeof(double));

      cout << endl;
      for(int i=0; i<nCells; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
		  
	      foreach(const typename Rapid::PartialDerivatives::value_type& ii, rMomentumCell[i][j]._dv)
		{
		  //First check whether the derivatives are with respect to X momentum (or U)
		  if (ii.first < 5*nCells)
		    {
		      nnzElementsRowdRdUVWP[5*i+j]++;
		      nCoeffUVWP++;
		    }
		  else
		    {			  
		      nnzElementsRowdRdBeta[5*i+j]++;
		      nCoeffBeta++;
		    }
		}
	    }
	  //cout << "P-Res = " << rContinuityCell[i] << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  nnzElementsRowdRdUVWP[5*i+3]++;
		  nCoeffUVWP++;
		}
	      else
		{
		  nnzElementsRowdRdBeta[5*i+3]++;
		  nCoeffBeta++;
		}
	    }
	  //Adding temperature residual
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  nnzElementsRowdRdUVWP[5*i+4]++;
		  nCoeffUVWP++;
		}
	      else
		{
		  nnzElementsRowdRdBeta[5*i+4]++;
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

      //Number of rows for dRdUVWP = 5N
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);	  

      //Number of columns for dRdUVWP = 5N
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	  
      //Number of rows for dRdBeta = 5N
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
      for(int i=0; i<5*nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdUVWP[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
	}
      for(int i=0; i<5*nCells; i++)
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
		  if (ii.first < 5*nCells)
		    {
		      bigEndianInteger = ii.first;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		    }
		  else
		    {
		      bigEndianInteger = ii.first-5*nCells;
		      LowEndian2BigEndian(bigEndianInteger);
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		    }
		}

	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-5*nCells;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-5*nCells;
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
		  if (ii.first < 5*nCells)
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
		  else
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }	    
	}

      fclose(dRdUVWPFileBinary);
      fclose(dRdBetaFileBinary);
    }
}

#endif
#endif
