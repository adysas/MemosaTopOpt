
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
      _filter.applyFilter(mesh, _thermalFields.sensitivity, _thermalFields.beta);
    }
}

void applyFilterConstraintFunction()
{
    
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      _filter.applyFilter(mesh, _thermalFields.sensitivityConstraint, _thermalFields.beta);
    }
}


#endif


void initResidualFieldThermal()
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

/////////Copying fields and values ///////////////////

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


//void interpolateThermalConductivityFieldSIMP() Replaced the function with SIMP on k, rho, cp
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
      TArray& cond = dynamic_cast<TArray&>(_thermalFields.conductivity[cells]);
      const TArray& beta = dynamic_cast<const TArray&>(_thermalFields.beta[cells]);
	
      kHigh = vc["thermalConductivityHigh"];
      kLow = vc["thermalConductivityLow"];
	
      p = _options["p"];
	
      //cout << endl << "Conductivity" << endl; 
      for ( int c=0; c<nCells; c++)
	{
	  cond[c] = (kHigh - kLow)*pow(beta[c],p) + kLow; 
	}
    }
}


void interpolateMaterialThermalPropertiesFieldSIMP()
{
  const int numMeshes = _meshes.size();
  T kHigh(0.);
  T kLow(0.);
  T rhoHigh(0.);
  T rhoLow(0.);
  T cpHigh(0.);
  T cpLow(0.);  
  T p(0.);

  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const ThermalVC<T>& vc = *_vcMap[mesh.getID()];
      TArray& cond = dynamic_cast<TArray&>(_thermalFields.conductivity[cells]);
      TArray& cp = dynamic_cast<TArray&>(_thermalFields.specificHeat[cells]);
      TArray& rho = dynamic_cast<TArray&>(_thermalFields.density[cells]);
      const TArray& beta = dynamic_cast<const TArray&>(_thermalFields.beta[cells]);
	
      kHigh = vc["thermalConductivityHigh"];
      kLow = vc["thermalConductivityLow"];
      rhoHigh = vc["densityHigh"];
      rhoLow = vc["densityLow"];
      cpHigh = vc["specificHeatHigh"];
      cpLow = vc["specificHeatLow"];
	
      p = _options["p"];
	
      //cout << endl << "Conductivity" << endl; 
      for ( int c=0; c<nCells; c++)
	{
	  cond[c] = (kHigh - kLow)*pow(beta[c],p) + kLow; 
	  rho[c] = (rhoHigh - rhoLow)*beta[c] + rhoLow; 
	  cp[c] = (cpHigh - cpLow)*beta[c] + cpLow; 
	  //cout << c << ":" << cond[c] << endl;
	}
    }
}

void interpolateThermalConductivityFieldSIMPReverse()
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
	  cond[c] = (kHigh - kLow)*pow((T(1.0) - beta[c]),p) + kLow; 
	  //cout << c << ":" << cond[c] << endl;
	}
    }
}


void interpolateThermalConductivityFieldSIMPRamp()
{
  const int numMeshes = _meshes.size();
  T kSolid(0.);
  T kFluid(0.);
  T Ck(0.);
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
	
      kSolid = vc["thermalConductivityHigh"];
      kFluid = vc["thermalConductivityLow"];
	
      Ck = kFluid/kSolid;
      p = _options["p"];
	
      //cout << endl << "Conductivity" << endl; 
      for ( int c=0; c<nCells; c++)
	{
	  //cond[c] = (kLow + (kHigh - kLow)*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	  //cond[c] = (kHigh + (kLow - kHigh)*beta[c]*(1+_options["p"])/(beta[c]+_options["p"]));
	  //cond[c] = (kFluid + (kSolid - kFluid)*(1-beta[c])/(1 + beta[c]*_options["p"]));
	  //cond[c] = 1./(PeSolid + (PeFluid - PeSolid)*(1-beta[c])/(1 + beta[c]*_options["p"]));
	  //cond[c] = 1./(PeFluid + (PeSolid - PeFluid)*(1-beta[c])/(1 + beta[c]*_options["p"]));
	  cond[c] = kFluid*(beta[c]*(Ck*(1.+p)-1.) + 1.)/(Ck*(1.+p*beta[c]));
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


void interpolateHeatTransferCoefficient()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& htc = dynamic_cast<TArray&>(_thermalFields.heatTransferCoefficient[cells]);
      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const int nCells = cells.getCount();		
      for ( int c=0; c<nCells; c++)
	{
	  htc[c] = _options["heatTransferCoefficient"]*beta[c];
	}
    }   
}

void interpolateRadiallyVaryingHeatTransferCoefficient()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& htc = dynamic_cast<TArray&>(_thermalFields.heatTransferCoefficient[cells]);
      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      
      const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

      const int nCells = cells.getCount();		
      for ( int c=0; c<nCells; c++)
	{
	  VectorT3 ds = cellCentroid[c];
	  T_Scalar dsMag = mag(ds);
	  //htc[c] = _options["heatTransferCoefficient"]*dsMag*dsMag*beta[c];
	  htc[c] = _options["heatTransferCoefficient"]*dsMag*dsMag;
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

void setObjectiveFunctionThermal(const T objFunction)
{
  _objFunctionThermalModel = objFunction;
}

void setPressureDropComponent(const VectorT3 pDropComp, const int direction)
{
  _pressureDropComponent = pDropComp[direction];
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

  for ( unsigned int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
    
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& rTemperatureCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);

      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();

      int nCoeffPhi = 0;
      int nCoeffBeta = 0;

      double nnzElementsRowdRdPhi[nCells];
      memset(nnzElementsRowdRdPhi, 0, (nCells)*sizeof(double));
	
      double nnzElementsRowdRdBeta[nCells];
      memset(nnzElementsRowdRdBeta, 0, (nCells)*sizeof(double));

      for(int i=0; i<nCells; i++)
	{
	  {
	    foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	      {
		if (ii.first < nCells)
		  {
		    nnzElementsRowdRdPhi[i]++;
		    nCoeffPhi++;
		  }
		else
		  {
		    nnzElementsRowdRdBeta[i]++;
		    nCoeffBeta++;
		  }
	      }
	  }
	}

      //Write in binary format so that Petsc can read directly
      string dRdPhiFileNameBinary = "dRdPhi.binary";
      FILE *dRdPhiFileBinary = fopen(dRdPhiFileNameBinary.c_str(),"wb"); 
      string dRdBetaFileNameBinary = "dRdBeta.binary";
      FILE *dRdBetaFileBinary = fopen(dRdBetaFileNameBinary.c_str(),"wb"); 

      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      bigEndianInteger = 1211216;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
	
      //Number of rows for dRdPhi
      bigEndianInteger = nCells;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);

      //Number of rows for dRdBeta
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      //Number of columns for dRdPhi
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);

      //Number of columns for dRdBeta
      bigEndianInteger = nCellsWithoutGhosts;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      // Number of non-zero elements in dRdPhi
      bigEndianInteger = nCoeffPhi;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);
	
      // Number of non-zero elements in dRdBeta
      bigEndianInteger = nCoeffBeta;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      // Number of non-zeros in each row in both dRdPhi and dRdBeta

      for(int i=0; i<nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdPhi[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);

	  bigEndianInteger = nnzElementsRowdRdBeta[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
	}

      // Column indices of all non-zeros elements in both dRdPhi and dRdBeta (starting index is zero)
	
      for(int i=0; i<nCells; i++)
	{
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < nCells)
		{
		  bigEndianInteger = ii.first;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdPhiFileBinary);
		}
	      else
		{
		  bigEndianInteger = ii.first-nCells;
		  LowEndian2BigEndian(bigEndianInteger);
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);
		}
	    }
	}

      // Values of all non-zeros for both dRdPhi and dRdBeta
      for(int i=0; i<nCells; i++)
	{
	  //cout << rCell[i]._v << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdPhiFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }
	}
	
      fclose(dRdPhiFileBinary);
      fclose(dRdBetaFileBinary);
    }
}

void WriteObjectiveFunctionGradientToPetscThermal()
{
  //double objFunction(0.);
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();


      double temporaryArray[nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (nCells+_nDesignVariables)*sizeof(double));
    
      double scaleThermalModel = _options["scaleThermalModelObjFunction"]._v;
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionThermalModel._dv)
	{
	  temporaryArray[ii.first] = scaleThermalModel*ii.second;
	}

      //objFunction = scaleThermalModel*_objFunctionThermalModel._v;

      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing dcdPhi
      string dcdTFileNameBinary = "dcdT.binary";
      FILE *dcdTFileBinary = fopen(dcdTFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdTFileBinary);
                
      bigEndianInteger = nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdTFileBinary);
      
      for(int i=0; i<nCells; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdTFileBinary);
	}
      fclose(dcdTFileBinary);

      //Writing dcdBeta
      string dcdBetaFileNameBinary = "dcdBeta.binary";
      FILE *dcdBetaFileBinary = fopen(dcdBetaFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdBetaFileBinary);
      
      bigEndianInteger = _nDesignVariables;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdBetaFileBinary);
            
      for(int i=nCells; i<nCells+_nDesignVariables; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdBetaFileBinary);
	}
      fclose(dcdBetaFileBinary);
    }
  //return objFunction;
}
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
void dumpBetaSensistivityFields()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
	
      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      TArray& sensCell = dynamic_cast<TArray&>(_thermalFields.sensitivity[cells]);
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
	}
      fclose(sensitivityFileBinary);
      fclose(betaFileBinary);
	
    }
}

void dumpVolumeField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);	

      const int nCellsWithoutGhosts = cells.getSelfCount();
	
      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing Sensitivity
      string volumeFileNameBinary = "CellVolume.bin";
      FILE *volumeFileBinary = fopen(volumeFileNameBinary.c_str(),"wb"); 

      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, volumeFileBinary);

      bigEndianInteger = nCellsWithoutGhosts;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, volumeFileBinary);

      for(int i=0; i<nCellsWithoutGhosts; i++)
	{
	  bigEndianDouble = cellVolume[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, volumeFileBinary);
	}
      fclose(volumeFileBinary);	
    }
}

void dumpSensistivityConstraintFields()
{

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
	
      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      TArray& sensCell = dynamic_cast<TArray&>(_thermalFields.sensitivity[cells]);
      const int nCellsWithoutGhosts = cells.getSelfCount();
	
      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing Sensitivity
      string sensitivityFileNameBinary = "SensitivityConstraint.bin";
      FILE *sensitivityFileBinary = fopen(sensitivityFileNameBinary.c_str(),"wb"); 

      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, sensitivityFileBinary);

      bigEndianInteger = nCellsWithoutGhosts;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, sensitivityFileBinary);

      for(int i=0; i<nCellsWithoutGhosts; i++)
	{
	  bigEndianDouble = sensCell[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, sensitivityFileBinary);
	}
      fclose(sensitivityFileBinary);
	
    }
}

#endif



//Functions used to developing code

void printTemperatureField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& tCell = dynamic_cast<TArray&>(_thermalFields.temperature[cells]);
      const int nCells = cells.getCount();      
      cout << "\n";
      for(int n=0; n<nCells; n++)
	{
	  cout << "Nth Cell temperature" << tCell[n] << endl;
	}
    }    
}


void printBetaField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);
      TArray& betaCell = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const int nCells = cells.getCount();      
      cout << "\n";
      for(int n=0; n<nCells; n++)
	{
	  cout << n << "th cell beta" << betaCell[n] << endl;
	}
    }    
}


void initializeBetaField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);

      const int nCellsWithoutGhosts = cells.getSelfCount();
		
      for ( int c=0; c<nCellsWithoutGhosts; c++)
	{
	  beta[c] = ((double) rand() / (RAND_MAX));
	}

	
      beta[0] = 0.3456;
      beta[1] = 0.6789;
      beta[2] = 0.1243;
      beta[3] = 0.9087;
      beta[4] = 0.6745;
      beta[5] = 0.2341;
      beta[6] = 0.8967;
      beta[7] = 0.8123;
      beta[8] = 0.2367;
	

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
	      //beta[c1] = 0;
	    }
	}

    }    
}

void printConductivityField()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();

      MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

      TArray& condCell = dynamic_cast<TArray&>(_thermalFields.conductivity[cells]);

      const int nCells = cells.getCount();  
     
      cout << "\n";
      for(int n=0; n<nCells; n++)
	{
	  cout << n << "th Cell temperature" << condCell[n] << endl;
	}
    }    
}

void printBoundaryCellNumbers(const Mesh& mesh, const int faceGroupId)
{
  const StorageSite& cells = mesh.getCells();
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);
	      cout << "c0: " << c0 << ", c1: " << c1 << endl; 
	    }
	}
    }   
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void copyMassFlux(const Field& massFlux)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];

      const StorageSite& faces = mesh.getFaces();
      const StorageSite& cells = mesh.getCells();
      const CRConnectivity& faceCells = mesh.getFaceCells(faces);

      const int nFaces = faces.getCount();

      TArray& convectionFlux = dynamic_cast<TArray&>(_thermalFields.convectionFlux[faces]);
      TArray& cp = dynamic_cast<TArray&>(_thermalFields.specificHeat[cells]);
      
      const TArray& massFluxFlowModel = dynamic_cast<const TArray&>(massFlux[faces]);
	    
      for(int f=0; f<nFaces; f++)
	{
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);
	  const T cpF = 0.5*(cp[c0]+cp[c1]);
	  convectionFlux[f] = cpF*massFluxFlowModel[f];
	}
    }
}

T getDeltaEnthalpyIntegralThermoFluid(const Mesh& mesh, const int faceGroupId1, const int faceGroupId2)
{
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  const int nCellsWithoutGhosts = cells.getSelfCount();
  T r1(0.0);
  T r2(0.0);
  T r(0.0);
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
  r = r2 - r1;
  return (r);
}

T getBulkTemperature(const Mesh& mesh, const int faceGroupId)
{
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  const int nCellsWithoutGhosts = cells.getSelfCount();

  T Tb(0.0);
  T massFluxIntegral(0.0);
  const TArray& temperature = dynamic_cast<const TArray&>(_thermalFields.temperature[cells]);
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& convectionFlux = dynamic_cast<const TArray&>(_thermalFields.convectionFlux[faces]);	  
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c1 = faceCells(f,1);
	      //cout << convectionFlux[f] << endl;
	      Tb += temperature[c1]*convectionFlux[f];
	      massFluxIntegral += convectionFlux[f];
	    }	 
	  Tb /= massFluxIntegral;
	}      
    }
  return (Tb);
} 


T getRootMeanSquareTemperatureOutletTemperatureIdeal(const Mesh& mesh, const int faceGroupIdIdeal, const int faceGroupIdOut)
{
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCountLevel1();
  const int nCellsWithoutGhosts = cells.getSelfCount();

  T Tb = getBulkTemperature(mesh, faceGroupIdIdeal);

  T rms(0.0);
  const TArray& temperature = dynamic_cast<const TArray&>(_thermalFields.temperature[cells]);
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupIdOut)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const TArray& convectionFlux = dynamic_cast<const TArray&>(_thermalFields.convectionFlux[faces]);	  
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c1 = faceCells(f,1);
	      //rms += pow((temperature[c1] - Tb), 2);
	      rms += convectionFlux[f]*pow((temperature[c1] - Tb), 2);
	    }
	  rms = pow(rms, 0.5);
	}      
    }
  return (rms);
} 
 

void setObjectiveFunctionFlow(const VectorT3 objFunction, const int component)
{
  _objFunctionFlowModel = objFunction[component];
}

#if ( defined(USING_ATYPE_RAPID) )


void setIndependentStateVariablesThermoFluidLaminar()
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

void setIndependentDesignVariablesThermoFluidLaminar()
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


void WriteResidualGradientsToMatlabThermoFluidLaminar()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();    
	
      VectorT3Array& rMomentumCell = dynamic_cast<VectorT3Array&>(_thermalFields.momentumResidual[cells]);
      TArray& rContinuityCell = dynamic_cast<TArray&>(_thermalFields.continuityResidual[cells]);
      TArray& rTemperatureCell = dynamic_cast<TArray&>(_thermalFields.temperatureResidual[cells]);

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
		  if (ii.first < 5*nCells)	       
		    fprintf(DRDSVFile, "%d %d %22.15le \n", 5*i+j+1, ii.first+1, ii.second); 
		  else
		    fprintf(DRDDVFile, "%d %d %22.15le \n", 5*i+j+1, ii.first-5*nCells+1, ii.second); 
		}
	    }
	  //cout << rContinuityCell[i]._v << endl;
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  fprintf(DRDSVFile, "%d %d %22.15le \n", 5*i+3+1, ii.first+1, ii.second); 
		}
	      else
		{
		  fprintf(DRDDVFile, "%d %d %22.15le \n", 5*i+3+1, ii.first-5*nCells+1, ii.second);
		}
	    }	    
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      if (ii.first < 5*nCells)
		{
		  fprintf(DRDSVFile, "%d %d %22.15le \n", 5*i+4+1, ii.first+1, ii.second); 
		}
	      else
		{
		  fprintf(DRDDVFile, "%d %d %22.15le \n", 5*i+4+1, ii.first-5*nCells+1, ii.second);
		}
	    }	 
	}
      fclose(DRDSVFile);
      fclose(DRDDVFile);
    }
}


void WriteResidualGradientsToPetscThermoFluidLaminar()
{
  for ( int id = 0; id < _meshes.size(); id++ )
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
      double nnzElementsRowdRdUVWPT[5*nCells];
      memset(nnzElementsRowdRdUVWPT, 0, (5*nCells)*sizeof(double));
	
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
		      nnzElementsRowdRdUVWPT[5*i+j]++;
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
		  nnzElementsRowdRdUVWPT[5*i+3]++;
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
		  nnzElementsRowdRdUVWPT[5*i+4]++;
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
      string dRdUVWPTFileNameBinary = "dRdUVWPT.binary";
      FILE *dRdUVWPTFileBinary = fopen(dRdUVWPTFileNameBinary.c_str(),"wb"); 
      string dRdBetaFileNameBinary = "dRdBeta.binary";
      FILE *dRdBetaFileBinary = fopen(dRdBetaFileNameBinary.c_str(),"wb"); 

      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;

      //some id number for petsc
      bigEndianInteger = 1211216;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      //Number of rows for dRdUVWP = 5N
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);	  

      //Number of columns for dRdUVWP = 5N
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
	  
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
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);

      // Number of non-zero elements in dRdBeta
      bigEndianInteger = nCoeffBeta;
      LowEndian2BigEndian(bigEndianInteger);	
      fwrite(&bigEndianInteger, sizeof(int), 1, dRdBetaFileBinary);

      // Number of non-zeros in each row in both dRdUVWP and dRdBeta
      for(int i=0; i<5*nCells; i++)
	{
	  bigEndianInteger = nnzElementsRowdRdUVWPT[i];
	  LowEndian2BigEndian(bigEndianInteger);	
	  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
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
		      fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
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
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
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
		  fwrite(&bigEndianInteger, sizeof(int), 1, dRdUVWPTFileBinary);
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
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPTFileBinary);
		  else
		    fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
		}
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rContinuityCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPTFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }
	  foreach(const typename Rapid::PartialDerivatives::value_type& ii, rTemperatureCell[i]._dv)
	    {
	      bigEndianDouble = ii.second;
	      LowEndian2BigEndian(bigEndianDouble);
	      if (ii.first < 5*nCells)
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdUVWPTFileBinary);
	      else
		fwrite(&bigEndianDouble, sizeof(double), 1, dRdBetaFileBinary);
	    }	    
	}

      fclose(dRdUVWPTFileBinary);
      fclose(dRdBetaFileBinary);
    }
}


void WriteObjectiveFunctionGradientToPetscThermoFluidLaminar()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();


      double temporaryArray[5*nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (5*nCells+_nDesignVariables)*sizeof(double));
    
      double scaleFlowModel = _options["scaleFlowModelObjFunction"]._v;            
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionFlowModel._dv)
	{
	  temporaryArray[ii.first] = scaleFlowModel*ii.second;
	}

      double scaleThermalModel = _options["scaleThermalModelObjFunction"]._v;
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionThermalModel._dv)
	{
	  temporaryArray[ii.first] += scaleThermalModel*ii.second;
	}

      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      //Writing dcdPhi
      string dcdUVWPTFileNameBinary = "dcdUVWPT.binary";
      FILE *dcdUVWPTFileBinary = fopen(dcdUVWPTFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdUVWPTFileBinary);
                
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dcdUVWPTFileBinary);
      
      for(int i=0; i<5*nCells; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dcdUVWPTFileBinary);
	}
      fclose(dcdUVWPTFileBinary);

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
}

void WritePressureDropGradientToPetscThermoFluidLaminar()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const int nCellsWithoutGhosts = cells.getSelfCount();


      double temporaryArray[5*nCells+_nDesignVariables];
      memset(temporaryArray, 0.0, (5*nCells+_nDesignVariables)*sizeof(double));
    
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _pressureDropComponent._dv)
	{
	  temporaryArray[ii.first] = ii.second;
	}

      //Write in binary format so that Petsc can read directly
      //Temporary variable for big endian version of the variable
      int bigEndianInteger;
      double bigEndianDouble;
      
      
      //Writing dPDdPhi
      string dPDdUVWPTFileNameBinary = "dcdUVWPT.binary";
      FILE *dPDdUVWPTFileBinary = fopen(dPDdUVWPTFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdUVWPTFileBinary);
                
      bigEndianInteger = 5*nCells;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdUVWPTFileBinary);
      
      for(int i=0; i<5*nCells; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dPDdUVWPTFileBinary);
	}
      fclose(dPDdUVWPTFileBinary);
      

      //Writing dcdBeta
      string dPDdBetaFileNameBinary = "dcdBeta.binary";
      FILE *dPDdBetaFileBinary = fopen(dPDdBetaFileNameBinary.c_str(),"wb"); 
      bigEndianInteger = 1211214;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdBetaFileBinary);
      
      bigEndianInteger = _nDesignVariables;
      LowEndian2BigEndian(bigEndianInteger);
      fwrite(&bigEndianInteger, sizeof(int), 1, dPDdBetaFileBinary);

      //for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
      for(int i=5*nCells; i<5*nCells+_nDesignVariables; i++)
	{
	  bigEndianDouble = temporaryArray[i];
	  LowEndian2BigEndian(bigEndianDouble);
	  fwrite(&bigEndianDouble, sizeof(double), 1, dPDdBetaFileBinary);
	}
      fclose(dPDdBetaFileBinary);
    }
}

//This is not fully implemented
void WritePressureDropGradientToMatlabThermoFluidLaminar()
{
  for ( int id = 0; id < _meshes.size(); id++ )
    {
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      
      double temporaryArray[5*nCells+_nDesignVariables];
      memset(temporaryArray, 0, (5*nCells+_nDesignVariables)*sizeof(double));    
       
      foreach(const typename Rapid::PartialDerivatives::value_type& ii, _objFunctionFlowModel._dv)
	{
	  temporaryArray[ii.first] = ii.second;
	}

      string DV = "DPressureDropDDV.mat";
      FILE *DPrDropDDVFile = fopen(DV.c_str(),"wb");


      //for(int i=4*nCells; i<4*nCells+nCellsWithoutGhosts; i++)
      for(int i=5*nCells; i<5*nCells+_nDesignVariables; i++)
	{
	  fprintf(DPrDropDDVFile, "%22.15le \n", temporaryArray[i]); 
	}
      fclose(DPrDropDDVFile);
    }
}

#endif


void copyBetaField(const Field& betaField)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
			
      const int nCells = cells.getCountLevel1();
      const int nCellsSelfCount = cells.getSelfCount();

      TArray& beta = dynamic_cast<TArray&>(_thermalFields.beta[cells]);
      const TArray& betaFieldFrom = dynamic_cast<const TArray&>(betaField[cells]);
	    
      for(int i=0; i<nCells; i++)
	{
	  beta[i] = betaFieldFrom[i];
	}
      _nDesignVariables = nCellsSelfCount; 
    }
}


void copyTemperatureBoundaryFace(const Field& temperatureDoubleField, Field& temperatureField, const StorageSite& facesDouble, const StorageSite& faces)
{
  //We are copying only one component of velocity
  const int nFaces = faces.getCount();
  const Array<double>& temperatureFaceDouble = dynamic_cast<const Array<double>&>(temperatureDoubleField[facesDouble]);
  TArray& temperatureFace = dynamic_cast<TArray&>(temperatureField[faces]);
  for(int f=0; f<nFaces; f++)
    {
      temperatureFace[f] = temperatureFaceDouble[f];
    }
}
