#ifndef _FLOWMODELTURBULENCE_H_
#define _FLOWMODELTURBULENCE_H_



void ComputeTotalViscosity()
{
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_PC))
applySpaldingsWallFunction();
#endif
const int numMeshes = _meshes.size();
 for (int n=0; n<numMeshes; n++)
   {
     const Mesh& mesh = *_meshes[n];
     const StorageSite& cells = mesh.getCells();
     const int nCells = cells.getCountLevel1();
	
     TArray& totalviscosity = dynamic_cast<TArray&>(_flowFields.totalviscosity[cells]);
     const TArray& eddyviscosity = dynamic_cast<const TArray&>(_flowFields.eddyviscosity[cells]);
     const TArray& viscosity = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);
      
     for (int i=0; i<nCells; i++)
       {
	 totalviscosity[i] = eddyviscosity[i] + viscosity[i];
       }
   }
}

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_RAPID) || defined(USING_ATYPE_PC))

void applySpaldingsWallFunction()
{

  _velocityGradientModel.compute();

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      TArray& eddyViscosity = dynamic_cast<TArray&>(_flowFields.eddyviscosity[cells]);
      const TArray& viscosity = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);
      const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);
      const TArray& wallDistance = dynamic_cast<const TArray&>(_flowFields.wallDistance[cells]);
      const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
      const VGradArray& vGrad = dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);



      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();

	  const FlowBC<T>& bc = *_bcMap[fg.id];
	  //const TArray& wallDistance = dynamic_cast<const TArray&>(_flowFields.wallDistance[faces]);
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	    

	  if (bc.bcType == "NoSlipWall")
	    {
	      for(int f=0; f<nFaces; f++)
		{
		  const int c0 = faceCells(f,0);
		  const int c1 = faceCells(f,1);      
		  const VectorT3& Af = faceArea[f];
		  const T nuWall = viscosity[c0]/density[c0];
		  const T nutWall = eddyViscosity[c0]/density[c0];
		  const VectorT3 gradVelocityWall = vGrad[c0]*Af;
		  const T magGradVelocityWall = mag(gradVelocityWall);
		  const T y = wallDistance[c0];
		  const T magDeltaVelocity = mag(velocity[c0] - velocity[c1]);
		  //Init uTau
		  T uTau = sqrt((nutWall + nuWall)*magGradVelocityWall);
		  //Correct uTau
		  uTau = ComputeUTau(magDeltaVelocity, y, nuWall, uTau);
		  eddyViscosity[c1] = density[c1]*max(T(0.0), (pow(uTau,T(2.))/(magGradVelocityWall + T(1.0e-150)) - nuWall));
		  //eddyViscosity[c1] = eddyViscosity[c0];
		}
	    }
	}
    }
}
  
#endif

T ComputeUTau(const T magUp, const T y, const T nuw, T uTau)
{
  //Obtained this from courtesy OpenFOAM.
  const T kappa(0.41);
  const T E(9.8);
  int iter(0);
  T ROOTVSMALL(1.e-150);  
  T err(1.0e15);
  do 
    {
      T kUu = min(kappa*magUp/uTau, T(50.0));
      T fkUu = exp(kUu) - 1. - kUu*(1. + 0.5*kUu);

      T f = -uTau*y/nuw + magUp/uTau + 1/E*(fkUu - 1.0/6.0*kUu*pow(kUu,T(2.)));
      T df = y/nuw + magUp/pow(uTau,T(2.)) + 1/E*kUu*fkUu/uTau;
      T uTauNew = uTau + f/df;
      err = abs((uTau - uTauNew)/uTau);
      uTau = uTauNew;

    } while(uTau > ROOTVSMALL && err > 0.01 && ++iter < 10);
  return(max(T(0.0),uTau));
}

void ComputeWallStress(const Mesh& mesh, const int faceGroupId)
{

  VectorT3 r(VectorT3::getZero());
  bool found = false;
  const StorageSite& cells = mesh.getCells();
  _velocityGradientModel.compute();

  const TArray& eddyViscosity = dynamic_cast<const TArray&>(_flowFields.eddyviscosity[cells]);
  const TArray& viscosity = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);
  const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);
  const TArray& wallDistance = dynamic_cast<const TArray&>(_flowFields.wallDistance[cells]);
  const VectorT3Array& velocity = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
  const VGradArray& vGrad = dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

	  shared_ptr<TArray> tauWallCell(new TArray(faces.getCount()));
	  tauWallCell->zero();
	  _flowFields.tauWall.addArray(faces,tauWallCell);

	  TArray& wallStress = dynamic_cast<TArray&>(_flowFields.tauWall[faces]);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c0 = faceCells(f,0);
	      const int c1 = faceCells(f,1);      
	      const VectorT3& Af = faceArea[f];
	      const T nuWall = viscosity[c0]/density[c0];
	      const T nutWall = eddyViscosity[c0]/density[c0];
	      const VectorT3 gradVelocityWall = vGrad[c0]*Af;
	      const T magGradVelocityWall = mag(gradVelocityWall);
	      const T y = wallDistance[c0];
	      const T magDeltaVelocity = mag(velocity[c0] - velocity[c1]);
	      //Init uTau
	      T uTau = sqrt((nutWall + nuWall)*magGradVelocityWall);
	      //Correct uTau
	      uTau = ComputeUTau(magDeltaVelocity, y, nuWall, uTau);
	      wallStress[f] = uTau*uTau;
	    }
	  found=true;
	}
    }
  if (!found)
    throw CException("Invalid faceGroupID");
  //return r;

}

void ComputeVorticityMagnitude()
{
    
  const int numMeshes = _meshes.size();
 
  for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();
      _velocityGradientModel.compute();

      const VGradArray& vGrad = dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);
      TArray& vorticityMag = dynamic_cast<TArray&>(_flowFields.vorticityMagnitude[cells]);

      for(int n=0; n<nCells; n++)
	{
	  const VGradType& vg = vGrad[n];
	  VGradType vgMinusTranspose = vGrad[n];

	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++)
	      vgMinusTranspose[i][j] -= vg[j][i];       

	  vorticityMag[n] = sqrt((pow((vgMinusTranspose[0][1]),2) +  
				  pow((vgMinusTranspose[1][2]),2) + 
				  pow((vgMinusTranspose[2][0]),2)));
	  vorticityMag[n] = abs(vg[1][0] - vg[0][1]);
	}

    }
}


void ComputeShearStress(const Mesh& mesh, const int faceGroupId)
{
  VectorT3 r(VectorT3::getZero());
  bool found = false;
  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCount();
  _velocityGradientModel.compute();

  const VGradArray& vGrad = dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);
  const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);
  const TArray& viscosity = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);

  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupId)
	{
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	  shared_ptr<TArray> tauWallCell(new TArray(faces.getCount()));
	  tauWallCell->zero();
	  _flowFields.tauWall.addArray(faces,tauWallCell);

	  TArray& shearStressWall = dynamic_cast<TArray&>(_flowFields.tauWall[faces]);
	  for(int f=0; f<nFaces; f++)
	    {
	      const int c1 = faceCells(f,1);
	      const VGradType& vg = vGrad[c1];
	      shearStressWall[f] = -viscosity[c1]/density[c1]*vg[1][0];
	    }
	  found=true;
	}
    }
  if (!found)
    throw CException("Invalid faceGroupID");
  //return r;
}

#endif
