// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Mesh.h"
#include <sstream>

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
//#include "FieldSet.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GenericBCS.h"
#include "Vector.h"
#include "DiffusionDiscretization.h"
//#include "DiffusionDiscretizationHamiltonJacobiSpecific.h"
#include "ConvectionDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"
#include "TimeDerivativeDiscretization.h"
#include "SourceDiscretizationSP.h"

template<class T>
class HamiltonJacobiModel<T>::Impl
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
 Impl(const GeomFields& geomFields,
      HamiltonJacobiFields& hamiltonJacobiFields,
      const MeshList& meshes) :
  _meshes(meshes),
    _geomFields(geomFields),
    _hamiltonJacobiFields(hamiltonJacobiFields),
    _hamiltonJacobiGradientModel(_meshes,_hamiltonJacobiFields.phi,
				 _hamiltonJacobiFields.phiGradient,_geomFields),
    _initialNorm(),
    _niters(0),
    _filter(_geomFields)
    {
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  HamiltonJacobiVC<T> *vc(new HamiltonJacobiVC<T>());
	  vc->vcType = "flow";
	  _vcMap[mesh.getID()] = vc;
        
	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      HamiltonJacobiBC<T> *bc(new HamiltonJacobiBC<T>());
            
	      _bcMap[fg.id] = bc;

	      if ((fg.groupType == "wall") ||
		  (fg.groupType == "symmetry"))
		{
		  bc->bcType = "SpecifiedPhiFlux";
		}
	      else if ((fg.groupType == "velocity-inlet") ||
		       (fg.groupType == "pressure-outlet"))
		{
		  bc->bcType = "SpecifiedPhi";
		}
	      else
		throw CException("HamiltonJacobiModel: unknown face group type "
				 + fg.groupType);
	    }
	}
    }

  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        const HamiltonJacobiVC<T>& vc = *_vcMap[mesh.getID()];

	//phi
        shared_ptr<TArray> tCell(new TArray(cells.getCountLevel1()));
        *tCell = _options["initialPhi"];
        _hamiltonJacobiFields.phi.addArray(cells,tCell);
	
	//conductivity
        shared_ptr<TArray> gammaCell(new TArray(cells.getCountLevel1()));
        *gammaCell = vc["gamma"];
        _hamiltonJacobiFields.gamma.addArray(cells,gammaCell);
	
	//source 
	shared_ptr<TArray> sCell(new TArray(cells.getCountLevel1()));
	*sCell = T(-1.);
	_hamiltonJacobiFields.source.addArray(cells,sCell);

	//source solid
	shared_ptr<TArray> ssolidCell(new TArray(cells.getCountLevel1()));
	*ssolidCell = T(0.);
	_hamiltonJacobiFields.sourceSolid.addArray(cells,ssolidCell);

	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	*zeroCell = T(0.0);
	_hamiltonJacobiFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	*oneCell = T(1.0);
	_hamiltonJacobiFields.one.addArray(cells,oneCell);

	//create continuity field
	shared_ptr<TArray> contCell(new TArray(cells.getCountLevel1()));
	*contCell = T(0.0);
	_hamiltonJacobiFields.continuityResidual.addArray(cells,contCell);

	//create wall distance
	shared_ptr<TArray> wdCell(new TArray(cells.getCountLevel1()));
	*wdCell = T(0.0);
	_hamiltonJacobiFields.wallDistance.addArray(cells,wdCell);

	//beta
        shared_ptr<TArray> betaCell(new TArray(cells.getCountLevel1()));
        *betaCell = _options["volumeFraction"];
        _hamiltonJacobiFields.beta.addArray(cells,betaCell);

	//initial temparature gradient array
	shared_ptr<TGradArray> gradPhi(new TGradArray(cells.getCountLevel1()));
	gradPhi->zero();
	_hamiltonJacobiFields.phiGradient.addArray(cells,gradPhi);
        
	//inital convection flux at faces

	
	const StorageSite& allFaces = mesh.getFaces();
	shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	convFlux->zero();
	_hamiltonJacobiFields.convectionFlux.addArray(allFaces,convFlux);
	
	
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _hamiltonJacobiFields.phiFlux.addArray(faces,fluxFace);
          
	  }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _hamiltonJacobiFields.phiFlux.addArray(faces,fluxFace);
          
	  }

	
      }
    _hamiltonJacobiFields.gamma.syncLocal();
    _niters  = 0;
    _initialNorm = MFRPtr();
  }
  


  
  HamiltonJacobiBCMap& getBCMap() {return _bcMap;}
  HamiltonJacobiVCMap& getVCMap() {return _vcMap;}

  HamiltonJacobiBC<T>& getBC(const int id) {return *_bcMap[id];}

  HamiltonJacobiModelOptions<T>& getOptions() {return _options;}

  void initLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_hamiltonJacobiFields.phi,&cells);

        ls.getX().addArray(tIndex,_hamiltonJacobiFields.phi.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_hamiltonJacobiFields.phiFlux,&faces);
            ls.getX().addArray(fIndex,_hamiltonJacobiFields.phiFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	  }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_hamiltonJacobiFields.phiFlux,&faces);
            ls.getX().addArray(fIndex,_hamiltonJacobiFields.phiFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	  }

      }
  }

  void ComputeConvectionFlux()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];

	const HamiltonJacobiVC<T>& vc = *_vcMap[mesh.getID()];
            
	const StorageSite& cells = mesh.getCells();
	const StorageSite& faces = mesh.getFaces();
	const int nFaces = faces.getCount();
	const CRConnectivity& faceCells = mesh.getAllFaceCells();

	const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

	const TGradArray& phiGrad = dynamic_cast<const TGradArray&>(_hamiltonJacobiFields.phiGradient[cells]);
	TArray& convectingFlux = dynamic_cast<TArray&>(_hamiltonJacobiFields.convectionFlux[faces]);

	for(int f=0; f<nFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    //convectingFlux[f] = 0.5*( dot(phiGrad[c0],faceArea[f]) + dot(phiGrad[c1],faceArea[f]) );
	    convectingFlux[f] = 0.5*(phiGrad[c0][0]*faceArea[f][0]+phiGrad[c0][1]*faceArea[f][1]+phiGrad[c0][2]*faceArea[f][2] + phiGrad[c1][0]*faceArea[f][0]+phiGrad[c1][1]*faceArea[f][1]+phiGrad[c1][2]*faceArea[f][2]);
	  }
      }

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    const CRConnectivity& faceCells = mesh.getAllFaceCells();
	    TArray& convectingFlux = dynamic_cast<TArray&>(_hamiltonJacobiFields.convectionFlux[faces]);
                
	    const VectorT3Array& faceArea =
	      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
                
	    const TGradArray& phiGrad = dynamic_cast<const TGradArray&>(_hamiltonJacobiFields.phiGradient[cells]);
	    
	    for(int f=0; f<nFaces; f++)
	      {
		const int c1 = faceCells(f,1);
		//convectingFlux[f] = dot(phiGrad[c1],faceArea[f]);
		convectingFlux[f] = phiGrad[c1][0]*faceArea[f][0]+phiGrad[c1][1]*faceArea[f][1]+phiGrad[c1][2]*faceArea[f][2];
	      }

	  }
      }
  }

  void ComputeContinuityResidual()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();
	const StorageSite& faces = mesh.getFaces();

	TArray& r = dynamic_cast<TArray&>(_hamiltonJacobiFields.continuityResidual[cells]);
	const TArray& convectingFlux = dynamic_cast<const TArray&>(_hamiltonJacobiFields.convectionFlux[faces]);

	const CRConnectivity& faceCells = mesh.getAllFaceCells();
	const int nFaces = faces.getCount();

	r.zero();
	for(int f=0; f<nFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);

	    r[c0] += convectingFlux[f];
	    r[c1] -= convectingFlux[f];
	  }
      }
  }

  void linearize(LinearSystem& ls)
  {
    _hamiltonJacobiGradientModel.compute();
    
    ComputeConvectionFlux();
    ComputeContinuityResidual();

    DiscrList discretizations;   

    /*
      shared_ptr<Discretization>
      dd(new DiffusionDiscretizationHamiltonJacobi<T,T,T>
      (_meshes,_geomFields,
      _hamiltonJacobiFields.phi,
      _hamiltonJacobiFields.gamma,
      _hamiltonJacobiFields.phiGradient,
      _options["epsilon"]));
      discretizations.push_back(dd);
    */

    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _hamiltonJacobiFields.phi,
	  _hamiltonJacobiFields.one,
	  _hamiltonJacobiFields.phiGradient));

    discretizations.push_back(dd);

    /*    
	  shared_ptr<Discretization>
	  cd(new ConvectionDiscretization<T,T,T>
	  (_meshes,_geomFields,
	  _hamiltonJacobiFields.phi,
	  _hamiltonJacobiFields.convectionFlux,
	  _hamiltonJacobiFields.continuityResidual,
	  _hamiltonJacobiFields.phiGradient,
	  _options.useCentralDifference));

	  discretizations.push_back(cd);
    */

    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _hamiltonJacobiFields.phi,
	  _hamiltonJacobiFields.source));

    discretizations.push_back(sd);
    


    shared_ptr<Discretization>
      ssolid(new SourceDiscretizationSP<T, T, T>
	     (_meshes,
	      _geomFields, 
	      _hamiltonJacobiFields.phi,
	      _hamiltonJacobiFields.sourceSolid));

    discretizations.push_back(ssolid);
          
    Linearizer linearizer;

    linearizer.linearize(discretizations, _meshes, ls.getMatrix(), ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const HamiltonJacobiBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _hamiltonJacobiFields.phi,
                                  _hamiltonJacobiFields.phiFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedPhi")
	      {
                FloatValEvaluator<T>
                  bPhi(bc.getVal("specifiedPhi"),faces);
                if (_hamiltonJacobiFields.convectionFlux.hasArray(faces))
		  {
                    const TArray& convectingFlux =
                      dynamic_cast<const TArray&>
                      (_hamiltonJacobiFields.convectionFlux[faces]);
                    const int nFaces = faces.getCount();
                                
                    for(int f=0; f<nFaces; f++)
		      {
                        if (convectingFlux[f] > 0.)
			  {
                            gbc.applyExtrapolationBC(f);
			  }
                        else
			  {
                            gbc.applyDirichletBC(f,bPhi[f]);
			  }
		      }
		  }
                else
		  gbc.applyDirichletBC(bPhi);
	      }
            else if (bc.bcType == "SpecifiedPhiFlux")
	      {
                FloatValEvaluator<T>
		  bPhiFlux(bc.getVal("specifiedPhiFlux"),faces);
                    
                const int nFaces = faces.getCount();
                                
                for(int f=0; f<nFaces; f++)
		  {                        
		    gbc.applyNeumannBC(f, bPhiFlux[f]);
		  }                              
	      }
            else if (bc.bcType == "Symmetry")
	      {
		T zeroFlux(NumTypeTraits<T>::getZero());
		gbc.applyNeumannBC(zeroFlux);
	      }	 
            else
              throw CException(bc.bcType + " not implemented for HamiltonJacobiModel");
	  }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _hamiltonJacobiFields.phi,
                                  _hamiltonJacobiFields.phiFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
	  }
      }

    /*
      DiscrList discretizations2;
      shared_ptr<Discretization>
      ud(new Underrelaxer<T,T,T>
      (_meshes,_hamiltonJacobiFields.phi,
      _options["phiURF"]));

      discretizations2.push_back(ud);

      linearizer.linearize(discretizations2,_meshes,ls.getMatrix(),
      ls.getX(), ls.getB());*/

  }
  

  T getPhiFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
	  {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& phiFlux =
              dynamic_cast<const TArray&>(_hamiltonJacobiFields.phiFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += phiFlux[f];
            found=true;
	  }
      }
    if (!found)
      throw CException("getPhiFluxIntegral: invalid faceGroupID");
    return r;
  }

  void ComputeWallDistance()
  {
    _hamiltonJacobiGradientModel.compute();
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	TArray& wallDistance = dynamic_cast<TArray&>(_hamiltonJacobiFields.wallDistance[cells]);
	const TArray& phi = dynamic_cast<const TArray&>(_hamiltonJacobiFields.phi[cells]);
	const TGradArray& phiGrad = dynamic_cast<const TGradArray&>(_hamiltonJacobiFields.phiGradient[cells]);
	
	for (int i=0; i<nCells; i++)
	  {
	    wallDistance[i] = sqrt(pow(phiGrad[i][0],2) + pow(phiGrad[i][1],2) + pow(phiGrad[i][2],2));
	  }
      }
  }


  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)
      { 
	LinearSystem ls;
	initLinearization(ls);
        
	ls.initAssembly();

	linearize(ls);

	ls.initSolve();

	MFRPtr rNorm(_options.getLinearSolver().solve(ls));

	if (!_initialNorm) _initialNorm = rNorm;
        
	MFRPtr normRatio((*rNorm)/(*_initialNorm));

	cout << _niters << ": " << *rNorm << endl;

        
	_options.getLinearSolver().cleanup();

	ls.postSolve();
	ls.updateSolution();

	_niters++;

	if (*rNorm < _options.absoluteTolerance || *normRatio < _options.relativeTolerance)
	  break;
      }
  }


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  
  void dumpMatrix(const string fileBase)
  {
    LinearSystem ls;
    initLinearization(ls);
    
    ls.initAssembly();
    
    linearize(ls);
    
    ls.initSolve();

    MultiFieldMatrix& matrix = ls.getMatrix();
    MultiField& b = ls.getB();
    for ( unsigned int id = 0; id < _meshes.size(); id++ ){
      const Mesh& mesh = *_meshes[id];
      const StorageSite& cells = mesh.getCells();
    
      MultiField::ArrayIndex tIndex(&_hamiltonJacobiFields.phi,&cells);

      T_Matrix& tMatrix =
	dynamic_cast<T_Matrix&>(matrix.getMatrix(tIndex,tIndex));

      TArray& tDiag = tMatrix.getDiag();
      TArray& tCoeff = tMatrix.getOffDiag();

      TArray& rCell = dynamic_cast<TArray&>(b[tIndex]);

      const CRConnectivity& cr = tMatrix.getConnectivity();

      const Array<int>& row = cr.getRow();
      const Array<int>& col = cr.getCol();
    
      const int nCells = cells.getSelfCount();
      int nCoeffs = nCells;

      for(int i=0; i<nCells; i++)
	for(int jp=row[i]; jp<row[i+1]; jp++)
	  {
	    const int j = col[jp];
	    if (j<nCells) nCoeffs++;
	  }
      stringstream ss;
      ss << id;
      string matFileName = fileBase + "_mesh" + ss.str() +  ".mat";


      FILE *matFile = fopen(matFileName.c_str(),"wb");
    
      fprintf(matFile,"%%%%MatrixMarket matrix coordinate real general\n");
      fprintf(matFile,"%d %d %d\n", nCells,nCells,nCoeffs);

      for(int i=0; i<nCells; i++)
	{
	  fprintf(matFile,"%d %d %le\n", i+1, i+1, tDiag[i]);
	  for(int jp=row[i]; jp<row[i+1]; jp++)
	    {
	      const int j = col[jp];
	      if (j<nCells)
		fprintf(matFile,"%d %d %le\n", i+1, j+1, tCoeff[jp]);
	    }
	}

      fclose(matFile);

      string rhsFileName = fileBase + ".rhs";
      FILE *rhsFile = fopen(rhsFileName.c_str(),"wb");
    
      for(int i=0; i<nCells; i++)
	fprintf(rhsFile,"%lf\n",-rCell[i]);

      fclose(rhsFile);

    }
  }

#endif


  /////////////////////////////////////////////////////////////////////////////////////////////////////////


 private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  HamiltonJacobiFields& _hamiltonJacobiFields;

  HamiltonJacobiBCMap _bcMap;
  HamiltonJacobiVCMap _vcMap;
  HamiltonJacobiModelOptions<T> _options;
  GradientModel<T> _hamiltonJacobiGradientModel;
  
  MFRPtr _initialNorm;
  int _niters;
  Filter<T> _filter; 

  int _nDesignVariables = 0;
  
  
};

template<class T>
HamiltonJacobiModel<T>::HamiltonJacobiModel(const GeomFields& geomFields,
					    HamiltonJacobiFields& hamiltonJacobiFields,
					    const MeshList& meshes) :
Model(meshes),
  _impl(new Impl(geomFields,hamiltonJacobiFields,meshes))
{
  logCtor();
}


template<class T>
HamiltonJacobiModel<T>::~HamiltonJacobiModel()
{
  logDtor();
}

template<class T>
void
HamiltonJacobiModel<T>::init()
{
  _impl->init();
}
  


template<class T>
typename HamiltonJacobiModel<T>::HamiltonJacobiBCMap&
HamiltonJacobiModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename HamiltonJacobiModel<T>::HamiltonJacobiVCMap&
HamiltonJacobiModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
HamiltonJacobiBC<T>&
HamiltonJacobiModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
HamiltonJacobiModelOptions<T>&
HamiltonJacobiModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
HamiltonJacobiModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}



#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
template<class T>
void
HamiltonJacobiModel<T>::dumpMatrix(const string fileBase)
{
  _impl->dumpMatrix(fileBase);
}
#endif


template<class T>
T
HamiltonJacobiModel<T>::getPhiFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getPhiFluxIntegral(mesh, faceGroupId);
}

template<class T>
void
HamiltonJacobiModel<T>::ComputeWallDistance()
{
  return _impl->ComputeWallDistance();
}
