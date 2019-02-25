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
#include "ConvectionDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"
#include "SourceDiscretizationSpalartAllmarasTerm1.h"
#include "SourceDiscretizationSpalartAllmarasTerm2.h"
#include "TimeDerivativeDiscretization.h"

template<class T>
class SpalartAllmarasModel<T>::Impl
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  
 Impl(const GeomFields& geomFields,
      SpalartAllmarasFields& spalartAllmarasFields,
      const MeshList& meshes) :
  _meshes(meshes),
    _geomFields(geomFields),
    _spalartAllmarasFields(spalartAllmarasFields),
    _nuTildeGradientModel(_meshes,_spalartAllmarasFields.nuTilde,
			  _spalartAllmarasFields.nuTildeGradient,_geomFields),
    _initialNorm(),
    _niters(0)
      {
	const int numMeshes = _meshes.size();
	for (int n=0; n<numMeshes; n++)
	  {
	    const Mesh& mesh = *_meshes[n];
	    SpalartAllmarasVC<T> *vc(new SpalartAllmarasVC<T>());
	    vc->vcType = "flow";
	    _vcMap[mesh.getID()] = vc;
        
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		SpalartAllmarasBC<T> *bc(new SpalartAllmarasBC<T>());
            
		_bcMap[fg.id] = bc;

		if ((fg.groupType == "wall") ||
		    (fg.groupType == "symmetry"))
		  {
		    bc->bcType = "SpecifiedNuTildeFlux";
		  }
		else if ((fg.groupType == "velocity-inlet") ||
			 (fg.groupType == "pressure-outlet"))
		  {
		    bc->bcType = "SpecifiedNuTilde";
		  }
		else
		  throw CException("SpalartAllmarasModel: unknown face group type "
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
        const SpalartAllmarasVC<T>& vc = *_vcMap[mesh.getID()];

	//nuTilde
        shared_ptr<TArray> nuTildeCell(new TArray(cells.getCountLevel1()));
        *nuTildeCell = _options["initialNuTilde"];
        _spalartAllmarasFields.nuTilde.addArray(cells,nuTildeCell);
	
	if(_options.transient)
	  {
	    _spalartAllmarasFields.nuTildeN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(nuTildeCell->newCopy()));
	    if (_options.timeDiscretizationOrder > 1)
	      _spalartAllmarasFields.nuTildeN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(nuTildeCell->newCopy()));
	  }
	
	//density
        shared_ptr<TArray> densityCell(new TArray(cells.getCountLevel1()));
        *densityCell = vc["density"];
        _spalartAllmarasFields.density.addArray(cells,densityCell);
	
	//nu
        shared_ptr<TArray> nuCell(new TArray(cells.getCountLevel1()));
        *nuCell = vc["viscosity"]/vc["density"];
        _spalartAllmarasFields.nu.addArray(cells,nuCell);

	//initial nuTilde gradient array
	shared_ptr<TGradArray> gradnuTilde(new TGradArray(cells.getCountLevel1()));
	gradnuTilde->zero();
	_spalartAllmarasFields.nuTildeGradient.addArray(cells,gradnuTilde);

	//sourceTerm1
	shared_ptr<TArray> sTerm1Cell(new TArray(cells.getSelfCount()));
	*sTerm1Cell = T(0.);
	_spalartAllmarasFields.sourceTerm1.addArray(cells,sTerm1Cell);

	//sourceTerm2
	shared_ptr<TArray> sTerm2Cell(new TArray(cells.getSelfCount()));
	*sTerm2Cell = T(0.);
	_spalartAllmarasFields.sourceTerm2.addArray(cells,sTerm2Cell);

	//inital convection flux at faces

	const StorageSite& allFaces = mesh.getFaces();
	shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	convFlux->zero();
	_spalartAllmarasFields.convectionFlux.addArray(allFaces,convFlux);


	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	*zeroCell = T(0.0);
	_spalartAllmarasFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	*oneCell = T(1.0);
	_spalartAllmaraslFields.one.addArray(cells,oneCell);



	//diffusivity for term 1
	shared_ptr<TArray> diffusivityTerm1Cell(new TArray(cells.getCountLevel1()));
	*diffusivityTerm1Cell = T(0.0);
	_spalartAllmaraslFields.diffusivityTerm1.addArray(cells,diffusivityTerm1Cell);

	//diffusivity for term 2
	shared_ptr<TArray> diffusivityTerm2Cell(new TArray(cells.getCountLevel1()));
	*diffusivityTerm2Cell = T(-1.0);
	_spalartAllmaraslFields.diffusivityTerm2.addArray(cells,diffusivityTerm2Cell);

	//vorticity
	shared_ptr<TArray> vorticityMagCell(new TArray(cells.getCountLevel1()));
	*vorticityMagCell = T(0.0);
	_spalartAllmarasFields.vorticityMagnitude.addArray(cells,vorticityMagCell);

	// chi
	shared_ptr<TArray> chiCell(new TArray(cells.getSelfCount()));
	*chiCell = T(0.);
	_spalartAllmarasFields.chi.addArray(cells,chiCell);

	// fv1
	shared_ptr<TArray> fv1Cell(new TArray(cells.getSelfCount()));
	*fv1Cell = T(0.);
	_spalartAllmarasFields.fv1.addArray(cells,fv1Cell);

	// fv2
	shared_ptr<TArray> fv2Cell(new TArray(cells.getSelfCount()));
	*fv2Cell = T(0.);
	_spalartAllmarasFields.fv2.addArray(cells,fv2Cell);

	// r
	shared_ptr<TArray> rSACell(new TArray(cells.getSelfCount()));
	*rSACell = T(0.);
	_spalartAllmarasFields.r.addArray(cells,rSACell);

	// g
	shared_ptr<TArray> rSACell(new TArray(cells.getSelfCount()));
	*gCell = T(0.);
	_spalartAllmarasFields.g.addArray(cells,gCell);

	// fw
	shared_ptr<TArray> fwCell(new TArray(cells.getSelfCount()));
	*fwCell = T(0.);
	_spalartAllmarasFields.fw.addArray(cells,fwCell);

	// ft2
	shared_ptr<TArray> ft2Cell(new TArray(cells.getSelfCount()));
	*ft2Cell = T(0.);
	_spalartAllmarasFields.ft2.addArray(cells,ft2Cell);

	//vorticityTurbulent
	shared_ptr<TArray> vorticityTurbulentCell(new TArray(cells.getCountLevel1()));
	*vorticityTurbulentCell = T(0.0);
	_spalartAllmarasFields.vorticityTurbulent.addArray(cells,vorticityTurbulentCell);

	//effective viscosity
	shared_ptr<TArray> effectiveViscosityCell(new TArray(cells.getCountLevel1()));
	*effectiveViscosityCell = T(0.0);
	_spalartAllmarasFields.effectiveViscosity.addArray(cells,effectiveViscosityCell);

	// wall distance
	shared_ptr<TArray> wdCell(new TArray(cells.getSelfCount()));
	*wdCell = T(0.);
	_spalartAllmarasFields.wallDistance.addArray(cells,wdCell);
	
	//beta
        shared_ptr<TArray> betaCell(new TArray(cells.getCountLevel1()));
        *betaCell = _options["volumeFraction"];
        _spalartAllmarasFields.beta.addArray(cells,betaCell);

	//heat flux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _spalartAllmarasFields.nuTildeFlux.addArray(faces,fluxFace);
          
	  }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _spalartAllmarasFields.nuTildeFlux.addArray(faces,fluxFace);
          
	  }
	
	shared_ptr<VectorT3Array> resMomCell(new VectorT3Array(cells.getCountLevel1()));
	resMomCell->zero();
	_spalartAllmarasFields.momentumResidual.addArray(cells,resMomCell);
	
	shared_ptr<TArray> resContCell(new TArray(cells.getCountLevel1()));
	resContCell->zero();
	_spalartAllmarasFields.continuityResidual.addArray(cells, resContCell);
	
	shared_ptr<TArray> resTempCell(new TArray(cells.getCountLevel1()));
	resSACell->zero();
	_spalartAllmarasFields.spalartAllmarasResidual.addArray(cells, resSACell);


      }
    _spalartAllmarasFields.conductivity.syncLocal();
    _niters  =0;
    _initialNorm = MFRPtr();
  }
  
  SpalartAllmarasBCMap& getBCMap() {return _bcMap;}
  SpalartAllmarasVCMap& getVCMap() {return _vcMap;}

  SpalartAllmarasBC<T>& getBC(const int id) {return *_bcMap[id];}

  SpalartAllmarasModelOptions<T>& getOptions() {return _options;}

  void initLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_spalartAllmarasFields.nuTilde,&cells);

        ls.getX().addArray(tIndex,_spalartAllmarasFields.nuTilde.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_spalartAllmarasFields.nuTildeFlux,&faces);
            ls.getX().addArray(fIndex,_spalartAllmarasFields.nuTildeFlux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_spalartAllmarasFields.heatFlux,&faces);
            ls.getX().addArray(fIndex,_spalartAllmarasFields.nuTildeFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	  }

      }
  }

  void ComputeDiffusivityTerm1()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {
     
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
	
	TArray& diffusivityTerm1 = dynamic_cast<TArray&>(_spalartAllmarasFields.diffusivityTerm1[cells]);
	const TArray& nu = dynamic_cast<const TArray&>(_spalartAllmarasFields.nu[cells]);
	const TArray& density = dynamic_cast<const TArray&>(_spalartAllmarasFields.density[cells]);
	const TArray& nuTilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.nuTilde[cells]);
      
	T Cb2 = _options["Cb2"];
	T sigma = _options["sigma"];

	for (int i=0; i<nCells; i++)
	  {
	    diffusivity[i] = (1+Cb2)/sigma[i]*density[i]*(nu[i] + nuTilde[i]);
	  }
      }
  }

  void ComputeDiffusivityTerm2()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {
     
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
	
	TArray& diffusivityTerm2 = dynamic_cast<TArray&>(_spalartAllmarasFields.diffusivityTerm2[cells]);
	const TArray& nu = dynamic_cast<const TArray&>(_spalartAllmarasFields.nu[cells]);
	const TArray& density = dynamic_cast<const TArray&>(_spalartAllmarasFields.density[cells]);
	const TArray& nuTilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.nuTilde[cells]);
      
	T Cb2 = _options["Cb2"];
	T sigma = _options["sigma"];

	for (int i=0; i<nCells; i++)
	  {
	    diffusivityTerm2[i] = Cb2/sigma[i]*density[i]*(nu[i] + nuTilde[i]);
	  }
      }
  }
   
  void ComputeSpalartAllmarasClosureFunctions()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {
     
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	
	const TArray& nu = dynamic_cast<const TArray&>(_spalartAllmarasFields.nu[cells]);
	const TArray& nuTilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.nuTilde[cells]);
	const TArray& vorticityMagnitude = dynamic_cast<const TArray&>(_spalartAllmarasFields.vorticityMagnitude[cells]);
	const TArray& wallDistance = dynamic_cast<const TArray&>(_spalartAllmarasFields.wallDistance[cells]);

	TArray& chi = dynamic_cast<TArray&>(_spalartAllmarasFields.chi[cells]);
	TArray& fv1 = dynamic_cast<TArray&>(_spalartAllmarasFields.fv1[cells]);
	TArray& fv2 = dynamic_cast<TArray&>(_spalartAllmarasFields.fv2[cells]);
	TArray& STilde = dynamic_cast<TArray&>(_spalartAllmarasFields.STilde[cells]);
	TArray& r = dynamic_cast<TArray&>(_spalartAllmarasFields.r[cells]);
	TArray& g = dynamic_cast<TArray&>(_spalartAllmarasFields.g[cells]);
	TArray& fw = dynamic_cast<TArray&>(_spalartAllmarasFields.fw[cells]);
	TArray& ft2 = dynamic_cast<TArray&>(_spalartAllmarasFields.ft2[cells]);


	const T kappa = _options["kappa"];
	const T Cv1 = _options["Cv1"];
	const T Cw2 = _options["Cw2"];
	const T Cw3 = _options["Cw3"];
	const T Ct3 = _options["Ct3"];
	const T Ct4 = _options["Ct4"];
	
	for (int i=0; i<nCells; i++)
	  {
	    chi[i] = nuTilde[i]/(nu[i]);
	    fv1[i] = pow(chi[i],3)/(pow(chi[i],3) + pow(Cv1,3));
	    fv2[i] = 1 - chi[i]/(1+chi[i]*fv1[i]);
	    Stilde[i] = vorticityMagnitude[i] + nuTilde[i]/(pow(kappa,2)*pow(wallDistance[i],2))*fv2[i];
	    r[i] = min(nuTilde[i]/(pow(kappa,2)*Stilde[i]*pow(wallDistance[i],2)),10);
	    g[i] = r[i] + Cw2*(pow(r[i],6) - r[i]);
	    fw[i] = g[i]*pow(((1+pow(Cw3,6))/(pow(g,6) + pow(Cw3,6))),1./6.);
	    ft2[i] = Ct3*exp(-Ct4*pow(chi[i],2));
	  }
      }
    
  }

  void ComputeWallDistance()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	TArray& wallDistance = dynamic_cast<TArray&>(_spalartAllmarasFields.wallDistance[cells]);
	const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
	
	for (int i=0; i<nCells; i++)
	  {
	    wallDistance[i] = cellCentroid[i][1]; 
	  }
      }
  }

  void ComputeSourceTerm1()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {

	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();

	TArray& sourceTerm1 = dynamic_cast<TArray&>(_spalartAllmarasFields.sourceTerm1[cells]);
	const TArray& nuTilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.nuTilde[cells]);
	const TArray& density = dynamic_cast<const TArray&>(_spalartAllmarasFields.density[cells]);
	const TArray& STilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.STilde[cells]);
	const TArray& ft2 = dynamic_cast<const TArray&>(_spalartAllmarasFields.ft2[cells]);
       
	const T Cb1 = _options["Cb1"];

	for (int i=0; i<nCells; i++)
	  {
	    sourceTerm1[i] = Cb1*density[i]*(1-ft2[i])*STilde[i]*nuTilde[i];
	  }
      }
  }

  void ComputeSourceTerm2()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    
      {

	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();

	TArray& sourceTerm2 = dynamic_cast<TArray&>(_spalartAllmarasFields.sourceTerm2[cells]);
	const TArray& nuTilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.nuTilde[cells]);
	const TArray& density = dynamic_cast<const TArray&>(_spalartAllmarasFields.density[cells]);
	const TArray& STilde = dynamic_cast<const TArray&>(_spalartAllmarasFields.STilde[cells]);
	const TArray& ft2 = dynamic_cast<const TArray&>(_spalartAllmarasFields.ft2[cells]);
	const TArray& fw = dynamic_cast<const TArray&>(_spalartAllmarasFields.fw[cells]);
	const TArray& wallDistance = dynamic_cast<const TArray&>(_spalartAllmarasFields.wallDistance[cells]);

	const T Cb1 = _options["Cb1"];
	const T Cw1 = _options["Cw1"];

	for (int i=0; i<nCells; i++)
	  {
	    sourceTerm2[i] = -density[i]*((Cw1*fw[i] - Cb1/pow(kappa,2)*ft2)*nuTilde[i]*pow(1/wallDistance[i],2)); 
	  }
      }
  }

  void linearize(LinearSystem& ls)
  {
    _temperatureGradientModel.compute();
    ComputeDiffusivityTerm1();
    ComputeDiffusivityTerm2();
    ComputeSpalartAllmarasClosureFunctions();
    ComputeSourceTerm1();
    ComputeSourceTerm2();

    DiscrList discretizations;
   
    shared_ptr<Discretization>
      ddTerm1(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _spalartAllmarasFields.nuTilde,
	  _spalartAllmarasFields.diffusivityTerm1,
	  _spalartAllmarasFields.nuTildeGradient));
    discretizations.push_back(ddTerm1);

    shared_ptr<Discretization>
      ddTerm2(new DiffusionDiscretizationSpalartAllmarasSpecific<T,T,T>
	 (_meshes,_geomFields,
	  _spalartAllmarasFields.nuTilde,
	  _spalartAllmarasFields.diffusivityTerm2,
	  _spalartAllmarasFields.nuTildeGradient));
    discretizations.push_back(ddTerm2);

   
    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _spalartAllmarasFields.nuTilde,
	  _spalartAllmarasFields.convectionFlux,
	  _spalartAllmarasFields.zero,
	  _spalartAllmarasFields.nuTildeGradient));
    discretizations.push_back(cd);
    
    
    shared_ptr<Discretization>
      sdTerm1(new SourceDiscretizationSpalartAllmarasTerm1<T>
	      (_meshes, 
	       _geomFields, 
	       _spalartAllmarasFields.nuTilde,
	       _spalartAllmarasFields.sourceTerm1));
    discretizations.push_back(sdTerm1);
    
    shared_ptr<Discretization>
      sdTerm2(new SourceDiscretizationSpalartAllmarasTerm2<T>
	      (_meshes, 
	       _geomFields, 
	       _spalartAllmarasFields.nuTilde,
	       _spalartAllmarasFields.sourceTerm2));
    discretizations.push_back(sdTerm2);


    /*
    if (_options.transient)
      {
	shared_ptr<Discretization>
	  td(new TimeDerivativeDiscretization<T, T, T>
	     (_meshes, _geomFields, 
	      _thermalFields.temperature, 
	      _thermalFields.temperatureN1,
	      _thermalFields.temperatureN2,
	      _thermalFields.specificHeat,
	      _options["timeStep"]));
	discretizations.push_back(td);
      }
       
    */


    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
			 ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];

	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;

	    const SpalartAllmarasBC<T>& bc = *_bcMap[fg.id];
            

	    GenericBCS<T,T,T> gbc(faces,mesh,
				  _geomFields,
				  _spalartAllmarasFields.nuTilde,
				  _spalartAllmarasFields.nuTildeFlux,
				  ls.getMatrix(), ls.getX(), ls.getB());

	    if (bc.bcType == "SpecifiedNuTilde")
	      {
		FloatValEvaluator<T>
		  bNu(bc.getVal("specifiedNuTilde"),faces);
		if (_spalartAllmarasFields.convectionFlux.hasArray(faces))
		  {
		    const TArray& convectingFlux =
		      dynamic_cast<const TArray&>
		      (_spalartAllmarasFields.convectionFlux[faces]);
		    const int nFaces = faces.getCount();
                                
		    for(int f=0; f<nFaces; f++)
		      {
			if (convectingFlux[f] > 0.)
			  {
			    gbc.applyExtrapolationBC(f);
			  }
			else
			  {
			    gbc.applyDirichletBC(f,bNu[f]);
			  }
		      }
		  }
		else
		  gbc.applyDirichletBC(bNu);
	      }
	    else if (bc.bcType == "SpecifiedNuTildeFlux")
	      {
		FloatValEvaluator<T>
		  bNuTildeFlux(bc.getVal("specifiedNuTildeFlux"),faces);
                    
		const int nFaces = faces.getCount();
                                
		for(int f=0; f<nFaces; f++)
		  {                        
		    gbc.applyNeumannBC(f, bNuTildeFlux[f]);
		  }                              
	      }
	    else if (bc.bcType == "Symmetry")
	      {
		T zeroFlux(NumTypeTraits<T>::getZero());
		gbc.applyNeumannBC(zeroFlux);
	      }
	    else
	      throw CException(bc.bcType + " not implemented for Spalart Allmaras Model");
	  }

	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    GenericBCS<T,T,T> gbc(faces,mesh,
				  _geomFields,
				  _spalartAllmarasFields.nuTilde,
				  _spalartAllmarasFields.nuTildeFlux,
				  ls.getMatrix(), ls.getX(), ls.getB());
	    gbc.applyInterfaceBC();
	  }
      }

    DiscrList discretizations2;
    shared_ptr<Discretization>
      ud(new Underrelaxer<T,T,T>
	 (_meshes,_spalartAllmarasFields.nuTilde,_options["nuTildeURF"]));

    discretizations2.push_back(ud);

    linearizer.linearize(discretizations2,_meshes,ls.getMatrix(),
			 ls.getX(), ls.getB());

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
	if (*rNorm < _options.absoluteTolerance ||
	    *normRatio < _options.relativeTolerance)
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
    
      MultiField::ArrayIndex tIndex(&_spalartAllmarasFields.nuTilde,&cells);

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
	  fprintf(matFile,"%d %d %lf\n", i+1, i+1, tDiag[i]);
	  for(int jp=row[i]; jp<row[i+1]; jp++)
	    {
	      const int j = col[jp];
	      if (j<nCells)
		fprintf(matFile,"%d %d %lf\n", i+1, j+1, tCoeff[jp]);
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

#include "SpalartAllmarasTopologyOptimization.h"
  

 private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  SpalartAllmarasFields _spalartAllmarasFields;

  SpalartAllmarasBCMap _bcMap;
  SpalartAllmarasVCMap _vcMap;
  SpalartAllmarasModelOptions<T> _options;
  GradientModel<T> _spalartAllmarasGradientModel;
  
  MFRPtr _initialNorm;
  int _niters;
  T _objFunctionSAModel;
  T _objFunctionFlowModel;
  int _nDesignVariables = 0;
};

template<class T>
SpalartAllmarasModel<T>::SpalartAllmarasModel(const GeomFields& geomFields,
					      SpalartAllmarasFields& spalartAllmarasFields,	      
					      const MeshList& meshes) :
Model(meshes),
  _impl(new Impl(geomFields,spalartAllmarasFields,meshes))
{
  logCtor();
}


template<class T>
SpalartAllmarasModel<T>::~SpalartAllmarasModel()
{
  logDtor();
}

template<class T>
void
SpalartAllmarasModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename SpalartAllmarasModel<T>::SpalartAllmarasBCMap&
SpalartAllmarasModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename SpalartAllmarasModel<T>::SpalartAllmarasVCMap&
SpalartAllmarasModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
SpalartAllmarasBC<T>&
SpalartAllmarasModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
SpalartAllmarasModelOptions<T>&
SpalartAllmarasModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
SpalartAllmarasModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))

template<class T>
void
SpalartAllmarasModel<T>::dumpMatrix(const string fileBase)
{
  _impl->dumpMatrix(fileBase);
}
#endif

#include "SpalartAllmarasTopologyOptimizationImpl.h"
