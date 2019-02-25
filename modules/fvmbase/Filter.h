#ifndef _FILTER_H_
#define _FILTER_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include <vector>
#include <math.h>

template<class X>
class Filter
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef map<int, T_Scalar> NeighborListWeight;
  //typedef map<int, double>::iterator Iter;
  //typedef map<int,double>::const_iterator CIter;
  
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef std::vector<NeighborListWeight> NLWArray;
  typedef std::vector<T_Scalar> AvgnlwArray;
  
 Filter(const GeomFields& geomFields):
  _geomFields(geomFields)
    {}
  
  void buildNeighborListWeight(const Mesh& mesh, const T_Scalar rMin)
  {
 
    const StorageSite& cells = mesh.getCells();
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
    
    const int nCells = cells.getSelfCount();
    //cout << "\n Number of cells in Filter class - " << nCells << endl; 

    //shared_ptr<TArray> abcd(new TArray(nCells));
    //*abcd = X(0.);
    
    _nlwArray.clear();
    _avgnlwArray.clear();

    for (int c=0; c<nCells; c++)
      {
	NeighborListWeight neighborWeight;
	T_Scalar avgnlw(0.0);
	for (int i=0; i<nCells; i++)
	  {
	    VectorT3 ds=cellCentroid[i]-cellCentroid[c];
	    T_Scalar dsMag = mag(ds);
	    if (dsMag < rMin)
	      {		
		T_Scalar delta = rMin - dsMag;
		neighborWeight[i] = max((T_Scalar)0.0, delta);
		avgnlw += delta;
	      }
	  }
	//(*abcd)[c] = 3.14;
	_nlwArray.push_back(neighborWeight);

	//avgnlw = avgnlw/_nlwArray.size();
	_avgnlwArray.push_back(avgnlw);
      }
    
    /*
    for (int c=0; c<nCells; c++)
      {
	NeighborListWeight neighborWeight;
	neighborWeight = _nlwArray[c];
	cout << endl;
	cout << "Cell# " << c <<"->  WeightSum = " << _avgnlwArray[c] <<"; ";
	foreach(const typename NeighborListWeight::value_type& ii, neighborWeight)
	  {
	    cout << ii.first <<":" << ii.second <<", ";
	  }
	cout << endl;
	}*/
    
  }
  
        
  void applyFilter(const Mesh& mesh, Field& sensitivity, Field& beta)
  {
    const StorageSite& cells = mesh.getCells();
    XArray& sensitivityArray = dynamic_cast<XArray&>(sensitivity[cells]);
    XArray& betaArray = dynamic_cast<XArray&>(beta[cells]);

    const XArray& cellVolume = dynamic_cast<const XArray&>(_geomFields.volume[cells]);

    Field sensitivityCopy("Trial");
    shared_ptr<TArray> zeroCell(new TArray(cells.getSelfCount()));
    *zeroCell = X(0.0);
    sensitivityCopy.addArray(cells,zeroCell);    
    XArray& sensitivityCopyArray = dynamic_cast<XArray&>(sensitivityCopy[cells]);
   
    const int nCells = cells.getSelfCount();
    for ( int c=0; c<nCells; c++)
      {
	sensitivityCopyArray[c] = sensitivityArray[c];
      }

    
    for (int c=0; c<nCells; c++)
      {
	NeighborListWeight neighborWeight;
	neighborWeight = _nlwArray[c];
	
	T_Scalar tempVar(0.0);
	foreach(const typename NeighborListWeight::value_type& ii, neighborWeight)
	  {
	    tempVar += ii.second*betaArray[ii.first]*sensitivityCopyArray[ii.first];
	  }
	sensitivityArray[c] = tempVar/_avgnlwArray[c]/(max(0.001,betaArray[c]));
      }
    

    /*
    for (int c=0; c<nCells; c++)
      {
	NeighborListWeight neighborWeight;
	neighborWeight = _nlwArray[c];
	T_Scalar tempVar(0.0);
	foreach(const typename NeighborListWeight::value_type& ii, neighborWeight)
	  {
	    tempVar += ii.second*betaArray[ii.first]*cellVolume[ii.first]*sensitivityCopyArray[ii.first];
	  }
	sensitivityArray[c] = tempVar/_avgnlwArray[c]/(max(0.001,betaArray[c]*cellVolume[c]));
      }
    */
    /*
    for ( int c=0; c<nCells; c++)
      {
	cout << c << "th cell filtered sensitivity = " <<  sensitivityArray[c] << endl;
      }
    */
  }
  
 private:
  const GeomFields& _geomFields;
  AvgnlwArray _avgnlwArray;
  NLWArray _nlwArray;
  
};

#endif


/*
  void applyFilter(const Mesh& mesh, Field &varField, Field& sensitivity, Field& beta)
  {
    const StorageSite& cells = mesh.getCells();
    //const int nCells = cells.getSelfCount();

    Field& sensitivityCopy = new Field("Trial");
    //const MultiField::ArrayIndex cVarIndex(&varField, &cells);
    shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
    *zeroCell = X(0.0);
    sensitivityCopy.addArray(cells,zeroCell);
    
    // sensitivityCopy.addArray(cVarIndex, varField.getArrayPtr(cells));

    //create a zero field
    //shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
    //*zeroCell = X(0.0);
    //sensitivityCopy.addArray(cells,zeroCell);
    
    XArray& sensCopy = dynamic_cast<XArray&>(sensitivity[cells]);
    const int nCells = cells.getCountLevel1();
    for ( int c=0; c<nCells; c++)
      {
	sensCopy[c] = 3.14;
      }
    
  }
 */

/*
  template<class X>
  class Filter
  {
  public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef map<int, double> NeighborWeight;
  typedef map<int, double>::iterator Iter;
  typedef map<int,double>::const_iterator CIter;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  
  Filter(const GeomFields& geomFields,
  const double rMin):
  _geomFields(geomFields), 
  _rMin(rMin)
  {}
  
  void buildNeighborListWeight(const Mesh& mesh)
  {
 
  const StorageSite& cells = mesh.getCells();
  const VectorT3Array& cellCentroid =
  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
    
  const int nCells = cells.getSelfCount();
  cout << "\n Number of cells in Filter class - " << nCells << endl;  
  for (int c=0; c<nCells; c++)
  {
	
  //VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
  //T_Scalar dsMag = mag(ds);
	
  }
  }
  
  private:
  const GeomFields& _geomFields;
  const double _rMin;
  int _cellId;
  NeighborWeight _neighborWeight;
  
  };
*/
