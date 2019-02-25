// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLUENTREADER_H_
#define _FLUENTREADER_H_

#include "misc.h"
#include "StorageSite.h"
#include "SchemeReader.h"
#include "Array.h"
#include "Vector.h"
#include "CRConnectivity.h"
#include "Mesh.h"

class OneToOneIndexMap;

typedef map<int, shared_ptr<OneToOneIndexMap> > GhostCellMapsMap;

struct FluentZone 
{
  int ID;
  int iBeg;
  int iEnd;
  int partnerId;
  int threadType;
  string zoneName;
  string zoneType;
  string zoneVars;
};

struct FluentFaceZone : public FluentZone
{
  int leftCellZoneId;
  int rightCellZoneId;
};

struct FluentCellZone : public FluentZone
{
  vector<int> boundaryZoneIds;
  vector<int> interiorZoneIds;
  vector<int> interfaceZoneIds;

  Mesh *mesh;
  GhostCellMapsMap ghostCellMaps;
  shared_ptr<Array<int> > globalToLocalNodeMap;
};

struct FluentFacePairs
{
  FluentFacePairs(const int count_, const int leftID_, const int rightID_,
                  shared_ptr<Array<int> > leftFaces_,
                  shared_ptr<Array<int> > rightFaces_) :
    count(count_),
    leftID(leftID_),
    rightID(rightID_),
    leftFaces(leftFaces_),
    rightFaces(rightFaces_)
  {}

  const int count;
  const int leftID;
  const int rightID;
  shared_ptr<Array<int> > leftFaces;
  shared_ptr<Array<int> > rightFaces;
  
};


typedef map<int,FluentFaceZone*> FaceZonesMap;
typedef map<int,FluentCellZone*> CellZonesMap;
typedef map<int,shared_ptr<FluentFacePairs> > FacePairsMap;

class FluentReader : public SchemeReader
{
public:

  typedef Vector<double,3> Vec3;
  

  
  FluentReader(const string& fileName);
  virtual ~FluentReader();

  void readMesh();
  //void orderCellFacesAndNodes();


  MeshList getMeshList();

  int getNumCells() {return _numCells;}

  string getVars() {return _rpVars;}

  FaceZonesMap& getFaceZones() {return _faceZones;}
  CellZonesMap& getCellZones() {return _cellZones;}
  
protected:
  int _dimension;
  int _numCells;
  int _numFaces;
  int _numNodes;
  int _numBoundaryFaces;

  StorageSite _cells;
  StorageSite _faces;
  StorageSite _nodes;
  
  shared_ptr<CRConnectivity> _faceNodes;
  shared_ptr<CRConnectivity> _faceCells;
  shared_ptr<CRConnectivity> _cellFaces;
  shared_ptr<CRConnectivity> _cellNodes;
  shared_ptr<CRConnectivity> _nodeCells;
  
  FaceZonesMap _faceZones;
  CellZonesMap _cellZones;
  FacePairsMap _facePairs;
  
  Array<Vec3> _coords;
  int _rpVarStringLength;
  string _rpVars;
  map<int,int> _zoneVarStringLength;

  void read(const int pass);
  void readNodes(const int pass, const bool isBinary,
                 const bool isDP, const int id);
  void readCells(const int pass, const bool isBinary, const int id);
  void readFaces(const int pass, const bool isBinary, const int id);
  void readFacePairs(const int pass, const bool isBinary, const int id);

  
  void readVectorData(Array<Vec3>& a,
                      const int iBeg, const int iEnd, const bool isBinary,
                      const bool isDP);

  void buildZones();

  const CRConnectivity& getCellFaces();
  const CRConnectivity& getCellNodes();
  const CRConnectivity& getNodeCells();
  
  // get the cell zone ID for a particular cell
  int getCellZoneID(const int c) const;

  shared_ptr<OneToOneIndexMap>
  getGhostCellMap(const FluentCellZone& cz, const Array<int>& indices);

  shared_ptr<OneToOneIndexMap>
  getCommonNodeMap(const FluentCellZone& cz0, const FluentCellZone& cz1);

  Mesh* createMesh(const int cellZoneID, Array<int>&);

};
#endif
