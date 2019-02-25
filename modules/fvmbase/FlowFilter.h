#ifndef _FLOWFILTER_H_
#define _FLOWFILTER_H_

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


#endif
