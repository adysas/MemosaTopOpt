#ifndef _THERMALFILTERFUNCTIONPROTOTYPES_H_
#define _THERMALFILTERFUNCTIONPROTOTYPES_H_

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
  void buildNeighborsList(const T rMin); 
  void applyFilter(); 
#endif

#endif
