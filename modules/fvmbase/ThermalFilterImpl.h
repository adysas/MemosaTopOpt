#ifndef _THERMALFILTERIMPL_H_
#define _THERMALFILTERIMPL_H_


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
template<class T>
void
ThermalModel<T>::buildNeighborsList(const T rMin)
{
   _impl->buildNeighborsList(rMin);
}

template<class T>
void
ThermalModel<T>::applyFilter()
{
   _impl->applyFilter();
}
#endif

#endif
