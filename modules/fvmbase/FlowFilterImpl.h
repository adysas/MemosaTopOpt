#ifndef _FLOWFILTERIMPL_H_
#define _FLOWFILTERIMPL_H_

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
template<class T>
void
FlowModel<T>::buildNeighborsList(const T rMin)
{
   _impl->buildNeighborsList(rMin);
}

template<class T>
void
FlowModel<T>::applyFilter()
{
   _impl->applyFilter();
}
#endif

#endif
