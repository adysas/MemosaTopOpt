#ifndef _FLOWMODELTURBULENCEIMPL_H_
#define _FLOWMODELTURBULENCEIMPL_H_

template<class T>
void
FlowModel<T>::ComputeVorticityMagnitude()
{
  _impl->ComputeVorticityMagnitude();
}

template<class T>
void
FlowModel<T>::ComputeTotalViscosity()
{
  _impl->ComputeTotalViscosity();
}

template<class T>
void
FlowModel<T>::ComputeShearStress(const Mesh& mesh, const int faceGroupId)
{
  _impl->ComputeShearStress(mesh, faceGroupId);
}

template<class T>
void
FlowModel<T>::ComputeWallStress(const Mesh& mesh, const int faceGroupId)
{
  _impl->ComputeWallStress(mesh, faceGroupId);
}

#endif
