
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
template<class T>
void
SpalartAllmarasModel<T>::buildNeighborsList(const T rMin)
{
  _impl->buildNeighborsList(rMin);
}
template<class T>
void
SpalartAllmarasModel<T>::applyFilter()
{
  _impl->applyFilter();
}

#endif
template<class T>
void
SpalartAllmarasModel<T>::initResidualFieldSpalartAllmaras()
{
  _impl->initResidualFieldSpalartAllmaras();
}


/////////Copying fields and values ///////////////////
template<class T>
void
SpalartAllmarasModel<T>::copyNuTildeField(const Field& nuTildeDoubleField, const MeshList& meshesDouble)
{
  _impl->copyNuTildeField(nuTildeDoubleField, meshesDouble);
}

template<class T>
void
SpalartAllmarasModel<T>::copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble)
{
  _impl->copyBetaField(betaDoubleField, meshesDouble);
}

template<class T>
void
SpalartAllmarasModel<T>::interpolateSourceFieldSIMP()
{
  _impl->interpolateSourceFieldSIMP();
}

template<class T>
void
SpalartAllmarasModel<T>::updateBetaFieldGhostCells()
{
  _impl->updateBetaFieldGhostCells();
}

template<class T>
void
SpalartAllmarasModel<T>::ComputeSpalartAllmarasResidual()
{
  _impl->ComputeSpalartAllmarasResidual();
}

template<class T>
void
SpalartAllmarasModel<T>::setObjectiveFunctionSA(const T objFunction)
{
  _impl->setObjectiveFunctionSA(objFunction);
}


template<class T>
void
SpalartAllmarasModel<T>::printResidualField()
{
  _impl->printResidualField();
}


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
template<class T>
void
SpalartAllmarasModel<T>::dumpBetaSensistivityFields()
{
  _impl->dumpBetaSensistivityFields();
}

#endif

//Functions used to developing code
template<class T>
void
SpalartAllmarasModel<T>::printNuTildeField()
{
  _impl->printNuTildeField();
}

template<class T>
void
SpalartAllmarasModel<T>::printBetaField()
{
  _impl->printBetaField();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class T>
void
SpalartAllmarasModel<T>::copyMomentumResidual(const Field& momentumRes)
{
  _impl->copyMomentumResidual(momentumRes);
}

template<class T>
void
SpalartAllmarasModel<T>::copyContinuityResidual(const Field& continuityRes)
{
  _impl->copyContinuityResidual(continuityRes);
}

template<class T>
void
SpalartAllmarasModel<T>::copyWallDistanceModelResidual(const Field& wallDistanceRes)
{
  _impl->copyWallDistanceModelResidual(wallDistanceRes);
}

template<class T>
void
SpalartAllmarasModel<T>::copyMassFlux(const Field& massFlux)
{
  _impl->copyMassFlux(massFlux);
}


template<class T>
void
SpalartAllmarasModel<T>::copyWallDistanceField(const Field& wallDistanceField)
{
  _impl->copyWallDistanceField(wallDistanceField);
}


template<class T>
void
SpalartAllmarasModel<T>::copyVorticityMagnitudeField(const Field& vorticityMagnitudeFlowField)
{
  _impl->copyVorticityMagnitudeField(vorticityMagnitudeFlowField);
}


template<class T>
void
SpalartAllmarasModel<T>::setObjectiveFunctionFlow(const VectorT3 objFunction, const int component)
{
  _impl->setObjectiveFunctionFlow(objFunction, component);
}


#if ( defined(USING_ATYPE_RAPID) )

template<class T>
void
SpalartAllmarasModel<T>::setIndependentStateVariablesFlowTurbulent()
{
  _impl->setIndependentStateVariablesFlowTurbulent();
}
/*
Don't need this, since we copy from flow model
template<class T>
void
SpalartAllmarasModel<T>::setIndependentDesignVariablesFlowTurbulent()
{
  _impl->setIndependentDesignVariablesFlowTurbulent();
}
*/
template<class T>
void
SpalartAllmarasModel<T>::WriteResidualGradientsToMatlabFlowTurbulent()
{
  _impl->WriteResidualGradientsToMatlabFlowTurbulent();
}

template<class T>
void
SpalartAllmarasModel<T>::WriteResidualGradientsToPetscFlowTurbulent()
{
  _impl->WriteResidualGradientsToPetscFlowTurbulent();
}


template<class T>
void
SpalartAllmarasModel<T>::WriteObjectiveFunctionGradientToPetscFlowTurbulent()
{
  _impl->WriteObjectiveFunctionGradientToPetscFlowTurbulent();
}

template<class T>
void
SpalartAllmarasModel<T>::WriteObjectiveFunctionGradientToMatlabFlowTurbulent()
{
  _impl->WriteObjectiveFunctionGradientToMatlabFlowTurbulent();
}

#endif

template<class T>
void
SpalartAllmarasModel<T>::copyBetaField(const Field& betaField)
{
  _impl->copyBetaField(betaField);
}
