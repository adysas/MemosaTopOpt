
template<class T>
void
FlowModel<T>::interpolateMassFluxApproximate()
{
  return _impl->interpolateMassFluxApproximate();
}

template<class T>
void
FlowModel<T>::interpolateMassFluxInteriorFacesSetFaceCellPressureOnPressureBoundarySetMassFluxOtherBoundaries()
{
  return _impl->interpolateMassFluxInteriorFacesSetFaceCellPressureOnPressureBoundarySetMassFluxOtherBoundaries();
}

template<class T>
void
FlowModel<T>::updateBetaFieldGhostCells()
{
  return _impl->updateBetaFieldGhostCells();
}

template<class T>
void
FlowModel<T>::initMomApField()
{
  return _impl->initMomApField();
}

template<class T>
void
FlowModel<T>::initpreviousVelocityField()
{
  return _impl->initpreviousVelocityField();
}

template<class T>
void
FlowModel<T>::interpolateSourceFieldSIMP()
{
  return _impl->interpolateSourceFieldSIMP();
}

template<class T>
void
FlowModel<T>::interpolateSourceFieldSIMPRamp()
{
  return _impl->interpolateSourceFieldSIMPRamp();
}


template<class T>
void
FlowModel<T>::interpolateViscosityFieldSIMP()
{
  return _impl->interpolateViscosityFieldSIMP();
}

template<class T>
void
FlowModel<T>::initResidualFieldsFlow()
{
  return _impl->initResidualFieldsFlow();
}

template<class T>
void
FlowModel<T>::copyVelocityPressureField(const Field& velocityDoubleField, const Field& pressureDoubleField, const MeshList& meshesDouble)
{
  return _impl->copyVelocityPressureField(velocityDoubleField, pressureDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copyVelocityBoundaryFace(const Field& velocityComponentDoubleField, Field& velocityComponentField, const StorageSite& facesDouble, const StorageSite& faces)
{
  return _impl->copyVelocityBoundaryFace(velocityComponentDoubleField, velocityComponentField, facesDouble, faces);
}

template<class T>
void
FlowModel<T>::copyMomApField(const Field& momApDoubleField, const MeshList& meshesDouble)
{
  return _impl->copyMomApField(momApDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copyMomentumResidual(const Field& momentumRes, const MeshList& meshesFrom)
{
  return _impl->copyMomentumResidual(momentumRes, meshesFrom);
}

template<class T>
void
FlowModel<T>::copyContinuityResidual(const Field& continuityRes, const MeshList& meshesFrom)
{
  return _impl->copyContinuityResidual(continuityRes, meshesFrom);
}

template<class T>
void
FlowModel<T>::copypreviousVelocityField(const Field& previousVelocityDoubleField, const MeshList& meshesDouble)
{
  return _impl->copypreviousVelocityField(previousVelocityDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copypressureGradientField(const Field& pressureGradientDoubleField, const MeshList& meshesDouble)
{
  return _impl->copypressureGradientField(pressureGradientDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copyfacePressureField(const Field& facePressureDoubleField, const MeshList& meshesDouble)
{
  return _impl->copyfacePressureField(facePressureDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copyMassFluxField(const Field& massFluxFaceDoubleField, const MeshList& meshesDouble)
{
  return _impl->copyMassFluxField(massFluxFaceDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble)
{
  return _impl->copyBetaField(betaDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::copySourceField(const Field& sourceDoubleField, const MeshList& meshesDouble)
{
  return _impl->copySourceField(sourceDoubleField, meshesDouble);
}

template<class T>
void
FlowModel<T>::ComputeMomApCoefficients()
{
  return _impl->ComputeMomApCoefficients();
}

template<class T>
void
FlowModel<T>::updateFacePressure()
{
  return _impl->updateFacePressure();
}

template<class T>
void
FlowModel<T>::updatePressureGradient()
{
  return _impl->updatePressureGradient();
}

template<class T>
void
FlowModel<T>::ComputeContinuityResidual()
{
  return _impl->ComputeContinuityResidual();
}

template<class T>
void
FlowModel<T>::ComputeMomentumResidual()
{
  return _impl->ComputeMomentumResidual();
}

template<class T>
void
FlowModel<T>::ComputeMomentumResidualSaveMomAp()
{
  return _impl->ComputeMomentumResidualSaveMomAp();
}

template<class T>
void
FlowModel<T>::setObjectiveFunctionFlow(const Vector<T,3> objFunction, const int component)
{
  return _impl->setObjectiveFunctionFlow(objFunction, component); 
}

template<class T>
void
FlowModel<T>::setObjectiveFunctionFlow(const T objFunction)
{
  return _impl->setObjectiveFunctionFlow(objFunction); 
}

#if ( defined(USING_ATYPE_RAPID) )
template<class T>
void
FlowModel<T>::setIndependentStateVariablesFlow()
{
  return _impl->setIndependentStateVariablesFlow();
}

template<class T>
void
FlowModel<T>::setIndependentDesignVariablesFlow()
{
  return _impl->setIndependentDesignVariablesFlow();
}

template<class T>
void
FlowModel<T>::setViscosityAsIndependentVariable()
{
  return _impl->setViscosityAsIndependentVariable();
}

template<class T>
void
FlowModel<T>::setSourceAsIndependentVariable()
{
  return _impl->setSourceAsIndependentVariable();
}

template<class T>
void
FlowModel<T>::setSourceAsIndependentVariable(const int cellId)
{
  return _impl->setSourceAsIndependentVariable(cellId);
}

template<class T>
void
FlowModel<T>::setSourceFieldAsIndependentVariables()
{
  return _impl->setSourceFieldAsIndependentVariables();
}

template<class T>
void
FlowModel<T>::setDensityAsIndependentVariable()
{
  return _impl->setDensityAsIndependentVariable();
}

template<class T>
void
FlowModel<T>::WriteResidualGradientsToPetscFlow()
{
  return _impl->WriteResidualGradientsToPetscFlow();
}

template<class T>
void
FlowModel<T>::WriteResidualGradientsToMatlabFlow()
{
  return _impl->WriteResidualGradientsToMatlabFlow();
}


template<class T>
void
FlowModel<T>::WriteObjectiveFunctionGradientToPetscFlow()
{
  return _impl->WriteObjectiveFunctionGradientToPetscFlow();
}

template<class T>
void
FlowModel<T>::WritePressureDropGradientToPetsc(const int faceId, const int direction)
{
  return _impl->WritePressureDropGradientToPetsc(faceId, direction);
}

template<class T>
void
FlowModel<T>::WriteObjectiveFunctionGradientToMatlabFlow()
{
  return _impl->WriteObjectiveFunctionGradientToMatlabFlow();
}
#endif

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID) )

template<class T>
void
FlowModel<T>::dumpBetaSensistivityFields()
{
  return _impl->dumpBetaSensistivityFields();
}

template<class T>
void
FlowModel<T>::dumpVelocityVector(const string fileBase)
{
  return _impl->dumpVelocityVector(fileBase);
}

template<class T>
void
FlowModel<T>::ComputeVelocityGradient()
{
  return _impl->ComputeVelocityGradient();
}
#endif

// Functions for developing the code

template<class T>
void
FlowModel<T>::printPressureIntegrals()
{
  return _impl->printPressureIntegrals();
}

template<class T>
void
FlowModel<T>::printMomentumFluxIntegrals()
{
  return _impl->printMomentumFluxIntegrals();
}

template<class T>
void
FlowModel<T>::printMassFluxIntegrals()
{
  return _impl->printMassFluxIntegrals();
}

template<class T>
void
FlowModel<T>::printCellPressure(string fileName)
{
  return _impl->printCellPressure(fileName);
}

template<class T>
void
FlowModel<T>::printCellVelocity(string fileName)
{
  return _impl->printCellVelocity(fileName);
}

template<class T>
void
FlowModel<T>::printMomApField()
{
  return _impl->printMomApField();
}

template<class T>
void
FlowModel<T>::printMomApField(string fileName)
{
  return _impl->printMomApField(fileName);
}

template<class T>
void
FlowModel<T>::printpreviousVelocityField(string fileName)
{
  return _impl->printpreviousVelocityField(fileName);
}

template<class T>
void
FlowModel<T>::printpressureGradientField()
{
  return _impl->printpressureGradientField();
}

template<class T>
void
FlowModel<T>::printpressureGradientField(string fileName)
{
  return _impl->printpressureGradientField(fileName);
}

template<class T>
void
FlowModel<T>::printfacePressureField()
{
  return _impl->printfacePressureField();
}

template<class T>
void
FlowModel<T>::printfacePressureField(string fileName)
{
  return _impl->printfacePressureField(fileName);
}

template<class T>
void
FlowModel<T>::printMassFlux()
{
  return _impl->printMassFlux();
}

template<class T>
void
FlowModel<T>::printMassFlux(string fileName)
{
  return _impl->printMassFlux(fileName);
}

template<class T>
void
FlowModel<T>::printMomentumResidual()
{
  return _impl->printMomentumResidual();
}

template<class T>
void
FlowModel<T>::printMomentumResidual(string fileName)
{
  return _impl->printMomentumResidual(fileName);
}

template<class T>
void
FlowModel<T>::printContinuityResidual()
{
  return _impl->printContinuityResidual();
}

template<class T>
void
FlowModel<T>::printContinuityResidual(string fileName)
{
  return _impl->printContinuityResidual(fileName);
}

template<class T>
void
FlowModel<T>::printSourceField()
{
  return _impl->printSourceField();
}

#if ( defined(USING_ATYPE_RAPID) )

template<class T>
void
FlowModel<T>::printIndependentStateVariables()
{
  return _impl->printIndependentStateVariables();
}

template<class T>
void
FlowModel<T>::printIndependentDesignVariables()
{
  return _impl->printIndependentDesignVariables();
}

#endif


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


//////////////////////////////
#if ( defined(USING_ATYPE_RAPID) )
template<class T>
void
FlowModel<T>::setIndependentStateVariablesThermoFluidLaminar()
{
  return _impl->setIndependentStateVariablesThermoFluidLaminar();
}
template<class T>
void
FlowModel<T>::setIndependentDesignVariablesThermoFluidLaminar()
{
  return _impl->setIndependentDesignVariablesThermoFluidLaminar();
}
template<class T>
void
FlowModel<T>::setIndependentStateVariablesTurbulentFlow()
{
  return _impl->setIndependentStateVariablesTurbulentFlow();
}
template<class T>
void
FlowModel<T>::setIndependentDesignVariablesTurbulentFlow()
{
  return _impl->setIndependentDesignVariablesTurbulentFlow();
}
#endif

template<class T>
void
FlowModel<T>::copyEddyViscosityField(const Field& eddyViscosityField)
{
  return _impl->copyEddyViscosityField(eddyViscosityField);
}

template<class T>
Vector<T,3>
  FlowModel<T>::getPressureGradientIntegral()
{
  return _impl->getPressureGradientIntegral();
}

template<class T>
Vector<T,3> 
  FlowModel<T>::getTotalPressureLossesBetweenInletOutlet(const Mesh& mesh, const int inlet, const int outlet)
{
  return _impl->getTotalPressureLossesBetweenInletOutlet(mesh, inlet, outlet);
}

template<class T>
Vector<T,3> 
  FlowModel<T>::getDynamicPressureIntegral(const Mesh& mesh, const int faceId)
{
  return _impl->getDynamicPressureIntegral(mesh, faceId);
}


/*
template<class T>
Vector<T,3> 
  FlowModel<T>::getDynamicPressureIntegralOutlet(const Mesh& mesh, const int outlet)
{
  return _impl->getDynamicPressureIntegralOutlet(mesh, outlet);
}

template<class T>
Vector<T,3> 
  FlowModel<T>::getDynamicPressureIntegralInlet(const Mesh& mesh, const int inlet)
{
  return _impl->getDynamicPressureIntegralInlet(mesh, inlet);
}
*/
