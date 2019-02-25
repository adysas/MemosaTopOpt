//////////////////////////Thermal topology optimization///////////////////////////////////////////////////////////

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
void buildNeighborsList(const T rMin); 
void applyFilter(); 
#endif


void initResidualFieldThermal();
void updateBetaFieldGhostCells();
void copyTemperatureField(const Field& tempDoubleField, const MeshList& meshesDouble);
void copyBetaField(const Field& betaDoubleField, const MeshList& meshesDouble);
void updateSourceFieldType1();
void updateSourceFieldType2();
void updateSourceFieldType3();
void interpolateThermalConductivityFieldSIMP();
void interpolateThermalConductivityFieldSIMPRamp();
void ComputeTemperatureResidual();
T getCellTemperature(const Mesh& mesh, const int cellId);
T getAverageTemperature(const Mesh& mesh);
void setObjectiveFunctionThermal(const T objFunction);

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC) || defined(USING_ATYPE_RAPID))
void dumpBetaSensistivityFields();
#endif

#if (defined(USING_ATYPE_RAPID))
void setIndependentStateVariablesThermal();
void setIndependentDesignVariablesThermal();
void WriteResidualGradientsToPetscThermal();
void WriteObjectiveFunctionGradientToPetscThermal();
#endif

void initializeBetaField();
void printResidualField();
void printTemperatureField();
void printBetaField();
void printConductivityField();
void printBoundaryCellNumbers(const Mesh& mesh, const int faceGroupId);


/////////////////////////////////////////////////////////////////////////////////

T getDeltaEnthalpyIntegralThermoFluid(const Mesh& mesh, const int faceGroupId1, const int faceGroupId2);
T getBulkTemperature(const Mesh& mesh, const int faceGroupId);
void copyMomentumResidual(const Field& momentumRes);
void copyContinuityResidual(const Field& continuityRes);
void copyMassFlux(const Field& massFlux);
void setObjectiveFunctionFlow(const VectorT3 objFunction, const int direction);

#if ( defined(USING_ATYPE_RAPID) )
void setIndependentStateVariablesThermoFluidLaminar();
void setIndependentDesignVariablesThermoFluidLaminar();
void WriteResidualGradientsToPetscThermoFluidLaminar();
void WriteResidualGradientsToMatlabThermoFluidLaminar();
void WriteObjectiveFunctionGradientToPetscThermoFluidLaminar();
#endif
////////////////////////////////////////////////////////////////////////////////
