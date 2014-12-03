double StgFEM_Vrms( FeVariable* velsq, Swarm* swarm );
void StgFEM_InterpolateValue_WithNi( void* _feVariable, Element_LocalIndex lElement_I, double* Ni, double* value );
void StgFEM_InterpolateDerivatives_WithGNx( void* _feVariable, Element_LocalIndex lElement_I, double** GNx, double* value );
