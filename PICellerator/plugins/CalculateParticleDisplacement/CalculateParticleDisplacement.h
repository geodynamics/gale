
	typedef struct ParticleDisplacementInfo {
		Coord     originalCoord;
		XYZ       displacement;
	} ParticleDisplacementInfo;

	/** The CalculateParticleDisplacementPlugin extends upon the Codelet class */
	#define __CalculateParticleDisplacementPlugin \
		__Codelet \
		\
		ExtensionInfo_Index     particleDisplacementInfo_Handle; \
		SwarmVariable*          particleOriginalCoordSwarmVariable; \
		SwarmVariable*          particleDisplacementSwarmVariable;  \
		OperatorSwarmVariable*  particleDisplacementMagSwarmVariable;

	typedef struct CalculateParticleDisplacementPlugin 
		{ __CalculateParticleDisplacementPlugin } CalculateParticleDisplacementPlugin;

	extern const Type PICellerator_CalculateParticleDisplacement_Type;

	void _PICellerator_CalculateParticleDisplacement_StoreOriginalPos( PICelleratorContext* context );
	void _PICellerator_CalculateParticleDisplacement_UpdateDisplacement( PICelleratorContext* context );
	
	
