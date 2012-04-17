/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Victorian Partnership for Advanced Computing (VPAC) Ltd, Australia
** (C) 2003-2005 All Rights Reserved
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**
** Role:
**	Defines any header info, such as new structures, needed by this plugin
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_plugins_EulerDeform_EulerDeform_h__
#define __Underworld_plugins_EulerDeform_EulerDeform_h__

	extern Name	EULERDEFORM_PLUGIN_TAG;


	Index EulerDeform_Register( PluginsManager* pluginsMgr );

	void* EulerDeform_DefaultNew( Name name );

	void EulerDeform_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data );

	void EulerDeform_Build( void* component, void* data );

	void EulerDeform_Destroy( void* component, void* data );

	Variable* EulerDeform_RegisterLocalNodeCoordsAsVariables( EulerDeform_System* sys, void* _variable_Register, 
								  Variable** variableList );

	void EulerDeform_IntegrationSetup( void* _timeIntegrator, void* _velField );

	Bool EulerDeform_TimeDeriv( void* crdAdvector, Index arrayInd, double* timeDeriv );

	void EulerDeform_Remesh( TimeIntegrand* crdAdvector, EulerDeform_Context* edCtx );

	void EulerDeform_WrapTopSurface( EulerDeform_System* sys, double** oldCrds );

	void EulerDeform_WrapBottomSurface( EulerDeform_System* sys, double** oldCrds );

	void EulerDeform_WrapLeftSurface( EulerDeform_System* sys, double** oldCrds );
        void EulerDeform_Advection_Correction(void* sle, void* data);
        IndexSet* EulerDeform_CreateStaticSet(EulerDeform_System* sys);
        void EulerDeform_Advection_Correction(void* sle, void* data);
        void EulerDeform_InternalLoop(EulerDeform_System* sys, Grid* grm,
                                      double** oldCrds, unsigned* ijk,
                                      unsigned curDim, const bool &top);
        Bool _EulerDeform_QuadYInterp(double** crds, const double* pnt, double* val);
        Bool _EulerDeform_LineInterp(double** crds, const double &pnt,
                                     unsigned fromDim, unsigned toDim, double* val);
        void EulerDeform_Remesh(TimeIntegrand* crdAdvector, EulerDeform_Context* edCtx);
        void EulerDeform_Remesh_Corner(Mesh *mesh, const int &corner, const int &inside,
                                       const double &side_coord, const int &boundary_dim,
                                       const int &height_dim, const int &tangent_dim);
        void EulerDeform_WrapSurface(EulerDeform_System* sys,double** oldCrds,
                                     const bool &top);
        void EulerDeform_FloatLeftTop(EulerDeform_System* sys, Grid *grid,
                                      double** crds);
        void EulerDeform_FloatRightTop(EulerDeform_System* sys, Grid *grid,
                                       double** crds);
        Bool EulerDeform_TimeDeriv( void* crdAdvector, Index arrayInd, double* timeDeriv );
        void EulerDeform_IntegrationSetup( void* _timeIntegrator, void* context );

        const Type Underworld_EulerDeform_Type = "EulerDeform";
#endif
