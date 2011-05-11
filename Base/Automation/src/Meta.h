/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
*/
/** \file
**  Role:
**	Funtions to obtain a components meta information from a dictionary
**
** Assumptions:
**
** Comments:
**
** $Id: Stg_ComponentMeta.h 3367 2005-12-09 07:39:53Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Automation_Stg_ComponentMeta_h__
#define __StGermain_Base_Automation_Stg_ComponentMeta_h__
	
	char* Stg_Meta_GetType( Dictionary* dictionary );
	char* Stg_Meta_GetCreator( Dictionary* dictionary );
	char* Stg_Meta_GetPublisher( Dictionary* dictionary );
	char* Stg_Meta_GetRights( Dictionary* dictionary );
	char* Stg_Meta_GetSource( Dictionary* dictionary );
	char* Stg_Meta_GetSubject( Dictionary* dictionary );
	char* Stg_Meta_GetDescription( Dictionary* dictionary );

	char* Stg_Meta_GetExampleDocumentation( Dictionary* dictionary );
	char* Stg_Meta_GetExampleCode( Dictionary* dictionary );
	char* Stg_Meta_GetInherits( Dictionary* dictionary );

	char* Stg_Meta_GetReference( Dictionary* dictionary );
	char* Stg_Meta_GetEquation( Dictionary* dictionary );

	Index Stg_Meta_GetParameterCount( Dictionary* dictionary );
	char* Stg_Meta_GetParameterName( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetParameterType( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetParameterDefault( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetParameterDocumentation( Dictionary* dictionary, Index i );

	Index Stg_Meta_GetAssociationCount( Dictionary* dictionary );
	char* Stg_Meta_GetAssociationName( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetAssociationType( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetAssociationNillable( Dictionary* dictionary, Index i );
	char* Stg_Meta_GetAssociationDocumentation( Dictionary* dictionary, Index i );
        void Stg_Meta_Print( Dictionary* dictionary, Stream* stream );
#endif /* __StGermain_Base_Automation_Stg_ComponentMeta_h__ */
