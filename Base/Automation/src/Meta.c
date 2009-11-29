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

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"

#include "Automation.h"

/* Info parts --------------------------------------------------------------------------------------------------------------------*/

char* Stg_Meta_GetType( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "title" ) );
}

char* Stg_Meta_GetCreator( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "creator" ) );
}

char* Stg_Meta_GetPublisher( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "publisher" ) );
}

char* Stg_Meta_GetRights( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "rights" ) );
}

char* Stg_Meta_GetSource( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "source" ) );
}

char* Stg_Meta_GetSubject( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "subject" ) );
}

char* Stg_Meta_GetDescription( Dictionary* dictionary ) {
	Dictionary* info = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "info" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( info, "description" ) );
}


/* Code parts --------------------------------------------------------------------------------------------------------------------*/

char* Stg_Meta_GetExampleDocumentation( Dictionary* dictionary ) {
	Dictionary* code = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "code" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( code, "example-documentation" ) );
}

char* Stg_Meta_GetExampleCode( Dictionary* dictionary ) {
	Dictionary* code = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "code" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( code, "example-code" ) );
}

char* Stg_Meta_GetInherits( Dictionary* dictionary ) {
	Dictionary* code = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "code" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( code, "inherits" ) );
}


/* Implements parts --------------------------------------------------------------------------------------------------------------*/

char* Stg_Meta_GetReference( Dictionary* dictionary ) {
	Dictionary* implements = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "implements" ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( implements, "reference" ) );
}

char* Stg_Meta_GetEquation( Dictionary* dictionary ) {
	Dictionary* implements = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "implements"  ) );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( implements, "equation" ) );
}


/* Parameter parts ---------------------------------------------------------------------------------------------------------------*/

Index Stg_Meta_GetParameterCount( Dictionary* dictionary ) {
	Dictionary* parameters = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "parameters" ) );
	return Dictionary_GetCount( parameters );
}

char* Stg_Meta_GetParameterName( Dictionary* dictionary, Index i ) {
	Dictionary* parameters = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "parameters" ) );
	Dictionary* parameter = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( parameters, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( parameter, "name" ) );
}

char* Stg_Meta_GetParameterType( Dictionary* dictionary, Index i ) {
	Dictionary* parameters = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "parameters" ) );
	Dictionary* parameter = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( parameters, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( parameter, "type" ) );
}

char* Stg_Meta_GetParameterDefault( Dictionary* dictionary, Index i ) {
	Dictionary* parameters = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "parameters" ) );
	Dictionary* parameter = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( parameters, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( parameter, "default" ) );
}

char* Stg_Meta_GetParameterDocumentation( Dictionary* dictionary, Index i ) {
	Dictionary* parameters = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "parameters" ) );
	Dictionary* parameter = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( parameters, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( parameter, "documentation" ) );
}


/* Association parts -------------------------------------------------------------------------------------------------------------*/

Index Stg_Meta_GetAssociationCount( Dictionary* dictionary ) {
	Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "associations" ) );
	return Dictionary_GetCount( associations );
}

char* Stg_Meta_GetAssociationName( Dictionary* dictionary, Index i ) {
	Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "associations" ) );
	Dictionary* association = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( associations, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( association, "name" ) );
}

char* Stg_Meta_GetAssociationType( Dictionary* dictionary, Index i ) {
	Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "associations" ) );
	Dictionary* association = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( associations, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( association, "type" ) );
}

char* Stg_Meta_GetAssociationNillable( Dictionary* dictionary, Index i ) {
	Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "associations" ) );
	Dictionary* association = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( associations, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( association, "nillable" ) );
}

char* Stg_Meta_GetAssociationDocumentation( Dictionary* dictionary, Index i ) {
	Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "associations" ) );
	Dictionary* association = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( associations, i )  );
	return Dictionary_Entry_Value_AsString( Dictionary_Get( association, "documentation" ) );
}

/* Print function ----------------------------------------------------------------------------------------------------------------*/
void Stg_Meta_Print( Dictionary* dictionary, Stream* stream ) {
	Index i;
	char* str;

	/* Info parts */
	Journal_Printf( stream, "Type: %s\n", (str = Stg_Meta_GetType( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Creator: %s\n", (str = Stg_Meta_GetCreator( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Publisher: %s\n", (str = Stg_Meta_GetPublisher( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Rights: %s\n", (str = Stg_Meta_GetRights( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Source: %s\n", (str = Stg_Meta_GetSource( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Subject: %s\n", (str = Stg_Meta_GetSubject( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Description: %s\n", (str = Stg_Meta_GetDescription( dictionary )) ? str : "(not provided)" );

	/* Code parts */
	Journal_Printf( stream, "Example documentation: %s\n", (str = Stg_Meta_GetExampleDocumentation( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Example code: %s\n", (str = Stg_Meta_GetExampleCode( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Inherits: %s\n", (str = Stg_Meta_GetInherits( dictionary )) ? str : "Stg_Component (assumed)" );

	/* Implements parts */
	Journal_Printf( stream, "Reference: %s\n", (str = Stg_Meta_GetReference( dictionary )) ? str : "(not provided)" );
	Journal_Printf( stream, "Equation: %s\n", (str = Stg_Meta_GetEquation( dictionary )) ? str : "(not provided)" );

	/* Parameter parts */
	Journal_Printf( stream, "Parameters:\n" );
	Stream_Indent( stream );
	for ( i = 0; i < Stg_Meta_GetParameterCount( dictionary ); i++ ) {
		Journal_Printf( stream, "{\n" );
		Stream_Indent( stream );
		Journal_Printf( stream, "Type: %s\n", (str = Stg_Meta_GetParameterType( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Name: %s\n", (str = Stg_Meta_GetParameterName( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Default: %s\n", (str = Stg_Meta_GetParameterDefault( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Documentation: %s\n", (str = Stg_Meta_GetParameterDocumentation( dictionary, i )) ? str : "(not provided)" );
		Stream_UnIndent( stream );
		Journal_Printf( stream, "}\n" );
	}
	Stream_UnIndent( stream );
	
	/* Association parts */
	Journal_Printf( stream, "Associations:\n" );
	Stream_Indent( stream );
	for ( i = 0; i < Stg_Meta_GetAssociationCount( dictionary ); i++ ) {
		Journal_Printf( stream, "{\n" );
		Stream_Indent( stream );
		Journal_Printf( stream, "Type: %s\n", (str = Stg_Meta_GetAssociationType( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Name: %s\n", (str = Stg_Meta_GetAssociationName( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Nillable: %s\n", (str = Stg_Meta_GetAssociationNillable( dictionary, i )) ? str : "(not provided)" );
		Journal_Printf( stream, "Documentation: %s\n", (str = Stg_Meta_GetAssociationDocumentation( dictionary, i )) ? str : "(not provided)" );
		Stream_UnIndent( stream );
		Journal_Printf( stream, "}\n" );
	}
	Stream_UnIndent( stream );
}



