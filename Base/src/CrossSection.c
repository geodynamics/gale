#include <StGermain/StGermain.h>

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "CrossSection.h"

/* Returns a cross section struct parsed from XML string */
lucCrossSection* lucCrossSection_Read( Stg_ComponentFactory* cf, Name component)
{
   Name crossSectionStr;
   char axisChar;
   char crossSectionVal[20];
   char modifierChar = ' ';
   lucCrossSection* self = (lucCrossSection*)Memory_Alloc_Bytes_Unnamed(sizeof(lucCrossSection), "lucCrossSection");
   self->value = 0.0;
   self->axis = 0;
   self->interpolate = False;

   /* Read the cross section string specification from xml */
   crossSectionStr = Stg_ComponentFactory_GetString( cf, component, (Dictionary_Entry_Key)"crossSection", "z=min" );

   /* axis=value    : draw at this exact value on axis
    * axis=min      : draw at minimum of range on axis
    * axis=max      : draw at maximum of range on axis
    * axis=value%   : draw at interpolated percentage value of range on axis
    * Axis is a single character, one of [xyzXYZ] */

   /* Parse the input string */
   if ( sscanf( crossSectionStr, "%c=%s", &axisChar, crossSectionVal ) == 2 ) 
   {
      /* Axis X/Y/Z */
      if ( toupper( axisChar ) >= 'X' )
		   self->axis = toupper( axisChar ) - 'X';   /* x=0 y=1 z=2 */

    	if (sscanf( crossSectionVal, "%lf%c", &self->value, &modifierChar) >= 1)
      {
         /* Found a numeric value  + optional modifier character */
         //fprintf(stderr, "CROSS SECTION VALUE %lf on Axis %c\n",self->value, axisChar);

         /* Interpolate cross section using percentage value */
         if (modifierChar == '%')
         {
            /* Interpolate between max and min value using provided value as percentage */
            self->interpolate = True;
            //fprintf(stderr, "PERCENTAGE %lf %% CROSS SECTION on Axis %c\n", self->value, axisChar);
            self->value *= 0.01;
         }
      }
      /* Max or Min specified? */
      else if (strcmp(crossSectionVal, "min") == 0) 
      {
         self->value = 0.0;
         self->interpolate = True;
         //fprintf(stderr, "MIN CROSS SECTION AT %lf on Axis %c\n", self->value, axisChar);
      }
      else if (strcmp(crossSectionVal, "max") == 0) 
      {
         self->value = 1.0;
         self->interpolate = True;
         //fprintf(stderr, "MAX CROSS SECTION AT %lf on Axis %c\n", self->value, axisChar);
      }
	}

   /* Return cross section data */
   return self;
}

/* Setup cross section values from passed parameters 
 * Returns pointer to passed in cross section so can be used in function calls 
 * If input object is NULL a new one is created */
lucCrossSection* lucCrossSection_Set(lucCrossSection* self, double value, Axis axis, Bool interpolate)
{
   /* Allocate if necessary */
   if (self == NULL)
      self = (lucCrossSection*)Memory_Alloc_Bytes_Unnamed(sizeof(lucCrossSection), "lucCrossSection");

   /* Copy values */
   self->value = value;
   self->axis = axis;
   self->interpolate = interpolate;

   /* Return pointer */
   return self; 
}

/* Returns the cross section value, interpolating where necessary */
double lucCrossSection_GetValue(lucCrossSection* self, double min, double max)
{
   if (self->interpolate)
      /* Interpolation factor 0-1 provided to determine cross-section value */
	   return min + self->value * (max - min);
   else
      /* Exact value provided */
	   return self->value;
}

/* Free cross-section memory */
void lucCrossSection_Delete(lucCrossSection* self)
{
   if (self != NULL) Memory_Free(self);
}


