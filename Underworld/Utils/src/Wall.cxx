#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "types.h"

Name WallEnumToStr[Wall_Size] = {
       "back",
       "left",
       "bottom",
       "right",
       "top",
       "front" };
