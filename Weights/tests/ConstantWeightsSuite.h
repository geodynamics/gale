#ifndef StGermain_ConstantWeights_h
#define StGermain_ConstantWeights_h

void WeightsSuite_TestElementIntegral(
      DomainContext* context,
      Name           funcName,
      Index          count,
      double         meanTolerance,
      double         stdDevTolerance,
      double         expectedMean,
      double         expectedStdDev );

void ConstantWeightsSuite( pcu_suite_t* suite );

#endif
