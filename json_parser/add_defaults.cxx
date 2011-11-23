#include <json_spirit_reader.h>
#include <stdexcept>

void add_components(json_spirit::Object &root,
                    const json_spirit::Object &components);

void add_defaults(json_spirit::Value &root)
{
  std::string toolbox_string
    ("{"
    "\"import\":"
    "{"
    "    \"toolbox\": \"Underworld\""
    "}}");

  std::string plugins_string
    ("{"
    "\"plugins\":"
    "["
    "    {"
    "        \"Type\": \"StgFEM_StandardConditionFunctions\","
    "        \"Context\": \"context\""
    "    }"
    "]}");


  std::string basic_components
    ("{"
    "\"context\":"
    "{"
    "    \"Type\": \"UnderworldContext\""
    "},"
    "\"v-mesh\":"
    "{"
    "    \"Type\": \"FeMesh\","
    "    \"elementType\": \"quadratic\""
    "},"
    "\"v-mesh-generator\":"
    "{"
    "    \"Type\": \"C2Generator\","
    "    \"mesh\": \"v-mesh\","
    "    \"dim\": \"dim\","
    "    \"shadowDepth\": \"shadowDepth\","
    "    \"regular\": \"False\","
    "    \"size\": ["
    "        \"nx\","
    "        \"ny\","
    "        \"nz\""
    "    ],"
    "    \"minCoord\": ["
    "        \"minX\","
    "        \"minY\","
    "        \"minZ\""
    "    ],"
    "    \"maxCoord\": ["
    "        \"maxX\","
    "        \"maxY\","
    "        \"maxZ\""
    "    ]"
    "},"
    "\"velocity\":"
    "{"
    "    \"Type\": \"MeshVariable\","
    "    \"mesh\": \"v-mesh\","
    "    \"Rank\": \"Vector\","
    "    \"DataType\": \"Double\","
    "    \"VectorComponentCount\": \"dim\","
    "    \"names\": ["
    "        \"vx\","
    "        \"vy\","
    "        \"vz\""
    "    ]"
    "},"
    "\"velocityBCs\":"
    "{"
    "    \"Type\": \"CompositeVC\","
    "    \"Data\": \"v-mesh\""
    "},"
    "\"velocityICs\":"
    "{"
    "    \"Type\": \"CompositeVC\","
    "    \"Data\": \"v-mesh\""
    "},"
    "\"velocityDofLayout\":"
    "{"
    "    \"Type\": \"DofLayout\","
    "    \"mesh\": \"v-mesh\","
    "    \"BaseVariableCount\": \"dim\","
    "    \"BaseVariables\": ["
    "        \"vx\","
    "        \"vy\","
    "        \"vz\""
    "    ]"
    "},"
    "\"VelocityField\":"
    "{"
    "    \"Type\": \"FeVariable\","
    "    \"FEMesh\": \"v-mesh\","
    "    \"DofLayout\": \"velocityDofLayout\","
    "    \"BC\": \"velocityBCs\","
    "    \"IC\": \"velocityICs\","
    "    \"LinkedDofInfo\": \"velocityLinkedDofs\""
    "},"
    "\"cellLayout\":"
    "{"
    "    \"Type\": \"ElementCellLayout\","
    "    \"Mesh\": \"v-mesh\""
    "},"
    "\"particleLayout\":"
    "{"
    "    \"Type\": \"GaussParticleLayout\","
    "    \"gaussParticles\": \"3\""
    "},"
    "\"gaussSwarm\":"
    "{"
    "    \"Type\": \"IntegrationPointsSwarm\","
    "    \"CellLayout\": \"cellLayout\","
    "    \"ParticleLayout\": \"particleLayout\","
    "    \"FeMesh\": \"v-mesh\","
    "    \"TimeIntegrator\": \"timeIntegrator\","
    "    \"IntegrationPointMapper\": \"oneToManyMapper\""
    "},"
    "\"oneToManyMapper\":"
    "{"
    "    \"Type\": \"OneToManyMapper\","
    "    \"IntegrationPointsSwarm\": \"gaussSwarm\","
    "    \"MappedSwarm\": \"picIntegrationPoints\""
    "},"
    "\"timeIntegrator\":"
    "{"
    "    \"Type\": \"TimeIntegrator\","
    "    \"order\": \"timeIntegratorOrder\","
    "    \"simultaneous\": \"f\","
    "    \"Context\": \"context\""
    "},"
    "\"weights\":"
    "{"
    "    \"Type\": \"PCDVC\","
    "    \"resolutionX\": \"20\","
    "    \"resolutionY\": \"20\","
    "    \"resolutionZ\": \"20\","
    "    \"lowerT\": \"0.6\","
    "    \"upperT\": \"25\","
    "    \"maxDeletions\": \"3\","
    "    \"maxSplits\": \"3\","
    "    \"MaterialPointsSwarm\": \"materialSwarm\","
    "    \"Inflow\": \"True\""
    "},"
    "\"localLayout\":"
    "{"
    "    \"Type\": \"MappedParticleLayout\""
    "},"
    "\"picIntegrationPoints\":"
    "{"
    "    \"Type\": \"IntegrationPointsSwarm\","
    "    \"CellLayout\": \"cellLayout\","
    "    \"ParticleLayout\": \"localLayout\","
    "    \"FeMesh\": \"v-mesh\","
    "    \"WeightsCalculator\": \"weights\","
    "    \"TimeIntegrator\": \"timeIntegrator\","
    "    \"IntegrationPointMapper\": \"mapper\""
    "},"
    "\"mapper\":"
    "{"
    "    \"Type\": \"CoincidentMapper\","
    "    \"IntegrationPointsSwarm\": \"picIntegrationPoints\","
    "    \"MaterialPointsSwarm\": \"materialSwarm\""
    "},"
    "\"materialSwarmParticleLayout\":"
    "{"
    "    \"Type\": \"MeshParticleLayout\","
    "    \"cellParticleCount\": \"particlesPerCell\","
    "    \"mesh\": \"v-mesh\""
    "},"
    "\"pMovementHandler\":"
    "{"
    "    \"Type\": \"ParticleMovementHandler\""
    "},"
    "\"pShadowSync\":"
    "{"
    "    \"Type\": \"ParticleShadowSync\""
    "},"
    "\"materialSwarm\":"
    "{"
    "    \"Type\": \"MaterialPointsSwarm\","
    "    \"CellLayout\": \"cellLayout\","
    "    \"ParticleLayout\": \"materialSwarmParticleLayout\","
    "    \"FeMesh\": \"v-mesh\","
    "    \"ParticleCommHandlers\": ["
    "        \"pMovementHandler\","
    "        \"pShadowSync\""
    "    ],"
    "    \"SplittingRoutine\": \"splittingRoutine\","
    "    \"RemovalRoutine\": \"removalRoutine\","
    "    \"EscapedRoutine\": \"escapedRoutine\""
    "},"
    "\"materialSwarmAdvector\":"
    "{"
    "    \"Type\": \"SwarmAdvector\","
    "    \"Swarm\": \"materialSwarm\","
    "    \"TimeIntegrator\": \"timeIntegrator\","
    "    \"VelocityField\": \"VelocityField\","
    "    \"PeriodicBCsManager\": \"periodicBCsManager\","
    "    \"allowFallbackToFirstOrder\": \"True\""
    "},"
    "\"escapedRoutine\":"
    "{"
    "    \"Type\": \"EscapedRoutine\","
    "    \"idealParticleCount\": \"0\""
    "},"
    "\"velocityRemesher\":"
    "{"
    "    \"Type\": \"RegularRemesherCmpt\","
    "    \"mesh\": \"v-mesh\","
    "    \"remeshDims\": ["
    "        \"0\","
    "        \"1\","
    "        \"2\""
    "    ]"
    "}"
    "}");

  std::string stokes_components
    ("{"
    "\"pressure-mesh\":"
    "{"
    "    \"Type\": \"FeMesh\","
    "    \"elementType\": \"linear-inner\""
    "},"
    "\"pressure-mesh-generator\":"
    "{"
    "    \"Type\": \"InnerGenerator\","
    "    \"mesh\": \"pressure-mesh\","
    "    \"elementMesh\": \"v-mesh\","
    "    \"dim\": \"dim\","
    "    \"shadowDepth\": \"shadowDepth\","
    "    \"regular\": \"False\","
    "    \"size\": ["
    "        \"nx\","
    "        \"ny\","
    "        \"nz\""
    "    ],"
    "    \"minCoord\": ["
    "        \"minX\","
    "        \"minY\","
    "        \"minZ\""
    "    ],"
    "    \"maxCoord\": ["
    "        \"maxX\","
    "        \"maxY\","
    "        \"maxZ\""
    "    ]"
    "},"
    "\"VelocityGradientsField\":"
    "{"
    "    \"Type\": \"OperatorFeVariable\","
    "    \"Operator\": \"Gradient\","
    "    \"FeVariable\": \"VelocityField\""
    "},"
    "\"StrainRateField\":"
    "{"
    "    \"Type\": \"OperatorFeVariable\","
    "    \"Operator\": \"TensorSymmetricPart\","
    "    \"FeVariable\": \"VelocityGradientsField\""
    "},"
    "\"StrainRateInvariantField\":"
    "{"
    "    \"Type\": \"OperatorFeVariable\","
    "    \"Operator\": \"SymmetricTensor_Invariant\","
    "    \"FeVariable\": \"StrainRateField\""
    "},"
    "\"pressure\":"
    "{"
    "    \"Type\": \"MeshVariable\","
    "    \"mesh\": \"pressure-mesh\","
    "    \"Rank\": \"Scalar\","
    "    \"DataType\": \"Double\""
    "},"
    "\"pressureDofLayout\":"
    "{"
    "    \"Type\": \"DofLayout\","
    "    \"mesh\": \"pressure-mesh\","
    "    \"BaseVariables\": ["
    "        \"pressure\""
    "    ]"
    "},"
    "\"PressureField\":"
    "{"
    "    \"Type\": \"FeVariable\","
    "    \"FEMesh\": \"pressure-mesh\","
    "    \"DofLayout\": \"pressureDofLayout\","
    "    \"LinkedDofInfo\": \"pressureLinkedDofs\""
    "},"
    "\"solutionVelocity\":"
    "{"
    "    \"Type\": \"SolutionVector\","
    "    \"FeVariable\": \"VelocityField\""
    "},"
    "\"solutionPressure\":"
    "{"
    "    \"Type\": \"SolutionVector\","
    "    \"FeVariable\": \"PressureField\""
    "},"
    "\"mom_force\":"
    "{"
    "    \"Type\": \"ForceVector\","
    "    \"FeVariable\": \"VelocityField\","
    "    \"ExtraInfo\": \"context\""
    "},"
    "\"cont_force\":"
    "{"
    "    \"Type\": \"ForceVector\","
    "    \"FeVariable\": \"PressureField\","
    "    \"ExtraInfo\": \"context\""
    "},"
    "\"k_matrix\":"
    "{"
    "    \"Type\": \"StiffnessMatrix\","
    "    \"RowVariable\": \"VelocityField\","
    "    \"ColumnVariable\": \"VelocityField\","
    "    \"RHS\": \"mom_force\","
    "    \"allowZeroElementContributions\": \"False\""
    "},"
    "\"constitutiveMatrix\":"
    "{"
    "    \"Type\": \"ConstitutiveMatrixCartesian\","
    "    \"Swarm\": \"gaussSwarm\","
    "    \"StiffnessMatrix\": \"k_matrix\""
    "},"
    "\"g_matrix\":"
    "{"
    "    \"Type\": \"StiffnessMatrix\","
    "    \"RowVariable\": \"VelocityField\","
    "    \"ColumnVariable\": \"PressureField\","
    "    \"RHS\": \"mom_force\","
    "    \"transposeRHS\": \"cont_force\","
    "    \"allowZeroElementContributions\": \"False\""
    "},"
    "\"gradientStiffnessMatrixTerm\":"
    "{"
    "    \"Type\": \"GradientStiffnessMatrixTerm\","
    "    \"Swarm\": \"gaussSwarm\","
    "    \"StiffnessMatrix\": \"g_matrix\""
    "},"
    "\"preconditioner\":"
    "{"
    "    \"Type\": \"StiffnessMatrix\","
    "    \"RowVariable\": \"PressureField\","
    "    \"ColumnVariable\": \"PressureField\","
    "    \"RHS\": \"cont_force\","
    "    \"allowZeroElementContributions\": \"True\""
    "},"
    "\"preconditionerTerm\":"
    "{"
    "    \"Type\": \"UzawaPreconditionerTerm\","
    "    \"Swarm\": \"gaussSwarm\","
    "    \"StiffnessMatrix\": \"preconditioner\""
    "},"
    "\"uzawa\":"
    "{"
    "    \"Type\": \"Stokes_SLE_UzawaSolver\","
    "    \"velocitySolver\": \"matrixSolver\","
    "    \"Preconditioner\": \"preconditioner\","
    "    \"tolerance\": \"linearTolerance\","
    "    \"monitor\": \"false\","
    "    \"maxIterations\": \"5000\","
    "    \"minIterations\": \"1\""
    "},"
    "\"stokesEqn\":"
    "{"
    "    \"Type\": \"Stokes_SLE\","
    "    \"SLE_Solver\": \"uzawa\","
    "    \"Context\": \"context\","
    "    \"StressTensorMatrix\": \"k_matrix\","
    "    \"GradientMatrix\": \"g_matrix\","
    "    \"DivergenceMatrix\": \"\","
    "    \"VelocityVector\": \"solutionVelocity\","
    "    \"PressureVector\": \"solutionPressure\","
    "    \"ForceVector\": \"mom_force\","
    "    \"ContinuityForceVector\": \"cont_force\","
    "    \"nonLinearMinIterations\": \"nonLinearMinIterations\","
    "    \"nonLinearMaxIterations\": \"nonLinearMaxIterations\","
    "    \"nonLinearTolerance\": \"nonLinearTolerance\","
    "    \"makeConvergenceFile\": \"false\""
    "},"
    "\"storeViscosity\":"
    "{"
    "    \"Type\": \"StoreVisc\","
    "    \"MaterialPointsSwarm\": \"materialSwarm\""
    "},"
    "\"storeStress\":"
    "{"
    "    \"Type\": \"StoreStress\","
    "    \"MaterialPointsSwarm\": \"materialSwarm\""
    "}"
    "}");

  std::string thermal_components
    ("{"
     "\"T-mesh\":"
     "{"
     "    \"Type\": \"FeMesh\","
     "    \"elementType\": \"linear\""
     "},"
     "\"T-mesh-generator\":"
     "{"
     "    \"Type\": \"CartesianGenerator\","
     "    \"mesh\": \"T-mesh\","
     "    \"dim\": \"dim\","
     "    \"shadowDepth\": \"shadowDepth\","
     "    \"regular\": \"False\","
     "    \"size\": ["
     "        \"nx\","
     "        \"ny\","
     "        \"nz\""
     "    ],"
     "    \"minCoord\": ["
     "        \"minX\","
     "        \"minY\","
     "        \"minZ\""
     "    ],"
     "    \"maxCoord\": ["
     "        \"maxX\","
     "        \"maxY\","
     "        \"maxZ\""
     "    ]"
     "},"
     "\"temperature\":"
     "{"
     "\"Type\": \"MeshVariable\","
     "\"Rank\": \"Scalar\","
     "\"DataType\": \"Double\","
     "\"mesh\": \"T-mesh\""
     "},"
     "\"temperatureBCs\":"
     "{"
     "\"Type\": \"CompositeVC\","
     "\"Data\": \"T-mesh\""
     "},"
     "\"temperatureICs\":"
     "{"
     "\"Type\": \"CompositeVC\","
     "\"Data\": \"T-mesh\""
     "},"
     "\"temperatureDofLayout\":"
     "{"
     "\"Type\": \"DofLayout\","
     "\"mesh\": \"T-mesh\","
     "\"BaseVariables\":"
     "["
     "\"temperature\""
     "]"
     "},"
     "\"TemperatureField\":"
     "{"
     "\"Type\": \"FeVariable\","
     "\"FEMesh\": \"T-mesh\","
     "\"DofLayout\": \"temperatureDofLayout\","
     "\"BC\": \"temperatureBCs\","
     "\"IC\": \"temperatureICs\","
     "\"LinkedDofInfo\": \"temperatureLinkedDofs\""
     "},"
     "\"TemperatureGradientsField\":"
     "{"
     "\"Type\": \"OperatorFeVariable\","
     "\"Operator\": \"Gradient\","
     "\"FeVariable\": \"TemperatureField\""
     "},"
     "\"displacement\":"
     "{"
     "\"Type\": \"MeshVariable\","
     "\"Rank\": \"Vector\","
     "\"DataType\": \"Double\","
     "\"mesh\": \"v-mesh\","
     "\"VectorComponentCount\": \"dim\","
     "\"names\" : ["
     "\"dx\","
     "\"dy\","
     "\"dz\""
     "]"
     "},"
     "\"displacementBCs\":"
     "{"
     "\"Type\": \"CompositeVC\","
     "\"Data\": \"v-mesh\""
     "},"
     "\"displacementICs\":"
     "{"
     "\"Type\": \"CompositeVC\","
     "\"Data\": \"v-mesh\""
     "},"
     "\"displacementDofLayout\":"
     "{"
     "\"Type\": \"DofLayout\","
     "\"mesh\": \"v-mesh\","
     "\"BaseVariableCount\": \"dim\","
     "\"BaseVariables\":"
     "["
     "\"dx\","
     "\"dy\","
     "\"dz\""
     "]"
     "},"
     "\"DisplacementField\":"
     "{"
     "\"Type\": \"FeVariable\","
     "\"FEMesh\": \"v-mesh\","
     "\"DofLayout\": \"displacementDofLayout\","
     "\"BC\": \"displacementBCs\","
     "\"IC\": \"displacementICs\""
     "},"

     "\"residual\":"
     "{"
     "\"Type\": \"ForceVector\","
     "\"FeVariable\": \"TemperatureField\""
     "},"
     "\"massMatrix\":"
     "{"
     "\"Type\": \"ForceVector\","
     "\"FeVariable\": \"TemperatureField\""
     "},"
     "\"predictorMulticorrector\":"
     "{"
     "\"Type\": \"AdvDiffMulticorrector\""
     "},"
     "\"EnergyEqn\":"
     "{"
     "\"Type\": \"AdvectionDiffusionSLE\","
     "\"SLE_Solver\": \"predictorMulticorrector\","
     "\"Context\": \"context\","
     "\"PhiField\": \"TemperatureField\","
     "\"Residual\": \"residual\","
     "\"MassMatrix\": \"massMatrix\","
     "\"courantFactor\": \"0.25\""
     "},"
     "\"lumpedMassMatrixForceTerm\":"
     "{"
     "\"Type\": \"LumpedMassMatrixForceTerm\","
     "\"Swarm\": \"gaussSwarm\","
     "\"ForceVector\": \"massMatrix\""
     "},"
     "\"defaultResidualForceTerm\":"
     "{"
     "\"Type\": \"AdvDiffResidualForceTerm\","
     "\"Swarm\": \"gaussSwarm\","
     "\"ForceVector\": \"residual\","
     "\"ExtraInfo\": \"EnergyEqn\","
     "\"VelocityField\": \"VelocityField\","
     "\"defaultDiffusivity\": \"defaultDiffusivity\","
     "\"UpwindXiFunction\": \"Exact\""
     "},"
     "\"internalHeatingTerm\":"
     "{"
     "\"Type\": \"RadiogenicHeatingTerm\","
     "\"ForceVector\": \"residual\","
     "\"Swarm\": \"gaussSwarm\""
     "}"
     "}");

  std::string tracer_components
    ("{"
     "\"passiveSwarmMovementHandler\":"
     "{"
     "\"Type\": \"ParticleMovementHandler\""
     "},"
     "\"passiveTracerSwarm\":"
     "{"
     "\"Type\": \"MaterialPointsSwarm\","
     "\"CellLayout\": \"cellLayout\","
     "\"ParticleLayout\": \"pLayout\","
     "\"FeMesh\": \"v-mesh\","
     "\"ParticleCommHandlers\":"
     "["
     "\"passiveSwarmMovementHandler\""
     "],"
     "\"EscapedRoutine\": \"escapedRoutine\""
     "},"
     ""
     "\"passiveTracerAdvect\":"
     "{"
     "\"Type\": \"SwarmAdvector\","
     "\"Swarm\": \"passiveTracerSwarm\","
     "\"TimeIntegrator\": \"timeIntegrator\","
     "\"VelocityField\": \"VelocityField\""
     "}"
     "}");

  json_spirit::Value toolbox, plugins, basic, stokes, thermal, tracers;
  json_spirit::read_or_throw(toolbox_string,toolbox);
  json_spirit::read_or_throw(plugins_string,plugins);
  json_spirit::read_or_throw(basic_components,basic);
  json_spirit::read_or_throw(stokes_components,stokes);
  json_spirit::read_or_throw(thermal_components,thermal);
  json_spirit::read_or_throw(tracer_components,tracers);

  json_spirit::Object &o(root.get_obj());

  /* Add the toolbox */
  bool found=false;
  for(json_spirit::Object::const_iterator i=o.begin(); i!=o.end(); ++i)
    {
      if(i->name_=="import")
        {
          found=true;
          break;
        }
    }
  if(!found)
    o.insert(o.end(),toolbox.get_obj().begin(),toolbox.get_obj().end());

  /* Add the plugins */
  found=false;
  for(json_spirit::Object::iterator i=o.begin(); i!=o.end(); ++i)
    {
      if(i->name_=="plugins")
        {
          if(i->value_.type()!=json_spirit::array_type)
            throw std::runtime_error("'plugins' must be an array'");
          json_spirit::Array
            &default_plugins=(*(plugins.get_obj().begin())).value_.get_array();
          json_spirit::Array &a(i->value_.get_array());

          a.insert(a.end(),default_plugins.begin(),default_plugins.end());
          found=true;
          break;
        }
    }
  if(!found)
    o.insert(o.end(),plugins.get_obj().begin(),plugins.get_obj().end());

  /* Add the standard components */
  /* First find the components */
  json_spirit::Object::iterator components(o.end());
  bool enable_thermal(false), enable_stokes(true), enable_tracers(false);
  for(json_spirit::Object::iterator i=o.begin(); i!=o.end(); ++i)
    {
      if(i->name_=="components")
        components=i;
      else if(i->name_=="enable-stokes")
        {
          if(i->value_.type()!=json_spirit::bool_type)
            throw std::runtime_error("The field enable-stokes must be either"
                                     "true or false (no quotes)");
          enable_stokes=i->value_.get_bool();
        }
      else if(i->name_=="enable-thermal")
        {
          if(i->value_.type()!=json_spirit::bool_type)
            throw std::runtime_error("The field enable-thermal must be either"
                                     "true or false (no quotes)");
          enable_thermal=i->value_.get_bool();
        }
      else if(i->name_=="enable-tracers")
        {
          if(i->value_.type()!=json_spirit::bool_type)
            throw std::runtime_error("The field enable-tracers must be either"
                                     "true or false (no quotes)");
          enable_tracers=i->value_.get_bool();
        }
    }

  /* New components get prepended, so add in reverse order */
  if(enable_tracers)
    add_components(components->value_.get_obj(),tracers.get_obj());
  if(enable_thermal)
    add_components(components->value_.get_obj(),thermal.get_obj());
  if(enable_stokes)
    add_components(components->value_.get_obj(),stokes.get_obj());
  add_components(components->value_.get_obj(),basic.get_obj());
}
