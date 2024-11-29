within ExternalMedia.Carnot_Battery_Models.MediaTests;

model ReferenceAir_vs_EM
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  package EM = ExternalMedia.Media.CoolPropMedium(mediumName = "Air", substanceNames = {"Air"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
  package RefA = Modelica.Media.Air.ReferenceAir.Air_pT;
  package IdAir = Modelica.Media.Air.DryAirNasa;
  parameter SI.Pressure P_set = 559304*1;
  parameter SI.Pressure beta_set = 5.445;
  //STATE 3 discharge
  SI.Pressure p_3 = P_set "Pressure after recup";
  SI.Temperature T_3 = from_degC(556) "outlet temperature after HEX";
  EM.ThermodynamicState state_3 "thermodynamic state of HEX outlet";
  EM.SpecificEnthalpy h_3(start = 579480) "HEX outlet enthalpy";
  EM.SpecificEntropy s_3 "HEX outlet spec. entropy";
  //STATE 4 isentropic
  SI.Temperature T_4_is(start = from_degC(250)) "isentropic outlet temperature of turbine";
  EM.ThermodynamicState state_4_is "isentropic state of turbine outlet";
  EM.SpecificEnthalpy h_4_is "turbine outlet enthalpy";
  EM.SpecificEntropy s_4_is "isentropic turbine outlet spec. entropy";
  //STATE 3 discharge
  SI.Pressure p_3_RefA = P_set "Pressure after recup";
  SI.Temperature T_3_RefA = from_degC(556) "outlet temperature after HEX";
  RefA.ThermodynamicState state_3_RefA "thermodynamic state of HEX outlet";
  RefA.BaseProperties state_3_RefA_BP "base properties of HEX outlet";
  RefA.SpecificEnthalpy h_3_RefA(start = 579480) "HEX outlet enthalpy";
  RefA.SpecificEntropy s_3_RefA "HEX outlet spec. entropy";
  //STATE 4 isentropic
  SI.Temperature T_4_is_RefA(start = from_degC(250)) "isentropic outlet temperature of turbine";
  RefA.ThermodynamicState state_4_is_RefA "isentropic state of turbine outlet";
  RefA.SpecificEntropy s_4_is_RefA "isentropic turbine outlet spec. entropy";
  /*  
    RefA.BaseProperties state_4_is_RefA_BP "base properties of HEX outlet";      
    RefA.SpecificEnthalpy h_4_is_RefA "turbine outlet enthalpy";  
  
  */
  //STATE 3 discharge
  SI.Pressure p_3_id = P_set "Pressure after recup";
  SI.Temperature T_3_id = from_degC(556) "outlet temperature after HEX";
  IdAir.ThermodynamicState state_3_id "thermodynamic state of HEX outlet";
  IdAir.SpecificEntropy s_3_id "HEX outlet spec. entropy";
  //STATE 4 isentropic
  SI.Temperature T_4_is_id(start = from_degC(250)) "isentropic outlet temperature of turbine";
  IdAir.ThermodynamicState state_4_is_id "isentropic state of turbine outlet";
  IdAir.SpecificEntropy s_4_is_id "isentropic turbine outlet spec. entropy";
  //STATE 4 discharge
  SI.Pressure p_4 = P_set/beta_set " pressure at turb outlet";
  /*  
    WorkingFluid.SpecificEnthalpy h_4 "turbine outlet enthalpy";
   SI.Temperature T_4(start=from_degC(270)) "outlet temperature of turbine";
    WorkingFluid.ThermodynamicState state_4(p(start=105039), T(start=from_degC(270))) "thermodynamic state of turbine outlet";
    WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";    
       */
  //STATE 4 discharge
  SI.Pressure p_4_RefA = P_set/beta_set " pressure at turb outlet";
  /*  
    WorkingFluid.SpecificEnthalpy h_4 "turbine outlet enthalpy";
   SI.Temperature T_4(start=from_degC(270)) "outlet temperature of turbine";
    WorkingFluid.ThermodynamicState state_4(p(start=105039), T(start=from_degC(270))) "thermodynamic state of turbine outlet";
    WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";    
       */
  //STATE 4 discharge
  SI.Pressure p_4_id = P_set/beta_set " pressure at turb outlet";
  /*  
    WorkingFluid.SpecificEnthalpy h_4 "turbine outlet enthalpy";
   SI.Temperature T_4(start=from_degC(270)) "outlet temperature of turbine";
    WorkingFluid.ThermodynamicState state_4(p(start=105039), T(start=from_degC(270))) "thermodynamic state of turbine outlet";
    WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";    
       */
       EM.BaseProperties air_bp "Medium properties of port_a";
equation
//state 3 DISCHARGE
  state_3 = EM.setState_pT(p_3, T_3);
  h_3 = state_3.h;
  s_3 = EM.specificEntropy(state_3);
//state 4 isentropic
  s_4_is = s_3;
  state_4_is = EM.setState_ps(p_4, s_4_is) "isentropic state of turbine outlet";
  T_4_is = state_4_is.T "isentropic turbine outlet temperature";
  h_4_is = state_4_is.h;
/*
  //state 4
  state_4 = WorkingFluid.setState_ph(p_4, h_4);
  T_4 = state_4.T " turbine outlet spec. entropy";
  s_4 = WorkingFluid.specificEntropy(state_4);
  */
//state 3 DISCHARGE
  state_3_RefA = RefA.setState_pT(p_3_RefA, T_3_RefA);
  h_3_RefA = state_3_RefA.h;
  s_3_RefA = RefA.specificEntropy(state_3_RefA);
  state_3_RefA_BP.p = p_3_RefA;
  state_3_RefA_BP.T = T_3_RefA;
//state 4 isentropic
  s_4_is_RefA = s_3_RefA;
  state_4_is_RefA = RefA.setState_ps(p_4_RefA, s_4_is_RefA) "isentropic state of turbine outlet";
  T_4_is_RefA = state_4_is_RefA.T "isentropic turbine outlet temperature";
/*
  h_4_is_RefA=state_4_is_RefA.h;
    */
//state 3 DISCHARGE
  state_3_id = IdAir.setState_pT(p_3_id, T_3_id);
  s_3_id = IdAir.specificEntropy(state_3_id);
//state 4 isentropic
  s_4_is_id = s_3_id;
  state_4_is_id = IdAir.setState_ps(p_4_id, s_4_is_id) "isentropic state of turbine outlet";
  T_4_is_id = state_4_is_id.T "isentropic turbine outlet temperature";
/*
  h_4_is_RefA=state_4_is_RefA.h;
    */
    air_bp.p=101325;
    air_bp.T=293.15;
  annotation(
    Documentation(info = "<html><head></head><body>Fazit: Reference Air and External Media AIr have the same equations of State</body></html>"),
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
  experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 1),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,evaluateAllParameters,NLSanalyticJacobian");
end ReferenceAir_vs_EM;
