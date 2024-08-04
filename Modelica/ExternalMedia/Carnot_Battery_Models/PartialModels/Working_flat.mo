within ExternalMedia.Carnot_Battery_Models.PartialModels;

model Working_flat
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  //import Constants
  import Modelica.Constants.pi;
  import Modelica.Constants.g_n;

  package NaK "NaK properties from CoolProp"
    extends ExternalMedia.Media.IncompressibleCoolPropMedium(mediumName = "NaK", substanceNames = {"NaK"});
  end NaK;

  replaceable package HotTESLiquid = NaK constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
  /*
    package WorkingFluid = ExternalMedia.Media.CoolPropMedium(mediumName = "Nitrogen", substanceNames = {"N2"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
   */
  package WorkingFluid = Modelica.Media.Air.ReferenceAir.Air_pT;

  package Methanol "NaK properties from CoolProp"
    extends ExternalMedia.Media.IncompressibleCoolPropMedium(mediumName = "MMA", substanceNames = {"MMA[0.6]"});
  end Methanol;

  replaceable package ColdTESLiquid = Methanol constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
  //--------------------------PARAMETERS & VARIABLES-----------------------------//
  //ambient
  parameter SI.Pressure p_amb = 101325;
  //geometry for both solar salt tanks
  parameter SI.Diameter D_tank_solsalt = 37.34;
  parameter SI.Height h_tank_solsalt = 12.44;
  parameter SI.Area A_cross_tank_solsalt = pi*(D_tank_solsalt/2)^2 "Cross-section Tank";
  parameter SI.Volume V_tank_solsalt = A_cross_tank_solsalt*h_tank_solsalt "Volume Tank";
  parameter SI.Thickness d_insulation_tank_solsalt = 0.4 "thickness of insulation wall layer";
  parameter SI.Radius r_tank_solsalt = D_tank_solsalt/2 "radius at start of insulation";
  parameter SI.Radius r_outer = D_tank_solsalt/2 + d_insulation_tank_solsalt "outer radius";
  parameter SI.Area A_W_tank_solsalt = (pi*D_tank_solsalt*h_tank_solsalt) + A_cross_tank_solsalt*2 "Total Surface Wall";
  //geometry for both methanol tanks
  parameter SI.Diameter D_tank_coldliq = 37.50;
  parameter SI.Height h_tank_coldliq = 12.50;
  parameter SI.Area A_cross_tank_coldliq = pi*(D_tank_coldliq/2)^2 "Cross-section Tank";
  parameter SI.Volume V_tank_coldliq = A_cross_tank_coldliq*h_tank_coldliq "Volume Tank";
  parameter SI.Thickness d_insulation_tank_coldliq = 0.4 "thickness of insulation wall layer";
  parameter SI.Radius r_tank_coldliq = D_tank_coldliq/2 "radius at start of insulation";
  parameter SI.Radius r_outer_coldliq = D_tank_coldliq/2 + d_insulation_tank_coldliq "outer radius";
  parameter SI.Area A_W_tank_coldliq = (pi*D_tank_coldliq*h_tank_coldliq) + A_cross_tank_coldliq*2 "Total Surface Wall";
  //--------------------------TANK 1
  SI.HeatFlowRate Q_dot_to_amb_tank1(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank1;
  //nominal
  parameter SI.Temperature T_nom_tank1 = from_degC(560);
  parameter SI.Pressure p_nom_tank1 = p_amb "unpressurized tank";
  HotTESLiquid.ThermodynamicState solsalt_tank1_nom(p(start = p_nom_tank1), T(start = T_nom_tank1)) "Nominal thermodynamic state";
  HotTESLiquid.ThermodynamicState solsalt_tank1_state(p(start = p_nom_tank1), T(start = T_nom_tank1)) "Nominal thermodynamic state";
  SI.Energy Energy_nom_tank1(displayUnit = "MWh");
  parameter SI.Mass m_nom_tank1 = V_tank_solsalt*rho_start_tank1;
  //start values
  parameter SI.Temperature T_start_tank1 = from_degC(560);
  parameter SI.Mass m_start_tank1 = m_nom_tank1;
  parameter SI.Density rho_start_tank1 = 1734;
  //state variables
  SI.Mass m_tank1(start = m_nom_tank1);
  Real SOC_tank1;
  Real SOE_tank1;
  HotTESLiquid.BaseProperties solsalt_tank1(u(start = 298675)) "Medium properties of port_a";
  SI.Temperature T_tank1(start = T_start_tank1);
  //other
  SI.Energy Energy_content_tank1(displayUnit = "MWh");
  SI.PathLength x_tank1 "salt-level";
  //outflow
  SI.MassFlowRate m_out_tank1;
  parameter SI.MassFlowRate m_in_tank1 = 0;
  parameter SI.SpecificEnthalpy h_in_tank1 = 1;
  //dummy value TO_DO
  SI.HeatFlowRate Q_port_tank1(displayUnit = "MW") "Heat Flow from tank 2";
  parameter SI.Height x_min_tank1 = 0.4;
  //--------------------------TANK 2
  SI.HeatFlowRate Q_dot_to_amb_tank2(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank2;
  //nominal
  parameter SI.Temperature T_nom_tank2 = from_degC(271.5);
  parameter SI.Pressure p_nom_tank2 = p_amb "unpressurized tank";
  HotTESLiquid.ThermodynamicState solsalt_tank2_nom(p(start = p_nom_tank2), T(start = T_nom_tank2)) "Nominal thermodynamic state";
  HotTESLiquid.ThermodynamicState solsalt_tank2_state(p(start = p_nom_tank2), T(start = T_nom_tank2)) "Nominal thermodynamic state";
  SI.Energy Energy_nom_tank2(displayUnit = "MWh");
  parameter SI.Mass m_nom_tank2 = V_tank_solsalt*rho_start_tank2;
  parameter SI.Mass m_empty_tank2 = x_min_tank2*A_cross_tank_solsalt*rho_start_tank2;
  //start values
  parameter SI.Temperature T_start_tank2 = from_degC(271.5);
  parameter SI.Mass m_start_tank2 = m_nom_tank2;
  parameter SI.Density rho_start_tank2 = 1896.02;
  //state variables
  SI.Mass m_tank2(start = m_nom_tank2);
  Real SOC_tank2;
  Real SOE_tank2;
  HotTESLiquid.BaseProperties solsalt_tank2(u(start = 298675)) "Medium properties of port_a";
  SI.Temperature T_tank2(start = T_start_tank2);
  SI.Temperature T_in_tank2(start = T_start_tank2);
  //other
  SI.Energy Energy_content_tank2(displayUnit = "MWh");
  SI.PathLength x_tank2 "salt-level";
  //inflow
  SI.SpecificEnthalpy h_in_tank2;
  //outflow
  parameter SI.MassFlowRate m_out_tank2 = 0;
  SI.MassFlowRate m_in_tank2;
  parameter SI.SpecificEnthalpy h_out_tank2 = 1;
  //dummy value TO_DO
  SI.HeatFlowRate Q_port_tank2(displayUnit = "MW") "Heat Flow into tank 2";
  parameter SI.Height x_min_tank2 = 0.4;
  //--------------------------TANK 3
  //nominal
  parameter SI.Temperature T_nom_tank3 = from_degC(25);
  parameter SI.Pressure p_nom_tank3 = p_amb "unpressurized tank";
  parameter SI.Mass m_nom_tank3 = 13116801.563929563;
  ColdTESLiquid.ThermodynamicState coldliq_tank3_nom(p(start = p_nom_tank3), T(start = T_nom_tank3)) "Nominal thermodynamic state";
  //start values
  parameter SI.Temperature T_start_tank3 = from_degC(23);
  parameter SI.Mass m_start_tank3 = x_min_tank3*rho_start_tank3*A_cross_tank_coldliq;
  parameter SI.Density rho_start_tank3 = 890.69;
  //actual
  SI.Temperature T_tank3(start = T_start_tank3);
  SI.Temperature T_in_tank3(start = T_start_tank3);
  ColdTESLiquid.ThermodynamicState coldliq_tank3_state(p(start = p_nom_tank3), T(start = T_nom_tank3)) "Nominal thermodynamic state";
  SI.HeatFlowRate Q_dot_to_amb_tank3(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank3;
  //state variables
  SI.Mass m_tank3(start = m_start_tank3);
  Real SOC_tank3;
  //Real SOE_tank3;
  ColdTESLiquid.BaseProperties coldliq_tank3_base "Base Medium properties tank 3";
  SI.PathLength x_tank3 "cold liquid level";
  //port
  parameter SI.MassFlowRate m_out_tank3 = 0;
  SI.MassFlowRate m_in_tank3;
  parameter SI.SpecificEnthalpy h_in_tank3 = 16974;
  //dummy value TO_DO
  //SI.HeatFlowRate Q_port_tank1(displayUnit = "MW") "Heat Flow from tank 2";
  parameter SI.Height x_min_tank3 = 0.4;
  //--------------------------TANK 4
  //nominal
  parameter SI.Temperature T_nom_tank4 = from_degC(-60);
  parameter SI.Pressure p_nom_tank4 = 101325 "unpressurized tank";
  ColdTESLiquid.ThermodynamicState coldliq_tank4_nom(p(start = p_nom_tank4), T(start = T_nom_tank4)) "Nominal thermodynamic state";
  //SI.Energy Energy_nom_tank4(displayUnit = "MWh");
  parameter SI.Mass m_nom_tank4 = V_tank_coldliq*rho_start_tank4;
  //start values
  parameter SI.Temperature T_start_tank4 = from_degC(-60);
  parameter SI.Mass m_start_tank4 = m_nom_tank4;
  parameter SI.Density rho_start_tank4 = 900;
  //actual
  SI.Temperature T_tank4(start = T_start_tank4);
  ColdTESLiquid.ThermodynamicState coldliq_tank4_state(p(start = p_nom_tank4), T(start = T_nom_tank4)) "Nominal thermodynamic state";
  SI.HeatFlowRate Q_dot_to_amb_tank4(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank4;
  //state variables
  SI.Mass m_tank4(start = m_nom_tank4);
  Real SOC_tank4;
  // Real SOE_tank4;
  ColdTESLiquid.BaseProperties coldliq_tank4_base "Base Medium properties tank 4";
  //other
  // SI.Energy Energy_content_tank1(displayUnit = "MWh");
  SI.PathLength x_tank4 "salt-level";
  //outflow
  SI.MassFlowRate m_out_tank4;
  parameter SI.MassFlowRate m_in_tank4 = 0;
  parameter SI.SpecificEnthalpy h_in_tank4 = 1;
  //dummy value TO_DO
  //SI.HeatFlowRate Q_port_tank4(displayUnit = "MW") "Heat Flow from tank 2";
  parameter SI.Height x_min_tank4 = 0.4;
  //-------------ALL HEX//
  parameter Real f_p = 0.01625;
  //pressure loss factor
  //-------------HEX HOT//
  parameter Real UA_HEX_Hot = 17456614;
  Real C_WF_HEX_Hot "Heat capacity rate of cold side of hot HEX";
  Real C_solsalt_HEX_Hot "Heat capacity rate of hot side of hot HEX";
  SI.MassFlowRate m_solsalt_HEX_Hot "mass flow rate required for balanced HEX";
  Real NTU_HEX_Hot "Number of transfer units of hot HEX";
  Real eff_HEX_Hot "effectiveness of hot HEX";
  HotTESLiquid.BaseProperties solsalt_HEX1_tank2(T(start = 570)) "Medium properties of HEX port at interface to tank 2";
  //-------------HEX RECUP//
  parameter Real UA_HEX_recup = 25443162;
  Real C_4_HEX_recup "Heat capacity rate of hot side of recuperation HEX (after turbine)";
  Real C_2_HEX_recup "Heat capacity rate of cold side of recuperation HEX (after compressor)";
  Real C_min_HEX_recup;
  Real NTU_HEX_recup "Number of transfer units of recuperation HEX";
  Real eff_HEX_recup "effectiveness of hot HEX";
  Real C_max_HEX_recup;
  //-------------HEAT REJECTION//
  SI.HeatFlowRate Q_4_a_1_a "heat flow rate rejected to district heating or such";
  //-------------HEX COLD//
  parameter Real UA_HEX_Cold = 9500395;
  Real C_T_1_a_HEX_Cold "Heat capacity rate of hot side of cold HEX (coming from heat rejection)";
  Real C_tank4_HEX_Cold "Heat capacity rate of cold side of cold HEX (coming from cold tank 4)";
  Real NTU_HEX_Cold "Number of transfer units of cold HEX";
  Real eff_HEX_Cold "effectiveness of cold HEX";
  SI.MassFlowRate m_coldliq_HEX_Cold "mass flow rate required for balanced HEX";
  ColdTESLiquid.ThermodynamicState coldliq_HEX_in_tank3(p(start = p_nom_tank3), T(start = T_nom_tank3)) " thermodynamic state";
  //SI.Temperature T_1_;
  //-------------COMPRESSOR//
  parameter Real k = 1.4 "isentropic expansion factor";
  parameter Real p = 1.8 "compressor map factor";
  parameter Real m = 1.4 "compressor map factor";
  parameter Real c4 = 0.3 "factor 4, see Zhang2002";
  Real c1 "factor 1";
  Real c2 "factor 2";
  Real c3 "factor 3";
  //design
  parameter Real beta_CO_nom = 5.911 "design compression ratio";
  parameter Real n_CO_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom = 0.88 "design isentropic efficiency";
  //actual
  Real beta_CO(start = beta_CO_nom) "absolute compression ratio";
  parameter Real n_CO = 3000 "actual speed";
  SI.Efficiency eta_is_CO "absolute isentropic efficiency";
  //reduced
  Real beta_CO_red(start = 1) "reduced compression ratio";
  Real n_CO_red(start = 1) "reduced speed";
  SI.Efficiency eta_is_CO_red "reduced isentropic efficiency";
  //-------------TURBINE//
  parameter Real t = 0.3 "parameter, see Zhang2002";
  Real alpha "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom = 5.445 "design expansion ratio";
  parameter Real n_TU_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom = 0.93 "design isentropic efficiency";
  //actual
  Real beta_TU(start = beta_TU_nom) "absolute expansion ratio";
  parameter Real n_TU = 3000 "actual speed";
  SI.Efficiency eta_is_TU "absolute isentropic turbine efficiency";
  //reduced
  Real beta_TU_red(start = 1) "reduced expansion ratio";
  Real n_TU_red(start = 1) "reduced speed";
  SI.Efficiency eta_is_TU_red "reduced isentropic turbine efficiency";
  //-------------STATES//
  // INITIAL
  parameter SI.Pressure p_1_initial = 100000 "fixed pressure point";
  parameter SI.Temperature T_1_initial = from_degC(-52.760) "ambient temperature";
  parameter SI.Temperature T_2_initial = from_degC(113.6) "ambient temperature";
  parameter SI.SpecificEnthalpy h_1_initial = -53000;
  parameter SI.SpecificEntropy s_1_initial = -50;
  WorkingFluid.ThermodynamicState state_initial "thermodynamic state initially";
  parameter SI.Temperature T_3_a_initial = from_degC(263.38);
  parameter SI.Temperature T_4_a_initial = from_degC(121.64);
  parameter SI.Temperature T_4_initial = from_degC(272);
  parameter SI.Pressure p_4_initial = 105039;
  parameter SI.Pressure p_4_a_initial = 103276;
  parameter WorkingFluid.ThermodynamicState state_4_initial(p = p_4_initial, T = T_4_initial) "thermodynamic state of turbine outlet";
  SI.Temperature T_1(start = T_1_initial) "ambient temperature";
  SI.Pressure p_1(start = p_1_initial) "ambient pressure";
  WorkingFluid.ThermodynamicState state_1 "thermodynamic state of inlet";
  WorkingFluid.SpecificEnthalpy h_1(start = h_1_initial) "inlet enthalpy";
  WorkingFluid.SpecificEntropy s_1(start = s_1_initial) "inlet spec. entropy";
  //STATE 2 isentropic
  SI.Temperature T_2_is(start = from_degC(100)) "isentropic outlet temperature of compressor";
  WorkingFluid.ThermodynamicState state_2_is "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEntropy s_2_is "compressor outlet spec. entropy";
  //STATE 2
  SI.Temperature T_2(start = from_degC(114)) "actual outlet temperature of compressor";
  SI.Pressure p_2(start = 591100) "pressure coming out of compressor";
  WorkingFluid.ThermodynamicState state_2 "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_2(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_2 "compressor outlet spec. entropy";
  //STATE 3a
  WorkingFluid.ThermodynamicState state_3_a "thermodynamic state of compressor outlet";
  SI.Temperature T_3_a(start = from_degC(260)) "outlet temperature after recuperation";
  WorkingFluid.SpecificEnthalpy h_3_a(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a "compressor outlet spec. entropy";
  SI.Pressure p_3_a(start = 581494) "Pressure after recup";
  //STATE 3
  WorkingFluid.ThermodynamicState state_3 "thermodynamic state of HEX outlet";
  SI.Temperature T_3(start = from_degC(560)) "outlet temperature after HEX";
  SI.Pressure p_3(start = 572045) "Pressure after HEX";
  WorkingFluid.SpecificEnthalpy h_3(start = 579480) "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3 "HEX outlet spec. entropy";
  //STATE 4 isentropic
  SI.Temperature T_4_is(start = from_degC(250)) "isentropic outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4_is "isentropic state of turbine outlet";
  WorkingFluid.SpecificEntropy s_4_is "isentropic turbine outlet spec. entropy";
  //STATE 4
  SI.Pressure p_4(start = 105039) " pressure at turb outlet";
  SI.Temperature T_4 "outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4 "thermodynamic state of turbine outlet";
  WorkingFluid.SpecificEnthalpy h_4 "turbine outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";
  //STATE 4 a
  SI.Pressure p_4_a(start = 103332) " pressure at recuperation outlet";
  SI.Temperature T_4_a(start = from_degC(270)) "outlet temperature after recuperation";
  WorkingFluid.ThermodynamicState state_4_a(p(start = p_4_initial), T(start = T_4_initial)) "thermodynamic state of turbine-side recuperation outlet";
  WorkingFluid.SpecificEnthalpy h_4_a "turbine-side recuperation outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_a "turbine-side recuperation outlet spec. entropy";
  //STATE 1 a
  SI.Pressure p_1_a "pressure after Heat rejection ";
  SI.Temperature T_1_a "outlet temperature after Heat rejection ";
  WorkingFluid.ThermodynamicState state_1_a "thermodynamic state after Heat rejection ";
  WorkingFluid.SpecificEnthalpy h_1_a "turbine-side recuperation after Heat rejection ";
  WorkingFluid.SpecificEntropy s_1_a "turbine-side recuperation after Heat rejection ";
  //-------------CYCLE//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  parameter SI.Temperature T_3_nom = from_degC(555) "design turbine inlet temperature";
  parameter SI.Pressure p_3_nom = 572336 "design turbine inlet pressure";
  parameter SI.Temperature T_1_nom = from_degC(-52.760) "design compressor inlet temperature";
  parameter SI.Pressure p_1_nom = 100000 "design compressor inlet pressure";
  //actual
  parameter SI.MassFlowRate m_dot_WF = 762 "mass flow rate";
  Real G_CO_red(start = 1) "reduced mass flow rate compressor";
  Real G_TU_red(start = 1) "reduced mass flow rate turbine";
  SI.Power P_12(displayUnit = "MW") "Driving Power for compressor";
  SI.Power P_34(displayUnit = "MW") "Power delivered by turbine";
  SI.HeatFlowRate Q_dot_3a3(displayUnit = "MW") "Heat flow rate added to working fluid by Hot-TES-HEX";
  SI.Power P_shaft(displayUnit = "MW") "Shaft Power available to electrical machinery";
  SI.Efficiency eta_cycle "cycle efficiency Q_therm/P_shaft";
  Real back_work_ratio "ratio of compressor power to turbine power [0.4...0.8 open cycle]";
  Real work_ratio "ratio of turbine power to compressor power ";
  SI.Energy E_dis(displayUnit = "MW.h", start = 0, fixed = true);
  //--------------------------EQUATIONS-----------------------------//
initial equation
  m_tank1 = m_nom_tank1;
  T_tank1 = T_start_tank1;
  p_1 = p_1_initial;
  m_tank2 = m_empty_tank2;
  T_tank2 = T_start_tank2;
  T_tank3 = T_start_tank3;
  T_tank4 = T_start_tank4;
  m_tank3 = m_start_tank3;
equation
//-------------TANK 1//
  solsalt_tank1_nom = HotTESLiquid.setState_pT(p_nom_tank1, T_nom_tank1);
  solsalt_tank1_state = HotTESLiquid.setState_pT(p_nom_tank1, T_tank1);
  Energy_nom_tank1 = m_nom_tank1*solsalt_tank1_nom.cp*T_nom_tank1;
//solar salt properties
  solsalt_tank1.p = p_nom_tank1;
  solsalt_tank1.T = T_tank1;
//state description
  m_tank1 = solsalt_tank1.d*A_cross_tank_solsalt*x_tank1;
// get tank level
  SOC_tank1 = m_tank1/m_nom_tank1;
  SOE_tank1 = Energy_content_tank1/Energy_nom_tank1;
  Energy_content_tank1 = m_tank1*solsalt_tank1.u;
//mass balance
  der(m_tank1) = m_in_tank1 + m_out_tank1;
//energy balance
  der(m_tank1*solsalt_tank1.u) = m_in_tank1*h_in_tank1 + m_out_tank1*solsalt_tank1.h + Q_dot_to_amb_tank1;
  Q_div_A_tank1 = (0.00017*(T_tank1 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank1 = -Q_div_A_tank1*A_W_tank_solsalt;
  Q_port_tank1 = m_in_tank1*h_in_tank1 + m_out_tank1*solsalt_tank1.h;
//-------------TANK 2//
  solsalt_tank2_nom = HotTESLiquid.setState_pT(p_nom_tank2, T_nom_tank2);
  solsalt_tank2_state = HotTESLiquid.setState_pT(p_nom_tank2, T_tank2);
  Energy_nom_tank2 = m_nom_tank2*solsalt_tank2_nom.cp*T_nom_tank2;
//solar salt properties
  solsalt_tank2.p = p_nom_tank2;
  solsalt_tank2.T = T_tank2;
//state description
  m_tank2 = solsalt_tank2.d*A_cross_tank_solsalt*x_tank2;
// get tank level
  SOC_tank2 = m_tank2/m_nom_tank2;
  SOE_tank2 = Energy_content_tank2/Energy_nom_tank2;
  Energy_content_tank2 = m_tank1*solsalt_tank2.u;
//mass balance
  der(m_tank2) = m_in_tank2 + m_out_tank2;
//energy balance
//der(m_tank2*cp_tank2*T_tank2) = m_in_tank2*cp_tank2*T_in_tank2 + m_out_tank2*h_tank2 + Q_dot_to_amb_tank2;
  der(m_tank2*solsalt_tank2.u) = m_in_tank2*h_in_tank2 + m_out_tank2*solsalt_tank2.h + Q_dot_to_amb_tank2;
  Q_div_A_tank2 = (0.00017*(T_tank2 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank2 = -Q_div_A_tank2*x_tank2*pi*D_tank_solsalt;
  Q_port_tank2 = m_in_tank2*h_in_tank2 + m_out_tank2*solsalt_tank2.h;
//-------------TANK 3//
  coldliq_tank3_nom = ColdTESLiquid.setState_pT(p_nom_tank3, T_nom_tank3);
  coldliq_tank3_state = ColdTESLiquid.setState_pT(p_nom_tank3, T_tank3);
//Energy_nom_tank1 = m_nom_tank1*solsalt_tank1_nom.cp*T_nom_tank1;
//cold liquid properties
  coldliq_tank3_base.p = p_nom_tank3;
  coldliq_tank3_base.T = T_tank3;
//state description
  m_tank3 = coldliq_tank3_base.d*A_cross_tank_coldliq*x_tank3;
// get tank level
  SOC_tank3 = m_tank3/m_nom_tank3;
// SOE_tank1 = Energy_content_tank1/Energy_nom_tank1;
//Energy_content_tank1 = m_tank1*solsalt_tank1.u;
//mass balance
  der(m_tank3) = m_in_tank3 + m_out_tank3;
//energy balance
  der(m_tank3*coldliq_tank3_base.u) = m_in_tank3*coldliq_HEX_in_tank3.h + m_out_tank3*coldliq_tank3_base.h + Q_dot_to_amb_tank3;
  Q_div_A_tank3 = 1;
//dummyvalue
  Q_dot_to_amb_tank3 = 0;
//Q_port_tank1 = m_in_tank1*h_in_tank1 + m_out_tank1*solsalt_tank1.h;
//-------------TANK 4//
// T_tank4=from_degC(-60);
  coldliq_tank4_nom = ColdTESLiquid.setState_pT(p_nom_tank4, T_nom_tank4);
  coldliq_tank4_state = ColdTESLiquid.setState_pT(p_nom_tank4, T_tank4);
//Energy_nom_tank1 = m_nom_tank1*solsalt_tank1_nom.cp*T_nom_tank1;
//cold liquid properties
  coldliq_tank4_base.p = p_nom_tank4;
  coldliq_tank4_base.T = T_tank4;
//state description
  m_tank4 = coldliq_tank4_base.d*A_cross_tank_coldliq*x_tank4;
// get tank level
  SOC_tank4 = m_tank4/m_nom_tank4;
// SOE_tank1 = Energy_content_tank1/Energy_nom_tank1;
//Energy_content_tank1 = m_tank1*solsalt_tank1.u;
//mass balance
  der(m_tank4) = m_in_tank4 + m_out_tank4;
//energy balance
  der(m_tank4*coldliq_tank4_base.u) = m_in_tank4*h_in_tank4 + m_out_tank4*coldliq_tank4_base.h + Q_dot_to_amb_tank4;
  Q_div_A_tank4 = 1;
//dummy value
  Q_dot_to_amb_tank4 = 0;
//Q_port_tank1 = m_in_tank1*h_in_tank1 + m_out_tank1*solsalt_tank1.h;
//-------------HEX HOT//
  C_WF_HEX_Hot = m_dot_WF*WorkingFluid.specificHeatCapacityCp(state_3_a);
  C_solsalt_HEX_Hot = C_WF_HEX_Hot;
  m_solsalt_HEX_Hot = C_solsalt_HEX_Hot/solsalt_tank1_state.cp;
  m_out_tank1 + m_solsalt_HEX_Hot = 0;
  m_solsalt_HEX_Hot = m_in_tank2;
  NTU_HEX_Hot = UA_HEX_Hot/C_WF_HEX_Hot;
  eff_HEX_Hot = NTU_HEX_Hot/(1 + NTU_HEX_Hot);
  eff_HEX_Hot = (T_tank1 - T_in_tank2)/(T_tank1 - T_3_a);
//calc T_tank2
  (T_in_tank2 - T_tank1) = (T_3_a - T_3);
//calc T_3 //because C_min=C_max (balanced HEX)
//-------------HEX RECUP//
  C_4_HEX_recup = m_dot_WF*WorkingFluid.specificHeatCapacityCp(state_4);
//hotstate_4.cp
  C_2_HEX_recup = m_dot_WF*WorkingFluid.specificHeatCapacityCp(state_2);
//coldstate_2.cp
  C_min_HEX_recup = min(C_4_HEX_recup, C_2_HEX_recup);
  C_max_HEX_recup = max(C_4_HEX_recup, C_2_HEX_recup);
  NTU_HEX_recup = UA_HEX_recup/C_min_HEX_recup;
  eff_HEX_recup = NTU_HEX_recup/(1 + NTU_HEX_recup);
  eff_HEX_recup = (C_4_HEX_recup/C_min_HEX_recup)*(T_4 - T_4_a)/(T_4 - T_2);
//calc T_4_a
  C_min_HEX_recup/C_max_HEX_recup = (T_3_a - T_2)/(T_4 - T_4_a);
  solsalt_HEX1_tank2.T = T_in_tank2;
//To-Do unneccessary state
  solsalt_HEX1_tank2.p = p_amb;
//To-Do unneccessary state
  h_in_tank2 = solsalt_HEX1_tank2.h;
//-------------HEAT REJECTION//
  Q_4_a_1_a = m_dot_WF*(h_4_a - h_1_a);
  T_1_a = T_4_a - 89;
//-------------HEX COLD//
  C_T_1_a_HEX_Cold = m_dot_WF*WorkingFluid.specificHeatCapacityCp(state_1_a);
//state_1_a.cp
  C_tank4_HEX_Cold = C_T_1_a_HEX_Cold;
//balanced HEX
  NTU_HEX_Cold = UA_HEX_Cold/C_T_1_a_HEX_Cold;
  eff_HEX_Cold = NTU_HEX_Cold/(1 + NTU_HEX_Cold);
  eff_HEX_Cold = (T_1_a - T_1)/(T_1_a - T_tank4);
//calc T_1_
  (T_in_tank3 - T_tank4) = (T_1_a - T_1);
//calc T_in_tank3 //because C_min=C_max (balanced HEX)
  m_coldliq_HEX_Cold = C_T_1_a_HEX_Cold/coldliq_tank4_state.cp;
  m_out_tank4 + m_coldliq_HEX_Cold = 0;
  m_coldliq_HEX_Cold = m_in_tank3;
  coldliq_HEX_in_tank3 = ColdTESLiquid.setState_pT(p_nom_tank3, T_in_tank3);
//-------------STATES//
  state_initial = WorkingFluid.setState_pT(p_1_initial, T_1_initial);
//state 1
  p_1 = p_1_initial;
//fixed pressure point
//T_1 = T_1_initial;
  state_1 = WorkingFluid.setState_pT(p_1, T_1);
  h_1 = WorkingFluid.specificEnthalpy(state_1);
  s_1 = WorkingFluid.specificEntropy(state_1);
//state 2 isentropic
  s_2_is = s_1 "isentropic compressor outlet spec. entropy";
  state_2_is = WorkingFluid.setState_ps(p_2, s_2_is) "isentropic state of compressor outlet";
  T_2_is = state_2_is.T;
//state 2
  state_2 = WorkingFluid.setState_pT(p_2, T_2);
  h_2 = WorkingFluid.specificEnthalpy(state_2);
  s_2 = WorkingFluid.specificEntropy(state_2);
//state 3a
  p_3_a = p_2*(1 - f_p);
  state_3_a = WorkingFluid.setState_pT(p_2, T_3_a);
  h_3_a = WorkingFluid.specificEnthalpy(state_3_a);
  s_3_a = WorkingFluid.specificEntropy(state_3_a);
//state 3
  state_3 = WorkingFluid.setState_pT(p_3, T_3);
  p_3 = p_3_a*(1 - f_p);
  h_3 = state_3.h;
  s_3 = WorkingFluid.specificEntropy(state_3);
//state 4 isentropic
  s_4_is = s_3;
  state_4_is = WorkingFluid.setState_ps(p_4, s_4_is) "isentropic state of turbine outlet";
  T_4_is = state_4_is.T "isentropic turbine outlet spec. entropy";
//state 4
  p_4 = p_4_a*(1 + f_p);
  state_4 = WorkingFluid.setState_pT(p_4, T_4);
  h_4 = WorkingFluid.specificEnthalpy(state_4);
  s_4 = WorkingFluid.specificEntropy(state_4);
//state 4_a
  p_4_a = p_1_a*(1 + f_p);
  state_4_a = WorkingFluid.setState_pT(p_4_a, T_4_a);
  h_4_a = WorkingFluid.specificEnthalpy(state_4_a);
  s_4_a = WorkingFluid.specificEntropy(state_4_a);
//state 1_a
  p_1_a = p_1*(1 + f_p);
  state_1_a = WorkingFluid.setState_pT(p_1_a, T_1_a);
  h_1_a = WorkingFluid.specificEnthalpy(state_1_a);
  s_1_a = WorkingFluid.specificEntropy(state_1_a);
//-------------COMPRESSOR//
//reduced values compressor
  eta_is_CO_red = eta_is_CO/eta_is_CO_nom;
  beta_CO_red = beta_CO/beta_CO_nom;
  G_CO_red = m_dot_WF*sqrt(T_1)/p_1/(m_dot_WF_nom*sqrt(T_1_nom)/p_1_nom);
  eta_is_CO_red = (1 - c4*(1 - n_CO_red)^2)*(n_CO_red/G_CO_red)*(2 - n_CO_red/G_CO_red);
  n_CO_red = n_CO/sqrt(T_1)/(n_CO_nom/sqrt(T_1_nom));
//other compressor equations
  beta_CO = p_2/p_1;
  eta_is_CO = (T_2_is - T_1)/(T_2 - T_1);
  P_12 = m_dot_WF*(h_2 - h_1);
  beta_CO_red = c1*G_CO_red^2 + c2*G_CO_red + c3;
  c1 = n_CO_red/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c2 = (p - 2*m*n_CO_red^2)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c3 = -1*(p*m*n_CO_red - m^2*n_CO_red^3)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
//-------------TURBINE//
//reduced values turbine
  eta_is_TU_red = eta_is_TU/eta_is_TU_nom;
  n_TU_red = n_TU/sqrt(T_3)/(n_TU_nom/sqrt(T_3_nom));
//other turbine equations
  alpha = sqrt(1.4 - 0.4*n_TU_red);
  beta_TU = p_3/p_4;
  beta_TU_red = beta_TU/beta_TU_nom;
  G_TU_red = alpha*sqrt(T_3_nom/T_3)*sqrt((beta_TU^2 - 1)/(beta_TU_nom^2 - 1));
  eta_is_TU_red = (1 - t*(1 - n_TU_red)^2)*(n_TU_red/G_TU_red)*(2 - ((n_TU_red/G_TU_red)));
  eta_is_TU = (T_3 - T_4)/(T_3 - T_4_is);
  P_34 = m_dot_WF*(h_3 - h_4);
//-------------CYCLE//
  P_shaft = P_34 - P_12;
  eta_cycle = P_shaft/Q_dot_3a3;
  Q_dot_3a3 = m_dot_WF*(h_3 - h_3_a);
  back_work_ratio = P_12/P_34;
  work_ratio = P_34/P_12;
//work_per_unit_mass = P_shaft/(10^6)/m_dot_WF;
  der(E_dis) = P_shaft;
//--------------------------COMPONENT MANAGEMENT--------------------------//
  if x_tank1 < x_min_tank1 then
    terminate("Minimum fill level of tank reached");
  end if;
end Working_flat;
