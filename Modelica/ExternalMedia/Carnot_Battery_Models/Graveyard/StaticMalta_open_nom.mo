within ExternalMedia.Carnot_Battery_Models.Graveyard;

model StaticMalta_open_nom
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
  //package WorkingFluid = Modelica.Media.Air.ReferenceAir.Air_pT;

  package Methanol "NaK properties from CoolProp"
    extends ExternalMedia.Media.IncompressibleCoolPropMedium(mediumName = "MMA", substanceNames = {"MMA[0.6]"});
  end Methanol;

  replaceable package ColdTESLiquid = Methanol constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
  package WorkingFluid = ExternalMedia.Media.CoolPropMedium(mediumName = "Air", substanceNames = {"Air"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
  //--------------------------PARAMETERS & VARIABLES-----------------------------//
  //-------------CYCLE//
  //-------------Charge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom_charge = 766 "design mass flow rate";
  //actual
  parameter SI.MassFlowRate m_dot_WF_charge = 766 "mass flow rate";
  parameter Real f_p_charge = 0.01075 "pressure_loss_factor percent";
  parameter SI.Pressure p_fix_charge = 97889 "fixed pressure point through expansion vessel at p_4";
  //-------------Discharge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  //actual
  parameter SI.MassFlowRate m_dot_WF = 762 "mass flow rate";
  parameter Real f_p = 0.01625 "pressure_loss_factor percent";
  parameter SI.Pressure p_fix = 100000 "fixed pressure point through expansion vessel";
  //--------------------------TANK 1
  //nominal
  parameter SI.Temperature T_nom_tank1 = from_degC(565);
  parameter SI.Temperature T_tank1 = from_degC(565);
  parameter SI.Pressure p_nom_tank1 = 101325 "unpressurized tank";
  HotTESLiquid.ThermodynamicState solsalt_tank1_nom(p(start = p_nom_tank1), T(start = T_nom_tank1)) "Nominal thermodynamic state";
  HotTESLiquid.ThermodynamicState solsalt_tank1_state(p(start = p_nom_tank1), T(start = T_nom_tank1)) "Nominal thermodynamic state";
  HotTESLiquid.BaseProperties solsalt_tank1(p(start = p_nom_tank1), T(start = T_nom_tank1)) "Medium properties of port_a";
  //--------------------------TANK 2
  //nominal
  parameter SI.Temperature T_nom_tank2 = from_degC(271.5);
  parameter SI.Temperature T_tank2 = from_degC(271.5);
  parameter SI.Pressure p_nom_tank2 = 101325 "unpressurized tank";
  HotTESLiquid.ThermodynamicState solsalt_tank2_nom(p(start = p_nom_tank2), T(start = T_nom_tank2)) "Nominal thermodynamic state";
  HotTESLiquid.ThermodynamicState solsalt_tank2_state(p(start = p_nom_tank2), T(start = T_nom_tank2)) "Nominal thermodynamic state";
  HotTESLiquid.BaseProperties solsalt_tank2(p(start = p_nom_tank2), T(start = T_nom_tank2)) "Medium properties of port_a";
  //--------------------------TANK 3
  //nominal
  parameter SI.Temperature T_nom_tank3 = from_degC(25.195);
  parameter SI.Temperature T_tank3 = from_degC(25.195);
  parameter SI.Pressure p_nom_tank3 = 101325 "unpressurized tank";
  ColdTESLiquid.ThermodynamicState coldliq_tank3_nom(p(start = p_nom_tank3), T(start = T_nom_tank3)) "Nominal thermodynamic state";
  ColdTESLiquid.ThermodynamicState coldliq_tank3_state(p(start = p_nom_tank3), T(start = T_nom_tank3)) "Nominal thermodynamic state";
  //--------------------------TANK 4
  //nominal
  parameter SI.Temperature T_nom_tank4 = from_degC(-59.85);
  parameter SI.Temperature T_tank4 = from_degC(-59.85);
  parameter SI.Pressure p_nom_tank4 = 101325 "unpressurized tank";
  ColdTESLiquid.ThermodynamicState coldliq_tank4_nom(p(start = p_nom_tank4), T(start = T_nom_tank4)) "Nominal thermodynamic state";
  ColdTESLiquid.ThermodynamicState coldliq_tank4_state(p(start = p_nom_tank4), T(start = T_nom_tank4)) "Nominal thermodynamic state";
  //-------------HEX 1//
  //design
  parameter Real UA_HEX1_nom = 25906701;
  SI.MassFlowRate m_dot_solsalt_HEX1_nom = 553.021;
  //actual
  Real UA_HEX1_charge;
  SI.MassFlowRate m_dot_solsalt_HEX1_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  HotTESLiquid.ThermodynamicState outlet_hotside_guess_HEX1_charge(T(start = T_nom_tank2)) "Medium properties of HEX port at interface to tank 2";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX1_charge(T(start = T_nom_tank1)) "Medium properties of HEX port at interface to tank 2";
  SI.SpecificHeatCapacity cp_hot_ave_HEX1_charge;
  SI.SpecificHeatCapacity cp_cold_ave_HEX1_charge;
  Real C_cold_HEX1_charge "Heat capacity rate of cold side of hot HEX";
  Real C_hot_HEX1_charge "Heat capacity rate of hot side of hot HEX";
  Real C_min_HEX1_charge;
  Real C_max_HEX1_charge;
  Real C_r_HEX1_charge;
  Real NTU_HEX1_charge "Number of transfer units of hot HEX";
  Real eff_HEX1_charge "effectiveness of hot HEX";
  SI.HeatFlowRate Q_dot_HEX1_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX1_charge(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX1_charge;
  SI.SpecificEnthalpy h_cold_out_HEX1_charge;
  WorkingFluid.ThermodynamicState outlet_hotside_HEX1_charge(T(start = T_nom_tank2)) "Medium properties of HEX port at interface to tank 2";
  HotTESLiquid.ThermodynamicState outlet_coldside_HEX1_charge(T(start = T_nom_tank1)) "Medium properties of HEX port at interface to tank 1";
  //-------------HEX 2//
  //design
  parameter Real UA_HEX2_nom = 13279998;
  //actual
  Real UA_HEX2_charge;
  //cold-side
  Real C_cold_HEX2_charge "Heat capacity rate of cold side of recuperation HEX ";
  SI.SpecificHeatCapacity cp_cold_HEX2_charge;
  //hot-side
  Real C_hot_HEX2_charge "Heat capacity rate of cold side of recuperation HEX (after compressor)";
  SI.SpecificHeatCapacity cp_hot_HEX2_charge;
  //variables for effectiveness
  Real C_min_HEX2_charge;
  Real C_max_HEX2_charge;
  Real C_r_HEX2_charge;
  Real NTU_HEX2_charge "Number of transfer units of recuperation HEX";
  Real eff_HEX2_charge "effectiveness of hot HEX";
  //heat transferred
  SI.HeatFlowRate Q_dot_HEX2_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX2_charge(displayUnit = "MW");
  //energy balance and outlet states
  SI.SpecificEnthalpy h_hot_out_HEX2_charge;
  SI.SpecificEnthalpy h_cold_out_HEX2_charge;
  WorkingFluid.ThermodynamicState outlet_hotside_HEX2_charge(T(start = 309)) "Medium properties";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX2_charge(T(start = 540)) "Medium properties";
  //-------------HEX 3//
  //design
  parameter Real UA_HEX3_nom = 9505824;
  parameter SI.MassFlowRate m_dot_methanol_HEX3_nom = 311.770;
  //actual
  Real UA_HEX3_charge;
  SI.MassFlowRate m_dot_methanol_HEX3_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  ColdTESLiquid.ThermodynamicState outlet_hotside_guess_HEX3_charge(T(start = T_nom_tank1)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX3_charge(T(start = T_nom_tank3)) "Medium properties ";
  SI.SpecificHeatCapacity cp_hot_ave_HEX3_charge;
  SI.SpecificHeatCapacity cp_cold_ave_HEX3_charge;
  Real C_cold_HEX3_charge "Heat capacity rate of cold side of hot HEX";
  Real C_hot_HEX3_charge "Heat capacity rate of hot side of hot HEX";
  Real C_min_HEX3_charge;
  Real C_max_HEX3_charge;
  Real C_r_HEX3_charge;
  Real NTU_HEX3_charge "Number of transfer units of hot HEX";
  Real eff_HEX3_charge "effectiveness of hot HEX";
  SI.HeatFlowRate Q_dot_HEX3_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX3_charge(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX3_charge;
  SI.SpecificEnthalpy h_cold_out_HEX3_charge;
  ColdTESLiquid.ThermodynamicState outlet_hotside_HEX3_charge(T(start = T_nom_tank4)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX3_charge(T(start = T_nom_tank3)) "Medium properties ";
  //-------------COMPRESSOR CHARGE//
  parameter Real k = 1.4 "isentropic expansion factor";
  parameter Real p = 1.8 "compressor map factor";
  parameter Real m = 1.4 "compressor map factor";
  parameter Real c4 = 0.3 "factor 4, see Zhang2002";
  Real c1_charge "factor 1";
  Real c2_charge "factor 2";
  Real c3_charge "factor 3";
  //design
  parameter Real beta_CO_nom_charge = 4.592 "design compression ratio";
  parameter Real n_CO_nom_charge = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom_charge = 0.88 "design isentropic efficiency";
  SI.Temperature T_4_nom_charge = from_degC(267.533) "state 4 temperature";
  SI.Pressure p_4_nom_charge = p_fix_charge "state 4 pressure";
  //actual
  Real beta_CO_charge(start = beta_CO_nom_charge) "absolute compression ratio";
  parameter Real n_CO_charge = 3000 "actual speed";
  SI.Efficiency eta_is_CO_charge "absolute isentropic efficiency";
  //reduced
  Real beta_CO_red_charge(start = 1) "reduced compression ratio";
  Real n_CO_red_charge(start = 1) "reduced speed";
  SI.Efficiency eta_is_CO_red_charge "reduced isentropic efficiency";
  Real G_CO_red_charge(start = 1) "reduced mass flow rate compressor";
  //-------------EXPANDER CHARGE//
  parameter Real t = 0.3 "parameter, see Zhang2002";
  Real alpha_charge "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom_charge = 4.331 "design expansion ratio";
  parameter Real n_TU_nom_charge = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom_charge = 0.93 "design isentropic efficiency";
  SI.Temperature T_2_nom_charge = 309.3 "state temperature";
  SI.Pressure p_2_nom_charge = 433100 "state pressure";
  //actual
  parameter Real n_TU_charge = 3000 "actual speed";
  Real beta_TU_charge(start = beta_TU_nom_charge) "absolute expansion ratio";
  SI.Efficiency eta_is_TU_charge "absolute isentropic turbine efficiency";
  //reduced
  Real n_TU_red_charge(start = 1) "reduced speed";
  Real G_TU_red_charge(start = 1) "reduced mass flow rate turbine";
  Real beta_TU_red_charge(start = 1) "reduced expansion ratio";
  SI.Efficiency eta_is_TU_red_charge "reduced isentropic turbine efficiency";
  /*
    //-------------COMPRESSOR//
    parameter Real k = 1.4 "isentropic expansion factor";
    parameter Real p = 1.8 "compressor map factor";
    parameter Real m = 1.4 "compressor map factor";
    parameter Real c4 = 0.3 "factor 4, see Zhang2002";
    Real c1 "factor 1";
    Real c2 "factor 2";
    Real c3 "factor 3";
    //design
    parameter Real beta_CO_nom = 6 "design compression ratio";
    parameter Real n_CO_nom = 1500 "design speed";
    parameter SI.Efficiency eta_is_CO_nom = 0.9 "design isentropic efficiency";
    //actual
    Real beta_CO(start = beta_CO_nom) "absolute compression ratio";
    parameter Real n_CO = 1500 "actual speed";
    SI.Efficiency eta_is_CO "absolute isentropic efficiency";
    //reduced
    Real beta_CO_red(start = 1) "reduced compression ratio";
    Real n_CO_red(start = 1) "reduced speed";
    SI.Efficiency eta_is_CO_red "reduced isentropic efficiency";
    //-------------TURBINE//
    parameter Real t = 0.3 "parameter, see Zhang2002";
    Real alpha "factor, see Zhang2002";
    //design
    parameter Real beta_TU_nom = beta_CO_nom "design expansion ratio";
    parameter Real n_TU_nom = 1500 "design speed";
    parameter SI.Efficiency eta_is_TU_nom = 0.9 "design isentropic efficiency";
    //actual
    Real beta_TU(start = beta_TU_nom) "absolute expansion ratio";
    parameter Real n_TU = 1500 "actual speed";
    SI.Efficiency eta_is_TU "absolute isentropic turbine efficiency";
    //reduced
    Real beta_TU_red(start = 1) "reduced expansion ratio";
    Real n_TU_red(start = 1) "reduced speed";
    SI.Efficiency eta_is_TU_red "reduced isentropic turbine efficiency"; 
    */
  //-------------STATES//
  //state 4 charge
  SI.Temperature T_4_charge = from_degC(267.533) "state 4 temperature";
  //SI.Temperature T_4_charge(start=from_degC(267.533)) "state 4 temperature";
  SI.Pressure p_4_charge = p_fix_charge "state 4 pressure";
  WorkingFluid.ThermodynamicState state_4_charge "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_charge "enthalpy";
  WorkingFluid.SpecificEntropy s_4_charge "spec. entropy";
  //STATE 3 isentropic charge
  SI.Temperature T_3_is_charge "isentropic outlet temperature of compressor";
  WorkingFluid.ThermodynamicState state_3_is_charge "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEntropy s_3_is_charge "compressor outlet spec. entropy";
  //STATE 3 charge
  SI.Temperature T_3_charge "actual outlet temperature of compressor";
  SI.Pressure p_3_charge "pressure coming out of compressor";
  WorkingFluid.ThermodynamicState state_3_charge "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_charge "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_charge "compressor outlet spec. entropy";
  //STATE 3a
  WorkingFluid.ThermodynamicState state_3_a_charge "thermodynamic state of HEX outlet";
  SI.Temperature T_3_a_charge "outlet temperature after HEX";
  WorkingFluid.SpecificEnthalpy h_3_a_charge "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a_charge "HEX outlet spec. entropy";
  SI.Pressure p_3_a_charge "Pressure after HEX";
  //state 4a charge  guess
  SI.Temperature T_4_a_charge_guess = 291 "state  temperature";
  SI.Pressure p_4_a_charge_guess = 98942 "state  pressure";
  WorkingFluid.ThermodynamicState state_4_a_charge_guess "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_a_charge_guess "enthalpy";
  WorkingFluid.SpecificEntropy s_4_a_charge_guess "spec. entropy";
  //state 2 charge
  SI.Pressure p_2_charge "state  pressure";
  SI.Temperature T_2_charge(start = 309.3) "temperature ";
  WorkingFluid.ThermodynamicState state_2_charge "thermodynamic state ";
  WorkingFluid.SpecificEnthalpy h_2_charge "enthalpy";
  WorkingFluid.SpecificEntropy s_2_charge "spec. entropy";
  //state 1 charge
  SI.Pressure p_1_charge(start = 100000) " pressure ";
  SI.Temperature T_1_charge(start = from_degC(-66.388)) " temperature";
  WorkingFluid.ThermodynamicState state_1_charge "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_1_charge "spec enthalpy";
  WorkingFluid.SpecificEntropy s_1_charge " spec. entropy";
  //STATE 1 isentropic
  SI.Temperature T_1_is_charge "isentropic outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_1_is_charge "isentropic state of turbine outlet";
  WorkingFluid.SpecificEntropy s_1_is_charge "isentropic turbine outlet spec. entropy";
  //state 4a charge
  SI.Pressure p_4_a_charge "state  pressure";
  SI.Temperature T_4_a_charge "state  temperature";
  WorkingFluid.ThermodynamicState state_4_a_charge "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_a_charge "enthalpy";
  WorkingFluid.SpecificEntropy s_4_a_charge "spec. entropy";
  /*
    parameter SI.Temperature T_1=from_degC(-52.760) "ambient temperature";
    parameter SI.Pressure p_1=100000 "ambient pressure";
    WorkingFluid.BaseProperties state_1 "thermodynamic state of inlet";
    WorkingFluid.SpecificEnthalpy h_1 "inlet enthalpy";
    WorkingFluid.SpecificEntropy s_1 "inlet spec. entropy";
   
    //STATE 2 isentropic
    SI.Temperature T_2_is "isentropic outlet temperature of compressor";
    WorkingFluid.ThermodynamicState state_2_is "thermodynamic state of compressor outlet";
    WorkingFluid.SpecificEntropy s_2_is "compressor outlet spec. entropy";
    //STATE 2
    SI.Temperature T_2 "actual outlet temperature of compressor";
    SI.Pressure p_2 "pressure coming out of compressor";
    WorkingFluid.ThermodynamicState state_2 "thermodynamic state of compressor outlet";
    WorkingFluid.SpecificEnthalpy h_2 "compressor outlet enthalpy";
    WorkingFluid.SpecificEntropy s_2 "compressor outlet spec. entropy";  
    */
equation
//-------------TANK 1//
  solsalt_tank1_nom = HotTESLiquid.setState_pT(p_nom_tank1, T_nom_tank1);
  solsalt_tank1_state = HotTESLiquid.setState_pT(p_nom_tank1, T_tank1);
//solar salt properties
  solsalt_tank1.p = p_nom_tank1;
  solsalt_tank1.T = T_tank1;
//-------------TANK 2//
  solsalt_tank2_nom = HotTESLiquid.setState_pT(p_nom_tank2, T_nom_tank2);
  solsalt_tank2_state = HotTESLiquid.setState_pT(p_nom_tank2, T_tank2);
//solar salt properties
  solsalt_tank2.p = p_nom_tank2;
  solsalt_tank2.T = T_tank2;
//-------------TANK 3//
  coldliq_tank3_nom = ColdTESLiquid.setState_pT(p_nom_tank3, T_nom_tank3);
  coldliq_tank3_state = ColdTESLiquid.setState_pT(p_nom_tank3, T_tank3);
//-------------TANK 4//
  coldliq_tank4_nom = ColdTESLiquid.setState_pT(p_nom_tank4, T_nom_tank4);
  coldliq_tank4_state = ColdTESLiquid.setState_pT(p_nom_tank4, T_tank4);
//-------------COMPRESSOR CHARGE//
//P_12=compressor.P_12;
//reduced values compressor
  eta_is_CO_red_charge = eta_is_CO_charge/eta_is_CO_nom_charge;
  beta_CO_red_charge = beta_CO_charge/beta_CO_nom_charge;
  G_CO_red_charge = m_dot_WF_charge*sqrt(T_4_charge)/p_4_charge/(m_dot_WF_nom_charge*sqrt(T_4_nom_charge)/p_4_nom_charge);
  eta_is_CO_red_charge = (1 - c4*(1 - n_CO_red_charge)^2)*(n_CO_red_charge/G_CO_red_charge)*(2 - n_CO_red_charge/G_CO_red_charge);
  n_CO_red_charge = n_CO_charge/sqrt(T_4_charge)/(n_CO_nom_charge/sqrt(T_4_nom_charge));
//other compressor equations
  beta_CO_charge = p_3_charge/p_4_charge;
  eta_is_CO_charge = (T_3_is_charge - T_4_charge)/(T_3_charge - T_4_charge);
//P_12 = m_dot_WF*(h_2 - h_1);
  beta_CO_red_charge = c1_charge*G_CO_red_charge^2 + c2_charge*G_CO_red_charge + c3_charge;
  c1_charge = n_CO_red_charge/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c2_charge = (p - 2*m*n_CO_red_charge^2)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c3_charge = -1*(p*m*n_CO_red_charge - m^2*n_CO_red_charge^3)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
//-------------EXPANDER CHARGE//
//parameters
  alpha_charge = sqrt(1.4 - 0.4*n_TU_red_charge);
//reduced values turbine
  n_TU_red_charge = n_TU_charge/sqrt(T_2_charge)/(n_TU_nom_charge/sqrt(T_2_nom_charge));
  G_TU_red_charge = (m_dot_WF_charge*sqrt(T_2_charge)/p_2_charge)/(m_dot_WF_nom_charge*sqrt(T_2_nom_charge)/p_2_nom_charge);
  G_TU_red_charge = alpha_charge*sqrt(T_2_nom_charge/T_2_charge)*sqrt((beta_TU_charge^2 - 1)/(beta_TU_nom_charge^2 - 1));
  beta_TU_red_charge = beta_TU_charge/beta_TU_nom_charge;
  eta_is_TU_red_charge = (1 - t*(1 - n_TU_red_charge)^2)*(n_TU_red_charge/G_TU_red_charge)*(2 - ((n_TU_red_charge/G_TU_red_charge)));
  eta_is_TU_red_charge = eta_is_TU_charge/eta_is_TU_nom_charge;
//other turbine equations
  beta_TU_charge = p_2_charge/p_1_charge;
  eta_is_TU_charge = (T_2_charge - T_1_charge)/(T_2_charge - T_1_is_charge);
//-------------HEX 1 CHARGE//
//pressure loss
  p_3_charge = p_3_a_charge*(1 - f_p);
//off-design
  UA_HEX1_charge/UA_HEX1_nom = ((m_dot_WF_charge^0.8*m_dot_solsalt_HEX1_charge^0.8)/(m_dot_WF_nom_charge^0.8*m_dot_solsalt_HEX1_nom^0.8))*((m_dot_WF_nom_charge^0.8 + m_dot_solsalt_HEX1_nom^0.8)/(m_dot_WF_charge^0.8 + m_dot_solsalt_HEX1_charge^0.8));
//charge: inlet cold side= tank 2 inlet hot side=state 3
  outlet_hotside_guess_HEX1_charge = WorkingFluid.setState_pT(p_3_a_charge, T_tank2);
//WF
  outlet_coldside_guess_HEX1_charge = HotTESLiquid.setState_pT(p_nom_tank1, T_3_charge);
//NaK
  cp_cold_ave_HEX1_charge = (solsalt_tank2_state.cp + outlet_coldside_guess_HEX1_charge.cp)/2;
  cp_hot_ave_HEX1_charge = (WorkingFluid.specificHeatCapacityCp(state_3_charge) + WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX1_charge))/2;
  C_hot_HEX1_charge = m_dot_WF_charge*cp_hot_ave_HEX1_charge;
  C_cold_HEX1_charge = C_hot_HEX1_charge;
//balanced operation
  m_dot_solsalt_HEX1_charge = C_cold_HEX1_charge/cp_cold_ave_HEX1_charge;
  C_min_HEX1_charge = min(C_cold_HEX1_charge, C_hot_HEX1_charge);
  C_max_HEX1_charge = max(C_cold_HEX1_charge, C_hot_HEX1_charge);
  C_r_HEX1_charge = C_min_HEX1_charge/C_max_HEX1_charge;
  NTU_HEX1_charge = UA_HEX1_charge/C_min_HEX1_charge;
  if C_r_HEX1_charge == 1 then
    eff_HEX1_charge = NTU_HEX1_charge/(1 + NTU_HEX1_charge);
  else
    eff_HEX1_charge = (1 - exp(-NTU_HEX1_charge*(1 - C_r_HEX1_charge)))/(1 - C_r_HEX1_charge*exp(-NTU_HEX1_charge*(1 - C_r_HEX1_charge)));
  end if;
  Q_dot_HEX1_charge = eff_HEX1_charge*Q_dot_max_HEX1_charge;
  Q_dot_max_HEX1_charge = C_min_HEX1_charge*(T_tank2 - T_3_charge);
//calc temperature of fluid going into T_tank1 from cold side energy balance
  Q_dot_HEX1_charge = m_dot_solsalt_HEX1_charge*(solsalt_tank2.h - h_cold_out_HEX1_charge);
//calc T_3_a from hot side energy balance
  Q_dot_HEX1_charge = m_dot_WF_charge*(h_hot_out_HEX1_charge - state_3_charge.h);
//set outlets
  outlet_coldside_HEX1_charge = HotTESLiquid.setState_ph(p_nom_tank1, h_cold_out_HEX1_charge);
//NaK
  outlet_hotside_HEX1_charge = WorkingFluid.setState_ph(p_3_a_charge, h_hot_out_HEX1_charge);
//WF
  T_3_a_charge = outlet_hotside_HEX1_charge.T;
//-------------HEX2 (Recuperation)//
//off-design
  UA_HEX2_charge/UA_HEX2_nom = (m_dot_WF_charge^0.8)/(m_dot_WF_nom_charge^0.8);
//cold side
  cp_cold_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_3_a_charge);
//just inlet cp, not average cp
  C_cold_HEX2_charge = m_dot_WF_charge*cp_cold_HEX2_charge;
// hot side
  cp_hot_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_4_a_charge_guess);
//just inlet cp, not average cp
  C_hot_HEX2_charge = m_dot_WF_charge*cp_hot_HEX2_charge;
//variables for effectiveness
  C_min_HEX2_charge = min(C_cold_HEX2_charge, C_hot_HEX2_charge);
  C_max_HEX2_charge = max(C_cold_HEX2_charge, C_hot_HEX2_charge);
  C_r_HEX2_charge = C_min_HEX2_charge/C_max_HEX2_charge;
  NTU_HEX2_charge = UA_HEX2_charge/C_min_HEX2_charge;
  if C_r_HEX2_charge == 1 then
    eff_HEX2_charge = NTU_HEX2_charge/(1 + NTU_HEX2_charge);
  else
    eff_HEX2_charge = (1 - exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)))/(1 - C_r_HEX2_charge*exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)));
  end if;
//heat transferred
  Q_dot_HEX2_charge = eff_HEX2_charge*Q_dot_max_HEX2_charge;
  Q_dot_max_HEX2_charge = C_min_HEX2_charge*(T_3_a_charge - T_4_a_charge_guess);
//should acutally be T_4
//energy balance and outlet states
// hot side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(state_3_a_charge.h - h_hot_out_HEX2_charge);
//cold side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(h_cold_out_HEX2_charge - state_4_a_charge.h);
//instead of huess
//set outlets
  outlet_hotside_HEX2_charge = WorkingFluid.setState_ph(p_2_charge, h_hot_out_HEX2_charge);
//
  outlet_coldside_HEX2_charge = WorkingFluid.setState_ph(p_4_charge, h_cold_out_HEX2_charge);
  T_2_charge = outlet_hotside_HEX2_charge.T;
//T_4_charge= outlet_coldside_HEX2_charge.T;
//-------------HEX 3 CHARGE//
//pressure loss
  p_4_a_charge = p_1_charge*(1 - f_p_charge);
//off-design
  UA_HEX3_charge/UA_HEX3_nom = ((m_dot_WF_charge^0.8*m_dot_methanol_HEX3_charge^0.8)/(m_dot_WF_nom_charge^0.8*m_dot_methanol_HEX3_nom^0.8))*((m_dot_WF_nom_charge^0.8 + m_dot_methanol_HEX3_nom^0.8)/(m_dot_WF_charge^0.8 + m_dot_methanol_HEX3_charge^0.8));
//charge: inlet cold side= state 1 inlet hot side=tank 3
  outlet_hotside_guess_HEX3_charge = ColdTESLiquid.setState_pT(p_nom_tank4, T_tank4);
  outlet_coldside_guess_HEX3_charge = WorkingFluid.setState_pT(p_4_a_charge, T_tank3);
  cp_cold_ave_HEX3_charge = (WorkingFluid.specificHeatCapacityCp(state_1_charge) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX3_charge))/2;
  cp_hot_ave_HEX3_charge = (coldliq_tank4_state.cp + coldliq_tank3_state.cp)/2;
//heat capacity rates
  C_cold_HEX3_charge = m_dot_WF_charge*cp_cold_ave_HEX3_charge;
  C_hot_HEX3_charge = C_cold_HEX3_charge;
//balanced operation
  m_dot_methanol_HEX3_charge = C_hot_HEX3_charge/cp_hot_ave_HEX3_charge;
  C_min_HEX3_charge = min(C_cold_HEX3_charge, C_hot_HEX3_charge);
  C_max_HEX3_charge = max(C_cold_HEX3_charge, C_hot_HEX3_charge);
  C_r_HEX3_charge = C_min_HEX3_charge/C_max_HEX3_charge;
  NTU_HEX3_charge = UA_HEX3_charge/C_min_HEX3_charge;
  if C_r_HEX3_charge == 1 then
    eff_HEX3_charge = NTU_HEX3_charge/(1 + NTU_HEX3_charge);
  else
    eff_HEX3_charge = (1 - exp(-NTU_HEX3_charge*(1 - C_r_HEX3_charge)))/(1 - C_r_HEX3_charge*exp(-NTU_HEX3_charge*(1 - C_r_HEX3_charge)));
  end if;
  Q_dot_HEX3_charge = eff_HEX3_charge*Q_dot_max_HEX3_charge;
  Q_dot_max_HEX3_charge = C_min_HEX3_charge*(T_tank3 - T_1_charge);
//calc temperature of fluid going into T_tank4 from hot side energy balance
  Q_dot_HEX3_charge = m_dot_methanol_HEX3_charge*(coldliq_tank3_state.h - h_hot_out_HEX3_charge);
//calc T_4_a_charge from cold side energy balance
  Q_dot_HEX3_charge = m_dot_WF_charge*(h_cold_out_HEX3_charge - state_1_charge.h);
//set outlets
  outlet_hotside_HEX3_charge = ColdTESLiquid.setState_ph(p_nom_tank4, h_hot_out_HEX3_charge);
  outlet_coldside_HEX3_charge = WorkingFluid.setState_ph(p_4_a_charge, h_cold_out_HEX3_charge);
  T_4_a_charge = outlet_coldside_HEX3_charge.T;
//-------------STATES//
//STATE 4 charge
  state_4_charge = WorkingFluid.setState_pT(p_4_charge, T_4_charge);
  h_4_charge = WorkingFluid.specificEnthalpy(state_4_charge);
  s_4_charge = WorkingFluid.specificEntropy(state_4_charge);
//STATE 3 isentropic charge
  s_3_is_charge = s_4_charge "isentropic compressor outlet spec. entropy";
  state_3_is_charge = WorkingFluid.setState_ps(p_3_charge, s_3_is_charge) "isentropic state of compressor outlet";
  T_3_is_charge = state_3_is_charge.T;
//STATE 3 charge
  state_3_charge = WorkingFluid.setState_pT(p_3_charge, T_3_charge);
  h_3_charge = WorkingFluid.specificEnthalpy(state_3_charge);
  s_3_charge = WorkingFluid.specificEntropy(state_3_charge);
//STATE 3a charge
  state_3_a_charge = WorkingFluid.setState_pT(p_3_a_charge, T_3_a_charge);
  h_3_a_charge = WorkingFluid.specificEnthalpy(state_3_a_charge);
  s_3_a_charge = WorkingFluid.specificEntropy(state_3_a_charge);
//state 4a charge  guess
  state_4_a_charge_guess = WorkingFluid.setState_pT(p_4_a_charge_guess, T_4_a_charge_guess);
  h_4_a_charge_guess = WorkingFluid.specificEnthalpy(state_4_a_charge_guess);
  s_4_a_charge_guess = WorkingFluid.specificEntropy(state_4_a_charge_guess);
//state 2 charge
  p_2_charge = p_3_a_charge*(1 - f_p_charge);
  state_2_charge = WorkingFluid.setState_pT(p_2_charge, T_2_charge);
  h_2_charge = WorkingFluid.specificEnthalpy(state_2_charge);
  s_2_charge = WorkingFluid.specificEntropy(state_2_charge);
//state 1 charge
  state_1_charge = WorkingFluid.setState_pT(p_1_charge, T_1_charge);
  h_1_charge = WorkingFluid.specificEnthalpy(state_1_charge);
  s_1_charge = WorkingFluid.specificEntropy(state_1_charge);
//state 1  isentropic charge
  s_1_is_charge = s_2_charge;
  state_1_is_charge = WorkingFluid.setState_ps(p_1_charge, s_1_is_charge) "isentropic state of turbine outlet";
  T_1_is_charge = state_1_is_charge.T "isentropic turbine outlet spec. entropy";
//STATE 4_a charge
  state_4_a_charge = WorkingFluid.setState_pT(p_4_a_charge, T_4_a_charge);
  h_4_a_charge = WorkingFluid.specificEnthalpy(state_4_a_charge);
  s_4_a_charge = WorkingFluid.specificEntropy(state_4_a_charge);
  annotation(
    Documentation(info = "<html><head></head><body>Static Malta charge cycle<div>Everything nominal, no adjustments yet</div><div>T_4 is set&nbsp;</div><div>and cycle is open (so T_4 result is not calculated)</div></body></html>"));
end StaticMalta_open_nom;
