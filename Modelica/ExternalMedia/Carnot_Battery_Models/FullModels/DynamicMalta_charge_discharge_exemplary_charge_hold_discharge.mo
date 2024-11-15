within ExternalMedia.Carnot_Battery_Models.FullModels;

model DynamicMalta_charge_discharge_exemplary_charge_hold_discharge
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.NonSI; 
  import Modelica.Units.Conversions.from_degC;
  import Modelica.ComplexMath;
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

  package EthyleneGlycol "Ethylene Glycol properties from CoolProp"
    extends ExternalMedia.Media.IncompressibleCoolPropMedium(mediumName = "ZM", substanceNames = {"ZM[0.6]"});
  end EthyleneGlycol;

  replaceable package RejectionHeatTransferFluid = EthyleneGlycol constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
  /*
                package WorkingFluid = ExternalMedia.Media.CoolPropMedium(mediumName = "Nitrogen", substanceNames = {"N2"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph); 
                */
  //--------------------------INPUTS
  //input Integer Mode(start = 1);
  //parameter Integer Mode = 0;
    Integer Mode(start = 1);
    
    parameter Real SOC_tank1_start = 0;
    parameter SI.Temperature T_tank1_start = from_degC(565);
    parameter Real SOC_tank2_start = 1;
    parameter SI.Temperature T_tank2_start = from_degC(279);
    parameter Real SOC_tank3_start = 1;
    parameter SI.Temperature T_tank3_start = from_degC(25.1);
    parameter Real SOC_tank4_start = 0;
    parameter SI.Temperature T_tank4_start = from_degC(-59.75);
  
  /*
  parameter Real SOC_tank1_start = 1;
  parameter SI.Temperature T_tank1_start = from_degC(565);
  parameter Real SOC_tank2_start = 0;
  parameter SI.Temperature T_tank2_start = from_degC(279);
  parameter Real SOC_tank3_start = 0;
  parameter SI.Temperature T_tank3_start = from_degC(25.1);
  parameter Real SOC_tank4_start = 1;
  parameter SI.Temperature T_tank4_start = from_degC(-59.75);
  */
  //--------------------------PARAMETERS & VARIABLES SYSTEM-----------------------------//
  parameter SI.Temperature T0 = T_amb;
  parameter SI.Temperature T_amb = from_degC(25);
  WorkingFluid.ThermodynamicState state_amb_air(p(start = 101315), T(start = from_degC(25))) "thermodynamic state of rejec outlet";
  //-------------Tanks//
  SI.Energy exergy_total_tanks(displayUnit = "MWh");
  SI.Energy exergy_hot_tanks(displayUnit = "MWh");
  SI.Energy exergy_cold_tanks(displayUnit = "MWh");
  SI.Energy int_energy_total_tanks(displayUnit = "MWh");
  //-------------Charge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom_charge = 766 "design mass flow rate";
  //actual
  SI.MassFlowRate m_dot_WF_charge(start = 766) "mass flow rate";
  parameter Real f_p_charge = 0.01075 "pressure_loss_factor percent";
  //not needed anymore
  parameter Real k_p_charge = 0.0041715 "pressure_loss_factor";
  parameter SI.Pressure p_fix_charge = 100000 "fixed pressure point through expansion vessel at p_4";
  SI.Energy Elec_energy_charge(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_elec_charge(displayUnit = "MW");
  Real m_dot_div_p_charge(start = 0.0076);
   SI.Energy exergy_total_loss_irr_charge(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_total_loss_irr_charge(displayUnit = "MW");
  SI.Energy E_mech_shaft_charge(displayUnit = "MWh", start = 0, fixed = true);
  parameter Real hot_to_cold_mass_flow_ratio_charge = 2.044;
    Real COP_system_charge;
  SI.Energy E_total_loss_irr_charge(displayUnit = "MWh", start = 0, fixed = true);
  //-------------Discharge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  //actual
  SI.MassFlowRate m_dot_WF(start = 762) "mass flow rate";
  parameter Real f_p = 0.01625 "pressure_loss_factor percent";
  parameter Real k_p = 0.0062084 "pressure_loss_factor";
  parameter SI.Pressure p_fix = 100000 "fixed pressure point through expansion vessel";
  SI.Energy Elec_energy_discharge(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_elec(displayUnit = "MW");
  //other
  SI.HeatFlowRate Q_dot_hightemp_res(displayUnit = "MW");
  Real m_dot_div_p(start = 0.0076);
   SI.Energy exergy_total_loss_irr(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_total_loss_irr(displayUnit = "MW");
  SI.Energy E_mech_shaft(displayUnit = "MWh", start = 0, fixed = true);
  parameter Real hot_to_cold_mass_flow_ratio = 2.044;
      Real COP_system;
  SI.Energy E_total_loss_irr(displayUnit = "MWh", start = 0, fixed = true);  
  //--------------------------PARAMETERS & VARIABLES TANKS-----------------------------//
  parameter SI.Mass m_working_solar_salt = 19386000;
  parameter SI.Mass m_working_methanol = 9486000;
  //geometry for both solar salt tanks
  parameter SI.Diameter D_tank_solsalt = 37.395499;
  parameter SI.Height h_tank_solsalt = 11.218649;
  parameter SI.Area A_cross_tank_solsalt = pi*(D_tank_solsalt/2)^2 "Cross-section Tank";
  parameter SI.Volume V_tank_solsalt = A_cross_tank_solsalt*h_tank_solsalt "Volume Tank";
  parameter SI.Thickness d_insulation_tank_solsalt = 0.4 "thickness of insulation wall layer";
  parameter SI.Radius r_tank_solsalt = D_tank_solsalt/2 "radius at start of insulation";
  parameter SI.Radius r_outer = D_tank_solsalt/2 + d_insulation_tank_solsalt "outer radius";
  parameter SI.Area A_W_tank_solsalt = (pi*D_tank_solsalt*h_tank_solsalt) + A_cross_tank_solsalt*2 "Total Surface Wall";
  parameter SI.Height x_tank_solsalt_min = 0.4;
  //geometry for both methanol tanks
  parameter SI.Diameter D_tank_coldliq = 36.773594;
  parameter SI.Height h_tank_coldliq = 11.032078;
  parameter SI.Area A_cross_tank_coldliq = pi*(D_tank_coldliq/2)^2 "Cross-section Tank";
  parameter SI.Volume V_tank_coldliq = A_cross_tank_coldliq*h_tank_coldliq "Volume Tank";
  parameter SI.Thickness d_insulation_tank_coldliq = 0.4 "thickness of insulation wall layer";
  parameter SI.Radius r_tank_coldliq = D_tank_coldliq/2 "radius at start of insulation";
  parameter SI.Radius r_outer_coldliq = D_tank_coldliq/2 + d_insulation_tank_coldliq "outer radius";
  parameter SI.Area A_W_tank_coldliq = (pi*D_tank_coldliq*h_tank_coldliq) + A_cross_tank_coldliq*2 "Total Surface Wall";
  parameter SI.Height x_tank_coldliq_min = 0.4;
  //--------------------------TANK 1
  //nominal
  parameter SI.Temperature T_tank1_nom = from_degC(565.6);
  parameter SI.Pressure p_tank1_nom = 101325 "unpressurized tank";
  parameter SI.Density rho_tank1_nom = 1730.66;
  parameter SI.Mass m_tank1_min = x_tank_solsalt_min*A_cross_tank_solsalt*rho_tank1_nom;
  HotTESLiquid.ThermodynamicState solsalt_tank1_nom(p(start = p_tank1_nom), T(start = T_tank1_nom)) "Nominal thermodynamic state";
  //start
  parameter SI.Mass m_tank1_start = m_tank1_min + (SOC_tank1_start*m_working_solar_salt);
  //state variables
  SI.Mass m_tank1;
  SI.Temperature T_tank1;
  Real SOC_tank1;
  SI.PathLength x_tank1 "salt-level";
  SI.Energy int_energy_tank1(displayUnit = "MWh");
  SI.Energy exergy_tank1(displayUnit = "MWh");
  //thermodynamic states
  HotTESLiquid.ThermodynamicState solsalt_tank1_state(p(start = p_tank1_nom), T(start = T_tank1_nom)) "Nominal thermodynamic state";
  HotTESLiquid.BaseProperties solsalt_tank1(p(start = p_tank1_nom), T(start = T_tank1_nom)) "Medium properties of port_a";
  //losses
  SI.HeatFlowRate Q_dot_to_amb_tank1(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank1;
  //flows
  SI.MassFlowRate m_out_tank1;
  SI.MassFlowRate m_in_tank1;
  //ports
  SI.SpecificEnthalpy h_in_tank1;
  //--------------------------TANK 2
  //nominal
  parameter SI.Temperature T_tank2_nom = from_degC(271.5);
  parameter SI.Pressure p_tank2_nom = 101325 "unpressurized tank";
  parameter SI.Density rho_tank2_nom = 1917.33;
  parameter SI.Mass m_tank2_min = x_tank_solsalt_min*A_cross_tank_solsalt*rho_tank2_nom;
  HotTESLiquid.ThermodynamicState solsalt_tank2_nom(p(start = p_tank2_nom), T(start = T_tank2_nom)) "Nominal thermodynamic state";
  //start
  parameter SI.Mass m_tank2_start = m_tank2_min + (SOC_tank2_start*m_working_solar_salt);
  //state variables
  SI.Mass m_tank2;
  SI.Temperature T_tank2;
  Real SOC_tank2;
  SI.PathLength x_tank2 "salt-level";
  SI.Energy int_energy_tank2(displayUnit = "MWh");
  SI.Energy exergy_tank2(displayUnit = "MWh");
  //thermodynamic states
  HotTESLiquid.ThermodynamicState solsalt_tank2_state(p(start = p_tank2_nom), T(start = T_tank2_nom)) "Nominal thermodynamic state";
  HotTESLiquid.BaseProperties solsalt_tank2(p(start = p_tank2_nom), T(start = T_tank2_nom)) "Medium properties of port_a";
  //losses
  SI.HeatFlowRate Q_dot_to_amb_tank2(displayUnit = "kW") "Heat Flow to ambient";
  Real Q_div_A_tank2;
  //flows
  SI.MassFlowRate m_out_tank2;
  SI.MassFlowRate m_in_tank2;
  //ports
  SI.SpecificEnthalpy h_in_tank2;
  //--------------------------TANK 3
  //nominal
  parameter SI.Temperature T_tank3_nom = from_degC(25.195);
  parameter SI.Pressure p_tank3_nom = 101325 "unpressurized tank";
  parameter SI.Density rho_tank3_nom = 890.546;
  parameter SI.Mass m_tank3_min = x_tank_coldliq_min*A_cross_tank_coldliq*rho_tank3_nom;
  ColdTESLiquid.ThermodynamicState coldliq_tank3_nom(p(start = p_tank3_nom), T(start = T_tank3_nom)) "Nominal thermodynamic state";
  //start
  parameter SI.Mass m_tank3_start = m_tank3_min + (SOC_tank3_start*m_working_methanol);
  //state variables
  SI.Mass m_tank3;
  SI.Temperature T_tank3;
  Real SOC_tank3;
  SI.PathLength x_tank3 "salt-level";
  SI.Energy int_energy_tank3(displayUnit = "MWh");
  SI.Energy exergy_tank3(displayUnit = "MWh");
  //thermodynamic states
  ColdTESLiquid.ThermodynamicState coldliq_tank3_state(p(start = p_tank3_nom), T(start = T_tank3_nom)) "Nominal thermodynamic state";
  ColdTESLiquid.BaseProperties coldliq_tank3(p(start = p_tank3_nom), T(start = T_tank3_nom)) "Medium properties of port_a";
  //flows
  SI.MassFlowRate m_out_tank3;
  SI.MassFlowRate m_in_tank3;
  //ports
  SI.SpecificEnthalpy h_in_tank3;
  //--------------------------TANK 4
  //nominal
  parameter SI.Temperature T_tank4_nom = from_degC(-59.85);
  parameter SI.Pressure p_tank4_nom = 101325 "unpressurized tank";
  parameter SI.Density rho_tank4_nom = 949.56;
  parameter SI.Mass m_tank4_min = x_tank_coldliq_min*A_cross_tank_coldliq*rho_tank4_nom;
  ColdTESLiquid.ThermodynamicState coldliq_tank4_nom(p(start = p_tank4_nom), T(start = T_tank4_nom)) "Nominal thermodynamic state";
  //start
  parameter SI.Mass m_tank4_start = m_tank4_min + (SOC_tank4_start*m_working_methanol);
  //state variables
  SI.Mass m_tank4;
  SI.Temperature T_tank4;
  Real SOC_tank4;
  SI.PathLength x_tank4 "salt-level";
  SI.Energy int_energy_tank4(displayUnit = "MWh");
  SI.Energy exergy_tank4(displayUnit = "MWh");
  //thermodynamic states
  ColdTESLiquid.ThermodynamicState coldliq_tank4_state(p(start = p_tank4_nom), T(start = T_tank4_nom)) "Nominal thermodynamic state";
  ColdTESLiquid.BaseProperties coldliq_tank4(p(start = p_tank4_nom), T(start = T_tank4_nom)) "Medium properties of port_a";
  //flows
  SI.MassFlowRate m_out_tank4;
  SI.MassFlowRate m_in_tank4;
  //ports
  SI.SpecificEnthalpy h_in_tank4;
  //--------------------------PARAMETERS & VARIABLES CHARGE-----------------------------//
  //-------------HEX 1 CHARGE//
  //design
  parameter Real UA_HEX1_nom = 25906701;
  SI.MassFlowRate m_dot_solsalt_HEX1_nom = 540;
  //actual
  SI.Pressure delta_P_HEX1_charge;
  Real UA_HEX1_charge;
  SI.MassFlowRate m_dot_solsalt_HEX1_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX1_charge(p(start = 456931), T(start = T_tank2_nom), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
  HotTESLiquid.ThermodynamicState outlet_coldside_guess_HEX1_charge(p(start = 101325), T(start = T_tank1_nom), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
  SI.SpecificHeatCapacity cp_hot_ave_HEX1_charge;
  SI.SpecificHeatCapacity cp_cold_ave_HEX1_charge;
  //variables for effectiveness
  Real C_cold_HEX1_charge "Heat capacity rate of cold side of hot HEX";
  Real C_hot_HEX1_charge "Heat capacity rate of hot side of hot HEX";
  Real C_min_HEX1_charge;
  Real C_max_HEX1_charge;
  Real C_r_HEX1_charge;
  Real NTU_HEX1_charge "Number of transfer units of hot HEX";
  Real eff_HEX1_charge(start = 0.97) "effectiveness of hot HEX";
  SI.HeatFlowRate Q_dot_HEX1_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX1_charge(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX1_charge;
  SI.SpecificEnthalpy h_cold_out_HEX1_charge;
  WorkingFluid.ThermodynamicState outlet_hotside_HEX1_charge(p(start = 456931), T(start = T_tank2_nom), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
  HotTESLiquid.ThermodynamicState outlet_coldside_HEX1_charge(p(start = 101325), T(start = T_tank1_nom), phase(start = 1)) "Medium properties of HEX port at interface to tank 1";
  SI.Power P_loss_irr_HEX1_charge(displayUnit = "MW");
   SI.Energy E_loss_irr_HEX1_charge(displayUnit = "MWh", start = 0, fixed = true);
  //-------------HEX 2 CHARGE//
  //design
  parameter Real UA_HEX2_nom = 13279998;
  //actual
  SI.Pressure delta_P_HEX2_charge;
  Real UA_HEX2_charge;
  //outlet guess states
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
  Real eff_HEX2_charge(start = 0.96) "effectiveness of hot HEX";
  SI.Power P_loss_irr_HEX2_charge(displayUnit = "MW");
   SI.Energy E_loss_irr_HEX2_charge(displayUnit = "MWh", start = 0, fixed = true);
  //heat transferred
  SI.HeatFlowRate Q_dot_HEX2_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX2_charge(displayUnit = "MW");
  //energy balance and outlet states
  SI.SpecificEnthalpy h_hot_out_HEX2_charge;
  SI.SpecificEnthalpy h_cold_out_HEX2_charge;
  WorkingFluid.ThermodynamicState outlet_hotside_HEX2_charge(p(start = 452019), T(start = 309), phase(start = 1)) "Medium properties";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX2_charge(p(start = 97889), T(start = 540), phase(start = 1)) "Medium properties";
  //-------------HEX 3 CHARGE//
  //design
  parameter Real UA_HEX3_nom = 9505824;
  parameter SI.MassFlowRate m_dot_methanol_HEX3_nom = 265;
  //actual
  SI.Pressure delta_P_HEX3_charge;
  Real UA_HEX3_charge;
  SI.MassFlowRate m_dot_methanol_HEX3_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  ColdTESLiquid.ThermodynamicState outlet_hotside_guess_HEX3_charge(p(start = 101325), T(start = T_tank4_nom), phase(start = 1)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX3_charge(p(start = 100000), T(start = T_tank3_nom), phase(start = 1)) "Medium properties ";
  SI.SpecificHeatCapacity cp_hot_ave_HEX3_charge;
  SI.SpecificHeatCapacity cp_cold_ave_HEX3_charge;
  //variables for effectiveness
  Real C_cold_HEX3_charge "Heat capacity rate of cold side of hot HEX";
  Real C_hot_HEX3_charge "Heat capacity rate of hot side of hot HEX";
  Real C_min_HEX3_charge;
  Real C_max_HEX3_charge;
  Real C_r_HEX3_charge;
  Real NTU_HEX3_charge "Number of transfer units of hot HEX";
  Real eff_HEX3_charge(start = 0.92) "effectiveness of hot HEX";
  SI.Power P_loss_irr_HEX3_charge(displayUnit = "MW");
   SI.Energy E_loss_irr_HEX3_charge(displayUnit = "MWh", start = 0, fixed = true);
   //heat transferred
  SI.HeatFlowRate Q_dot_HEX3_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX3_charge(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX3_charge;
  SI.SpecificEnthalpy h_cold_out_HEX3_charge;
  ColdTESLiquid.ThermodynamicState outlet_hotside_HEX3_charge(p(start = 101325), T(start = T_tank4_nom), phase(start = 1)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX3_charge(p(start = 100000), T(start = T_tank3_nom), phase(start = 1)) "Medium properties ";
  //-------------HEX rejection charge//
  SI.HeatFlowRate Q_dot_HEXrej_charge(displayUnit = "MW");
   SI.Heat Q_HEXrej_charge(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_loss_irr_HEXrej_charge(displayUnit = "MW");
    SI.Energy E_loss_irr_HEXrej_charge(displayUnit = "MWh", start = 0, fixed = true);
  //-------------COMPRESSOR CHARGE//
  parameter Real p = 0.42 "compressor map factor";
  parameter Real m = 1.06 "compressor map factor";
  parameter Real c4 = 0.3 "factor 4, see Zhang2002";
  Real c1_charge "factor 1";
  Real c2_charge "factor 2";
  Real c3_charge "factor 3";
  //design
  parameter Real beta_CO_nom_charge = 4.592 "design compression ratio";
  parameter Real n_CO_nom_charge = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom_charge = 0.90385  "design isentropic efficiency";
  SI.Temperature T_4_nom_charge = from_degC(267.533) "state 4 temperature";
  SI.Pressure p_4_nom_charge = 100000 "state 4 pressure";
  //actual
  Real beta_CO_charge(start = beta_CO_nom_charge) "absolute compression ratio";
  parameter Real n_CO_charge = 3000 "actual speed";
  SI.Efficiency eta_is_CO_charge "absolute isentropic efficiency";
  //reduced
  Real beta_CO_red_charge(start = 1) "reduced compression ratio";
  Real n_CO_red_charge(start = 1) "reduced speed";
  SI.Efficiency eta_is_CO_red_charge "reduced isentropic efficiency";
  Real G_CO_red_charge(start = 1) "reduced mass flow rate compressor";
  //other
  SI.Power P_mech_CO_charge(displayUnit = "MW");
  SI.Power P_loss_irr_CO_charge(displayUnit = "MW");
    SI.Energy E_loss_irr_CO_charge(displayUnit = "MWh", start = 0, fixed = true);
  //limits
  Real beta_CO_red_charge_min;
  Real beta_CO_red_charge_max;
  //-------------EXPANDER CHARGE//
  parameter Real t = 0.3 "parameter, see Zhang2002";
  Real alpha_charge "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom_charge = 4.331 "design expansion ratio";
  parameter Real n_TU_nom_charge = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom_charge = 0.92 "design isentropic efficiency";
  parameter SI.Temperature T_2_a_nom_charge = from_degC(25) "nominal inlet temperature";
  parameter SI.Pressure p_2_a_nom_charge = 442461 "state pressure";
  //actual
  Real n_TU_charge(start = 3000) "actual speed";
  Real beta_TU_charge(start = beta_TU_nom_charge) "absolute expansion ratio";
  SI.Efficiency eta_is_TU_charge "absolute isentropic turbine efficiency";
  //reduced
  Real n_TU_red_charge(start = 1) "reduced speed";
  Real G_TU_red_charge(start = 1) "reduced mass flow rate turbine";
  Real beta_TU_red_charge(start = 1) "reduced expansion ratio";
  SI.Efficiency eta_is_TU_red_charge "reduced isentropic turbine efficiency";
  //other
  SI.Power P_mech_TU_charge(displayUnit = "MW");
  SI.Power P_loss_irr_TU_charge(displayUnit = "MW");
   SI.Energy E_loss_irr_TU_charge(displayUnit = "MWh", start = 0, fixed = true);
  //limits
  parameter Real beta_TU_red_charge_min = 0.4;
  parameter Real beta_TU_red_charge_max = 1.4;
  //-------------SYSTEM CHARGE//
  SI.Power P_mech_shaft_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_pump_charge(displayUnit = "MW");
  Real COP(start = 1);
  Real work_ratio_charge;
  //-------------STATES//
  //-------------Charge//
  //fixed temperature point charge
  parameter SI.Temperature T_2_a_charge = from_degC(25) "outlet temperature after rejec";
  //fixed pressure point charge
  SI.Pressure p_4_charge = p_fix_charge "state 4 pressure";
  //state 1 charge
  SI.Pressure p_1_charge(start = 104896) " pressure ";
  SI.Temperature T_1_charge(start = from_degC(-66)) " temperature";
  WorkingFluid.ThermodynamicState state_1_charge(p(start = 102161), T(start = from_degC(-66))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_1_charge(start = 336932) "spec enthalpy";
  WorkingFluid.SpecificEntropy s_1_charge(start = 3510) " spec. entropy";
  //STATE 1 isentropic charge
  SI.Temperature T_1_is_charge(start = from_degC(-69.4)) "isentropic outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_1_is_charge(p(start = 102161), T(start = from_degC(-69.4))) "isentropic state of turbine outlet";
  WorkingFluid.SpecificEntropy s_1_is_charge(start = 3474) "isentropic turbine outlet spec. entropy";
  WorkingFluid.SpecificEnthalpy h_1_is_charge(start = 318577) "spec enthalpy";
  //state 4a charge
  SI.Pressure p_4_a_charge(start = 102448) "state  pressure";
  SI.Temperature T_4_a_charge(start = from_degC(18.14)) "state  temperature";
  WorkingFluid.ThermodynamicState state_4_a_charge(p(start = 100000), T(start = from_degC(18.14))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_a_charge(start = 417519) "enthalpy";
  WorkingFluid.SpecificEntropy s_4_a_charge(start = 3837) "spec. entropy";
  //state 4 a charge guess
  SI.Temperature T_4_a_charge_guess(start = from_degC(18.14)) "state 4 temperature";
  WorkingFluid.ThermodynamicState state_4_a_charge_guess(p(start = 456931), T(start = from_degC(18.1))) "thermodynamic state";
  //state 4 charge guess
  SI.Temperature T_4_charge_guess(start = from_degC(267.533)) "state 4 temperature";
  //WorkingFluid.ThermodynamicState state_4_charge_guess(p(start = p_fix_charge), T(start = from_degC(267.533))) "thermodynamic state";
  WorkingFluid.ThermodynamicState state_4_charge_guess(p(start = 100000), T(start = from_degC(267.533))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_charge_guess(start = 671458) "enthalpy";
  WorkingFluid.SpecificEntropy s_4_charge_guess(start = 4495) "spec. entropy";
  //state 4 charge
  SI.Temperature T_4_charge(start = from_degC(267.533)) "state 4 temperature";
  // WorkingFluid.ThermodynamicState state_4_charge(p(start = p_fix_charge), T(start = from_degC(267.533))) "thermodynamic state";
  WorkingFluid.ThermodynamicState state_4_charge(p(start = 100000), T(start = from_degC(267.533))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_charge(start = 671458) "enthalpy";
  WorkingFluid.SpecificEntropy s_4_charge(start = 4495) "spec. entropy";
  //STATE 3 isentropic charge
  SI.Temperature T_3_is_charge(start = from_degC(541.61)) "isentropic outlet temperature of compressor";
  WorkingFluid.ThermodynamicState state_3_is_charge(p(start = 449506), T(start = from_degC(541.61))) "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEntropy s_3_is_charge(start = 4495) "compressor outlet spec. entropy";
  WorkingFluid.SpecificEnthalpy h_3_is_charge(start = 964771);
  //STATE 3 charge
  SI.Temperature T_3_charge(start = from_degC(578.98)) "actual outlet temperature of compressor";
  SI.Pressure p_3_charge(start = 459203) "pressure coming out of compressor";
  WorkingFluid.ThermodynamicState state_3_charge(p(start = 449506), T(start = from_degC(578.98))) "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_charge(start = 1006139) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_charge(start = 4545) "compressor outlet spec. entropy";
  //STATE 3a charge
  WorkingFluid.ThermodynamicState state_3_a_charge(p(start = 456931), T(start = from_degC(281.29))) "thermodynamic state";
  SI.Temperature T_3_a_charge(start = from_degC(281.29)) "outlet temperature after HEX";
  SI.Pressure p_3_a_charge(start = 456754) "Pressure after HEX";
  WorkingFluid.SpecificEnthalpy h_3_a_charge(start = 685666) "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a_charge(start = 4079) "HEX outlet spec. entropy";
  //state 2 charge
  SI.Pressure p_2_charge(start = 454306) "state  pressure";
  SI.Temperature T_2_charge(start = 309.3) "temperature ";
  WorkingFluid.ThermodynamicState state_2_charge(p(start = 452019), T(start = from_degC(32.66))) "thermodynamic state ";
  WorkingFluid.SpecificEnthalpy h_2_charge(start = 431395) "enthalpy";
  WorkingFluid.SpecificEntropy s_2_charge(start = 3474) "spec. entropy";
  //STATE 2a charge
  WorkingFluid.ThermodynamicState state_2_a_charge(p(start = 442461), T(start = from_degC(25))) "thermodynamic state of rejec outlet";
  WorkingFluid.SpecificEnthalpy h_2_a_charge(start = 423634) "rejec outlet enthalpy";
  WorkingFluid.SpecificEntropy s_2_a_charge(start = 3448) "rejec outlet spec. entropy";
  SI.Pressure p_2_a_charge(start = 454306) "Pressure after rejec";
  //--------------------------PARAMETERS & VARIABLES DISCHARGE-----------------------------//
  //-------------HEX 1//
  //actual
  Real UA_HEX1;
  SI.Pressure delta_P_HEX1;
  //outlet guess states
  HotTESLiquid.ThermodynamicState outlet_hotside_guess_HEX1(p(start = 101315), T(start = 271.5)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX1(p(start = 572045), T(start = 556.099)) "Medium properties ";
  SI.SpecificHeatCapacity cp_hot_ave_HEX1;
  SI.SpecificHeatCapacity cp_cold_ave_HEX1;
  //variables for effectiveness
  Real C_cold_HEX1 "Heat capacity rate of cold side of hot HEX";
  Real C_hot_HEX1 "Heat capacity rate of hot side of hot HEX";
  SI.MassFlowRate m_dot_solsalt_HEX1 "mass flow rate required for balanced HEX";
  Real C_min_HEX1;
  Real C_max_HEX1;
  Real C_r_HEX1;
  Real NTU_HEX1 "Number of transfer units of hot HEX";
  Real eff_HEX1 "effectiveness of hot HEX";
  SI.HeatFlowRate Q_dot_HEX1(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX1(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX1;
  SI.SpecificEnthalpy h_cold_out_HEX1;
  HotTESLiquid.ThermodynamicState outlet_hotside_HEX1(T(start = 270)) "Medium properties of HEX port at interface to tank 2";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX1(T(start = 570)) "Medium properties of HEX port at interface to tank 2";
  SI.Power P_loss_irr_HEX1(displayUnit = "MW");
   SI.Energy E_loss_irr_HEX1(displayUnit = "MWh", start = 0, fixed = true);
  //-------------HEX 2 DISCHARGE//
  //actual
  SI.Pressure delta_P_HEX2;
  Real UA_HEX2;
  //outlet guess states
  //cold-side
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX2(p(start = 581494), T(start = 269.477)) "Medium properties of HEX port at interface to tank 2";
  SI.SpecificHeatCapacity cp_cold_ave_HEX2;
  //hot-side
  WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX2(p(start = 105039), T(start = 121.643)) "Medium properties of HEX port at interface to tank 2";
  SI.SpecificHeatCapacity cp_hot_ave_HEX2;
  //variables for effectiveness
  Real C_cold_HEX2 "Heat capacity rate of cold side of recuperation HEX (after compressor)";
  Real C_hot_HEX2 "Heat capacity rate of cold side of recuperation HEX (after compressor)";
  Real C_min_HEX2;
  Real C_max_HEX2;
  Real C_r_HEX2;
  Real NTU_HEX2 "Number of transfer units of recuperation HEX";
  Real eff_HEX2 "effectiveness of hot HEX";
  //heat transferred
  SI.HeatFlowRate Q_dot_HEX2(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX2(displayUnit = "MW");
  //energy balance and outlet states
  SI.SpecificEnthalpy h_hot_out_HEX2;
  SI.SpecificEnthalpy h_cold_out_HEX2;
  WorkingFluid.ThermodynamicState outlet_hotside_HEX2(p(start = 581494), T(start = 269.477)) "Medium properties of HEX port at interface to tank 2";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX2(p(start = 105039), T(start = 121.643)) "Medium properties of HEX port at interface to tank 2";
  SI.Power P_loss_irr_HEX2(displayUnit = "MW");
  SI.Energy E_loss_irr_HEX2(displayUnit = "MWh", start = 0, fixed = true);
  //-------------HEX 3 DISCHARGE//
  //actual
  Real UA_HEX3;
  SI.Pressure delta_P_HEX3;
  SI.MassFlowRate m_dot_methanol_HEX3 "mass flow rate required for balanced HEX";
  //outlet guess states
  ColdTESLiquid.ThermodynamicState outlet_coldside_guess_HEX3(p(start = 101325), T(start = T_tank3_nom)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX3(p(start = 100000), T(start = T_tank4_nom)) "Medium properties ";
  SI.SpecificHeatCapacity cp_hot_ave_HEX3;
  SI.SpecificHeatCapacity cp_cold_ave_HEX3;
  //variables for effectiveness
  Real C_cold_HEX3 "Heat capacity rate of cold side ";
  Real C_hot_HEX3 "Heat capacity rate of hot side";
  Real C_min_HEX3;
  Real C_max_HEX3;
  Real C_r_HEX3;
  Real NTU_HEX3 "Number of transfer units of hot HEX";
  Real eff_HEX3(start = 0.92) "effectiveness of hot HEX";
  //heat transferred
  SI.HeatFlowRate Q_dot_HEX3(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX3(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX3;
  SI.SpecificEnthalpy h_cold_out_HEX3;
  ColdTESLiquid.ThermodynamicState outlet_coldside_HEX3(p(start = 101325), T(start = T_tank4_nom)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_hotside_HEX3(p(start = 100000), T(start = T_tank3_nom)) "Medium properties ";
  SI.Power P_loss_irr_HEX3(displayUnit = "MW");
    SI.Energy E_loss_irr_HEX3(displayUnit = "MWh", start = 0, fixed = true);
  //-------------HEX 4 rejection discharge//
  SI.HeatFlowRate Q_dot_HEXrej(displayUnit = "MW");
   SI.Heat Q_HEXrej(displayUnit = "MWh", start = 0, fixed = true);
  SI.Power P_loss_irr_HEXrej(displayUnit = "MW");
   SI.Energy E_loss_irr_HEXrej(displayUnit = "MWh", start = 0, fixed = true);
  //-------------COMPRESSOR DISCHARGE/
  Real c1 "factor 1";
  Real c2 "factor 2";
  Real c3 "factor 3";
  //design
  parameter Real beta_CO_nom = 5.911 "design compression ratio";
  parameter Real n_CO_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom = 0.88 "design isentropic efficiency";
  SI.Temperature T_1_nom = from_degC(-53.76) "state 1 temperature";
  parameter SI.Pressure p_1_nom = 100000 "state 1 pressure";
  //actual
  Real beta_CO(start = beta_CO_nom) "absolute compression ratio";
  parameter Real n_CO = 3000 "actual speed";
  SI.Efficiency eta_is_CO "absolute isentropic efficiency";
  //reduced
  Real beta_CO_red(start = 1) "reduced compression ratio";
  Real n_CO_red(start = 1) "reduced speed";
  SI.Efficiency eta_is_CO_red "reduced isentropic efficiency";
  Real G_CO_red(start = 1) "reduced mass flow rate compressor";
  //other
  SI.Power P_mech_CO(displayUnit = "MW");
  SI.Power P_loss_irr_CO(displayUnit = "MW");
   SI.Energy E_loss_irr_CO(displayUnit = "MWh", start = 0, fixed = true);
  //limits
  Real beta_CO_red_min;
  Real beta_CO_red_max;
  //-------------EXPANDER DISCHARGE//
  Real alpha "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom = 5.445 "design expansion ratio";
  parameter Real n_TU_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom = 0.92 "design isentropic efficiency";
  parameter SI.Temperature T_3_nom = from_degC(555) "design turbine inlet temperature";
  parameter SI.Pressure p_3_nom = 572336 "design turbine inlet pressure";
  //actual
  Real beta_TU(start = beta_TU_nom) "absolute expansion ratio";
  Real n_TU(start = 3000) "actual speed";
  SI.Efficiency eta_is_TU "absolute isentropic turbine efficiency";
  //reduced
  Real beta_TU_red(start = 1) "reduced expansion ratio";
  Real n_TU_red(start = 1) "reduced speed";
  SI.Efficiency eta_is_TU_red "reduced isentropic turbine efficiency";
  Real G_TU_red(start = 1) "reduced mass flow rate turbine";
  //other
  SI.Power P_mech_TU(displayUnit = "MW");
  SI.Power P_loss_irr_TU(displayUnit = "MW");
    SI.Energy E_loss_irr_TU(displayUnit = "MWh", start = 0, fixed = true);
  //-------------SYSTEM DISCHARGE//
  SI.Power P_mech_shaft(displayUnit = "MW");
  SI.HeatFlowRate Q_pump(displayUnit = "MW");
  Real eta_heat_to_power(start = 1);
      Real eta_heat_to_power_sys(start = 1);  
  Real work_ratio;
     //limits
  parameter Real beta_TU_red_min = 0.4;
  parameter Real beta_TU_red_max = 1.4;
  //-------------STATES//
  //-------------Discharge//
  //fixed temperature point 
  SI.Temperature T_1_a = from_degC(34.105) "outlet temperature after rejec";
  //fixed pressure point 
  SI.Pressure p_1 = p_fix "state 1 pressure";
  //STATE 1 a discharge
  SI.Pressure p_1_a(start=103605) "pressure after Heat rejection ";
  WorkingFluid.ThermodynamicState state_1_a "thermodynamic state after Heat rejection ";
  WorkingFluid.SpecificEnthalpy h_1_a(start = 433595) "turbine-side recuperation after Heat rejection ";
  WorkingFluid.SpecificEntropy s_1_a "turbine-side recuperation after Heat rejection ";
  //state 1 discharge
  SI.Temperature T_1_guess(start = from_degC(-52.760)) " temperature";
  //state 1 discharge
  SI.Temperature T_1(start = from_degC(-52.760)) " temperature";
  WorkingFluid.ThermodynamicState state_1 "thermodynamic state of inlet";
  WorkingFluid.SpecificEnthalpy h_1(start = -53000) "inlet enthalpy";
  WorkingFluid.SpecificEntropy s_1(start = -50) "inlet spec. entropy";
  //STATE 2 isentropic discharge
  SI.Temperature T_2_is(start = from_degC(100)) "isentropic outlet temperature of compressor";
  WorkingFluid.ThermodynamicState state_2_is "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_2_is(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_2_is "compressor outlet spec. entropy";
  //STATE 2 discharge
  SI.Temperature T_2(start = from_degC(114)) "actual outlet temperature of compressor";
  SI.Pressure p_2(start = 591100) "pressure coming out of compressor";
  WorkingFluid.ThermodynamicState state_2 "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_2(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_2 "compressor outlet spec. entropy";
  //STATE 3a discharge
  SI.Pressure p_3_a(start = 581494) "Pressure after recup";
  SI.Temperature T_3_a(start = from_degC(260)) "outlet temperature after recuperation";
  WorkingFluid.ThermodynamicState state_3_a "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_a(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a "compressor outlet spec. entropy";
  //STATE 3 discharge
  SI.Pressure p_3(start = 572045) "Pressure after recup";
  SI.Temperature T_3(start = from_degC(556)) "outlet temperature after HEX";
  WorkingFluid.ThermodynamicState state_3 "thermodynamic state of HEX outlet";
  WorkingFluid.SpecificEnthalpy h_3(start = 579480) "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3 "HEX outlet spec. entropy";
  //STATE 3 discharge guess
  SI.Temperature T_3_guess(start = from_degC(556));
  //STATE 4 isentropic
  SI.Temperature T_4_is(start = from_degC(250)) "isentropic outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4_is "isentropic state of turbine outlet";
  WorkingFluid.SpecificEnthalpy h_4_is(start = 657274) "turbine outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_is "isentropic turbine outlet spec. entropy";
  //STATE 4 discharge
  SI.Pressure p_4(start = 105039) "pressure at turb outlet";
  WorkingFluid.SpecificEnthalpy h_4(start = 683188) "turbine outlet enthalpy";
  SI.Temperature T_4(start = from_degC(270)) "outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4(p(start = 105039), T(start = from_degC(270))) "thermodynamic state of turbine outlet";
  WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";
  //STATE 4 guess discharge
  SI.Temperature T_4_guess(start = from_degC(270)) "outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4_guess(p(start = 105039), T(start = from_degC(270))) "thermodynamic state of turbine outlet";
  WorkingFluid.SpecificEnthalpy h_4_guess(start = 683188) "turbine outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_guess "turbine outlet spec. entropy";
  //STATE 4 a discharge
  SI.Pressure p_4_a(start = 103332) "pressure at recuperation outlet";
  SI.Temperature T_4_a(start = from_degC(270)) "outlet temperature after recuperation";
  WorkingFluid.ThermodynamicState state_4_a "thermodynamic state of turbine-side recuperation outlet";
  WorkingFluid.SpecificEnthalpy h_4_a(start = 520945) "turbine-side recuperation outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_a "turbine-side recuperation outlet spec. entropy";
  //------------------ELECTRICAL MACHINERY
          //Grid connection point   
  Modelica.Units.SI.Power P_GC(displayUnit = "MW");  
  parameter Modelica.Units.SI.ReactivePower Q_GC_set(displayUnit = "Mvar") = -32.5*1000*1000;
  Modelica.Units.SI.ReactivePower S_GC;  
  Real cos_phi_GC(start=0.95);
  //connection point set parameters
  parameter SI.Voltage U_TR_set = 220*1000;
  parameter SI.Angle phi_TR_set = 0;
  SI.ReactivePower Q_TR_set(displayUnit = "MW",start = -50*1000*1000);  
 //------------------TRANSFORMER
  //transformer parameters
  parameter SI.ApparentPower S_TR_nom(displayUnit = "MVA") = 208*1000*1000;
  parameter SI.Voltage U_TR_HV_nom(displayUnit = "kV") = 220*1000 "RMS voltage of the high voltage side (fixed by upper grid), line-to-line";
  parameter SI.Voltage U_TR_LV_nom(displayUnit = "kV") = 15750 "RMS voltage of the low voltage side, line-to-line";
  parameter Real a = U_TR_HV_nom/U_TR_LV_nom "turns ratio";
  parameter Real u_TR_ohmic_OC = 0.33604996008458077;
  parameter Real i_TR_OC = 0.32272375619545934;
  parameter Real P_TR_loss_OC_div_S_TR_nom = 0.06121076111169815;
  parameter SI.Power P_TR_loss_OC(displayUnit = "kW") = (P_TR_loss_OC_div_S_TR_nom/100)*S_TR_nom;
  parameter Real u_TR_SC = 11.875135026230133;
  parameter SI.Current I_TR_max = 537.638 "maximum transformer current";
  SI.Angle phi_U_HV(start = 0) "Phase of voltage at high-voltage side of transformer";
  SI.Angle phi_U_HV_min_I_HV "Phase shift between voltage and current high-voltage side of transformer";
  SI.Angle phi_I_HV "Phase of current at high-voltage side of transformer";
  Real pf_TR_HV "Power factor at high-voltage side of transformer";
  //whole transformer
  Complex S_delta;
  SI.Power P_TR_loss(displayUnit = "kW");
  SI.Power P_TR_ohmic_loss_LV(displayUnit = "kW");
  SI.Power P_TR_ohmic_loss_HV(displayUnit = "kW");
  SI.Power P_TR_iron_loss(displayUnit = "kW");
  SI.Energy E_TR_loss(displayUnit = "MWh", start = 0, fixed = true);
  SI.ReactivePower Q_req(displayUnit = "Mvar", start = (S_TR_nom*u_TR_SC/100));
  Real eta_transformer(start = 0.995);
  //high voltage side
  SI.Impedance Z_TR_HV;
  SI.ReactivePower Q_TR_HV(displayUnit = "Mvar");
  SI.Power P_TR_HV(displayUnit = "MW");
  Complex S_TR_HV;
  SI.ApparentPower S_TR_HV_abs(displayUnit = "MVA", start = S_TR_nom);
  SI.Voltage U_TR_HV_phase(start = U_TR_HV_nom/sqrt(3)) "phase voltage of reference phase";
  SI.ComplexVoltage U_TR_HV(re(start = U_TR_HV_nom), im(start = 0)) "complex phase voltage of reference phase";
  //Current I_TR_HV_nom;
  SI.ComplexCurrent I_TR_HV;
  SI.ComplexCurrent I_mag;
  SI.ComplexVoltage U_TR_h_HV;
  Complex Z_TR_h;
  SI.Resistance R_TR_HV;
  SI.Reactance X_TR_HV;
  SI.Current I_TR_HV_abs(start = S_TR_nom/(U_TR_HV_nom*3));
  //low voltage side
  SI.ComplexCurrent I_TR_LV_transferred;
  SI.ComplexVoltage U_TR_LV_transferred;
  SI.ComplexCurrent I_TR_LV;
  SI.ComplexVoltage U_TR_LV(re(start = U_TR_LV_nom)) "phase voltage of reference phase";
  Complex S_TR_LV;
  SI.ApparentPower S_TR_LV_abs(displayUnit = "MVA");
  SI.Power P_TR_LV(displayUnit = "MW");
  SI.ReactivePower Q_TR_LV(displayUnit = "Mvar");
  SI.Resistance R_TR_LV;
  SI.Reactance X_TR_LV;
  SI.Resistance R_TR_LV_transferred;
  SI.Reactance X_TR_LV_transferred;
  SI.Impedance Z_TR_OC_HV;
  SI.Resistance R_TR_FE_HV;
  SI.Reactance X_TR_h_HV;
     //------------------SYNCHRONOUS MACHINE
  //parameters
  parameter SI.ApparentPower S_SM_nom(displayUnit = "MVA") = 208*1000*1000;
  parameter SI.Voltage U_SM_ST_nom(displayUnit = "kV") = 15.75*1000 "RMS voltage of the high voltage side";
  parameter SI.Frequency f = 50 " input frequency at stator";
  parameter Integer N_p = 1 "number of pole pairs";
  parameter NonSI.AngularVelocity_rpm n_SM_nom = 3000 "nominal synchronous speed in rpm";
  //geometry
  parameter SI.Diameter d_RO = 0.92 "Rotor diameter";
  parameter SI.Length l_FE = 4.715853752596681 "Rotor length";
  parameter SI.Mass m_ST_FE_teeth = 10875.044312382426 "iron mass of stator teeth";
  parameter SI.Mass m_ST_FE_yoke = 108758.707701592 "iron mass of stator yoke";
  parameter Real tau_p = 1.4451326206513049 "pole pitch";
  //Equivalent circuit parameters
  parameter SI.Inductance L_sigma_ST = 0.0011957997612770178;
  parameter SI.Inductance L_h = 0.006776198647236435;
  parameter SI.Inductance L_d = 0.007971998408513453;
  parameter Integer N_FD = 72 "number of exitation field windings";
  parameter Integer N_ST = 11 "number of stator windings";
  // losses
  parameter Real k_schuisky = 5 "experimental factor for correlation of Schuisky";
  parameter Real k_additional = 0.001 "factor for additional losses";
  parameter SI.Resistance R_FD = 0.13112806419264864 "resistance in excitation windings";
  parameter Real P10 = 1.3 "power loss factor in W/kg";
  parameter Real k_Fe_stator_yoke = 1.3 "correction factor for harmonics";
  parameter Real k_Fe_stator_teeth = 1.7 "correction factor for harmonics";
  parameter SI.MagneticFluxDensity B_ST_yoke = 1.35 "flux density in stator yoke";
  parameter SI.MagneticFluxDensity B_ST_teeth = 1.75 "flux density in stator teeth";
  //variables
  //energy
  Real eta_SM(start = 0.988);
  SI.Energy E_SM_loss(displayUnit = "MWh", start = 0, fixed = true);
  //stator-side variables
  SI.ReactivePower Q_SM_ST(displayUnit = "Mvar");
  SI.ComplexVoltage U_SM_ST(re(start = U_SM_ST_nom/sqrt(3)));
  //per phase value!
  SI.Power P_SM_ST(displayUnit = "MW");
  Complex S_SM_ST;
  SI.ApparentPower S_SM_ST_abs(displayUnit = "MVA",start=S_SM_nom);
  SI.ComplexCurrent I_SM_ST "stator side complex current";
  //air-gap
  SI.ComplexVoltage U_SM_h(re(start = U_SM_ST_nom/sqrt(3)));
  Complex S_SM_h;
  SI.Power P_airgap(displayUnit = "MW");
  //excitation
  SI.ComplexVoltage U_P(re(start = U_SM_ST_nom/sqrt(3)));
  SI.ComplexCurrent I_FD_ref "excitation current, referred to stator side";
  SI.Voltage U_FD "excitation voltage";
  SI.Power P_excitation(displayUnit = "MW") "power required for DC excitation";
  SI.Current I_FD(start = 900) "actual excitation current, DC";
  //equivalent circuit elements
  SI.Reactance X_sigma_ST;
  SI.Reactance X_h;
  SI.Reactance X_d;
  SI.Resistance R_SM_ST = 0.002981520432692308;
  //mechanical
  NonSI.AngularVelocity_rpm n_RO "rotor speed";
  SI.AngularVelocity omega_SM;
  SI.Power P_mech_RO(displayUnit = "MW");
  //losses
  SI.Power P_SM_loss(displayUnit = "MW") "total losses";
  SI.Velocity v_RO "rotor perimeter speed";
  SI.Power P_windage_ventilation(displayUnit = "MW") "losses due to windage and ventilation";
  SI.Power P_additional(displayUnit = "MW") "additional losses";
  SI.Power P_FE_ST_teeth(displayUnit = "MW") "iron losses in stator teeth";
  SI.Power P_FE_ST_yoke(displayUnit = "MW") "iron losses in stator yoke";
  SI.Power P_FE(displayUnit = "MW") "total iron losses";
  SI.Power P_Ohmic_Stator(displayUnit = "MW") "ohmic losses stator";

  Modelica.Blocks.Continuous.SecondOrder T4_a_guess_control(D = 0.4, w = 0.5) annotation(
    Placement(transformation(origin = {-66, 46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder T4_guess_control_discharge(w = 0.5, D = 0.4) annotation(
    Placement(transformation(origin = {56, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder T3_guess_control_discharge(w = 0.5, D = 0.4) annotation(
    Placement(transformation(origin = {58, -2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder T1_guess_control_discharge(w = 0.5, D = 0.4) annotation(
    Placement(transformation(origin = {58, -40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder T4_guess_control(w = 0.5, D = 0.4) annotation(
    Placement(transformation(origin = {-64, -38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Step step_down(height = -1, offset = 1, startTime = 36000)  annotation(
    Placement(transformation(origin = {-28, -76}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Step step_up(height = 2, offset = 0, startTime = 72000)  annotation(
    Placement(transformation(origin = {12, -76}, extent = {{-10, -10}, {10, 10}})));
initial equation
//--------------------------INITIAL EQUATIONS-----------------------------//
//control loops
//charge
  T_4_charge_guess = from_degC(274.5);
  T_4_a_charge_guess = from_degC(18.14);
//discharge
  T_4_guess = from_degC(270);
  T_1_guess = from_degC(-52.760);
  T_3_guess = from_degC(556);
//tanks
  m_tank1 = m_tank1_start;
  T_tank1 = T_tank1_start;
  m_tank2 = m_tank2_start;
  T_tank2 = T_tank2_start;
  m_tank3 = m_tank3_start;
  T_tank3 = T_tank3_start;
  m_tank4 = m_tank4_start;
  T_tank4 = T_tank4_start;
equation
  Mode=step_up.y+step_down.y;
//--------------------------EQUATIONS SYSTEM-----------------------------//
  state_amb_air = WorkingFluid.setState_pT(101315, T_amb);
  P_elec_charge = P_TR_HV;
  P_elec = P_TR_HV;
//-------------SYSTEM TANKS//
  exergy_total_tanks = exergy_tank1 + exergy_tank2 + exergy_tank3 + exergy_tank4;
    exergy_hot_tanks = exergy_tank1 + exergy_tank2;
  exergy_cold_tanks = exergy_tank3 + exergy_tank4;
  int_energy_total_tanks = int_energy_tank1 + int_energy_tank2 + int_energy_tank3 + int_energy_tank4;
  Q_dot_hightemp_res = der(int_energy_tank1 + int_energy_tank2);
//-------------SYSTEM CHARGE//
  P_mech_shaft_charge = P_mech_CO_charge + P_mech_TU_charge;
  Q_pump_charge = Q_dot_HEX1_charge;
  COP = Q_pump_charge/P_mech_shaft_charge;
  work_ratio_charge = P_mech_CO_charge/abs(P_mech_TU_charge);
  m_dot_div_p_charge = m_dot_WF_charge/p_4_charge;
    der(exergy_total_loss_irr_charge) = P_total_loss_irr_charge;
  P_total_loss_irr_charge = P_loss_irr_HEX1_charge + P_loss_irr_HEX2_charge + P_loss_irr_HEX3_charge + P_loss_irr_HEXrej_charge + P_loss_irr_CO_charge + P_loss_irr_TU_charge;
  der(E_mech_shaft_charge) = P_mech_shaft_charge;
  hot_to_cold_mass_flow_ratio_charge = m_dot_solsalt_HEX1_charge/m_dot_methanol_HEX3_charge;
  E_total_loss_irr_charge = E_loss_irr_CO_charge + E_loss_irr_TU_charge + E_loss_irr_HEX1_charge + E_loss_irr_HEX2_charge + E_loss_irr_HEX2_charge + E_loss_irr_HEXrej_charge + E_TR_loss + E_SM_loss;
//-------------SYSTEM DISCHARGE//
  P_mech_shaft = P_mech_CO + P_mech_TU;
  Q_pump = -Q_dot_HEX1;
  eta_heat_to_power = P_mech_shaft/Q_pump;
      eta_heat_to_power_sys = P_elec/Q_pump;  
  work_ratio = abs(P_mech_TU)/P_mech_CO;
  m_dot_div_p = m_dot_WF/p_1;
    der(exergy_total_loss_irr) = P_total_loss_irr;
  P_total_loss_irr = P_loss_irr_HEX1 + P_loss_irr_HEX2 + P_loss_irr_HEX3 + P_loss_irr_HEXrej + P_loss_irr_CO + P_loss_irr_TU;
  der(E_mech_shaft) = P_mech_shaft;
  hot_to_cold_mass_flow_ratio = m_dot_solsalt_HEX1/m_dot_methanol_HEX3;
  E_total_loss_irr = E_loss_irr_CO + E_loss_irr_TU + E_loss_irr_HEX1 + E_loss_irr_HEX2 + E_loss_irr_HEX3 + E_loss_irr_HEXrej + E_TR_loss + E_SM_loss;
//MODE 1 CHARGE
  if Mode == 1 then
    der(Elec_energy_charge) = P_elec_charge;
    der(Elec_energy_discharge) = 0;
    Q_GC_set=Q_TR_set;
    eta_transformer = abs(P_TR_LV)/(abs(P_TR_LV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;
    P_additional = k_additional*abs(P_mech_RO);
    P_SM_loss = abs(P_SM_ST) - abs(P_mech_RO);
    P_mech_RO = P_airgap - P_windage_ventilation - P_additional - P_FE - P_excitation;
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST +U_P ;
    eta_SM = abs(P_mech_RO)/abs(P_SM_ST);
    //connection to thermodynamic cycle
    P_mech_shaft_charge=P_mech_RO; 
        COP_system_charge=Q_pump_charge/P_elec_charge;     
    COP_system=Q_pump/P_elec;     
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    der(Elec_energy_charge) = 0;
    der(Elec_energy_discharge) = P_elec;
    Q_GC_set=Q_TR_set; 
    eta_transformer = abs(P_TR_HV)/(abs(P_TR_HV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;
    P_additional = k_additional*abs(P_mech_RO);
    P_SM_loss = abs(P_mech_RO) - abs(P_SM_ST);
    P_mech_RO = P_airgap - (P_windage_ventilation + P_additional + P_FE + P_excitation);
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST + U_P;
    eta_SM = abs(P_SM_ST)/abs(P_mech_RO);
    //connection to thermodynamic cycle
    P_mech_shaft=P_mech_RO;    
        COP_system_charge=Q_pump_charge/P_elec_charge;     
    COP_system=Q_pump/P_elec;
//MODE 0 HOLD
  else
    der(Elec_energy_charge) = 0;
    der(Elec_energy_discharge) = 0;
    Q_TR_set=0;
        P_TR_HV=0;
    eta_transformer = 0;
    I_TR_HV + I_TR_LV_transferred = Complex(0, 0);
    P_additional = 0;
    P_SM_loss = 0; 
    I_SM_ST = Complex(0, 0);
    U_P = Complex(0, 0);
    eta_SM = 0;
    P_mech_RO = 0;
    COP_system_charge=0;
    COP_system =0;   
    end if;
//--------------------------EQUATIONS TANKS-----------------------------//
//-------------TANK 1//
//nominal
  solsalt_tank1_nom = HotTESLiquid.setState_pT(p_tank1_nom, T_tank1_nom);
// connectors
  h_in_tank1 = outlet_coldside_HEX1_charge.h;
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank1 = 0;
    m_in_tank1 = m_dot_solsalt_HEX1_charge;
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    m_out_tank1 = -m_dot_solsalt_HEX1;
    m_in_tank1 = 0;
//MODE 0 HOLD
  else
    m_out_tank1 = 0;
    m_in_tank1 = 0;
  end if;
//mass balance
  der(m_tank1) = m_in_tank1 + m_out_tank1;
//energy balance
  der(m_tank1*solsalt_tank1.u) = m_out_tank1*solsalt_tank1.h + m_in_tank1*h_in_tank1 + Q_dot_to_amb_tank1;
//states
  SOC_tank1 = (m_tank1 - m_tank1_min)/m_working_solar_salt;
  m_tank1 = solsalt_tank1.d*A_cross_tank_solsalt*x_tank1;
  int_energy_tank1 = m_tank1*solsalt_tank1.u;
 exergy_tank1 = m_tank1*solsalt_tank1.u - T_amb*m_tank1*solsalt_tank1_state.s;
//thermodynamic states
  solsalt_tank1_state = HotTESLiquid.setState_pT(p_tank1_nom, T_tank1);
//solar salt properties
  solsalt_tank1.p = p_tank1_nom;
  solsalt_tank1.T = T_tank1;
//losses
  Q_div_A_tank1 = (0.00017*(T_tank1 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank1 = -Q_div_A_tank1*A_W_tank_solsalt;
// Q_dot_to_amb_tank1 = 0;
//-------------TANK 2//
//nominal
  solsalt_tank2_nom = HotTESLiquid.setState_pT(p_tank2_nom, T_tank2_nom);
// connectors
  h_in_tank2 = outlet_hotside_HEX1.h;
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank2 = -m_dot_solsalt_HEX1_charge;
    m_in_tank2 = 0;
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    m_out_tank2 = 0;
    m_in_tank2 = m_dot_solsalt_HEX1;
//MODE 0 HOLD
  else
    m_out_tank2 = 0;
    m_in_tank2 = 0;
  end if;
//mass balance
  der(m_tank2) = m_in_tank2 + m_out_tank2;
//energy balance
  der(m_tank2*solsalt_tank2.u) = m_out_tank2*solsalt_tank2.h + m_in_tank2*h_in_tank2 + Q_dot_to_amb_tank2;
//states
  SOC_tank2 = (m_tank2 - m_tank2_min)/m_working_solar_salt;
  m_tank2 = solsalt_tank2.d*A_cross_tank_solsalt*x_tank2;
  int_energy_tank2 = m_tank2*solsalt_tank2.u;
  exergy_tank2 = m_tank2*solsalt_tank2.u - T_amb*m_tank2*solsalt_tank2_state.s;
//thermodynamic states
  solsalt_tank2_state = HotTESLiquid.setState_pT(p_tank2_nom, T_tank2);
//solar salt properties
  solsalt_tank2.p = p_tank2_nom;
  solsalt_tank2.T = T_tank2;
//losses
  Q_div_A_tank2 = (0.00017*(T_tank2 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank2 = -Q_div_A_tank2*A_W_tank_solsalt;
// Q_dot_to_amb_tank2 = 0;
//-------------TANK 3//
//nominal
  coldliq_tank3_nom = ColdTESLiquid.setState_pT(p_tank3_nom, T_tank3_nom);
// connectors
  h_in_tank3 = outlet_coldside_HEX3.h;
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank3 = -m_dot_methanol_HEX3_charge;
    m_in_tank3 = 0;
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    m_out_tank3 = 0;
    m_in_tank3 = m_dot_methanol_HEX3;
//MODE 0 HOLD
  else
    m_out_tank3 = 0;
    m_in_tank3 = 0;
  end if;
//mass balance
  der(m_tank3) = m_in_tank3 + m_out_tank3;
//energy balance
  der(m_tank3*coldliq_tank3.u) = m_out_tank3*coldliq_tank3.h + m_in_tank3*h_in_tank3;
//states
  SOC_tank3 = (m_tank3 - m_tank3_min)/m_working_methanol;
  m_tank3 = coldliq_tank3.d*A_cross_tank_coldliq*x_tank3;
  int_energy_tank3 = m_tank3*coldliq_tank3.u;
 exergy_tank3 = m_tank3*coldliq_tank3.u - T_amb*m_tank3*coldliq_tank3_state.s;
//thermodynamic states
  coldliq_tank3_state = ColdTESLiquid.setState_pT(p_tank3_nom, T_tank3);
//solar salt properties
  coldliq_tank3.p = p_tank3_nom;
  coldliq_tank3.T = T_tank3;
//-------------TANK 4//
//nominal
  coldliq_tank4_nom = ColdTESLiquid.setState_pT(p_tank4_nom, T_tank4_nom);
// connectors
  h_in_tank4 = outlet_hotside_HEX3_charge.h;
//h_in_tank4=0; //dummy
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank4 = 0;
    m_in_tank4 = m_dot_methanol_HEX3_charge;
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    m_out_tank4 = -m_dot_methanol_HEX3;
    m_in_tank4 = 0;
//MODE 0 HOLD
  else
    m_out_tank4 = 0;
    m_in_tank4 = 0;
  end if;
//mass balance
  der(m_tank4) = m_in_tank4 + m_out_tank4;
//energy balance
  der(m_tank4*coldliq_tank4.u) = m_out_tank4*coldliq_tank4.h + m_in_tank4*h_in_tank4;
//states
  SOC_tank4 = (m_tank4 - m_tank4_min)/m_working_methanol;
  m_tank4 = coldliq_tank4.d*A_cross_tank_coldliq*x_tank4;
  int_energy_tank4 = m_tank4*coldliq_tank4.u;
  exergy_tank4 = m_tank4*(coldliq_tank4.u) - T_amb*m_tank4*coldliq_tank4_state.s;
//thermodynamic states
  coldliq_tank4_state = ColdTESLiquid.setState_pT(p_tank4_nom, T_tank4);
//solar salt properties
  coldliq_tank4.p = p_tank4_nom;
  coldliq_tank4.T = T_tank4;
//--------------------------EQUATIONS CHARGE-----------------------------//
//-------------COMPRESSOR CHARGE//
//reduced values compressor
  eta_is_CO_red_charge = eta_is_CO_charge/eta_is_CO_nom_charge;
  beta_CO_red_charge = beta_CO_charge/beta_CO_nom_charge;
  G_CO_red_charge = m_dot_WF_charge*sqrt(T_4_charge_guess)/p_4_charge/(m_dot_WF_nom_charge*sqrt(T_4_nom_charge)/p_4_nom_charge);
  eta_is_CO_red_charge = (1 - c4*(1 - n_CO_red_charge)^2)*(n_CO_red_charge/G_CO_red_charge)*(2 - n_CO_red_charge/G_CO_red_charge);
  n_CO_red_charge = n_CO_charge/sqrt(T_4_charge_guess)/(n_CO_nom_charge/sqrt(T_4_nom_charge));
//other compressor equations
  beta_CO_charge = p_3_charge/p_4_charge;
  eta_is_CO_charge = (h_3_is_charge - h_4_charge_guess)/(h_3_charge - h_4_charge_guess);
  beta_CO_red_charge = c1_charge*G_CO_red_charge^2 + c2_charge*G_CO_red_charge + c3_charge;
  c1_charge = n_CO_red_charge/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c2_charge = (p - 2*m*n_CO_red_charge^2)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c3_charge = -1*(p*m*n_CO_red_charge - m^2*n_CO_red_charge^3)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  P_mech_CO_charge = m_dot_WF_charge*(h_3_charge - h_4_charge);
  P_loss_irr_CO_charge = T0*m_dot_WF_charge*(s_3_charge - s_4_charge);
    der(E_loss_irr_CO_charge) = P_loss_irr_CO_charge;
//limits
  beta_CO_red_charge_min = 0.44008400621740806*(G_CO_red_charge)^2 + 0.15193846169229142*G_CO_red_charge + 0.1260521481189032;
  beta_CO_red_charge_max = -0.3505543472709579*(G_CO_red_charge)^2 + 1.6961020320713247*G_CO_red_charge - 0.17150246895420682;
//-------------EXPANDER CHARGE//
 n_TU_charge = n_CO_charge;
//parameters
  alpha_charge = sqrt(1.4 - 0.4*n_TU_charge/n_TU_nom_charge);
//reduced values expander
  n_TU_red_charge = n_TU_charge/sqrt(T_2_a_charge)/(n_TU_nom_charge/sqrt(T_2_a_nom_charge));
  G_TU_red_charge = m_dot_WF_charge*sqrt(T_2_a_charge)/p_2_a_charge/(m_dot_WF_nom_charge*sqrt(T_2_a_nom_charge)/p_2_a_nom_charge);
  G_TU_red_charge = alpha_charge*sqrt((1/(beta_TU_charge^2) - 1)/(1/(beta_TU_nom_charge^2) - 1));
  beta_TU_red_charge = beta_TU_charge/beta_TU_nom_charge;
  eta_is_TU_red_charge = (1 - t*(1 - n_TU_red_charge)^2)*(n_TU_red_charge/G_TU_red_charge)*(2 - ((n_TU_red_charge/G_TU_red_charge)));
  eta_is_TU_red_charge = eta_is_TU_charge/eta_is_TU_nom_charge;
//other turbine equations
  beta_TU_charge = p_2_a_charge/p_1_charge;
  eta_is_TU_charge = (h_2_a_charge - h_1_charge)/(h_2_a_charge - h_1_is_charge);
  P_mech_TU_charge = m_dot_WF_charge*(h_1_charge - h_2_a_charge);
  P_loss_irr_TU_charge = T0*m_dot_WF_charge*(s_1_charge - s_2_a_charge);
    der(E_loss_irr_TU_charge) = P_loss_irr_TU_charge;
//-------------HEX 1 CHARGE//
//pressure loss
  delta_P_HEX1_charge = k_p_charge*m_dot_WF_charge^2;
  p_3_charge = p_3_a_charge + delta_P_HEX1_charge;
//off-design
  UA_HEX1_charge/UA_HEX1_nom = ((m_dot_WF_charge^0.8*m_dot_solsalt_HEX1_charge^0.8)/(m_dot_WF_nom_charge^0.8*m_dot_solsalt_HEX1_nom^0.8))*((m_dot_WF_nom_charge^0.8 + m_dot_solsalt_HEX1_nom^0.8)/(m_dot_WF_charge^0.8 + m_dot_solsalt_HEX1_charge^0.8));
//charge: inlet cold side= tank 2 inlet hot side=state 3
  outlet_hotside_guess_HEX1_charge = WorkingFluid.setState_pT(p_3_a_charge, T_tank2 + 9);
//WF
  outlet_coldside_guess_HEX1_charge = HotTESLiquid.setState_pT(p_tank1_nom, T_3_charge - 9);
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
//eff_HEX1_charge = NTU_HEX1_charge/(1 + NTU_HEX1_charge);
  if C_r_HEX1_charge < 1 then
    eff_HEX1_charge = (1 - exp(-NTU_HEX1_charge*(1 - C_r_HEX1_charge)))/(1 - C_r_HEX1_charge*exp(-NTU_HEX1_charge*(1 - C_r_HEX1_charge)));
  else
    eff_HEX1_charge = NTU_HEX1_charge/(1 + NTU_HEX1_charge);
  end if;
  Q_dot_HEX1_charge = eff_HEX1_charge*Q_dot_max_HEX1_charge;
  Q_dot_max_HEX1_charge = C_min_HEX1_charge*(T_3_charge - T_tank2);
//calc temperature of fluid going into T_tank1 from cold side energy balance
  Q_dot_HEX1_charge = m_dot_solsalt_HEX1_charge*(h_cold_out_HEX1_charge - solsalt_tank2.h);
//calc T_3_a from hot side energy balance
  Q_dot_HEX1_charge = m_dot_WF_charge*(state_3_charge.h - h_hot_out_HEX1_charge);
//set outlets
  outlet_coldside_HEX1_charge = HotTESLiquid.setState_ph(p_tank1_nom, h_cold_out_HEX1_charge);
//NaK
  outlet_hotside_HEX1_charge = WorkingFluid.setState_ph(p_3_a_charge, h_hot_out_HEX1_charge);
//WF
  T_3_a_charge = outlet_hotside_HEX1_charge.T;
//irrev
  P_loss_irr_HEX1_charge = T0*((m_dot_WF_charge*(s_3_a_charge - s_3_charge)) + (m_dot_solsalt_HEX1_charge*(outlet_coldside_HEX1_charge.s - solsalt_tank2.s)));
  der(E_loss_irr_HEX1_charge) = P_loss_irr_HEX1_charge;
//-------------HEX2 (Recuperation) CHARGE//
//pressure loss
  delta_P_HEX2_charge = k_p_charge*m_dot_WF_charge^2;
  p_4_a_charge = p_4_charge + delta_P_HEX2_charge;
  p_2_charge = p_3_a_charge - delta_P_HEX2_charge;
//off-design
  UA_HEX2_charge/UA_HEX2_nom = (m_dot_WF_charge^0.8)/(m_dot_WF_nom_charge^0.8);
//cold side
  cp_cold_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_4_a_charge_guess);
//just inlet cp, not average cp
  C_cold_HEX2_charge = m_dot_WF_charge*cp_cold_HEX2_charge;
// hot side
  cp_hot_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_3_a_charge);
//just inlet cp, not average cp
  C_hot_HEX2_charge = m_dot_WF_charge*cp_hot_HEX2_charge;
//variables for effectiveness
  C_min_HEX2_charge = min(C_cold_HEX2_charge, C_hot_HEX2_charge);
  C_max_HEX2_charge = max(C_cold_HEX2_charge, C_hot_HEX2_charge);
  C_r_HEX2_charge = C_min_HEX2_charge/C_max_HEX2_charge;
  NTU_HEX2_charge = UA_HEX2_charge/C_min_HEX2_charge;
//eff_HEX2_charge = (1 - exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)))/(1 - C_r_HEX2_charge*exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)));
  if C_r_HEX2_charge < 1 then
    eff_HEX2_charge = (1 - exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)))/(1 - C_r_HEX2_charge*exp(-NTU_HEX2_charge*(1 - C_r_HEX2_charge)));
  else
    eff_HEX2_charge = NTU_HEX2_charge/(1 + NTU_HEX2_charge);
  end if;
//heat transferred
  Q_dot_HEX2_charge = eff_HEX2_charge*Q_dot_max_HEX2_charge;
  Q_dot_max_HEX2_charge = C_min_HEX2_charge*(T_3_a_charge - T_4_a_charge_guess);
//energy balance and outlet states
// hot side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(state_3_a_charge.h - h_hot_out_HEX2_charge);
//cold side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(h_cold_out_HEX2_charge - state_4_a_charge_guess.h);
//set outlets
  outlet_hotside_HEX2_charge = WorkingFluid.setState_ph(p_2_charge, h_hot_out_HEX2_charge);
  outlet_coldside_HEX2_charge = WorkingFluid.setState_ph(p_4_charge, h_cold_out_HEX2_charge);
  T_2_charge = outlet_hotside_HEX2_charge.T;
  T_4_charge = outlet_coldside_HEX2_charge.T;
//irrev
  P_loss_irr_HEX2_charge = T0*((m_dot_WF_charge*(s_2_charge - s_3_a_charge)) + (m_dot_WF_charge*(s_4_charge - s_4_a_charge)));
  der(E_loss_irr_HEX2_charge) = P_loss_irr_HEX2_charge;
//-------------HEX 3 CHARGE//
//pressure loss
  delta_P_HEX3_charge = k_p_charge*m_dot_WF_charge^2;
  p_1_charge = p_4_a_charge + delta_P_HEX3_charge;
//off-design
  UA_HEX3_charge/UA_HEX3_nom = ((m_dot_WF_charge^0.8*m_dot_methanol_HEX3_charge^0.8)/(m_dot_WF_nom_charge^0.8*m_dot_methanol_HEX3_nom^0.8))*((m_dot_WF_nom_charge^0.8 + m_dot_methanol_HEX3_nom^0.8)/(m_dot_WF_charge^0.8 + m_dot_methanol_HEX3_charge^0.8));
//charge: inlet cold side= state 1 inlet hot side=tank 3
  outlet_hotside_guess_HEX3_charge = ColdTESLiquid.setState_pT(p_tank4_nom, T_1_charge + 9);
//
//buffer to minimum temperature of methanol mixture through , which is fine, as it is a guess anyway, Tout, hot side will be higher
  outlet_coldside_guess_HEX3_charge = WorkingFluid.setState_pT(p_4_a_charge, T_tank3 - 7);
//
  cp_cold_ave_HEX3_charge = (WorkingFluid.specificHeatCapacityCp(state_1_charge) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX3_charge))/2;
  cp_hot_ave_HEX3_charge = (WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX3_charge) + coldliq_tank3_state.cp)/2;
//heat capacity rates
  C_cold_HEX3_charge = m_dot_WF_charge*cp_cold_ave_HEX3_charge;
  C_hot_HEX3_charge = C_cold_HEX3_charge;
//balanced operation
//m_dot_methanol_HEX3_charge = C_hot_HEX3_charge/cp_hot_ave_HEX3_charge;
  C_min_HEX3_charge = min(C_cold_HEX3_charge, C_hot_HEX3_charge);
  C_max_HEX3_charge = max(C_cold_HEX3_charge, C_hot_HEX3_charge);
  C_r_HEX3_charge = C_min_HEX3_charge/C_max_HEX3_charge;
  NTU_HEX3_charge = UA_HEX3_charge/C_min_HEX3_charge;
//eff_HEX3_charge = NTU_HEX3_charge/(1 + NTU_HEX3_charge);
//this eff_HEX3_charge formulation is computationally expensive
  if C_r_HEX3_charge < 1 then
    eff_HEX3_charge = (1 - exp(-NTU_HEX3_charge*(1 - C_r_HEX3_charge)))/(1 - C_r_HEX3_charge*exp(-NTU_HEX3_charge*(1 - C_r_HEX3_charge)));
  else
    eff_HEX3_charge = NTU_HEX3_charge/(1 + NTU_HEX3_charge);
  end if;
  Q_dot_HEX3_charge = eff_HEX3_charge*Q_dot_max_HEX3_charge;
  Q_dot_max_HEX3_charge = C_min_HEX3_charge*(T_tank3 - T_1_charge);
//calc temperature of fluid going into T_tank4 from hot side energy balance
  Q_dot_HEX3_charge = m_dot_methanol_HEX3_charge*(coldliq_tank3_state.h - h_hot_out_HEX3_charge);
//calc T_4_a_charge from cold side energy balance
  Q_dot_HEX3_charge = m_dot_WF_charge*(h_cold_out_HEX3_charge - state_1_charge.h);
//set outlets
  outlet_hotside_HEX3_charge = ColdTESLiquid.setState_ph(p_tank4_nom, h_hot_out_HEX3_charge);
  outlet_coldside_HEX3_charge = WorkingFluid.setState_ph(p_4_a_charge, h_cold_out_HEX3_charge);
  T_4_a_charge = outlet_coldside_HEX3_charge.T;
//irrev
  P_loss_irr_HEX3_charge = T0*((m_dot_WF_charge*(s_4_a_charge - s_1_charge)) + (m_dot_methanol_HEX3_charge*(outlet_hotside_HEX3_charge.s - coldliq_tank3_state.s)));
  der(E_loss_irr_HEX3_charge) = P_loss_irr_HEX3_charge;
  //-------------HEX 4 rejection charge//
  Q_dot_HEXrej_charge = (h_2_charge - h_2_a_charge)*m_dot_WF_charge;
    der(Q_HEXrej_charge) = Q_dot_HEXrej_charge;
//irrev
  P_loss_irr_HEXrej_charge = T0*((m_dot_WF_charge*(s_2_a_charge - s_2_charge - ((h_2_a_charge - h_2_charge)/T0))));
  der(E_loss_irr_HEXrej_charge) = P_loss_irr_HEXrej_charge;
//-------------STATES//
//state 1 charge
  state_1_charge = WorkingFluid.setState_ph(p_1_charge, h_1_charge);
  T_1_charge = state_1_charge.T;
  s_1_charge = WorkingFluid.specificEntropy(state_1_charge);
//state 1  isentropic charge
  s_1_is_charge = s_2_a_charge;
  state_1_is_charge = WorkingFluid.setState_ps(p_1_charge, s_1_is_charge) "isentropic state of turbine outlet";
  T_1_is_charge = state_1_is_charge.T "isentropic turbine outlet spec. entropy";
  h_1_is_charge = WorkingFluid.specificEnthalpy(state_1_is_charge);
//STATE 4_a charge guess
  T_4_a_charge = T4_a_guess_control.u;
  T_4_a_charge_guess = T4_a_guess_control.y;
  state_4_a_charge_guess = WorkingFluid.setState_pT(p_4_a_charge, T_4_a_charge_guess);
//STATE 4_a charge
  state_4_a_charge = WorkingFluid.setState_pT(p_4_a_charge, T_4_a_charge);
  h_4_a_charge = WorkingFluid.specificEnthalpy(state_4_a_charge);
  s_4_a_charge = WorkingFluid.specificEntropy(state_4_a_charge);
//STATE 4 charge guess
  T_4_charge = T4_guess_control.u;
  T_4_charge_guess = T4_guess_control.y;
  state_4_charge_guess = WorkingFluid.setState_pT(p_4_charge, T_4_charge_guess);
  h_4_charge_guess = WorkingFluid.specificEnthalpy(state_4_charge_guess);
  s_4_charge_guess = WorkingFluid.specificEntropy(state_4_charge_guess);
//STATE 4 charge
  state_4_charge = WorkingFluid.setState_pT(p_4_charge, T_4_charge);
  h_4_charge = WorkingFluid.specificEnthalpy(state_4_charge);
  s_4_charge = WorkingFluid.specificEntropy(state_4_charge);
//STATE 3 isentropic charge
  s_3_is_charge = s_4_charge_guess "isentropic compressor outlet spec. entropy";
  state_3_is_charge = WorkingFluid.setState_ps(p_3_charge, s_3_is_charge) "isentropic state of compressor outlet";
  T_3_is_charge = state_3_is_charge.T;
  h_3_is_charge = WorkingFluid.specificEnthalpy(state_3_is_charge);
//STATE 3 charge
  state_3_charge = WorkingFluid.setState_ph(p_3_charge, h_3_charge);
  T_3_charge = state_3_charge.T;
  s_3_charge = WorkingFluid.specificEntropy(state_3_charge);
//STATE 3a charge
  state_3_a_charge = WorkingFluid.setState_pT(p_3_a_charge, T_3_a_charge);
  h_3_a_charge = WorkingFluid.specificEnthalpy(state_3_a_charge);
  s_3_a_charge = WorkingFluid.specificEntropy(state_3_a_charge);
//state 2 charge
  state_2_charge = WorkingFluid.setState_pT(p_2_charge, T_2_charge);
  h_2_charge = WorkingFluid.specificEnthalpy(state_2_charge);
  s_2_charge = WorkingFluid.specificEntropy(state_2_charge);
//state 2_a charge
  p_2_a_charge = p_2_charge;
  state_2_a_charge = WorkingFluid.setState_pT(p_2_a_charge, T_2_a_charge);
  h_2_a_charge = WorkingFluid.specificEnthalpy(state_2_a_charge);
  s_2_a_charge = WorkingFluid.specificEntropy(state_2_a_charge);
//--------------------------EQUATIONS DISCHARGE-----------------------------//
//-------------COMPRESSOR DISCHARGE//
//reduced values compressor
  eta_is_CO_red = eta_is_CO/eta_is_CO_nom;
  beta_CO_red = beta_CO/beta_CO_nom;
  G_CO_red = m_dot_WF*sqrt(T_1_guess)/p_1/(m_dot_WF_nom*sqrt(T_1_nom)/p_1_nom);
  eta_is_CO_red = (1 - c4*(1 - n_CO_red)^2)*(n_CO_red/G_CO_red)*(2 - n_CO_red/G_CO_red);
  n_CO_red = n_CO/sqrt(T_1_guess)/(n_CO_nom/sqrt(T_1_nom));
//other compressor equations
  beta_CO = p_2/p_1;
  eta_is_CO = (h_2_is - h_1)/(h_2 - h_1);
  beta_CO_red = c1*G_CO_red^2 + c2*G_CO_red + c3;
  c1 = n_CO_red/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c2 = (p - 2*m*n_CO_red^2)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c3 = -1*(p*m*n_CO_red - m^2*n_CO_red^3)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  P_mech_CO = m_dot_WF*(h_2 - h_1);
  P_loss_irr_CO = T0*m_dot_WF*(s_2 - s_1);
  der(E_loss_irr_CO) = P_loss_irr_CO;
//limits
  beta_CO_red_min = 0.44008400621740806*(G_CO_red)^2 + 0.15193846169229142*G_CO_red + 0.1260521481189032;
  beta_CO_red_max = -0.3505543472709579*(G_CO_red)^2 + 1.6961020320713247*G_CO_red - 0.17150246895420682;
//-------------EXPANDER DISCHARGE//
  n_TU = n_CO;
//parameters
  alpha = sqrt(1.4 - 0.4*n_TU/n_TU_nom);
//reduced values turbine
  eta_is_TU_red = eta_is_TU/eta_is_TU_nom;
  n_TU_red = n_TU/sqrt(T_3_guess)/(n_TU_nom/sqrt(T_3_nom));
  beta_TU_red = beta_TU/beta_TU_nom;
   G_TU_red = m_dot_WF*sqrt(T_3_guess)/p_3/(m_dot_WF_nom*sqrt(T_3_nom)/p_3_nom);
  G_TU_red = alpha*sqrt((1/(beta_TU^2) - 1)/(1/(beta_TU_nom^2) - 1));
 eta_is_TU_red = (1 - t*(1 - n_TU_red)^2)*(n_TU_red/G_TU_red)*(2 - ((n_TU_red/G_TU_red)));
//other turbine equations
  beta_TU = p_3/p_4;
  eta_is_TU = (h_3 - h_4)/(h_3 - h_4_is);
  P_mech_TU = m_dot_WF*(h_4 - h_3);
  P_loss_irr_TU = T0*m_dot_WF*(s_4 - s_3);
    der(E_loss_irr_TU) = P_loss_irr_TU;
//-------------HEX 1 DISCHARGE//
//pressure loss
  delta_P_HEX1 = k_p*m_dot_WF^2;
  p_3 = p_3_a - delta_P_HEX1;
//off-design
  UA_HEX1/UA_HEX1_nom = ((m_dot_WF^0.8*m_dot_solsalt_HEX1^0.8)/(m_dot_WF_nom^0.8*m_dot_solsalt_HEX1_nom^0.8))*((m_dot_WF_nom^0.8 + m_dot_solsalt_HEX1_nom^0.8)/(m_dot_WF^0.8 + m_dot_solsalt_HEX1^0.8));
//discharge: inlet cold side= state_3_a inlet hot side=solsalt_tank1_state
  outlet_hotside_guess_HEX1 = HotTESLiquid.setState_pT(p_tank2_nom, T_3_a + 9);
//NaK
  outlet_coldside_guess_HEX1 = WorkingFluid.setState_pT(p_3, solsalt_tank1_state.T - 9);
//WF
  cp_hot_ave_HEX1 = (solsalt_tank1_state.cp + outlet_hotside_guess_HEX1.cp)/2;
  cp_cold_ave_HEX1 = (WorkingFluid.specificHeatCapacityCp(state_3_a) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX1))/2;
  C_cold_HEX1 = m_dot_WF*cp_cold_ave_HEX1;
  C_hot_HEX1 = C_cold_HEX1;
//balanced operation
//C_hot_HEX1 = m_dot_solsalt_HEX1*cp_hot_ave_HEX1;
//mass balance
  C_min_HEX1 = min(C_cold_HEX1, C_hot_HEX1);
  C_max_HEX1 = max(C_cold_HEX1, C_hot_HEX1);
  C_r_HEX1 = C_min_HEX1/C_max_HEX1;
  NTU_HEX1 = UA_HEX1/C_min_HEX1;
  if C_r_HEX1 == 1 then
    eff_HEX1 = NTU_HEX1/(1 + NTU_HEX1);
  else
    eff_HEX1 = (1 - exp(-NTU_HEX1*(1 - C_r_HEX1)))/(1 - C_r_HEX1*exp(-NTU_HEX1*(1 - C_r_HEX1)));
  end if;
  Q_dot_HEX1 = eff_HEX1*Q_dot_max_HEX1;
  Q_dot_max_HEX1 = C_min_HEX1*(T_tank1 - T_3_a);
//hot side energy balance
  Q_dot_HEX1 = m_dot_solsalt_HEX1*(solsalt_tank1.h - h_hot_out_HEX1);
// cold side energy balance
  Q_dot_HEX1 = m_dot_WF*(h_cold_out_HEX1 - state_3_a.h);
//set outlets
  outlet_hotside_HEX1 = HotTESLiquid.setState_ph(p_tank2_nom, h_hot_out_HEX1);
//NaK
  outlet_coldside_HEX1 = WorkingFluid.setState_ph(p_3, h_cold_out_HEX1);
//WF
  T_3 = outlet_coldside_HEX1.T;
//irrev
  P_loss_irr_HEX1 = T0*((m_dot_WF*(s_3 - s_3_a)) + (m_dot_solsalt_HEX1*(outlet_hotside_HEX1.s - solsalt_tank1.s)));
  der(E_loss_irr_HEX1) = P_loss_irr_HEX1;
//-------------HEX2 (Recuperation DISCHARGE)//
//off-design
  UA_HEX2/UA_HEX2_nom = (m_dot_WF^0.8)/(m_dot_WF_nom^0.8);
//pressure loss
  delta_P_HEX2 = k_p*m_dot_WF^2;
  p_3_a = p_2 - delta_P_HEX2;
  p_4_a = p_4 - delta_P_HEX2;
//cold side
  outlet_coldside_guess_HEX2 = WorkingFluid.setState_pT(p_3_a, T_4_guess);
  cp_cold_ave_HEX2 = (WorkingFluid.specificHeatCapacityCp(state_2) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX2))/2;
// hot side
  outlet_hotside_guess_HEX2 = WorkingFluid.setState_pT(p_4_a, T_2);
  cp_hot_ave_HEX2 = (WorkingFluid.specificHeatCapacityCp(state_4_guess) + WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX2))/2;
//variables for effectiveness
  C_cold_HEX2 = m_dot_WF*cp_cold_ave_HEX2;
  C_hot_HEX2 = m_dot_WF*cp_hot_ave_HEX2;
  C_min_HEX2 = min(C_cold_HEX2, C_hot_HEX2);
  C_max_HEX2 = max(C_cold_HEX2, C_hot_HEX2);
  C_r_HEX2 = C_min_HEX2/C_max_HEX2;
  NTU_HEX2 = UA_HEX2/C_min_HEX2;
  if C_r_HEX2 < 1 then
    eff_HEX2 = (1 - exp(-NTU_HEX2*(1 - C_r_HEX2)))/(1 - C_r_HEX2*exp(-NTU_HEX2*(1 - C_r_HEX2)));
  else
    eff_HEX2 = NTU_HEX2/(1 + NTU_HEX2);
  end if;
//heat transferred
  Q_dot_HEX2 = eff_HEX2*Q_dot_max_HEX2;
  Q_dot_max_HEX2 = C_min_HEX2*(T_4_guess - T_2);
//T_guess instead of T_4
//energy balance and outlet states
// hot side energy balance
  Q_dot_HEX2 = m_dot_WF*(state_4_guess.h - h_hot_out_HEX2);
//cold side energy balance
  Q_dot_HEX2 = m_dot_WF*(h_cold_out_HEX2 - state_2.h);
//set outlets
  outlet_hotside_HEX2 = WorkingFluid.setState_ph(p_4_a, h_hot_out_HEX2);
//WF
  outlet_coldside_HEX2 = WorkingFluid.setState_ph(p_3_a, h_cold_out_HEX2);
//WF
  T_3_a = outlet_coldside_HEX2.T;
  T_4_a = outlet_hotside_HEX2.T;
//irrev
  P_loss_irr_HEX2 = T0*((m_dot_WF*(s_4_a - s_4)) + (m_dot_WF*(s_3_a - s_2)));
  der(E_loss_irr_HEX2) = P_loss_irr_HEX2;
//-------------HEX 3 DISCHARGE//
//pressure loss
  p_1 = p_1_a - delta_P_HEX3;
  delta_P_HEX3 = k_p*m_dot_WF^2;
//discharge: inlet cold side=tank 4  inlet hot side=state 1a
  outlet_hotside_guess_HEX3 = WorkingFluid.setState_pT(p_1, T_tank4 + 7);
//
  outlet_coldside_guess_HEX3 = ColdTESLiquid.setState_pT(p_tank3_nom, T_1_a - 9);
//
//pinch point guess
  cp_cold_ave_HEX3 = (coldliq_tank4_state.cp + ColdTESLiquid.specificHeatCapacityCp(outlet_coldside_guess_HEX3))/2;
  cp_hot_ave_HEX3 = (WorkingFluid.specificHeatCapacityCp(state_1_a) + WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX3))/2;
//heat capacity rates
  C_hot_HEX3 = m_dot_WF*cp_hot_ave_HEX3;
  C_cold_HEX3 = C_hot_HEX3;
//balanced operation
  C_cold_HEX3 = m_dot_methanol_HEX3*cp_cold_ave_HEX3;
  C_min_HEX3 = min(C_cold_HEX3, C_hot_HEX3);
  C_max_HEX3 = max(C_cold_HEX3, C_hot_HEX3);
  C_r_HEX3 = C_min_HEX3/C_max_HEX3;
//off-design
  UA_HEX3/UA_HEX3_nom = ((m_dot_WF^0.8*m_dot_methanol_HEX3^0.8)/(m_dot_WF_nom^0.8*m_dot_methanol_HEX3_nom^0.8))*((m_dot_WF_nom^0.8 + m_dot_methanol_HEX3_nom^0.8)/(m_dot_WF^0.8 + m_dot_methanol_HEX3^0.8));
  NTU_HEX3 = UA_HEX3/C_min_HEX3;
//eff_HEX3_charge = NTU_HEX3_charge/(1 + NTU_HEX3_charge);
//this eff_HEX3_charge formulation is computationally expensive
  if C_r_HEX3 < 1 then
    eff_HEX3 = (1 - exp(-NTU_HEX3*(1 - C_r_HEX3)))/(1 - C_r_HEX3*exp(-NTU_HEX3*(1 - C_r_HEX3)));
  else
    eff_HEX3 = NTU_HEX3/(1 + NTU_HEX3);
  end if;
  Q_dot_HEX3 = eff_HEX3*Q_dot_max_HEX3;
  Q_dot_max_HEX3 = C_min_HEX3*(T_1_a - T_tank4);
// cold side energy balance
  Q_dot_HEX3 = m_dot_methanol_HEX3*(h_cold_out_HEX3 - coldliq_tank4_state.h);
//hot side energy balance
  Q_dot_HEX3 = m_dot_WF*(state_1_a.h - h_hot_out_HEX3);
//set outlets
  outlet_hotside_HEX3 = WorkingFluid.setState_ph(p_1, h_hot_out_HEX3);
  outlet_coldside_HEX3 = ColdTESLiquid.setState_ph(p_tank3_nom, h_cold_out_HEX3);
  T_1 = outlet_hotside_HEX3.T;
//irrev
  P_loss_irr_HEX3 = T0*((m_dot_WF*(s_1 - s_1_a))) + T0*(m_dot_methanol_HEX3*(outlet_coldside_HEX3.s - coldliq_tank4_state.s));
 der(E_loss_irr_HEX3) = P_loss_irr_HEX3;
//-------------HEX 4 rejection DISCHARGE//
  p_4_a = p_1_a;
  Q_dot_HEXrej = (h_4_a - h_1_a)*m_dot_WF;
  der(Q_HEXrej) = Q_dot_HEXrej;
//irrev
  P_loss_irr_HEXrej = T0*((m_dot_WF*(s_1_a - s_4_a - ((h_1_a - h_4_a)/T0))));
  der(E_loss_irr_HEXrej) = P_loss_irr_HEXrej;
//state 1_a DISCHARGE
  state_1_a = WorkingFluid.setState_pT(p_1_a, T_1_a);
  h_1_a = WorkingFluid.specificEnthalpy(state_1_a);
  s_1_a = WorkingFluid.specificEntropy(state_1_a);
//state 1 DISCHARGE guess
  T_1 = T1_guess_control_discharge.u;
  T_1_guess = T1_guess_control_discharge.y;
//state 1 DISCHARGE
  state_1 = WorkingFluid.setState_pT(p_1, T_1);
  h_1 = WorkingFluid.specificEnthalpy(state_1);
  s_1 = WorkingFluid.specificEntropy(state_1);
//state 2 isentropic DISCHARGE
  s_2_is = s_1 "isentropic compressor outlet spec. entropy";
  state_2_is = WorkingFluid.setState_ps(p_2, s_2_is) "isentropic state of compressor outlet";
  T_2_is = state_2_is.T;
  h_2_is = WorkingFluid.specificEnthalpy(state_2_is);
//state 2 DISCHARGE
  state_2 = WorkingFluid.setState_ph(p_2, h_2);
  T_2 = state_2.T;
  s_2 = WorkingFluid.specificEntropy(state_2);
//state 3a DISCHARGE
  state_3_a = WorkingFluid.setState_pT(p_3_a, T_3_a);
  h_3_a = WorkingFluid.specificEnthalpy(state_3_a);
  s_3_a = WorkingFluid.specificEntropy(state_3_a);
//state 3 DISCHARGE
  state_3 = WorkingFluid.setState_pT(p_3, T_3);
  h_3 = state_3.h;
  s_3 = WorkingFluid.specificEntropy(state_3);
//state 3 DISCHARGE guess
  T_3 = T3_guess_control_discharge.u;
  T_3_guess = T3_guess_control_discharge.y;
//state 4 isentropic
  s_4_is = s_3;
  state_4_is = WorkingFluid.setState_ps(p_4, s_4_is) "isentropic state of turbine outlet";
  T_4_is = state_4_is.T "isentropic turbine outlet temperature";
  h_4_is = state_4_is.h;
//state 4
  state_4 = WorkingFluid.setState_ph(p_4, h_4);
  T_4 = state_4.T " turbine outlet spec. entropy";
  s_4 = WorkingFluid.specificEntropy(state_4);
//state 4 guess DISCHARGE
  T_4 = T4_guess_control_discharge.u;
  T_4_guess = T4_guess_control_discharge.y;
  state_4_guess = WorkingFluid.setState_pT(p_4, T_4_guess);
  h_4_guess = WorkingFluid.specificEnthalpy(state_4_guess);
  s_4_guess = WorkingFluid.specificEntropy(state_4_guess);
//state 4_a DISCHARGE
  state_4_a = WorkingFluid.setState_pT(p_4_a, T_4_a);
  h_4_a = WorkingFluid.specificEnthalpy(state_4_a);
  s_4_a = WorkingFluid.specificEntropy(state_4_a);
//------------------ELECTRICAL MACHINERY
  S_GC = sqrt(P_GC^2 + Q_GC_set^2);
  cos_phi_GC = P_GC/S_GC;
  Q_TR_HV = Q_TR_set;
  P_TR_HV=P_GC; 
//------------------TRANSFORMER
//connection point to grid
  phi_U_HV = phi_TR_set;
  U_TR_HV_phase = U_TR_set/sqrt(3) "phase voltage of reference phase";
  U_TR_HV = Complex(U_TR_HV_phase*cos(phi_U_HV), U_TR_HV_phase*sin(phi_U_HV));
  phi_U_HV_min_I_HV = phi_I_HV - phi_U_HV;
  pf_TR_HV = cos(ComplexMath.arg(Complex(P_TR_HV, Q_TR_HV)));
  phi_U_HV_min_I_HV = ComplexMath.arg(Complex(P_TR_HV, Q_TR_HV));
//transformer parameters
  Z_TR_HV = (((u_TR_SC/100)*U_TR_HV_nom^2)/S_TR_nom)/2;
//OedingOswald eq 8.3
  R_TR_HV = ((u_TR_ohmic_OC/100)*U_TR_HV_nom^2/S_TR_nom)/2;
//OedingOswald eq 8.4
  X_TR_HV = sqrt(Z_TR_HV^2 - R_TR_HV^2);
  R_TR_LV_transferred = R_TR_HV;
  X_TR_LV_transferred = X_TR_HV;
  Z_TR_OC_HV = (100/(i_TR_OC))*U_TR_HV_nom^2/S_TR_nom;
//OedingOswald eq 8.1
  R_TR_FE_HV = U_TR_HV_nom^2/(P_TR_loss_OC);
//OedingOswald eq 8.2a
  X_TR_h_HV = (R_TR_FE_HV*Z_TR_OC_HV)/sqrt(R_TR_FE_HV^2 - Z_TR_OC_HV^2);
//OedingOswald eq 8.2b
//high voltage
  S_TR_HV = Complex(P_TR_HV, Q_TR_HV);
  S_TR_HV = 3*U_TR_HV*ComplexMath.conj(I_TR_HV);
  U_TR_HV = Complex(R_TR_HV, X_TR_HV)*I_TR_HV + U_TR_h_HV;
  U_TR_h_HV = ((Complex(R_TR_FE_HV, 0)*Complex(0, X_TR_h_HV))/(Complex(R_TR_FE_HV, 0) + Complex(0, X_TR_h_HV)))*I_mag;
  Z_TR_h = ((Complex(R_TR_FE_HV, 0)*Complex(0, X_TR_h_HV))/(Complex(R_TR_FE_HV, 0) + Complex(0, X_TR_h_HV)));
  I_TR_HV_abs = ComplexMath.abs(I_TR_HV);
//low voltage
  R_TR_LV_transferred = a^2*R_TR_LV;
  X_TR_LV_transferred = a^2*X_TR_LV;
  U_TR_LV_transferred = Complex(R_TR_LV_transferred, X_TR_LV_transferred)*I_TR_LV_transferred + U_TR_h_HV;
  I_TR_LV_transferred = I_TR_LV/a;
  U_TR_LV_transferred = a*U_TR_LV;
  S_TR_LV = 3*U_TR_LV*ComplexMath.conj(I_TR_LV);
  S_delta = S_TR_HV + S_TR_LV;
  S_delta.re = P_TR_loss;
  S_delta.im = Q_req;
  S_TR_LV.re = P_TR_LV;
  S_TR_LV.im = Q_TR_LV;
  S_TR_HV_abs = ComplexMath.abs(S_TR_HV);
  S_TR_LV_abs = ComplexMath.abs(S_TR_LV);
//losses
  P_TR_ohmic_loss_LV = 3*ComplexMath.abs(I_TR_LV)^2*R_TR_LV;
  P_TR_ohmic_loss_HV = 3*ComplexMath.abs(I_TR_HV)^2*R_TR_HV;
  P_TR_iron_loss = P_TR_loss - P_TR_ohmic_loss_LV - P_TR_ohmic_loss_HV;
  der(E_TR_loss) = P_TR_loss;
   //------------------SYNCHRONOUS MACHINE
//general
  omega_SM = 2*pi*f;
  N_p*n_RO = f*60;
//equivalent circuit elements
  X_h = omega_SM*L_h;
  X_d = omega_SM*L_d;
  X_sigma_ST = omega_SM*L_sigma_ST;
//Trafo side Interfaces
  0 = Q_TR_LV + Q_SM_ST;
  U_SM_ST = U_TR_LV;
  0 = P_TR_LV + P_SM_ST;
  S_SM_ST = Complex(P_SM_ST, Q_SM_ST);
    S_SM_ST_abs=ComplexMath.abs(S_SM_ST);
//airgap
  U_SM_ST = Complex(R_SM_ST, X_sigma_ST)*I_SM_ST + U_SM_h;
  S_SM_h = 3*U_SM_h*ComplexMath.conj(I_SM_ST);
  P_airgap = S_SM_h.re;
 P_Ohmic_Stator=P_SM_ST-P_airgap;
//excitation
  U_P = Complex(0, 1)*X_h*I_FD_ref;
//
  ComplexMath.abs(I_FD_ref) = I_FD*2*N_p*N_FD/(3*N_ST*0.9166);
//Binder p.5/41
  U_FD = R_FD*I_FD;
  P_excitation = U_FD*I_FD;
//losses
  v_RO = d_RO*pi*n_RO/60;
  P_windage_ventilation = k_schuisky*d_RO*(l_FE + 0.6*tau_p)*v_RO^2;
  P_FE_ST_teeth = k_Fe_stator_teeth*P10*B_ST_teeth^2*m_ST_FE_teeth;
  P_FE_ST_yoke = k_Fe_stator_yoke*P10*B_ST_yoke^2*m_ST_FE_yoke;
  P_FE = (f/50)^1.3*(P_FE_ST_teeth + P_FE_ST_yoke);
//energy
  der(E_SM_loss) = P_SM_loss;
//--------------------------COMPONENT MANAGEMENT--------------------------//
//MODE 1 CHARGE
  if Mode == 1 then
//tanks
    if SOC_tank2 < 0 then
      terminate("Minimum fill level of tank 2 reached");
    end if;
    if SOC_tank3 < 0 then
      terminate("Minimum fill level of tank 3 reached");
    end if;
    if SOC_tank1 > 1 then
      terminate("Maximum fill level of tank 1 reached");
    end if;
    if SOC_tank4 > 1 then
      terminate("Maximum fill level of tank 4 reached");
    end if;
    if T_tank1 > 873.150000 - 1 then
      terminate("Temperature in tank 1 too high");
    end if;
    if T_tank2 < 533.150000 + 1 then
      terminate("Temperature in tank 2 too low");
    end if;
//turbomachinery
    if beta_CO_red_charge > beta_CO_red_charge_max then
      terminate("surge line reached");
    end if;
    if beta_CO_red_charge < beta_CO_red_charge_min then
      terminate("choke line reached");
    end if;
        if beta_TU_red_charge > beta_TU_red_charge_max then
      terminate("maximum red pressure ratio reached");
    end if;
    if beta_TU_red_charge < beta_TU_red_charge_min then
      terminate("minimum red pressure ratio  reached");
    end if;
    if n_CO_red_charge < 0.8 then
      terminate("compressor reduced relative speed too low");
    end if;
    if n_CO_red_charge > 1.05 then
      terminate("compressor reduced relative speed too high");
      end if;
      //electrical machines      
    if S_TR_HV_abs > S_TR_nom then
    terminate("maximum transformer apparent power reached");
    end if;
   if S_SM_ST_abs > S_SM_nom then
    terminate("maximum synchronous machine apparent power reached");
    end if;  
//MODE 2 DISCHARGE
  elseif Mode == 2 then
//tanks
    if SOC_tank2 > 1 then
      terminate("Maximum fill level of tank 2 reached");
    end if;
    if SOC_tank3 > 1 then
      terminate("Maximum fill level of tank 3 reached");
    end if;
    if SOC_tank1 < 0 then
      terminate("Minimum fill level of tank 1 reached");
    end if;
    if SOC_tank4 < 0 then
      terminate("Minimum fill level of tank 4 reached");
    end if;
    if T_tank1 > 873.150000 - 1 then
      terminate("Temperature in tank 1 too high");
    end if;
    if T_tank2 < 533.150000 + 1 then
      terminate("Temperature in tank 2 too low");
    end if;
//turbomachinery
    if beta_CO_red > beta_CO_red_max then
      terminate("surge line reached");
    end if;
    if beta_CO_red < beta_CO_red_min then
      terminate("choke line reached");
    end if;
        if beta_TU_red > beta_TU_red_max then
      terminate("maximum red pressure ratio reached");
    end if;
    if beta_TU_red < beta_TU_red_min then
      terminate("minimum red pressure ratio  reached");
    end if;   
    if n_CO_red < 0.8 then
      terminate("compressor reduced relative speed too low");
    end if;
    if n_CO_red > 1.05 then
      terminate("compressor reduced relative speed too high");
    end if;
    //electrical machines    
    if S_TR_HV_abs > S_TR_nom then
    terminate("maximum transformer apparent power reached");
    end if;
   if S_SM_ST_abs > S_SM_nom then
    terminate("maximum synchronous machine apparent power reached");
    end if;
//MODE 0 HOLD
  else
    if T_tank1 > 873.150000 - 1 then
      terminate("Temperature in tank 1 too high");
    end if;
    if T_tank2 < 533.150000 + 1 then
      terminate("Temperature in tank 2 too low");
    end if;
  end if;
  annotation(
    Documentation(info = "<html><head></head><body>Dynamic Malta Charge &amp; discharge<div>Heat loss active</div><div><br></div><div><br></div></body></html>"));
end DynamicMalta_charge_discharge_exemplary_charge_hold_discharge;
