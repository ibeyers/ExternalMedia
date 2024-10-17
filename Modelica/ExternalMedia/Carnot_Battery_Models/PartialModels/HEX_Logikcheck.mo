within ExternalMedia.Carnot_Battery_Models.PartialModels;

model HEX_Logikcheck
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
  parameter Integer Mode=2;                
  parameter Real SOC_tank1_start = 1;
  parameter SI.Temperature T_tank1_start = from_degC(565);
  parameter Real SOC_tank2_start = 0;
  parameter SI.Temperature T_tank2_start = from_degC(279);
  
    //--------------------------PARAMETERS & VARIABLES SYSTEM-----------------------------//
  parameter SI.Temperature T0 = T_amb;
  parameter SI.Temperature T_amb = from_degC(25);
  WorkingFluid.ThermodynamicState state_amb_air(p(start = 101315), T(start = from_degC(25))) "thermodynamic state of rejec outlet";      
  //-------------Discharge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  //actual
  parameter Real relative_mass_flow=1;
  SI.MassFlowRate m_dot_WF(start = 762) "mass flow rate";
  parameter Real f_p = 0.01625 "pressure_loss_factor percent";
  parameter Real k_p = 0.0062084 "pressure_loss_factor";  
  //--------------------------PARAMETERS & VARIABLES TANKS-----------------------------//
  parameter SI.Mass m_working_solar_salt = 19386000;
  parameter SI.Mass m_working_methanol = 9486000;
  //geometry for both solar salt tanks
  parameter SI.Diameter D_tank_solsalt = 36.104957;
  parameter SI.Height h_tank_solsalt = 12.03498;
  parameter SI.Area A_cross_tank_solsalt = pi*(D_tank_solsalt/2)^2 "Cross-section Tank";
  parameter SI.Volume V_tank_solsalt = A_cross_tank_solsalt*h_tank_solsalt "Volume Tank";
  parameter SI.Thickness d_insulation_tank_solsalt = 0.4 "thickness of insulation wall layer";
  parameter SI.Radius r_tank_solsalt = D_tank_solsalt/2 "radius at start of insulation";
  parameter SI.Radius r_outer = D_tank_solsalt/2 + d_insulation_tank_solsalt "outer radius";
  parameter SI.Area A_W_tank_solsalt = (pi*D_tank_solsalt*h_tank_solsalt) + A_cross_tank_solsalt*2 "Total Surface Wall";
  parameter SI.Height x_tank_solsalt_min = 0.4;                
  //--------------------------TANK 1
  //nominal
  parameter SI.Temperature T_tank1_nom = from_degC(565.6);
  parameter SI.Pressure p_tank1_nom = 101325 "unpressurized tank";
  parameter SI.Density rho_tank1_nom = 1730.66;
  parameter SI.Mass m_tank1_min =  x_tank_solsalt_min*A_cross_tank_solsalt*rho_tank1_nom;
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
  parameter SI.Mass m_tank2_min =  x_tank_solsalt_min*A_cross_tank_solsalt*rho_tank2_nom;
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
  //-------------HEX 1//
  //design
  parameter Real UA_HEX1_nom = 25906701;
  SI.MassFlowRate m_dot_solsalt_HEX1_nom = 540;  
  //actual
  Real UA_HEX1;
  SI.Pressure delta_P_HEX1;
  //outlet guess states
  HotTESLiquid.ThermodynamicState outlet_hotside_guess_HEX1(p(start = 1001315), T(start = 271.5)) "Medium properties ";
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
//  SI.Power P_loss_irr_HEX1(displayUnit = "MW");
 // SI.Energy E_loss_irr_HEX1(displayUnit = "MWh", start = 0, fixed = true);
  //STATE 3 discharge
  SI.Pressure p_3(start = 572045) "Pressure after recup";
  //STATE 3a discharge
  SI.Pressure p_3_a= 581494 "Pressure after recup";
  SI.Temperature T_3_a= from_degC(260) "outlet temperature after recuperation";
  WorkingFluid.ThermodynamicState state_3_a "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_a(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a "compressor outlet spec. entropy";  
 initial equation          
  m_tank1 = m_tank1_start;
  T_tank1 = T_tank1_start;
  m_tank2 = m_tank2_start;
  T_tank2 = T_tank2_start;   
                 
equation
//--------------------------EQUATIONS SYSTEM-----------------------------//
  state_amb_air = WorkingFluid.setState_pT(101315, T_amb);
  m_dot_WF=relative_mass_flow*m_dot_WF_nom;
m_dot_solsalt_HEX1=relative_mass_flow*m_dot_solsalt_HEX1_nom;
//-------------TANK 1//
//nominal
  solsalt_tank1_nom = HotTESLiquid.setState_pT(p_tank1_nom, T_tank1_nom);
// connectors
  h_in_tank1 = 0;
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank1 = 0;
    m_in_tank1 = 0;
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
  SOC_tank1 = (m_tank1 - m_tank1_min)/m_working_solar_salt ;
  m_tank1 = solsalt_tank1.d*A_cross_tank_solsalt*x_tank1;
  int_energy_tank1 = m_tank1*solsalt_tank1.u;
  //exergy_tank1 = m_tank1*solsalt_tank1_state.cp*(T_tank1 - T_amb) - T_amb*m_tank1*solsalt_tank1_state.cp*log(T_tank1/T_amb);
  exergy_tank1 = m_tank1*solsalt_tank1.u - T_amb*m_tank1*solsalt_tank1_state.s;
//thermodynamic states
  solsalt_tank1_state = HotTESLiquid.setState_pT(p_tank1_nom, T_tank1);
//solar salt properties
  solsalt_tank1.p = p_tank1_nom;
  solsalt_tank1.T = T_tank1;
//losses
  //Q_div_A_tank1 = (0.00017*(T_tank1 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank1 = -Q_div_A_tank1*A_W_tank_solsalt;
  Q_dot_to_amb_tank1 = 0;
//-------------TANK 2//
//nominal
  solsalt_tank2_nom = HotTESLiquid.setState_pT(p_tank2_nom, T_tank2_nom);
// connectors
  h_in_tank2 = outlet_hotside_HEX1.h;
//MODE 1 CHARGE
  if Mode == 1 then
    m_out_tank2 = 0;
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
  //exergy_tank2 = m_tank2*solsalt_tank2_state.cp*(T_tank2 - T_amb) - T_amb*m_tank2*solsalt_tank2_state.cp*log(T_tank2/T_amb);
  exergy_tank2 = m_tank2*solsalt_tank2.u - T_amb*m_tank2*solsalt_tank2_state.s;
//thermodynamic states
  solsalt_tank2_state = HotTESLiquid.setState_pT(p_tank2_nom, T_tank2);
//solar salt properties
  solsalt_tank2.p = p_tank2_nom;
  solsalt_tank2.T = T_tank2;
//losses
  //Q_div_A_tank2 = (0.00017*(T_tank2 - 273.15) + 0.012)*1000;
  Q_dot_to_amb_tank2 = -Q_div_A_tank2*A_W_tank_solsalt;
  Q_dot_to_amb_tank2 = 0;
//-------------HEX 1 DISCHARGE//
//pressure loss
  delta_P_HEX1 = k_p*m_dot_WF^2;
  p_3 = p_3_a - delta_P_HEX1;
//off-design
  UA_HEX1/UA_HEX1_nom = ((m_dot_WF^0.8*m_dot_solsalt_HEX1^0.8)/(m_dot_WF_nom^0.8*m_dot_solsalt_HEX1_nom^0.8))*((m_dot_WF_nom^0.8 + m_dot_solsalt_HEX1_nom^0.8)/(m_dot_WF^0.8 + m_dot_solsalt_HEX1^0.8));
//discharge: inlet cold side= state_3_a inlet hot side=solsalt_tank1_state
    outlet_hotside_guess_HEX1 = HotTESLiquid.setState_pT(p_tank2_nom, T_tank2_nom);
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
 // T_3 = outlet_coldside_HEX1.T;
//irrev
// P_loss_irr_HEX1 = T0*((m_dot_WF*(s_3 - s_3_a)) + (m_dot_solsalt_HEX1*(outlet_hotside_HEX1.s - solsalt_tank1.s)));
 //P_loss_irr_HEX1 = T0*((m_dot_WF*(s_3 - s_3_a)) + (m_dot_solsalt_HEX1*((outlet_hotside_HEX1.cp*log(outlet_hotside_HEX1.T/T0)) - (solsalt_tank1.state.cp*log(solsalt_tank1.T/T0)))));
 //P_loss_irr_HEX1 = T0*((m_dot_WF*(s_3 - s_3_a)) + (m_dot_solsalt_HEX1*cp_hot_ave_HEX1*((log(outlet_hotside_HEX1.T/solsalt_tank1.T)))));
 // der(E_loss_irr_HEX1) = P_loss_irr_HEX1;
  
  //state 3a DISCHARGE
  state_3_a = WorkingFluid.setState_pT(p_3_a, T_3_a);
  h_3_a = WorkingFluid.specificEnthalpy(state_3_a);
  s_3_a = WorkingFluid.specificEntropy(state_3_a);

end HEX_Logikcheck;
