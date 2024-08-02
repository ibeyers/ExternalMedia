within ExternalMedia.Carnot_Battery_Models.StaticModels;

model StaticMalta_charge_discharge
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
  /*
    package WorkingFluid = ExternalMedia.Media.CoolPropMedium(mediumName = "Nitrogen", substanceNames = {"N2"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph); 
    */
  //--------------------------PARAMETERS & VARIABLES CYCLE-----------------------------//
  parameter SI.Temperature T0 = 273.15;  
  //-------------Charge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom_charge = 766 "design mass flow rate";
  //actual
  parameter SI.MassFlowRate m_dot_WF_charge = 766 "mass flow rate";
  parameter Real f_p_charge = 0.01075 "pressure_loss_factor percent";
  parameter SI.Pressure p_fix_charge = 100000 "fixed pressure point through expansion vessel at p_4";
  //-------------Discharge//
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  //actual
  parameter SI.MassFlowRate m_dot_WF = 762 "mass flow rate";
  parameter Real f_p = 0.01625 "pressure_loss_factor percent";
  parameter SI.Pressure p_fix = 100000 "fixed pressure point through expansion vessel";
  
  //--------------------------PARAMETERS & VARIABLES TANKS-----------------------------//
  //--------------------------TANK 1
  //nominal
  parameter SI.Temperature T_nom_tank1 = from_degC(565);
  //parameter SI.Temperature T_tank1 = from_degC(565);
  SI.Temperature T_tank1;
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

  //--------------------------PARAMETERS & VARIABLES CHARGE-----------------------------//

  //-------------HEX 1 CHARGE//
  //design
  parameter Real UA_HEX1_nom = 25906701;
  SI.MassFlowRate m_dot_solsalt_HEX1_nom = 553.021;
  //actual
  Real UA_HEX1_charge;
  SI.MassFlowRate m_dot_solsalt_HEX1_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX1_charge(p(start = 456931), T(start = T_nom_tank2), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
  HotTESLiquid.ThermodynamicState outlet_coldside_guess_HEX1_charge(p(start = 101325), T(start = T_nom_tank1), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
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
  WorkingFluid.ThermodynamicState outlet_hotside_HEX1_charge(p(start = 456931), T(start = T_nom_tank2), phase(start = 1)) "Medium properties of HEX port at interface to tank 2";
  HotTESLiquid.ThermodynamicState outlet_coldside_HEX1_charge(p(start = 101325), T(start = T_nom_tank1), phase(start = 1)) "Medium properties of HEX port at interface to tank 1";
  SI.Power P_loss_irr_HEX1_charge(displayUnit = "MW");
  //-------------HEX 2 CHARGE//
  //design
  parameter Real UA_HEX2_nom = 13279998;
  //actual
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
  parameter SI.MassFlowRate m_dot_methanol_HEX3_nom = 311.770;
  //actual
  Real UA_HEX3_charge;
  SI.MassFlowRate m_dot_methanol_HEX3_charge "mass flow rate required for balanced HEX";
  //outlet guess states
  ColdTESLiquid.ThermodynamicState outlet_hotside_guess_HEX3_charge(p(start = 101325), T(start = T_nom_tank4), phase(start = 1)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX3_charge(p(start = 100000), T(start = T_nom_tank3), phase(start = 1)) "Medium properties ";
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
  //heat transferred
  SI.HeatFlowRate Q_dot_HEX3_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEX3_charge(displayUnit = "MW");
  //outlet
  SI.SpecificEnthalpy h_hot_out_HEX3_charge;
  SI.SpecificEnthalpy h_cold_out_HEX3_charge;
  ColdTESLiquid.ThermodynamicState outlet_hotside_HEX3_charge(p(start = 101325), T(start = T_nom_tank4), phase(start = 1)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_HEX3_charge(p(start = 100000), T(start = T_nom_tank3), phase(start = 1)) "Medium properties ";

  //-------------HEX rejection charge//
  parameter SI.Temperature T_amb = from_degC(20);
  parameter Real NTU_HEXrej_charge = 7.8 "Number of transfer units of rejection HEX";
  WorkingFluid.ThermodynamicState state_amb_air(p(start = 1001315), T(start = from_degC(20))) "thermodynamic state of rejec outlet";
  //cold-side
  Real C_cold_HEXrej_charge "Heat capacity rate of cold side of rejection HEX ";
  SI.MassFlowRate m_dot_rej_charge "required mass flow rate to cool down to fixed T_2";
  //hot-side
  Real C_hot_HEXrej_charge "Heat capacity rate of hot side of rejection HEX";
  //variables for effectiveness
  Real C_min_HEXrej_charge;
  Real C_max_HEXrej_charge;
  Real C_r_HEXrej_charge;
  Real eff_HEXrej_charge(start = 0.99) "effectiveness of rej HEX";
  //source Farres-Artunez Dissertation
  SI.HeatFlowRate Q_dot_max_HEXrej_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_HEXrej_charge(displayUnit = "MW");
  SI.SpecificEnthalpy h_hot_out_HEXrej_charge(start = 431395);
  WorkingFluid.ThermodynamicState outlet_hotside_HEXrej_charge(p(start = p_2_a_nom_charge), T(start = from_degC(20)));
  SI.Power P_loss_irr_HEXrej_charge(displayUnit = "MW");

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
  //other
  SI.Power P_mech_CO_charge(displayUnit = "MW");
  SI.Power P_loss_irr_CO_charge(displayUnit = "MW");

 
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
  parameter Real n_TU_charge = 3000 "actual speed";
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
 
  //-------------SYSTEM CHARGE//
  SI.Power P_mech_shaft(displayUnit = "MW");
  SI.HeatFlowRate Q_pump_charge(displayUnit = "MW");
  Real COP(start = 1); 
 
  //-------------STATES//
  //-------------Charge//
  //fixed temperature point charge
  SI.Temperature T_2_a_charge = from_degC(20) "outlet temperature after rejec";  
//fixed pressure point charge
  SI.Pressure p_4_charge = p_fix_charge "state 4 pressure";    
  //state 1 charge
  SI.Pressure p_1_charge(start = 102161) " pressure ";
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
  SI.Pressure p_4_a_charge(start = 100000) "state  pressure";
  SI.Temperature T_4_a_charge(start = from_degC(18.14)) "state  temperature";
  WorkingFluid.ThermodynamicState state_4_a_charge(p(start = 100000), T(start = from_degC(18.14))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_a_charge(start = 417519) "enthalpy";
  WorkingFluid.SpecificEntropy s_4_a_charge(start = 3837) "spec. entropy";
  //state 4 charge
  SI.Temperature T_4_charge(start = from_degC(267.533)) "state 4 temperature";
  WorkingFluid.ThermodynamicState state_4_charge(p(start = p_fix_charge), T(start = from_degC(267.533))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_4_charge(start = 671458) "enthalpy";
  WorkingFluid.SpecificEntropy s_4_charge(start = 4495) "spec. entropy";
  //STATE 3 isentropic charge
  SI.Temperature T_3_is_charge(start = from_degC(541.61)) "isentropic outlet temperature of compressor";
  WorkingFluid.ThermodynamicState state_3_is_charge(p(start = 449506), T(start = from_degC(541.61))) "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEntropy s_3_is_charge(start = 4495) "compressor outlet spec. entropy";
  WorkingFluid.SpecificEnthalpy h_3_is_charge(start = 964771);
  //STATE 3 charge
  SI.Temperature T_3_charge(start = from_degC(578.98)) "actual outlet temperature of compressor";
  SI.Pressure p_3_charge(start = 449506) "pressure coming out of compressor";
  WorkingFluid.ThermodynamicState state_3_charge(p(start = 449506), T(start = from_degC(578.98))) "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_charge(start = 1006139) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_charge(start = 4545) "compressor outlet spec. entropy";
  //STATE 3a charge
  WorkingFluid.ThermodynamicState state_3_a_charge(p(start = 456931), T(start = from_degC(281.29))) "thermodynamic state";
  SI.Temperature T_3_a_charge(start = from_degC(281.29)) "outlet temperature after HEX";
  SI.Pressure p_3_a_charge(start = 456931) "Pressure after HEX";
  WorkingFluid.SpecificEnthalpy h_3_a_charge(start = 685666) "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a_charge(start = 4079) "HEX outlet spec. entropy";
  //state 3a charge  guess
  SI.Temperature T_3_a_charge_guess(start = from_degC(281.29)) "state  temperature";
  SI.Pressure p_3_a_charge_guess(start = 456931) "state  pressure";
  WorkingFluid.ThermodynamicState state_3_a_charge_guess(p(start = 456931), T(start = from_degC(281.29))) "thermodynamic state";
  WorkingFluid.SpecificEnthalpy h_3_a_charge_guess(start = 1006139) "enthalpy";
  WorkingFluid.SpecificEntropy s_3_a_charge_guess(start = 4545) "spec. entropy";
  //state 2 charge
  SI.Pressure p_2_charge(start = 452019) "state  pressure";
  SI.Temperature T_2_charge(start = 309.3) "temperature ";
  WorkingFluid.ThermodynamicState state_2_charge(p(start = 452019), T(start = from_degC(32.66))) "thermodynamic state ";
  WorkingFluid.SpecificEnthalpy h_2_charge(start = 431395) "enthalpy";
  WorkingFluid.SpecificEntropy s_2_charge(start = 3474) "spec. entropy";
  //STATE 2a charge
  WorkingFluid.ThermodynamicState state_2_a_charge(p(start = 442461), T(start = from_degC(25))) "thermodynamic state of rejec outlet";
  WorkingFluid.SpecificEnthalpy h_2_a_charge(start = 423634) "rejec outlet enthalpy";
  WorkingFluid.SpecificEntropy s_2_a_charge(start = 3448) "rejec outlet spec. entropy";
  SI.Pressure p_2_a_charge(start = 452019) "Pressure after rejec";
  SI.Pressure p_2_a_charge_guess(start = 452019) "Pressure after rejec";
  
//--------------------------PARAMETERS & VARIABLES DISCHARGE-----------------------------//
  //-------------HEX 1//
  //actual
  Real UA_HEX1;
  //outlet guess states
  HotTESLiquid.ThermodynamicState outlet_hotside_guess_HEX1(p(start=1001315),T(start = 271.5)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX1(p(start=572045),T(start = 556.099)) "Medium properties ";
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

  //-------------HEX 2 DISCHARGE//
  //actual
  //SI.Pressure delta_HEX2 "pressure loss over HEX";
   Real UA_HEX2;   
//outlet guess states
   //cold-side
   WorkingFluid.ThermodynamicState outlet_coldside_guess_HEX2(p(start=581494),T(start = 269.477)) "Medium properties of HEX port at interface to tank 2";   
   SI.SpecificHeatCapacity cp_cold_ave_HEX2;   
   //hot-side
   WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX2(p(start=105039),T(start = 121.643)) "Medium properties of HEX port at interface to tank 2";   
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
   WorkingFluid.ThermodynamicState outlet_hotside_HEX2(p(start=581494),T(start = 269.477)) "Medium properties of HEX port at interface to tank 2";
   WorkingFluid.ThermodynamicState outlet_coldside_HEX2(p(start=105039),T(start = 121.643)) "Medium properties of HEX port at interface to tank 2";       
   
/*   

  SI.Power P_loss_irr_HEX2(displayUnit = "MW");
  */
  
  //-------------HEX 3 DISCHARGE//
  //actual
  //SI.Pressure delta_HEX3 "pressure loss over HEX";
  Real UA_HEX3;
  SI.MassFlowRate m_dot_methanol_HEX3 "mass flow rate required for balanced HEX";
  //outlet guess states
  ColdTESLiquid.ThermodynamicState outlet_coldside_guess_HEX3(p(start = 101325), T(start = T_nom_tank3)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_hotside_guess_HEX3(p(start = 100000), T(start = T_nom_tank4)) "Medium properties ";
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
  ColdTESLiquid.ThermodynamicState outlet_coldside_HEX3(p(start = 101325), T(start = T_nom_tank4)) "Medium properties ";
  WorkingFluid.ThermodynamicState outlet_hotside_HEX3(p(start = 100000), T(start = T_nom_tank3)) "Medium properties ";   
  SI.Power P_loss_irr_HEX3(displayUnit = "MW");

    //-------------COMPRESSOR DISCHARGE/
  Real c1 "factor 1";
  Real c2 "factor 2";
  Real c3 "factor 3";
  //design
  parameter Real beta_CO_nom = 5.911 "design compression ratio";
  parameter Real n_CO_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom = 0.88 "design isentropic efficiency";
  SI.Temperature T_1_nom = from_degC(-53.76) "state 1 temperature";
  SI.Pressure p_1_nom = p_fix "state 1 pressure";
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
  parameter Real n_TU = 3000 "actual speed";
  SI.Efficiency eta_is_TU "absolute isentropic turbine efficiency";
  //reduced
  Real beta_TU_red(start = 1) "reduced expansion ratio";
  Real n_TU_red(start = 1) "reduced speed";
  SI.Efficiency eta_is_TU_red "reduced isentropic turbine efficiency";
  Real G_TU_red(start = 1) "reduced mass flow rate turbine";
  //other
  SI.Power P_mech_TU(displayUnit = "MW");
  SI.Power P_loss_irr_TU(displayUnit = "MW");

//-------------STATES//
  //-------------Discharge//
  //fixed temperature point charge
  SI.Temperature T_1_a = from_degC(32) "outlet temperature after rejec";  
//fixed pressure point charge
  SI.Pressure p_1 = p_fix "state 1 pressure";      

  //STATE 1 a discharge
  SI.Pressure p_1_a  "pressure after Heat rejection ";
  WorkingFluid.ThermodynamicState state_1_a "thermodynamic state after Heat rejection ";
  WorkingFluid.SpecificEnthalpy h_1_a "turbine-side recuperation after Heat rejection ";
  WorkingFluid.SpecificEntropy s_1_a "turbine-side recuperation after Heat rejection ";     
//state 1 discharge
  SI.Temperature T_1(start = from_degC(-52.760)) "ambient temperature";
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
  SI.Temperature T_3_a(start= from_degC(260)) "outlet temperature after recuperation";   
  WorkingFluid.ThermodynamicState state_3_a "thermodynamic state of compressor outlet";
  WorkingFluid.SpecificEnthalpy h_3_a(start = 264244) "compressor outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3_a "compressor outlet spec. entropy";
    //STATE 3 discharge
   SI.Pressure p_3(start = 572045) "Pressure after recup";    
  SI.Temperature T_3(start = from_degC(560)) "outlet temperature after HEX";
  WorkingFluid.ThermodynamicState state_3 "thermodynamic state of HEX outlet";  
  WorkingFluid.SpecificEnthalpy h_3(start = 579480) "HEX outlet enthalpy";
  WorkingFluid.SpecificEntropy s_3 "HEX outlet spec. entropy";   
    //STATE 4 isentropic
  SI.Temperature T_4_is(start=from_degC(250)) "isentropic outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4_is "isentropic state of turbine outlet";
  WorkingFluid.SpecificEnthalpy h_4_is "turbine outlet enthalpy";  
  WorkingFluid.SpecificEntropy s_4_is "isentropic turbine outlet spec. entropy";
  //STATE 4 discharge
  SI.Pressure p_4(start=105039) " pressure at turb outlet";
  WorkingFluid.SpecificEnthalpy h_4 "turbine outlet enthalpy";
 SI.Temperature T_4(start=from_degC(270)) "outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4(p(start=105039), T(start=from_degC(270))) "thermodynamic state of turbine outlet";
  WorkingFluid.SpecificEntropy s_4 "turbine outlet spec. entropy";
    //STATE 4 guess discharge
  SI.Temperature T_4_guess(start=from_degC(270)) "outlet temperature of turbine";
  WorkingFluid.ThermodynamicState state_4_guess(p(start=105039), T(start=from_degC(270))) "thermodynamic state of turbine outlet";
  WorkingFluid.SpecificEnthalpy h_4_guess "turbine outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_guess "turbine outlet spec. entropy";
    //STATE 4 a discharge
  SI.Pressure p_4_a(start=103332) " pressure at recuperation outlet";
  SI.Temperature T_4_a(start=from_degC(270) ) "outlet temperature after recuperation";  
  WorkingFluid.ThermodynamicState state_4_a "thermodynamic state of turbine-side recuperation outlet";
  WorkingFluid.SpecificEnthalpy h_4_a "turbine-side recuperation outlet enthalpy";
  WorkingFluid.SpecificEntropy s_4_a "turbine-side recuperation outlet spec. entropy"; 
  

  Modelica.Blocks.Continuous.SecondOrder T3_guess_control(D = 0.4, w = 0.5) annotation(
    Placement(transformation(origin = {0, 66}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder p3_guess_control(w = 0.5, D = 0.4) annotation(
    Placement(transformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.SecondOrder p_2a_guess_control(w = 0.5, D = 0.4)  annotation(
    Placement(transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}})));
 Modelica.Blocks.Continuous.SecondOrder T4_guess_control_discharge(w = 0.5, D = 0.4)  annotation(
    Placement(transformation(origin = {4, -56}, extent = {{-10, -10}, {10, 10}})));
initial equation
//--------------------------INITIAL EQUATIONS-----------------------------//
  T_3_a_charge_guess = from_degC(281.29);
  p_3_a_charge_guess = 456931;
  p_2_a_charge_guess = 442461;
  T_4_guess=from_degC(270);
  
equation
//--------------------------EQUATIONS TANKS-----------------------------//
//-------------TANK 1//
  solsalt_tank1_nom = HotTESLiquid.setState_pT(p_nom_tank1, T_nom_tank1);
  solsalt_tank1_state = HotTESLiquid.setState_pT(p_nom_tank1, T_tank1);
//solar salt properties
  solsalt_tank1.p = p_nom_tank1;
  solsalt_tank1.T = T_tank1;
//T_tank1=from_degC(565)-0.0006944444444444445 * time;
  T_tank1 = from_degC(565);
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
  //--------------------------EQUATIONS CHARGE-----------------------------//
//-------------COMPRESSOR CHARGE//
//reduced values compressor
  eta_is_CO_red_charge = eta_is_CO_charge/eta_is_CO_nom_charge;
  beta_CO_red_charge = beta_CO_charge/beta_CO_nom_charge;
  G_CO_red_charge = m_dot_WF_charge*sqrt(T_4_charge)/p_4_charge/(m_dot_WF_nom_charge*sqrt(T_4_nom_charge)/p_4_nom_charge);
  eta_is_CO_red_charge = (1 - c4*(1 - n_CO_red_charge)^2)*(n_CO_red_charge/G_CO_red_charge)*(2 - n_CO_red_charge/G_CO_red_charge);
  n_CO_red_charge = n_CO_charge/sqrt(T_4_charge)/(n_CO_nom_charge/sqrt(T_4_nom_charge));
//other compressor equations
  beta_CO_charge = p_3_charge/p_4_charge;
  eta_is_CO_charge = (h_3_is_charge - h_4_charge)/(h_3_charge - h_4_charge);
  beta_CO_red_charge = c1_charge*G_CO_red_charge^2 + c2_charge*G_CO_red_charge + c3_charge;
  c1_charge = n_CO_red_charge/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c2_charge = (p - 2*m*n_CO_red_charge^2)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c3_charge = -1*(p*m*n_CO_red_charge - m^2*n_CO_red_charge^3)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  P_mech_CO_charge = m_dot_WF_charge*(h_3_charge - h_4_charge);
  P_loss_irr_CO_charge = T0*m_dot_WF_charge*(s_3_charge - s_4_charge);

//-------------EXPANDER CHARGE//
//parameters
  alpha_charge = sqrt(1.4 - 0.4*n_TU_red_charge);
//reduced values expander
  n_TU_red_charge = n_TU_charge/sqrt(T_2_a_charge)/(n_TU_nom_charge/sqrt(T_2_a_nom_charge));
  G_TU_red_charge = alpha_charge*sqrt(T_2_a_nom_charge/T_2_a_charge)*sqrt((beta_TU_charge^2 - 1)/(beta_TU_nom_charge^2 - 1));
  beta_TU_red_charge = beta_TU_charge/beta_TU_nom_charge;
  eta_is_TU_red_charge = (1 - t*(1 - n_TU_red_charge)^2)*(n_TU_red_charge/G_TU_red_charge)*(2 - ((n_TU_red_charge/G_TU_red_charge)));
  eta_is_TU_red_charge = eta_is_TU_charge/eta_is_TU_nom_charge;
//other turbine equations
  beta_TU_charge = p_2_a_charge_guess/p_1_charge;
  eta_is_TU_charge = (h_2_a_charge - h_1_charge)/(h_2_a_charge - h_1_is_charge);
  P_mech_TU_charge = m_dot_WF_charge*(h_1_charge - h_2_a_charge);
  P_loss_irr_TU_charge = T0*m_dot_WF_charge*(s_1_charge - s_2_a_charge);

//-------------HEX 1 CHARGE//
//pressure loss
  p_3_charge = p_3_a_charge*(1 + f_p_charge);
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
  outlet_coldside_HEX1_charge = HotTESLiquid.setState_ph(p_nom_tank1, h_cold_out_HEX1_charge);
//NaK
  outlet_hotside_HEX1_charge = WorkingFluid.setState_ph(p_3_a_charge, h_hot_out_HEX1_charge);
//WF
  T_3_a_charge = outlet_hotside_HEX1_charge.T;
//irrev
  P_loss_irr_HEX1_charge = T0*((m_dot_WF_charge*(s_3_a_charge - s_3_charge)) + (m_dot_solsalt_HEX1_charge*(outlet_coldside_HEX1_charge.s - solsalt_tank2.s)));
  
//-------------HEX2 (Recuperation) CHARGE//
//pressure loss
  p_4_a_charge = p_4_charge*(1 + f_p_charge);
  p_2_charge = p_3_a_charge*(1 - f_p_charge);
//off-design
  UA_HEX2_charge/UA_HEX2_nom = (m_dot_WF_charge^0.8)/(m_dot_WF_nom_charge^0.8);
//cold side
  cp_cold_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_4_a_charge);
//just inlet cp, not average cp
  C_cold_HEX2_charge = m_dot_WF_charge*cp_cold_HEX2_charge;
// hot side
  cp_hot_HEX2_charge = WorkingFluid.specificHeatCapacityCp(state_3_a_charge_guess);
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
  Q_dot_max_HEX2_charge = C_min_HEX2_charge*(T_3_a_charge_guess - T_4_a_charge);
//energy balance and outlet states
// hot side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(state_3_a_charge.h - h_hot_out_HEX2_charge);
//cold side energy balance
  Q_dot_HEX2_charge = m_dot_WF_charge*(h_cold_out_HEX2_charge - state_4_a_charge.h);
//set outlets
  outlet_hotside_HEX2_charge = WorkingFluid.setState_ph(p_2_charge, h_hot_out_HEX2_charge);
  outlet_coldside_HEX2_charge = WorkingFluid.setState_ph(p_4_charge, h_cold_out_HEX2_charge);
  T_2_charge = outlet_hotside_HEX2_charge.T;
  T_4_charge = outlet_coldside_HEX2_charge.T;
//irrev
  P_loss_irr_HEX2_charge = T0*((m_dot_WF_charge*(s_2_a_charge - s_3_a_charge)) + (m_dot_WF_charge*(s_4_charge - s_4_a_charge)));
  
//-------------HEX 3 CHARGE//
//pressure loss
  p_1_charge = p_4_a_charge*(1 + f_p_charge);
//off-design
  UA_HEX3_charge/UA_HEX3_nom = ((m_dot_WF_charge^0.8*m_dot_methanol_HEX3_charge^0.8)/(m_dot_WF_nom_charge^0.8*m_dot_methanol_HEX3_nom^0.8))*((m_dot_WF_nom_charge^0.8 + m_dot_methanol_HEX3_nom^0.8)/(m_dot_WF_charge^0.8 + m_dot_methanol_HEX3_charge^0.8));
//charge: inlet cold side= state 1 inlet hot side=tank 3
  outlet_hotside_guess_HEX3_charge = ColdTESLiquid.setState_pT(p_nom_tank4, T_1_charge+6); //buffer to minimum temperature of methanol mixture through , which is fine, as it is a guess anyway, Tout, hot side will be higher
  outlet_coldside_guess_HEX3_charge = WorkingFluid.setState_pT(p_4_a_charge, T_tank3);
  cp_cold_ave_HEX3_charge = (WorkingFluid.specificHeatCapacityCp(state_1_charge) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX3_charge))/2;
  cp_hot_ave_HEX3_charge = (WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX3_charge) + coldliq_tank3_state.cp)/2;
//heat capacity rates
  C_cold_HEX3_charge = m_dot_WF_charge*cp_cold_ave_HEX3_charge;
  C_hot_HEX3_charge = C_cold_HEX3_charge;
//balanced operation
  m_dot_methanol_HEX3_charge = C_hot_HEX3_charge/cp_hot_ave_HEX3_charge;
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
  outlet_hotside_HEX3_charge = ColdTESLiquid.setState_ph(p_nom_tank4, h_hot_out_HEX3_charge);
  outlet_coldside_HEX3_charge = WorkingFluid.setState_ph(p_4_a_charge, h_cold_out_HEX3_charge);
  T_4_a_charge = outlet_coldside_HEX3_charge.T;
//irrev
  P_loss_irr_HEX3_charge = T0*((m_dot_WF_charge*(s_4_a_charge - s_1_charge)) + (m_dot_methanol_HEX3_charge*(outlet_hotside_HEX3_charge.s - coldliq_tank3_state.s)));


//-------------HEX 4 rejection charge//
  state_amb_air = WorkingFluid.setState_pT(1001315, T_amb);
  C_hot_HEXrej_charge = m_dot_WF_charge*WorkingFluid.specificHeatCapacityCp(state_2_charge);
  C_min_HEXrej_charge = C_hot_HEXrej_charge;
  C_cold_HEXrej_charge = C_max_HEXrej_charge;
  C_cold_HEXrej_charge = m_dot_rej_charge*WorkingFluid.specificHeatCapacityCp(state_amb_air);
  C_r_HEXrej_charge = C_min_HEXrej_charge/C_max_HEXrej_charge;
  Q_dot_max_HEXrej_charge = C_min_HEXrej_charge*(T_2_charge - T_amb);
  eff_HEXrej_charge = (1 - exp(-NTU_HEXrej_charge*(1 - C_r_HEXrej_charge)))/(1 - C_r_HEXrej_charge*exp(-NTU_HEXrej_charge*(1 - C_r_HEXrej_charge)));
  Q_dot_HEXrej_charge = eff_HEXrej_charge*Q_dot_max_HEXrej_charge;
//hot side energy balance
  outlet_hotside_HEXrej_charge = WorkingFluid.setState_pT(p_2_a_charge, T_2_a_charge);
  h_hot_out_HEXrej_charge = WorkingFluid.specificEnthalpy(outlet_hotside_HEXrej_charge);
  Q_dot_HEXrej_charge = (h_2_charge - h_hot_out_HEXrej_charge)*m_dot_WF_charge;
//irrev
  P_loss_irr_HEXrej_charge = T0*((m_dot_WF_charge*(s_2_a_charge - s_2_charge - ((h_2_a_charge - h_2_charge)/T0))));

//-------------SYSTEM CHARGE//
  P_mech_shaft = P_mech_CO_charge + P_mech_TU_charge;
  Q_pump_charge = Q_dot_HEX1_charge;
  COP = Q_pump_charge/P_mech_shaft;

  
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
//STATE 4_a charge
  state_4_a_charge = WorkingFluid.setState_pT(p_4_a_charge, T_4_a_charge);
  h_4_a_charge = WorkingFluid.specificEnthalpy(state_4_a_charge);
  s_4_a_charge = WorkingFluid.specificEntropy(state_4_a_charge);
//STATE 4 charge
  state_4_charge = WorkingFluid.setState_pT(p_4_charge, T_4_charge);
  h_4_charge = WorkingFluid.specificEnthalpy(state_4_charge);
  s_4_charge = WorkingFluid.specificEntropy(state_4_charge);
//STATE 3 a charge GUESS
  T_3_a_charge = T3_guess_control.u;
  T_3_a_charge_guess = T3_guess_control.y;
  p_3_a_charge = p3_guess_control.u;
  p_3_a_charge_guess = p3_guess_control.y;
  state_3_a_charge_guess = WorkingFluid.setState_pT(p_3_a_charge_guess, T_3_a_charge_guess);
  h_3_a_charge_guess = WorkingFluid.specificEnthalpy(state_3_a_charge_guess);
  s_3_a_charge_guess = WorkingFluid.specificEntropy(state_3_a_charge_guess);
//STATE 3 isentropic charge
  s_3_is_charge = s_4_charge "isentropic compressor outlet spec. entropy";
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
  state_2_a_charge = WorkingFluid.setState_pT(p_2_a_charge_guess, T_2_a_charge);
  p_2_a_charge=p_2a_guess_control.u;
  p_2a_guess_control.y=p_2_a_charge_guess;
  //state_2_a_charge = WorkingFluid.setState_pT(p_2_a_charge, T_2_a_charge);
  h_2_a_charge = WorkingFluid.specificEnthalpy(state_2_a_charge);
  s_2_a_charge = WorkingFluid.specificEntropy(state_2_a_charge);
//--------------------------EQUATIONS DISCHARGE-----------------------------//
//-------------COMPRESSOR DISCHARGE//
//reduced values compressor
  eta_is_CO_red = eta_is_CO/eta_is_CO_nom;
  beta_CO_red = beta_CO/beta_CO_nom;
  G_CO_red = m_dot_WF*sqrt(T_1)/p_1/(m_dot_WF_nom*sqrt(T_1_nom)/p_1_nom);
  eta_is_CO_red = (1 - c4*(1 - n_CO_red)^2)*(n_CO_red/G_CO_red)*(2 - n_CO_red/G_CO_red);
  n_CO_red = n_CO/sqrt(T_1)/(n_CO_nom/sqrt(T_1_nom));
//other compressor equations
  beta_CO = p_2/p_1;
  eta_is_CO = (h_2_is - h_1)/(h_2 - h_1);
  beta_CO_red = c1*G_CO_red^2 + c2*G_CO_red + c3;
  c1 = n_CO_red/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c2 = (p - 2*m*n_CO_red^2)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c3 = -1*(p*m*n_CO_red - m^2*n_CO_red^3)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  P_mech_CO = m_dot_WF*(h_2 - h_1);
  P_loss_irr_CO = T0*m_dot_WF*(s_2 - s_1);    

//-------------EXPANDER DISCHARGE//
//parameters
  alpha = sqrt(1.4 - 0.4*n_TU_red);
//reduced values turbine
  eta_is_TU_red = eta_is_TU/eta_is_TU_nom;
  n_TU_red = n_TU/sqrt(T_3)/(n_TU_nom/sqrt(T_3_nom));
  beta_TU_red = beta_TU/beta_TU_nom;
  G_TU_red = alpha*sqrt(T_3_nom/T_3)*sqrt((beta_TU^2 - 1)/(beta_TU_nom^2 - 1));
  eta_is_TU_red = (1 - t*(1 - n_TU_red)^2)*(n_TU_red/G_TU_red)*(2 - ((n_TU_red/G_TU_red)));
//other turbine equations
  beta_TU = p_3/p_4;
  eta_is_TU = (h_3 - h_4)/(h_3 -h_4_is);
  P_mech_TU = m_dot_WF*(h_4 - h_3);
  P_loss_irr_TU = T0*m_dot_WF*(s_4 - s_3);


//-------------HEX 1 DISCHARGE//
//pressure loss
  p_3 = p_3_a*(1 - f_p);
//off-design
  UA_HEX1/UA_HEX1_nom = ((m_dot_WF^0.8*m_dot_solsalt_HEX1^0.8)/(m_dot_WF_nom^0.8*m_dot_solsalt_HEX1_nom^0.8))*((m_dot_WF_nom^0.8 + m_dot_solsalt_HEX1_nom^0.8)/(m_dot_WF^0.8 + m_dot_solsalt_HEX1^0.8));
  //discharge: inlet cold side= state_3_a inlet hot side=solsalt_tank1_state
  outlet_hotside_guess_HEX1 = HotTESLiquid.setState_pT(p_nom_tank2, T_3_a);//NaK
  outlet_coldside_guess_HEX1 = WorkingFluid.setState_pT(p_3, solsalt_tank1_state.T);//WF
  cp_hot_ave_HEX1 = (solsalt_tank1_state.cp + outlet_hotside_guess_HEX1.cp)/2;
  cp_cold_ave_HEX1 = (WorkingFluid.specificHeatCapacityCp(state_3_a) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX1))/2;
  C_cold_HEX1 = m_dot_WF*cp_cold_ave_HEX1;
  C_hot_HEX1 = C_cold_HEX1;
//balanced operation
C_hot_HEX1 =  m_dot_solsalt_HEX1*cp_hot_ave_HEX1 ;
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
  outlet_hotside_HEX1 = HotTESLiquid.setState_ph(p_nom_tank2, h_hot_out_HEX1);//NaK
  outlet_coldside_HEX1 = WorkingFluid.setState_ph(p_3, h_cold_out_HEX1);//WF
  T_3 = outlet_coldside_HEX1.T;
//-------------HEX2 (Recuperation DISCHARGE)//
//off-design
  UA_HEX2/UA_HEX2_nom = (m_dot_WF^0.8)/(m_dot_WF_nom^0.8);
//pressure loss
  p_3_a = p_2*(1 - f_p);
  p_4 = p_4_a*(1 + f_p);
//cold side
  outlet_coldside_guess_HEX2 = WorkingFluid.setState_pT(p_3_a,T_4_guess);
  cp_cold_ave_HEX2=(WorkingFluid.specificHeatCapacityCp(state_2) + WorkingFluid.specificHeatCapacityCp(outlet_coldside_guess_HEX2))/2;
  // hot side
  outlet_hotside_guess_HEX2 = WorkingFluid.setState_pT(p_4_a,T_2);
  cp_hot_ave_HEX2=(WorkingFluid.specificHeatCapacityCp(state_4_guess) + WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX2))/2; //variables for effectiveness
  C_cold_HEX2 = m_dot_WF*cp_cold_ave_HEX2;   
  C_hot_HEX2 = m_dot_WF*cp_hot_ave_HEX2;     
   C_min_HEX2= min(C_cold_HEX2, C_hot_HEX2);
   C_max_HEX2 =max(C_cold_HEX2, C_hot_HEX2);
   C_r_HEX2=C_min_HEX2/C_max_HEX2;   
   NTU_HEX2 = UA_HEX2/C_min_HEX2;
   if C_r_HEX2 <1 then
    eff_HEX2 = (1 - exp(-NTU_HEX2*(1 - C_r_HEX2)))/(1 - C_r_HEX2*exp(-NTU_HEX2*(1 - C_r_HEX2)));   
  else
    eff_HEX2 = NTU_HEX2/(1 + NTU_HEX2);
  end if;
   //heat transferred
  Q_dot_HEX2=eff_HEX2*Q_dot_max_HEX2;
  Q_dot_max_HEX2= C_min_HEX2*(T_4_guess - T_2);  //T_guess instead of T_4
//energy balance and outlet states
// hot side energy balance
 Q_dot_HEX2 = m_dot_WF*(state_4_guess.h - h_hot_out_HEX2);
//cold side energy balance
  Q_dot_HEX2 = m_dot_WF*(h_cold_out_HEX2 - state_2.h);  //set outlets
 outlet_hotside_HEX2 =  WorkingFluid.setState_ph(p_4_a, h_hot_out_HEX2);//WF
 outlet_coldside_HEX2 = WorkingFluid.setState_ph(p_3_a, h_cold_out_HEX2);//WF
 T_3_a = outlet_coldside_HEX2.T;
 T_4_a = outlet_hotside_HEX2.T;  
 

  
  //-------------HEX 3 DISCHARGE//
  //pressure loss
//delta_HEX3 =
  p_1*(1 + f_p) = p_1_a;
//discharge: inlet cold side=tank 4  inlet hot side=state 1a
  outlet_hotside_guess_HEX3 = WorkingFluid.setState_pT(p_1, T_tank4);
  outlet_coldside_guess_HEX3 = ColdTESLiquid.setState_pT(p_nom_tank3, T_1_a);
  cp_cold_ave_HEX3 = (coldliq_tank4_state.cp + ColdTESLiquid.specificHeatCapacityCp(outlet_coldside_guess_HEX3))/2;  
  cp_hot_ave_HEX3 = (WorkingFluid.specificHeatCapacityCp(state_1_a) + WorkingFluid.specificHeatCapacityCp(outlet_hotside_guess_HEX3))/2;
//heat capacity rates
  C_hot_HEX3 = m_dot_WF*cp_hot_ave_HEX3;
  C_cold_HEX3 = C_hot_HEX3;//balanced operation
  C_cold_HEX3=m_dot_methanol_HEX3* cp_cold_ave_HEX3;
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
  Q_dot_max_HEX3 = C_min_HEX3*(T_1_a - T_tank4); // cold side energy balance
  Q_dot_HEX3 = m_dot_methanol_HEX3*(h_cold_out_HEX3-coldliq_tank4_state.h );   
//hot side energy balance
  Q_dot_HEX3 = m_dot_WF*(state_1_a.h-h_hot_out_HEX3);  //set outlets
  outlet_hotside_HEX3 = WorkingFluid.setState_ph(p_1, h_hot_out_HEX3);
  outlet_coldside_HEX3 = ColdTESLiquid.setState_ph(p_nom_tank3, h_cold_out_HEX3);  
  T_1 = outlet_hotside_HEX3.T;  
//irrev
  P_loss_irr_HEX3 = T0*((m_dot_WF*(s_1 - s_1_a)) + (m_dot_methanol_HEX3*(outlet_coldside_HEX3.s - coldliq_tank4_state.s)));    
  //-------------HEX  rejection DISCHARGE//
p_4_a = p_1_a;


  //state 1_a DISCHARGE
  state_1_a = WorkingFluid.setState_pT(p_1_a, T_1_a);
  h_1_a = WorkingFluid.specificEnthalpy(state_1_a);
  s_1_a = WorkingFluid.specificEntropy(state_1_a);
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
  T_2=state_2.T;
  s_2 = WorkingFluid.specificEntropy(state_2);
  //state 3a DISCHARGE
  state_3_a = WorkingFluid.setState_pT(p_3_a, T_3_a);
  h_3_a = WorkingFluid.specificEnthalpy(state_3_a);
  s_3_a = WorkingFluid.specificEntropy(state_3_a);
  //state 3 DISCHARGE
  state_3 = WorkingFluid.setState_pT(p_3, T_3);
  h_3 = state_3.h;
  s_3 = WorkingFluid.specificEntropy(state_3);
  //state 4 isentropic
  s_4_is = s_3;
  state_4_is = WorkingFluid.setState_ps(p_4, s_4_is) "isentropic state of turbine outlet";
  T_4_is = state_4_is.T "isentropic turbine outlet temperature";
  h_4_is=state_4_is.h;
  //state 4
  state_4 = WorkingFluid.setState_ph(p_4, h_4);
  T_4 = state_4.T "isentropic turbine outlet spec. entropy";
  s_4 = WorkingFluid.specificEntropy(state_4);
//state 4 guess DISCHARGE
  T_4 = T4_guess_control_discharge.u;
  T_4_guess= T4_guess_control_discharge.y;
  state_4_guess = WorkingFluid.setState_pT(p_4, T_4_guess);
  h_4_guess = WorkingFluid.specificEnthalpy(state_4_guess);
  s_4_guess = WorkingFluid.specificEntropy(state_4_guess);
//state 4_a DISCHARGE
  state_4_a = WorkingFluid.setState_pT(p_4_a, T_4_a);
  h_4_a = WorkingFluid.specificEnthalpy(state_4_a);
  s_4_a = WorkingFluid.specificEntropy(state_4_a);

  
  annotation(
    Documentation(info = "<html><head></head><body>Static Malta charge cycle<div>2nd order controller to guess T_3_a_charge</div><div>P_4_charge is set</div></body></html>"));
end StaticMalta_charge_discharge;