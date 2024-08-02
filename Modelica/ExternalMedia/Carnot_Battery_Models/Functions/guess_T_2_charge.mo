within ExternalMedia.Carnot_Battery_Models.Functions;


function guess_T_2_charge
package WorkingFluid = ExternalMedia.Media.CoolPropMedium(mediumName = "Air", substanceNames = {"Air"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph);  
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;  
  
  input SI.Temperature T_2_a_charge         "Independent variable";
  input SI.Pressure p_2_a_charge         "Independent variable";  
  input SI.MassFlowRate m_dot_WF_charge  "mass flow rate";  
  input SI.Temperature T_amb        "Independent variable";   
  output SI.Temperature T_2_charge        "Dependent variable";
protected
  Integer n_max=10;  //iteration
  Real NTU_HEXrej_charge=7.8 "Number of transfer units of rejection HEX";
  ThermodynamicState state_2_a_charge(p(start=452019),T(start=from_degC(25)))
 "thermodynamic state of rejec outlet";  
  ThermodynamicState state_amb_air(p(start=1001315),T(start=from_degC(20))) "thermodynamic state of rejec outlet";
  ThermodynamicState state_2_charge_guess(p(start=452019),T(start=from_degC(30)))
 "thermodynamic state of rejec outlet";    
    //hot-side
   Real C_hot_HEXrej_charge "Heat capacity rate of cold side of rejection HEX (after compressor)";   
   Real C_min_HEXrej_charge; 
  //cold-side
   Real C_cold_HEXrej_charge "Heat capacity rate of cold side of rejection HEX ";     
   Real C_max_HEXrej_charge;
   
  /*
  //Integer i; 
//-------------HEX rejection charge//
  //Real UA_HEX3_nom;
  parameter 
 
   //SI.SpecificHeatCapacity cp_cold_HEX2_charge;   
 
   //SI.SpecificHeatCapacity cp_hot_HEX2_charge;          
//variables for effectiveness

   Real C_r_HEXrej_charge;    
   parameter SI.Temperature T_amb=from_degC(20);
   
   Real eff_HEXrej_charge(start=0.92) "effectiveness of hot HEX";  
   //
  SI.HeatFlowRate Q_dot_HEXrej_charge(displayUnit = "MW");
  SI.HeatFlowRate Q_dot_max_HEXrej_charge(displayUnit = "MW");     
  SI.SpecificEnthalpy h_hot_in_HEXrej_charge(start=431395) ;
  WorkingFluid.ThermodynamicState inlet_hotside_HEXrej_charge(p(start=p_2_a_nom_charge),T(start=from_degC(36)));
  */  
algorithm
  T_2_charge := 300;
  state_amb_air:=WorkingFluid.setState_pT(1001315,T_amb);  
  state_2_a_charge := WorkingFluid.setState_pT(p_2_a_charge, T_2_a_charge);  
  C_hot_HEXrej_charge := m_dot_WF_charge * WorkingFluid.specificHeatCapacityCp(state_2_a_charge);  
  C_min_HEXrej_charge := C_hot_HEXrej_charge; 
  C_max_HEXrej_charge := C_min_HEXrej_charge*2;
  C_cold_HEXrej_charge :=  C_max_HEXrej_charge ;  
  /*
   //-------------HEX 4 rejection charge//



 C_r_HEXrej_charge = C_min_HEXrej_charge/C_max_HEXrej_charge;

 eff_HEXrej_charge = (1 - exp(-NTU_HEXrej_charge*(1 - C_r_HEXrej_charge)))/(1 - C_r_HEXrej_charge*exp(-NTU_HEXrej_charge*(1 - C_r_HEXrej_charge)));
 //
  Q_dot_HEXrej_charge = eff_HEXrej_charge*Q_dot_max_HEXrej_charge;
  Q_dot_max_HEXrej_charge =C_min_HEXrej_charge*(T_2_charge - T_amb);   
  // hot side energy balance
  Q_dot_HEXrej_charge = m_dot_WF_charge*(WorkingFluid.specificEnthalpy(inlet_hotside_HEXrej_charge)-h_hot_in_HEXrej_charge);    
  inlet_hotside_HEXrej_charge =WorkingFluid.setState_ph(p_2_charge, h_hot_in_HEXrej_charge);
  //T_2_charge=inlet_hotside_HEXrej_charge.T;
*/
  //while x>=ybar[i+1,1] loop
  //  i := i + 1;
  //end while;  
end guess_T_2_charge;
