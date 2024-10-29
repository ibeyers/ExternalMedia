within ExternalMedia.Carnot_Battery_Models.PartialModels;

model OnlyPressureComps_charge_const_HEX
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  //-------------charge//
  //design
  Real m_dot_div_p_charge(start = 0.00766);

  parameter SI.MassFlowRate m_dot_WF_nom_charge = 766 "design mass flow rate";  
  //actual
  SI.MassFlowRate m_dot_WF_charge(start= 766*scaling_factor) "mass flow rate";
  parameter Real f_p_charge = 0.01075 "pressure_loss_factor percent"; //not needed anymore
  parameter Real k_p_charge=0.0041715"pressure_loss_factor";
   SI.Pressure p_fix_charge(start = p_4_nom_charge*scaling_factor) "fixed pressure point through expansion vessel at p_4";   
  parameter Real scaling_factor=0.9;
//STATE 4 charge
  //fixed pressure point charge
  SI.Pressure p_4_charge(start=100000*scaling_factor) "state 4 pressure";
  SI.Pressure p_3_charge(start = 459203*scaling_factor) "pressure coming out of compressor";  
   SI.Pressure p_3_a_charge(start = 456754*scaling_factor) "Pressure after HEX";
       SI.Pressure p_2_charge(start = 454306*scaling_factor) "state  pressure";  
   SI.Pressure p_2_a_charge(start = 454306*scaling_factor) "Pressure after rejec";  
    SI.Pressure p_1_charge(start = 104896*scaling_factor) " pressure ";      
     SI.Pressure p_4_a_charge(start = 102448*scaling_factor) "state  pressure";


  SI.Pressure delta_P_HEX1_charge;
  SI.Pressure delta_P_HEX2_charge;
  SI.Pressure delta_P_HEX3_charge;

  SI.Temperature T_2_a_charge = from_degC(25) "outlet temperature after rejec";

  //state 4 charge guess
  SI.Temperature T_4_charge_guess(start = from_degC(267.533)) "state 4 temperature"; 
  //state 4 charge
 
  SI.Temperature T_4_charge(start = from_degC(267.533)) "state 4 temperature";

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
  parameter SI.Efficiency eta_is_CO_nom_charge = 0.90535  "design isentropic efficiency";
  SI.Temperature T_4_nom_charge = from_degC(267.533) "state 4 temperature";
  parameter SI.Pressure p_4_nom_charge = 100000 "state 4 pressure";
  //actual
  Real beta_CO_charge(start = beta_CO_nom_charge,min=0) "absolute compression ratio";
  parameter Real n_CO_charge = 3000 "actual speed";
  SI.Efficiency eta_is_CO_charge(start=0.9) "absolute isentropic efficiency";
  //reduced
  Real beta_CO_red_charge(start = 1,min=0) "reduced compression ratio";
  Real n_CO_red_charge(start = 1) "reduced speed";
  SI.Efficiency eta_is_CO_red_charge(start=1,min=0) "reduced isentropic efficiency";
  Real G_CO_red_charge(start = 1,min=0) "reduced mass flow rate compressor";
  
  //-------------EXPANDER CHARGE//
  parameter Real t = 0.3 "parameter, see Zhang2002";

  Real alpha_charge(start = 1) "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom_charge = 4.331 "design expansion ratio";
  parameter Real n_TU_nom_charge = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom_charge = 0.92 "design isentropic efficiency";
  parameter SI.Temperature T_2_a_nom_charge = from_degC(25) "nominal inlet temperature";
  parameter SI.Pressure p_2_a_nom_charge = 442461 "state pressure";
  //actual
  parameter Real n_TU_charge = 3000 "actual speed";
  Real beta_TU_charge(start = beta_TU_nom_charge) "absolute expansion ratio";
  SI.Efficiency eta_is_TU_charge(start=0.92) "absolute isentropic turbine efficiency";
 SI.Efficiency eta_is_TU_red_charge(start=1) "absolute isentropic turbine efficiency";
  //reduced
  Real n_TU_red_charge(start = 1) "reduced speed";
  Real G_TU_red_charge(start = 1) "reduced mass flow rate turbine";
  Real beta_TU_red_charge(start = 1) "reduced expansion ratio";
  
  
 Modelica.Blocks.Continuous.SecondOrder T4_guess_control(w = 0.5, D = 0.4)  annotation(
    Placement(transformation(origin = {-64, -38}, extent = {{-10, -10}, {10, 10}})));
initial equation

  T_4_charge_guess = from_degC(267.5);

equation
  p_4_charge = p_fix_charge;
  p_fix_charge = p_4_nom_charge*scaling_factor;
m_dot_div_p_charge = m_dot_WF_charge/p_4_charge;
//T_4_charge= 540.7 -0.0006944444444444445 * time;
  T_4_charge = from_degC(267.533);
//-------------COMPRESSOR CHARGE//
//reduced values compressor
eta_is_CO_red_charge = eta_is_CO_charge/eta_is_CO_nom_charge;
  beta_CO_red_charge = beta_CO_charge/beta_CO_nom_charge;
  G_CO_red_charge = m_dot_WF_charge*sqrt(T_4_charge_guess)/p_4_charge/(m_dot_WF_nom_charge*sqrt(T_4_nom_charge)/p_4_nom_charge);
eta_is_CO_red_charge = (1 - c4*(1 - n_CO_red_charge)^2)*(n_CO_red_charge/G_CO_red_charge)*(2 - n_CO_red_charge/G_CO_red_charge);
  n_CO_red_charge = n_CO_charge/sqrt(T_4_charge_guess)/(n_CO_nom_charge/sqrt(T_4_nom_charge));
//other compressor equations
  beta_CO_charge = p_3_charge/p_4_charge;
//  eta_is_CO_charge = (h_3_is_charge - h_4_charge)/(h_3_charge - h_4_charge);
  beta_CO_red_charge = c1_charge*G_CO_red_charge^2 + c2_charge*G_CO_red_charge + c3_charge;
  c1_charge = n_CO_red_charge/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c2_charge = (p - 2*m*n_CO_red_charge^2)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
  c3_charge = -1*(p*m*n_CO_red_charge - m^2*n_CO_red_charge^3)/(p*(1 - m/n_CO_red_charge) + n_CO_red_charge*(n_CO_red_charge - m)^2);
//-------------EXPANDER CHARGE//
//parameters
  alpha_charge = sqrt(1.4 - 0.4*n_TU_red_charge);
//reduced values expander
  n_TU_red_charge = n_TU_charge/sqrt(T_2_a_charge)/(n_TU_nom_charge/sqrt(T_2_a_nom_charge));
  G_TU_red_charge = m_dot_WF_charge*sqrt(T_2_a_charge)/p_2_a_charge/(m_dot_WF_nom_charge*sqrt(T_2_a_nom_charge)/p_2_a_nom_charge);
 m_dot_WF_charge/m_dot_WF_nom_charge = alpha_charge*sqrt(T_2_a_nom_charge/T_2_a_charge)*sqrt((beta_TU_charge^2 - 1)/(beta_TU_nom_charge^2 - 1));
  beta_TU_red_charge = beta_TU_charge/beta_TU_nom_charge;
eta_is_TU_red_charge = (1 - t*(1 - n_TU_red_charge)^2)*(n_TU_red_charge/G_TU_red_charge)*(2 - ((n_TU_red_charge/G_TU_red_charge)));
eta_is_TU_red_charge = eta_is_TU_charge/eta_is_TU_nom_charge;
//other turbine equations
  beta_TU_charge = p_2_a_charge/p_1_charge;
//pressure losses HEX
  delta_P_HEX1_charge = k_p_charge*m_dot_WF_nom_charge^2;
  delta_P_HEX2_charge = k_p_charge*m_dot_WF_nom_charge^2;
  delta_P_HEX3_charge = k_p_charge*m_dot_WF_nom_charge^2;
//-------------HEX 1 CHARGE//
//pressure loss
  p_3_charge = p_3_a_charge + delta_P_HEX1_charge;
//-------------HEX 2 DISCHARGE//
//pressure loss
  p_4_a_charge = p_4_charge + delta_P_HEX2_charge;
  p_2_charge = p_3_a_charge - delta_P_HEX2_charge;
//-------------HEX 3 DISCHARGE//
  p_1_charge = p_4_a_charge + delta_P_HEX3_charge;
//-------------HEX 4 rejection CHARGE//
  p_2_a_charge = p_2_charge;
/*
   T_3= from_degC(556)-0.0006944444444444445 * time;
  */
//STATE 4 charge guess
  T_4_charge = T4_guess_control.u;
  T_4_charge_guess = T4_guess_control.y;

end OnlyPressureComps_charge_const_HEX;
