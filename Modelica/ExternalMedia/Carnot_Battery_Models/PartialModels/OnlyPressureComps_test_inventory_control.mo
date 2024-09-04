within ExternalMedia.Carnot_Battery_Models.PartialModels;

model OnlyPressureComps_test_inventory_control
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  //-------------Discharge//
   Real m_dot_div_p(start = 0.0076);
  //design
  parameter SI.MassFlowRate m_dot_WF_nom = 762 "design mass flow rate";
  //actual
  SI.MassFlowRate m_dot_WF = 762 "mass flow rate";
  parameter Real f_p = 0.01625 "pressure_loss_factor percent";
  parameter Real k_p=0.0062084"pressure_loss_factor";
  
  parameter SI.Pressure p_fix = 100000 "fixed pressure point through expansion vessel";

  //fixed pressure point charge
  SI.Pressure p_1(start = p_fix) "state 1 pressure";
  //STATE 1 a discharge
  SI.Pressure p_1_a(start = 103606) "pressure after Heat rejection ";
  SI.Pressure p_2(start = 591100) "pressure coming out of compressor";
 //STATE 3a discharge
  SI.Pressure p_3_a(start = 587493) "Pressure after recup";
  //STATE 3 discharge
  SI.Pressure p_3(start = 583886) "Pressure after recup";
  //STATE 4 discharge
  SI.Pressure p_4(start = 107213) " pressure at turb outlet";
  //STATE 4 a discharge
  SI.Pressure p_4_a(start = 103606) " pressure at recuperation outlet";

  SI.Pressure delta_P_HEX1;
  SI.Pressure delta_P_HEX2;
  SI.Pressure delta_P_HEX3;


  //state 1 discharge
  SI.Temperature T_1= from_degC(-53.760) " temperature";
  SI.Temperature T_3 "outlet temperature after HEX";
  SI.Temperature T_3_guess;
  parameter Real p = 1.8 "compressor map factor";
  parameter Real m = 1.4 "compressor map factor";
  parameter Real c4 = 0.3 "factor 4, see Zhang2002";
  
  //-------------COMPRESSOR DISCHARGE/
  Real c1 "factor 1";
  Real c2 "factor 2";
  Real c3 "factor 3";
  
  //design
  parameter Real beta_CO_nom = 5.911 "design compression ratio";
  parameter Real n_CO_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_CO_nom = 0.87 "design isentropic efficiency";
  SI.Temperature T_1_nom = from_degC(-53.76) "state 1 temperature";
  SI.Pressure p_1_nom = p_fix "state 1 pressure";
  //actual
  Real beta_CO(start = beta_CO_nom) "absolute compression ratio";
  parameter Real n_CO = 3000 "actual speed";

  //reduced
  Real beta_CO_red(start = 1) "reduced compression ratio";
  Real n_CO_red(start = 1) "reduced speed";
  Real G_CO_red(start = 1) "reduced mass flow rate compressor";

  //-------------EXPANDER DISCHARGE//
  Real alpha "factor, see Zhang2002";
  //design
  parameter Real beta_TU_nom = 5.445 "design expansion ratio";
  parameter Real n_TU_nom = 3000 "design speed";
  parameter SI.Efficiency eta_is_TU_nom = 0.94 "design isentropic efficiency";
  parameter SI.Temperature T_3_nom = from_degC(555) "design turbine inlet temperature";
  parameter SI.Pressure p_3_nom = 583886 "design turbine inlet pressure";
  //actual
  Real beta_TU(start = beta_TU_nom) "absolute expansion ratio";
  Real n_TU(start = 3000) "actual speed";

  //reduced
  Real beta_TU_red(start = 1) "reduced expansion ratio";
  Real n_TU_red(start = 1) "reduced speed";

  Real G_TU_red(start = 1) "reduced mass flow rate turbine";
 Modelica.Blocks.Continuous.SecondOrder T3_guess_control_discharge(w = 0.5, D = 0.4)  annotation(
    Placement(transformation(extent = {{-10, -10}, {10, 10}})));
initial equation

  T_3_guess = from_degC(556);
equation
  m_dot_div_p = m_dot_WF/p_1;
  p_1=98000;
//-------------COMPRESSOR DISCHARGE//
//reduced values compressor
  beta_CO_red = beta_CO/beta_CO_nom;
  G_CO_red = m_dot_WF*sqrt(T_1)/p_1/(m_dot_WF_nom*sqrt(T_1_nom)/p_1_nom);
  n_CO_red = n_CO/sqrt(T_1)/(n_CO_nom/sqrt(T_1_nom));
//other compressor equations
  beta_CO = p_2/p_1;
  beta_CO_red = c1*G_CO_red^2 + c2*G_CO_red + c3;
  c1 = n_CO_red/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c2 = (p - 2*m*n_CO_red^2)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
  c3 = -1*(p*m*n_CO_red - m^2*n_CO_red^3)/(p*(1 - m/n_CO_red) + n_CO_red*(n_CO_red - m)^2);
//-------------EXPANDER DISCHARGE//
  n_TU=n_CO;
//parameters
  alpha = sqrt(1.4 - 0.4*n_TU_red);
//reduced values turbine
  n_TU_red = n_TU/sqrt(T_3_guess)/(n_TU_nom/sqrt(T_3_nom));
  beta_TU_red = beta_TU/beta_TU_nom;
  //m_dot_WF/m_dot_WF_nom = alpha*sqrt(T_3_nom/T_3_guess)*sqrt((beta_TU^2 - 1)/(beta_TU_nom^2 - 1));
  G_TU_red = m_dot_WF*sqrt(T_3_guess)/p_3/(m_dot_WF_nom*sqrt(T_3_nom)/p_3_nom);
//other turbine equations
  beta_TU = p_3/p_4;
//pressure losses HEX
  delta_P_HEX1 = k_p*m_dot_WF^2;
  delta_P_HEX2 = k_p*m_dot_WF^2;
  delta_P_HEX3 = k_p*m_dot_WF^2;
//-------------HEX 1 DISCHARGE//
//pressure loss
  p_3 = p_3_a - delta_P_HEX1;
//-------------HEX 2 DISCHARGE//
//pressure loss
  p_3_a = p_2 - delta_P_HEX2;
  p_4_a = p_4 - delta_P_HEX2;
//-------------HEX 3 DISCHARGE//
//delta_HEX3 =
  p_1 = p_1_a - delta_P_HEX3;
//-------------HEX 4 rejection DISCHARGE//
 p_4_a = p_1_a;
  
   T_3= from_degC(556)-0.0006944444444444445 * time;
  //state 3 DISCHARGE guess
  T_3 = T3_guess_control_discharge.u;
  T_3_guess = T3_guess_control_discharge.y;
  
end OnlyPressureComps_test_inventory_control;
