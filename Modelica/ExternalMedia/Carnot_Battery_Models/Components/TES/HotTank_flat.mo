within ExternalMedia.Carnot_Battery_Models.Components.TES;

model HotTank_flat
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
  parameter SI.Temperature T_start_tank1 = from_degC(563.149);
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
  parameter SI.SpecificEnthalpy h_in_tank1 = 1; //dummy value TO_DO
  SI.HeatFlowRate Q_port_tank1(displayUnit = "MW") "Heat Flow from tank 2";
  parameter SI.Height x_min_tank1 = 0.4;
initial equation
  m_tank1 = m_nom_tank1;
  T_tank1 = T_start_tank1;  
equation
//-------------TANK 1//
m_out_tank1=-100;
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
annotation(
    Icon(graphics = {Ellipse(origin = {0, 84}, extent = {{-100, 16}, {100, -16}}), Ellipse(origin = {0, -84}, extent = {{-100, 16}, {100, -16}}), Line(origin = {-83.8563, -0.968165}, points = {{-16.1437, 82.9682}, {-16.1437, -83.0318}, {-16.1437, -83.0318}}), Line(origin = {100, 0}, points = {{0, 84}, {0, -84}})}));
end HotTank_flat;
