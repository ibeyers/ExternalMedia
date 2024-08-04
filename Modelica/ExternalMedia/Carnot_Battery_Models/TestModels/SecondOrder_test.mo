within ExternalMedia.Carnot_Battery_Models.TestModels;

model SecondOrder_test
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  //import Constants
  import Modelica.Constants.pi;
  import Modelica.Constants.g_n;
  //import Functions
  //import ExternalMedia.Carnot_Battery_Models.Functions.guess_T_2_charge;

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
  
      //STATE 3a
  WorkingFluid.ThermodynamicState state_3_a_charge(p(start = 456931), T(start = from_degC(281.29))) "thermodynamic state";
  SI.Temperature T_3_a_charge(start = from_degC(281.29)) "outlet temperature after HEX";
SI.Pressure p_3_a_charge(start = 456931) "Pressure after HEX";
    WorkingFluid.SpecificEnthalpy h_3_a_charge(start=685666) "HEX outlet enthalpy";
    WorkingFluid.SpecificEntropy s_3_a_charge(start=4079) "HEX outlet spec. entropy";
     
    //state 3a charge  guess
    SI.Temperature T_3_a_charge_guess(start=from_degC(281.29))  "state  temperature";
    SI.Pressure p_3_a_charge_guess(start=456931)"state  pressure";
    WorkingFluid.ThermodynamicState state_3_a_charge_guess(p(start=456931),T(start=from_degC(281.29)));
    WorkingFluid.SpecificEnthalpy h_3_a_charge_guess(start=1006139) "enthalpy";
    WorkingFluid.SpecificEntropy s_3_a_charge_guess(start=4545) "spec. entropy";
  Modelica.Blocks.Continuous.SecondOrder secondOrder(w = 0.5, D = 0.4)  annotation(
    Placement(transformation(origin = {2, -4}, extent = {{-10, -10}, {10, 10}})));
initial equation
T_3_a_charge_guess=from_degC(281.29);
equation
 T_3_a_charge =from_degC(281.29)-0.06944444444444445 * time;
p_3_a_charge=456931-6.944444444444445 * time;
T_3_a_charge_guess=secondOrder.y;
T_3_a_charge=secondOrder.u;
p_3_a_charge_guess=456931;
//STATE 3 a charge GUESS
  state_3_a_charge_guess = WorkingFluid.setState_pT(p_3_a_charge_guess, T_3_a_charge_guess);
  h_3_a_charge_guess = WorkingFluid.specificEnthalpy(state_3_a_charge_guess);
  s_3_a_charge_guess = WorkingFluid.specificEntropy(state_3_a_charge_guess);    

 //STATE 3a charge
 state_3_a_charge = WorkingFluid.setState_pT(p_3_a_charge, T_3_a_charge);
  h_3_a_charge = WorkingFluid.specificEnthalpy(state_3_a_charge);
  s_3_a_charge = WorkingFluid.specificEntropy(state_3_a_charge);
    
end SecondOrder_test;
