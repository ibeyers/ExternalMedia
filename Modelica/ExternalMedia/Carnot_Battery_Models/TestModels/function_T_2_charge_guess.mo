within ExternalMedia.Carnot_Battery_Models.TestModels;
model function_T_2_charge_guess
  import Modelica.Units.SI;
  import Modelica.Units.Conversions.from_degC;
  import ExternalMedia.Carnot_Battery_Models.Functions.guess_T_2_charge;
 SI.Temperature T_2_a_charge=from_degC(21) "outlet temperature after rejec"; //fixed temperature point
  SI.Temperature T_2_charge(start = 309.3) "temperature ";
equation
  T_2_charge=guess_T_2_charge(T_2_a_charge);
end function_T_2_charge_guess;
