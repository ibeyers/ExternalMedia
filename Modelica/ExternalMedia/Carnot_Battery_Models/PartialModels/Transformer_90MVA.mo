within ExternalMedia.Carnot_Battery_Models.PartialModels;

model Transformer_90MVA
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.NonSI;
  import Modelica.Units.Conversions.from_degC;
  import Modelica.ComplexMath;
  //import Constants
  import Modelica.Constants.pi;
  import Modelica.Constants.g_n;
  //--------------------------INPUTS
  parameter Integer Mode = 2;
  parameter SI.Power P_set(displayUnit = "MW") = -81*1000*1000;
  parameter Modelica.Units.SI.ReactivePower Q_set(displayUnit = "Mvar") = -39.23*1000*1000;
  
  parameter SI.Power P_set_charge(displayUnit = "MW") = 150*1000*1000;
  //connection point set parameters
  parameter SI.Voltage U_TR_set = U_TR_HV_nom;
  parameter SI.Angle phi_TR_set = 0;
  SI.ReactivePower Q_TR_set(displayUnit = "MW",start = -39.23*1000*1000);

  //------------------TRANSFORMER
  //transformer parameters
  parameter SI.ApparentPower S_TR_nom(displayUnit = "MVA") = 90*1000*1000;
  parameter SI.Voltage U_TR_HV_nom(displayUnit = "kV") = 66*1000 "RMS voltage of the high voltage side (fixed by upper grid), line-to-line";
  parameter SI.Voltage U_TR_LV_nom(displayUnit = "kV") = 11500 "RMS voltage of the low voltage side, line-to-line";
  parameter Real a = U_TR_HV_nom/U_TR_LV_nom "turns ratio";
  parameter Real u_TR_ohmic_OC = 0.38119879894783126;
  parameter Real i_TR_OC = 0.41963432776291854;
  parameter Real P_TR_loss_OC_div_S_TR_nom = 0.07714640343180293;
  parameter SI.Power P_TR_loss_OC(displayUnit = "kW") = (P_TR_loss_OC_div_S_TR_nom/100)*S_TR_nom;
  parameter Real u_TR_SC = 10.680905579362;
  parameter SI.Current I_TR_max = 824.7860988423226 "maximum transformer current";
  
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
  
equation
//ensure that Q_GC_dis_set is negative, because OMPython has a problem with it
    Q_set = Q_TR_set;
    
//MODE 1 CHARGE
  if Mode == 1 then
    P_TR_HV = P_set_charge;
    Q_TR_HV = Q_TR_set;
    eta_transformer = abs(P_TR_LV)/(abs(P_TR_LV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;

//MODE 2 DISCHARGE
  elseif Mode == 2 then
       P_TR_HV = P_set;
       Q_TR_HV = Q_TR_set;
    eta_transformer = abs(P_TR_HV)/(abs(P_TR_HV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;

//MODE 0 HOLD
  else
    Q_TR_HV = 0;
    P_TR_HV = 0;
    eta_transformer = 0;
    I_TR_HV + I_TR_LV_transferred = Complex(0, 0);

  end if;
//------------------ELECTRICAL MACHINERY
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

annotation(
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,evaluateAllParameters,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end Transformer_90MVA;
