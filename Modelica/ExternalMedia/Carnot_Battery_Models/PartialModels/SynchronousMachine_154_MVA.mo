within ExternalMedia.Carnot_Battery_Models.PartialModels;

model SynchronousMachine_154_MVA
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
 
//parameter SI.Power P_set(displayUnit = "MW") = -153.95*1000*1000;
//parameter SI.ReactivePower Q_set(displayUnit = "Mvar") = 0*1000*1000;
    
parameter SI.Power P_set(displayUnit = "MW") = -130.9*1000*1000;
parameter SI.ReactivePower Q_set(displayUnit = "Mvar") = -81.12*1000*1000;
  
 // parameter SI.Power P_set(displayUnit = "MW") = -174.25*1000*1000;
 // parameter SI.ReactivePower Q_set(displayUnit = "Mvar") = -108*1000*1000;
      
  //------------------SYNCHRONOUS MACHINE
  //parameters
  parameter SI.ApparentPower S_SM_nom(displayUnit = "MVA") = 153.95*1000*1000;
  parameter SI.Voltage U_SM_ST_nom(displayUnit = "kV") = 15.75*1000 "RMS voltage of the high voltage side";
  parameter SI.Frequency f = 50 " input frequency at stator";
  parameter Integer N_p = 1 "number of pole pairs";
  parameter NonSI.AngularVelocity_rpm n_SM_nom = 3000 "nominal synchronous speed in rpm";
  //geometry
  parameter SI.Diameter d_RO = 0.92 "Rotor diameter";
  parameter SI.Length l_FE = 3.882622607396857 "Rotor length";
  parameter SI.Mass m_ST_FE_teeth = 8953.562836941044 "iron mass of stator teeth";
  parameter SI.Mass m_ST_FE_yoke = 89542.43270181026 "iron mass of stator yoke";
  parameter Real tau_p = 1.4451326206513049 "pole pitch";
  //Equivalent circuit parameters
  parameter SI.Inductance L_sigma_ST = 0.0016156307265061366;
  parameter SI.Inductance L_h = 0.009155240783534776;
  parameter SI.Inductance L_d = 0.010770871510040912;
  parameter Integer N_FD = 72 "number of exitation field windings";
  parameter Integer N_ST = 11 "number of stator windings";
  // losses
  parameter Real k_schuisky = 5 "experimental factor for correlation of Schuisky";
  parameter Real k_additional = 0.001 "factor for additional losses";
  parameter SI.Resistance R_FD = 0.11228917141027794 "resistance in excitation windings";
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
  SI.ComplexCurrent I_SM_ST "stator side complex current";
  SI.Current I_SM_ST_abs;
  //air-gap
  SI.ComplexVoltage U_SM_h(re(start = U_SM_ST_nom/sqrt(3)));
  Complex S_SM_h;
  SI.Power P_airgap(displayUnit = "MW");
  SI.Voltage U_SM_h_abs;  
  //excitation
  SI.ComplexVoltage U_P(re(start = U_SM_ST_nom/sqrt(3)));
  SI.Voltage U_P_abs;
  SI.ComplexCurrent I_FD_ref "excitation current, referred to stator side";
  SI.Voltage U_FD "excitation voltage";
  SI.Power P_excitation(displayUnit = "MW") "power required for DC excitation";
  SI.Current I_FD(start = 900) "actual excitation current, DC";
  //equivalent circuit elements
  SI.Reactance X_sigma_ST;
  SI.Reactance X_h;
  SI.Reactance X_d;
  SI.Resistance R_SM_ST =0.004833955829814876;
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
  SI.Power P_Ohmic_Stator(displayUnit = "MW") "total iron losses";
  
equation

//MODE 1 CHARGE
  if Mode == 1 then
    P_additional = k_additional*abs(P_mech_RO);
    P_SM_loss = abs(P_SM_ST) - abs(P_mech_RO);
    P_mech_RO = P_airgap - P_windage_ventilation - P_additional - P_FE - P_excitation;
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST +U_P ;
    eta_SM = abs(P_mech_RO)/abs(P_SM_ST);
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    P_additional = k_additional*abs(P_mech_RO);
//change
    P_SM_loss = abs(P_mech_RO) - abs(P_SM_ST);
    P_mech_RO = P_airgap - (P_windage_ventilation + P_additional + P_FE + P_excitation);
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST + U_P;
    eta_SM = abs(P_SM_ST)/abs(P_mech_RO);
//MODE 0 HOLD
  else
    P_mech_RO = 0;
    P_additional = 0;
    P_SM_loss = 0;
    I_SM_ST = Complex(0, 0);
    U_P = Complex(0, 0);
    eta_SM = 0;
  end if;
 U_P_abs=ComplexMath.abs(U_P);
//------------------ELECTRICAL MACHINERY

//------------------SYNCHRONOUS MACHINE
//general
  omega_SM = 2*pi*f;
  N_p*n_RO = f*60;
//equivalent circuit elements
  X_h = omega_SM*L_h;
  X_d = omega_SM*L_d;
  X_sigma_ST = omega_SM*L_sigma_ST;
//Stator side Interfaces
  Q_SM_ST = Q_set;
  U_SM_ST = Complex((U_SM_ST_nom/sqrt(3))*cos(0), (U_SM_ST_nom/sqrt(3))*sin(0));
  P_SM_ST = P_set;
  S_SM_ST = Complex(P_SM_ST, Q_SM_ST);
//airgap
  U_SM_ST = Complex(R_SM_ST, X_sigma_ST)*I_SM_ST + U_SM_h;
  I_SM_ST_abs=ComplexMath.abs(I_SM_ST);
  U_SM_h_abs=ComplexMath.abs(U_SM_h);
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
annotation(
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,evaluateAllParameters,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end SynchronousMachine_154_MVA;
