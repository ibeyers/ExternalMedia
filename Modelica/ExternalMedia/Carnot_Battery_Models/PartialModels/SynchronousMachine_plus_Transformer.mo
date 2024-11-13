within ExternalMedia.Carnot_Battery_Models.PartialModels;

model SynchronousMachine_plus_Transformer
  //--------------------------IMPORTS-----------------------------//
  import Modelica.Units.SI;
  import Modelica.Units.NonSI;
  import Modelica.Units.Conversions.from_degC;
  import Modelica.ComplexMath;
  //import Constants
  import Modelica.Constants.pi;
  import Modelica.Constants.g_n;
  //--------------------------INPUTS
  parameter Integer Mode = 1;
  Modelica.Units.SI.Power P_GC_set(displayUnit = "MW");  
  parameter SI.Power P_set(displayUnit = "MW") = -42*1000*1000;
  parameter SI.Power P_set_charge(displayUnit = "MW") = 150*1000*1000;
  //connection point set parameters
  parameter SI.Voltage U_TR_set = 220*1000;
  parameter SI.Angle phi_TR_set = 0;
  SI.ReactivePower Q_TR_set(displayUnit = "MW",start = -50*1000*1000);
  
  //------------------ELECTRICAL MACHINERY
    //Grid connection point
  Modelica.Units.SI.ReactivePower Q_GC_set(start = -50*1000*1000, displayUnit = "Mvar");
  Modelica.Units.SI.ReactivePower Q_GC(start=-50*1000*1000,displayUnit = "Mvar"); //helper variable
  Modelica.Units.SI.ReactivePower S_GC_set;  
  parameter Real cos_phi_GC=0.95;
  //------------------TRANSFORMER
  //transformer parameters
  parameter SI.ApparentPower S_TR_nom(displayUnit = "MVA") = 208*1000*1000;
  parameter SI.Voltage U_TR_HV_nom(displayUnit = "kV") = 220*1000 "RMS voltage of the high voltage side (fixed by upper grid), line-to-line";
  parameter SI.Voltage U_TR_LV_nom(displayUnit = "kV") = 15750 "RMS voltage of the low voltage side, line-to-line";
  parameter Real a = U_TR_HV_nom/U_TR_LV_nom "turns ratio";
  parameter Real u_TR_ohmic_OC = 0.33604996008458077;
  parameter Real i_TR_OC = 0.32272375619545934;
  parameter Real P_TR_loss_OC_div_S_TR_nom = 0.06121076111169815;
  parameter SI.Power P_TR_loss_OC(displayUnit = "kW") = (P_TR_loss_OC_div_S_TR_nom/100)*S_TR_nom;
  parameter Real u_TR_SC = 11.875135026230133;
  parameter SI.Current I_TR_max = 537.638 "maximum transformer current";


    
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
  //------------------SYNCHRONOUS MACHINE
  //parameters
  parameter SI.ApparentPower S_SM_nom(displayUnit = "MVA") = 208*1000*1000;
  parameter SI.Voltage U_SM_ST_nom(displayUnit = "kV") = 15.75*1000 "RMS voltage of the high voltage side";
  parameter SI.Frequency f = 50 " input frequency at stator";
  parameter Integer N_p = 1 "number of pole pairs";
  parameter NonSI.AngularVelocity_rpm n_SM_nom = 3000 "nominal synchronous speed in rpm";
  //geometry
  parameter SI.Diameter d_RO = 0.92 "Rotor diameter";
  parameter SI.Length l_FE = 4.715853752596681 "Rotor length";
  parameter SI.Mass m_ST_FE_teeth = 10875.044312382426 "iron mass of stator teeth";
  parameter SI.Mass m_ST_FE_yoke = 108758.707701592 "iron mass of stator yoke";
  parameter Real tau_p = 1.4451326206513049 "pole pitch";
  //Equivalent circuit parameters
  parameter SI.Inductance L_sigma_ST = 0.0011957997612770178;
  parameter SI.Resistance R_ST = 0.005963040865384616;
  parameter SI.Inductance L_h = 0.006776198647236435;
  parameter SI.Inductance L_d = 0.007971998408513453;
  parameter Integer N_FD = 72 "number of exitation field windings";
  parameter Integer N_ST = 11 "number of stator windings";
  // losses
  parameter Real k_schuisky = 5 "experimental factor for correlation of Schuisky";
  parameter Real k_additional = 0.001 "factor for additional losses";
  parameter SI.Resistance R_FD = 0.13112806419264864 "resistance in excitation windings";
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
  //air-gap
  SI.ComplexVoltage U_SM_h(re(start = U_SM_ST_nom/sqrt(3)));
  Complex S_SM_h;
  SI.Power P_airgap(displayUnit = "MW");
  //excitation
  SI.ComplexVoltage U_P(re(start = U_SM_ST_nom/sqrt(3)));
  SI.ComplexCurrent I_FD_ref "excitation current, referred to stator side";
  SI.Voltage U_FD "excitation voltage";
  SI.Power P_excitation(displayUnit = "MW") "power required for DC excitation";
  SI.Current I_FD(start = 900) "actual excitation current, DC";
  //equivalent circuit elements
  SI.Reactance X_sigma_ST;
  SI.Reactance X_h;
  SI.Reactance X_d;
  SI.Resistance R_SM_ST = 0.002981520432692308;
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
    S_GC_set = sqrt(P_GC_set^2 + Q_GC^2);
    cos_phi_GC = P_GC_set/S_GC_set;
    Q_GC_set = abs(Q_GC)*(0 - 1);
//ensure that Q_GC_dis_set is negative, because OMPython has a problem with it
    Q_GC_set = Q_TR_set;
    
//MODE 1 CHARGE
  if Mode == 1 then
    P_GC_set=P_set_charge;
    P_TR_HV = P_set_charge;
    Q_TR_HV = Q_TR_set;
    eta_transformer = abs(P_TR_LV)/(abs(P_TR_LV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;
    P_additional = k_additional*abs(P_mech_RO);
    P_SM_loss = abs(P_SM_ST) - abs(P_mech_RO);
    P_mech_RO = P_airgap - P_windage_ventilation - P_additional - P_FE - P_excitation;
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST +U_P ;
    eta_SM = abs(P_mech_RO)/abs(P_SM_ST);
//MODE 2 DISCHARGE
  elseif Mode == 2 then
    P_GC_set=P_set;  
    Q_TR_HV = Q_TR_set;
    P_TR_HV = P_set;
    eta_transformer = abs(P_TR_HV)/(abs(P_TR_HV) + P_TR_loss);
    I_TR_HV + I_TR_LV_transferred = I_mag;
    P_additional = k_additional*abs(P_mech_RO);
//change
    P_SM_loss = abs(P_mech_RO) - abs(P_SM_ST);
    P_mech_RO = P_airgap - (P_windage_ventilation + P_additional + P_FE + P_excitation);
    S_SM_ST = 3*U_SM_ST*ComplexMath.conj(I_SM_ST);
    U_SM_h=Complex(0,X_h)*I_SM_ST + U_P;
    eta_SM = abs(P_SM_ST)/abs(P_mech_RO);
//MODE 0 HOLD
  else
    P_GC_set=0;
    Q_TR_HV = 0;
    P_TR_HV = 0;
    eta_transformer = 0;
    I_TR_HV + I_TR_LV_transferred = Complex(0, 0);
    P_mech_RO = 0;
    P_additional = 0;
    P_SM_loss = 0;
    I_SM_ST = Complex(0, 0);
    U_P = Complex(0, 0);
    eta_SM = 0;
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
//------------------SYNCHRONOUS MACHINE
//general
  omega_SM = 2*pi*f;
  N_p*n_RO = f*60;
//equivalent circuit elements
  X_h = omega_SM*L_h;
  X_d = omega_SM*L_d;
  X_sigma_ST = omega_SM*L_sigma_ST;
//Trafo side Interfaces
  0 = Q_TR_LV + Q_SM_ST;
  U_SM_ST = U_TR_LV;
  0 = P_TR_LV + P_SM_ST;
  S_SM_ST = Complex(P_SM_ST, Q_SM_ST);
//airgap
  U_SM_ST = Complex(R_SM_ST, X_sigma_ST)*I_SM_ST + U_SM_h;
  S_SM_h = 3*U_SM_h*ComplexMath.conj(I_SM_ST);
  P_airgap = S_SM_h.re;
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
end SynchronousMachine_plus_Transformer;
