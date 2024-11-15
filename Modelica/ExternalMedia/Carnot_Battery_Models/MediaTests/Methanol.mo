within ExternalMedia.Carnot_Battery_Models.MediaTests;

model Methanol
  import Modelica.Units.SI.Temperature;
  import Modelica.Units.SI.Pressure;
  import Modelica.Units.SI.Density;  
  import Modelica.Units.Conversions.from_bar;
  import Modelica.Units.Conversions.from_degC;  
  package Methanol "NaK properties from CoolProp"
    extends ExternalMedia.Media.IncompressibleCoolPropMedium(mediumName = "MMA", substanceNames = {"MMA[0.6]"});
  end Methanol;

  replaceable package Medium =  Methanol constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model";  
    
    Medium.ThermodynamicState FluidState "Thermodynamic state of TES fluid"; 
 
   parameter Pressure p = from_bar(1) "TES fluid pressure";
   parameter Temperature T = from_degC(-70) "TES fluid temperature";
   Density d "density of TES fluid";   
   Medium.BaseProperties medium_a "Medium properties of port_a"; 
   Density d_basep;
equation
  FluidState = Medium.setState_pT(p,T); 
  d = FluidState.d;
  medium_a.p=p;
  medium_a.T=T;
  d_basep=medium_a.d;
end Methanol;
