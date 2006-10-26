/* *****************************************************************
 * Interface of the base class BaseTwoPhaseMedium.cpp
 *
 * The BaseTwoPhaseMedium class defines all the variables and member
 * functions which are needed to use external Modelica medium models
 * extending from PartialExternalTwoPhaseMedium.
 * 
 * Francesco Casella, Christoph Richter Sep 2006
 ********************************************************************/

#ifndef BASETWOPHASEMEDIUM_H_
#define BASETWOPHASEMEDIUM_H_

#include "include.h"

#include "basesolver.h"
#include "twophasemediumproperties.h"

class BaseTwoPhaseMedium{
public:
	BaseTwoPhaseMedium(const string &mediumName, const string &libraryName, 
					   const string &substanceName, BaseSolver *const solver,
					   const int &uniqueID);
	virtual ~BaseTwoPhaseMedium();

	int uniqueID() const;

	double beta() const;
	double cp() const;
	double cv() const;
	double d() const;
	double dd_dp_h() const;
	double dd_dh_p() const;
	double h() const;
	double kappa() const;
	double p() const;
	double s() const;
	double T() const;

	double ps() const;
	double Ts() const;

	double dl() const;
	double dv() const;
	double hl() const;
	double hv() const;
	double sl() const;
	double sv() const;

	double dc() const;
	double pc() const;
	double Tc() const;

	double MM() const;

	double eta() const;
	double lambda() const;
	double Pr() const;
	double sigma() const;

	virtual void setSat_p(const double &p);
	virtual void setSat_T(const double &T);

	virtual double saturationPressure(const double &T, const string &mediumName);
	virtual double saturationTemperature(const double &p, const string &mediumName);

	virtual void setState_dT(const double &d, const double &T, const int &phase);
	virtual void setState_ph(const double &p, const double &h, const int &phase);
	virtual void setState_ps(const double &p, const double &s, const int &phase);
	virtual void setState_pT(const double &p, const double &T);

protected:	
	// Pointer to medium property record
	TwoPhaseMediumProperties *_properties;

	// Pointer to solver
	BaseSolver *_solver;
};

#endif /*BASETWOPHASEMEDIUM_H_*/