#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"
#include "kiss/logger.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {


HadronicInteraction::HadronicInteraction(ref_ptr<Density> massDensity, double limit, int flag) {
	this->massDensity = massDensity;
	this->limit = limit;
	this->flag_Function = flag;
	setDescription("HadronicInteraction");
}

void HadronicInteraction::setMassDensity(ref_ptr<Density> dens) {
	this->massDensity = dens;
}

void HadronicInteraction::setLimit(double lim) {
	this->limit = lim;
}
// ------------------------------------ Kellner 2006 --------------------------------------------------------------------
double HadronicInteraction::KellnerSpectrumPion(double x, double ePrimary) const {
	// Kellner+ 2006, eqn. 12
	const double E0pi = 139.571 * MeV;
	if (x < E0pi / ePrimary)
		return 0.;
	const double L = std::log(ePrimary / TeV);
	const double a = 3.67 + 0.83 * L + 0.075 * L * L;
	const double Bpi = a + 0.25;
	const double alpha = 0.98 / std::sqrt(a);
	const double r = 2.6 / std::sqrt(a);
	const double xa = std::pow(x, alpha);
	const double term1 = 4. * alpha * Bpi * std::pow(x, alpha - 1.);
	const double term2 = std::pow((1. - xa) / (1. + r * xa * (1. - xa)), 4);
	const double term3 = 1. / (1. - xa) + r * (1. - 2. * xa) / (1. + r * xa * (1. - xa));
	const double term4 = std::sqrt(1. - E0pi  / (x * ePrimary));
	const double Fpi = term1 * term2 * term3 * term4;
	return Fpi;
}

double HadronicInteraction::KellnerSpectrumPhoton(double x, double ePrimary) const {
	// Kellner+ 2006, eqn. 58
	const double L = std::log(ePrimary / TeV);
	const double B = 1.30 + 0.14 * L + 0.011 * L * L;
	const double beta = 1. / (1.79 + 0.11 * L + 0.008 * L * L);
	const double k = 0.801 + 0.049 * L + 0.014 * L * L;
	const double xb = std::pow(x, beta);
	const double term1 = B * std::log(x) / x;
	const double term2 = (1. - xb) / (1 + k * xb * (1 - xb));
	const double term3 = 1 / std::log(x) - 4 * beta * xb / (1 - xb);
	const double term4 = 4 * k * beta * xb * (1 - 2 * xb) / (1 + k * xb * (1 - xb));
	const double Fgamma = term1 * std::pow(term2, 4) * (term3 - term4);
	return Fgamma;
}

double HadronicInteraction::KellnerXSection(double ePrimary) const {
	// Kellner+ 2006, eqn. 73
	const double L = std::log(ePrimary / TeV);
	const double A = 1 - std::pow(1.22 * 1e-3 * TeV / ePrimary, 4);
	return (34.3 + 1.88 * L + 0.25 * L * L) * A * A * 1e-31;
}

// ---------------------------------- Kamae 2006 ------------------------------------------------------------------------
double HadronicInteraction::KamaeXSection(double ePrimary) const {
	//KISS_LOG_WARNING << "Kamae 2006 does not provide a x-Section \n"
	//	<< "use Kellner x-Section instead\n";
	return KellnerXSection(ePrimary);
}

double HadronicInteraction::KamaeSpectrumPhoton(double Esec, double Tp) const {
	double tot = 0;
	tot += KamaeGammaND(Esec, Tp);
	tot += KamaeGammaDiff(Esec, Tp);
	tot += KamaeGamma1232(Esec, Tp);
	tot += KamaeGamma1600(Esec, Tp);
	return tot;
}

double HadronicInteraction::KamaeGammaND(double Esec, double Tp) const {
	Esec/GeV; // all in units of GeV
	Tp/GeV;
	double x, y, z;
    double Wl, Wh, Lmin, Lmax;
    double xa3, xa8;
    double pow1;
    double sigma;
    double cutoff, r_factor;
    double a[9];
    int i;             
    /* init some variables, given in table 2 */   
    x = std::log10(Esec);
   	y = log10(Tp*0.001);
    Lmin = -2.6;
    Lmax = L_MAX[0]*(y + 3.0);
    Wl = W_NDL[0];
    Wh = W_NDH[0];

     if ((Tp > 0.487) && (Tp < 512000.1))
    	{
        z = y + 3.3;
        a[0] = z*(-0.51187 + z*(7.6179 + z*(-2.1332 + 0.22184*z)));
        a[1] = -1.2592e-5 + 1.4439e-5*exp(-0.29360*(y + 3.4)) + 5.9363e-5/(y + 4.1485) + y*(2.2640e-6 - 3.3723e-7*y);
        a[2] = -174.83 + 152.78*log10(1.5682*(y + 3.4)) - 808.74/(y + 4.6157);
        a[3] = 0.81177 + y*(0.56385 + y*(0.0040031 + y*(-0.0057658 + 0.00012057*y)));
        z = y + 3.32;
        a[4] = z*(0.68631 + z*(10.145 + z*(-4.6176 + z*(0.86824 - 0.053741*z))));
        z = y + 4.7171;
        a[5] = 9.0466e-7 + 1.4539e-6*log10(0.015204*(y + 3.4)) + 0.00013253/(z*z) + y*(-4.1228e-7 + 2.2036e-7*y);
        a[6] = -339.45 + 618.73*log10(0.31595*(y + 3.9)) + 250.20/((y + 4.4395)*(y + 4.4395));
        a[7] = -35.105 + y*(36.167 + y*(-9.3575 + 0.33717*y));
        a[8] = 0.17554 + y*(0.37300 + y*(-0.014938 + y*(0.0032314 + 0.0025579*y)));
    	}
      else
      {
        for (i = 0; i < 9; i++)
      			a[i] = 1.e-10;
    	}
    /* calculate the flux due to non-diffractive process for given gamma-ray energy */
    xa3 = x - a[3] + a[2]*pow(x-a[3],2);
    xa8 = x - a[8] + a[6]*pow(x - a[8],2)  +  a[7] * pow(x-a[8],3.);
    sigma = a[0]*exp(-a[1]*xa3*xa3) + a[4]*exp(-a[5]*xa8*xa8) ;
    /* cutoff is the kinematic limit function as in the paper */
    cutoff = (1/(1 + exp(Wl*(Lmin - x))))*(1/(1 + exp(Wh*(x - Lmax))));
    sigma *= cutoff;
    		if (sigma < 0)
         sigma = 0;

            if (Tp <= 1.95)             /* renormalization  */
            {
                pow1 = (y + 3.25)/(1 + 8.08*(y + 3.25));
                r_factor = 3.05*exp(-107.0*pow1*pow1);
            }
            else
            {
                r_factor = 1.01;
            }
    sigma *= r_factor;
    return sigma;
}

double HadronicInteraction::KamaeGammaDiff(double Esec, double Tp) const {
	Esec = Esec/GeV;
	Tp = Tp/GeV;
    double x, y, z1, z2;
    double Wdiff, Lmax;
    double pow1,pow2;
    double sigma;
    double cutoff;
    double b[8];
    int i;             
    /* init some variables, given in table 2 */   
    x = log10(Esec);
    y = log10(Tp*0.001);
    if ((Tp > 1.94) && (Tp < 512000.1))
    	{
        if (Tp > 5.51)
        {
            z1 = y + 0.59913;
            z2 = y + 9.4773;
            b[0] = 60.142*tanh(-0.37555*(y + 2.2)) - 5.9564*z1*z1 + 0.0060162*z2*z2*z2*z2;
            z1 = y + 369.13;
            b[1] = 35.322 + 3.8026*tanh(-2.4979*(y + 1.9)) - 0.00021870*z1*z1;
            z1 = y + 252.43;
            b[2] = -15.732 - 0.082064*tanh(-1.9621*(y + 2.1)) + 0.00023355*z1*z1;
            pow1 = (y + 1.0444)/(1.0 + 0.27437*(y + 1.0444));
            b[3] = -0.086827 + 0.37646*exp(-0.53053*pow1*pow1);
        }
        else
        {
            b[0] = 1.0e-10;
            b[1] = 1.0e-10;
            b[2] = 1.0e-10;
            b[3] = 1.0e-10;
        }
        z1 = y + 2.95;
        pow1 = (y + 2.45) - 0.19717*(y + 2.45)*(y + 2.45);
        b[4] = 2.5982 + 0.39131*z1*z1 - 0.0049693*z1*z1*z1*z1 + 0.94131*exp(-24.347*pow1*pow1);
        z1 = (y - 0.83562)/(1.0 + 0.33933*(y - 0.83562));
        b[5] = 0.11198 + y*(-0.64582 + 0.16114*y) + 2.2853*exp(-0.0032432*z1*z1);
        b[6] = 1.7843 + y*(0.91914 + y*(0.050118 + y*(0.038096 + y*(-0.027334 + y*(-0.0035556 + 0.0025742*y)))));
        z1 = y + 1.8441;
        b[7] = -0.19870 + y*(-0.071003 + 0.019328*y) - 0.28321*exp(-6.0516*z1*z1);
      }
        else
        {
        for (i = 0; i < 8; i++)
            b[i] = 1.0e-10;
        }
              /* calculate sigma due to diffractive process for given gamma-ray energy */
    Lmax = y + 3.0;
    Wdiff = 75.0;

    pow1 = (x - b[2])/(1 + b[3]*(x - b[2]));
    pow2 = (x - b[6])/(1 + b[7]*(x - b[6]));
    sigma = b[0]*exp(-b[1]*pow1*pow1) + b[4]*exp(-b[5]*pow2*pow2);

    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma *= cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;

}

double HadronicInteraction::KamaeGamma1232(double Esec, double Tp) const {
	Esec = Esec/GeV;
	Tp = Tp/GeV;
    double pow1, pow2;
    double x, y, xc2, cutoff, sigma;
    double Wdiff, Lmax;
    int i;
    double c[5];

    x = log10(Esec);
    y = log10(Tp*0.001);
    Lmax = y + 3.0;
    Wdiff = 75.0;

    if ((Tp < 0.488) || (Tp > 1.95))
        {
        for (i = 0; i < 5; i++)
            c[i] = 1.0e-10;
        }
    else
    {
        pow1 = ((y + 3.1301)/(1.0 + 0.14921*(y + 3.1301)));
        c[0] = 2.4316*exp(-69.484*pow1*pow1) - (6.3003 + 9.5349/y - 0.38121*y*y);
        c[1] = 56.872 + y*(40.627 + 7.7528*y);
        c[2] = -5.4918 - 6.7872*tanh(4.7128*(y + 2.1)) + 0.68048*y;
        c[3] = -0.36414 + 0.039777*y;
        c[4] = -0.72807 + y*(-0.48828 - 0.092876*y);
    }

    xc2 = x - c[2];
    pow2 = xc2/(1 + xc2*(c[3] + c[4]*xc2));
    sigma = c[0]*exp(-c[1]*pow2*pow2);

    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma *= cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;
}

double HadronicInteraction::KamaeGamma1600(double Esec, double Tp) const {
	Esec = Esec/GeV;
	Tp = Tp/GeV;
    double pow1, pow2;
    double x, y, xd2, cutoff, sigma;
    double Wdiff, Lmax;
    int i;
    double d[5];

    x = log10(Esec);
    y = log10(Tp*0.001);
    Lmax = y + 3.0;
    Wdiff = 75.0;


    if ((Tp < 0.69) || (Tp > 2.76))
    {
        for (i = 0; i < 5; i++)
            d[i] = 1.e-10;
    }
    else
    {
        pow1 = ((y + 2.9507)/(1.0 + 1.2912*(y + 2.9507)));
        d[0] = 3.2433*exp(-57.133*pow1*pow1) - (1.0640 + 0.43925*y);
        d[1] = 16.901 + y*(5.9539 + y*(-2.1257 - 0.92057*y));
        d[2] = -6.6638 - 7.5010*tanh(30.322*(y + 2.1)) + 0.54662*y;
        d[3] = -1.50648 + y*(-0.87211 - 0.17097*y);
        d[4] = 0.42795 + y*(0.55136 + y*(0.20707 + 0.027552*y));
    }

    xd2 = x - d[2];
    pow2 = xd2/(1 + xd2*(d[3] + d[4]*xd2));
    sigma = d[0]*exp(-d[1]*pow2*pow2);

    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma *= cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;
}
// ------------------------------------ Dermer 1986 ---------------------------------------------------------------------
double HadronicInteraction::DermerXSection(double ePrimary) const {
	double momentum = sqrt(ePrimary*ePrimary/c_squared - mProton*mProton*c_squared) /(GeV/c_light);
	double Tp;
	double sinv;
  	double eta;
	double res;
	double pion_mass = mPion/GeV*c_squared; 
	
	if (momentum < 0.7765)
		{
			res = 1.e-20;
		}
	else if (momentum > 0.7765 && momentum < 0.954)
	{
		Tp = momentum*c_light/GeV;
		sinv = 2.*mProton*(Tp+2*mProton);
		eta =  sqrt( pow((  sinv - pion_mass*pion_mass
				- 4.0 * mProton*mProton),  2.0)
			-  16.0 * pion_mass*pion_mass * mProton*mProton );
		eta /= ( 2.0 * sqrt(sinv) * pion_mass);

		res =   0.032 * eta*eta
			+  0.04 * pow(eta, 6.0)
			+  0.047 * pow(eta, 8.0);
	}

	else if (momentum >= 0.954 && momentum < 1.294)
	{
		res =  32.6 * pow((momentum -  0.8), 3.21);
	}

	else if (momentum >= 1.294 && momentum < 8.0)
	{
		res =  5.4 * pow((momentum -  0.8), 0.81);
	}

	else if (momentum >= 8.0)
	res = 32.0*log(momentum)+(48.5/sqrt(momentum))-59.5;

	return(res*mb);
}

double HadronicInteraction::DermerSpectrumPhoton(double Esec, double TP) const {
  	double res;

      /* Calculates the gamma ray production due to pi0 decay in proton nuclei
 		collisions using the approach of Dermer 86 AA 157, 223.
 		Everything in GeV and per GeV */
	Esec = Esec/GeV;
	TP = TP/GeV;
	std::string log;

	int i;
	int n=200;
	int ipion, j, n6;
	double kinetic_energy_pion;
	double gamma_pion;
	double total_pion_energy;
	double momentum_pi;
	double a, b1, b2, b, c, s;
	double gr, gamma_cm, bs, esmx, u, x;
	double sum, fpi, mimina, mimaxa, arg1, arg2, eta, delta, md;
	double fnca, fncb, tot1, tot2, fac1, fac2, gds, bds, gdp, bdp;
	double gdm, bdm, gpip, g1,g2, g3, g4, fnc1, fnc2;
	double sig, fp1, fp2, dsig, fac, epi0, qgam, depi;
	double det, sn, tpn, gr1;
	double qpi;
	double bpip, coss;

	/* redefinition of constants for compability */
	double mpi = mPion/(GeV/c_squared);
	double mp = mProton/GeV*c_squared;
	double PI = M_PI;
	double mio = 1.232e0;
	double gamma = 0.0575e0;


	int NUMBER_OF_PION_BINS = 101;
	int GAUSS_LAGUERRE = 9;
	double pion_energy[101], qpion[101];
	/* calculate pion source function
		n is # of divisions in a Simpson's routine for Stecker's model
	//    n5 is the # of divisions per decade in proton KE integration */

				/*  Loop over pion energy */
	for (ipion=0; ipion<NUMBER_OF_PION_BINS; ipion++)            
	{
		pion_energy[ipion] = pow( 10.0, ipion / 10.0 -  3.0);   /* 0.001 < Tpi(GeV) < 1e7 */

			/* Calculate pion kinetic energy, gamma, and momentum */
		kinetic_energy_pion = pion_energy[ipion];

		gamma_pion = kinetic_energy_pion / mpi +  1.0;
		total_pion_energy = kinetic_energy_pion + mpi;
		momentum_pi = sqrt(total_pion_energy*total_pion_energy - mpi*mpi);

			/* Kinematics */
		a =   mp*mp + mpi*mpi - 2.0 * gamma_pion * mp * mpi;
		b1 =  4.0 * mp*mp - mpi*mpi;
		b2 =    2.0 * mp*mp * mpi*mpi * (gamma_pion*gamma_pion -  1.0);
		b =    -2.0 * ( (mp*mp - gamma_pion * mp*mpi) * b1  -  b2);
		c = mp*mp * b1*b1;

			det = sqrt(b*b -  4.0 * a * c);
		sn = (det + b) /  -2.0 / a;
		tpn = sn /  2.0 / mp -  2.0 * mp;
		gr1 = (tpn + mp)/mp;
		/* gr1 is minimum proton Lorentz factor to make a pion of energy
				kinetic_energy_pion */
		/*  Evaluate at TP, eqn.(1) Dermer (1986) */
			gr = (TP/mp) + 1.0;
		if (gr < gr1)
			{
				dsig = 1.e-40;
			}
		else
		{
			gamma_cm = sqrt((gr +  1.0) /  2.0);
			bs = sqrt( 1.0 -  1.0 / gamma_cm / gamma_cm);
			s =  2.0 * mp * (TP +  2.0 * mp);
			sum =  0.0;

		if (TP > 3.0)   /**********  TP > 3.0 **********/
		{
		/*
		Calculate dsig/de lab system pion energy spectrum for SB model
		Calculate lower integration limit in eqn. (11) Dermer (1986), sum=dN/dT_pi0 */
			esmx =   (s -  4.0 * mp*mp + mpi*mpi) / ( 2.0 * sqrt(s));
			u =   (   gamma_cm *  total_pion_energy  -  esmx) / (   bs *  gamma_cm *  momentum_pi);
				if (u <  1.0) {
				u = fmax(u,  -1.0);
				for (i=0; i<GAUSS_LAGUERRE; i++)
					{
						x =  c1[i];
						coss =  1.0-((  1.0 - u) * expl(-x));
						sum += ( 1.0 -  u ) *  d1[i] *  pionb_pi0(gr, total_pion_energy, coss);
					}

				/*        for (i = 0; i < 1000; ++i)
				{
				coss = u + i*(1.0-u)/1000.;
				sum += (( 1.0 -  u)/1000.)
					*  pionb_pi0( gr, total_pion_energy, coss, mpi);
				}               */


					/* sum = (dN/dT_pi) <zeta sigma_pi> according to SB */
				}
				sum *=  2.0 * M_PI * momentum_pi;
		}

		fpi =  0.0;
		/* Calculate dsig/de lab system pion energy spectrum for Stecker's model */
		if (TP < 7.0)             /**********  TP < 7.0 **********/
		{
			/* Dermer (1986) eqn. (8) integration limits	*/
			mimina = mp + mpi;
			mimaxa = sqrt(s) - mp;
				/* Calculate w_r(T_p) eqn. (10) Dermer (1986) */
			arg1 = (mimaxa - mio) / gamma;
			arg2 = (mimina - mio) / gamma;
			eta = PI / (atan(arg1) - atan(arg2));
				/* Initialize integration parameters */
			delta = (mimaxa - mimina) /  n;
			md = mimina;
			fnca = (fncb = (tot1 = (tot2 =  0.0)));
			fncb =  0.0;
			tot1 =  0.0;
			tot2 =  0.0;
			n6 = n-1;
				/*
				g=gamma, b=beta, ds=delta_star, pip=pion_prime, dp=delta_plus, dm=delta_minus
				*/
			for (j=0; j<n6; j++)
			{
				fac1 =  0.0;
				fac2 =  0.0;
				md = md + delta;
				gds =   (s + md*md - mp*mp)/ ( 2.0 * sqrt(s) * md);
				bds = sqrt( 1.0 -  1.0 / gds / gds);
				gpip =   (md*md + mpi*mpi - mp*mp)/ ( 2.0 * md * mpi);
				bpip = sqrt( 1.0 -  1.0 / gpip / gpip);
				gdp = gamma_cm * gds * (  1.0 + bs * bds);
				bdp = sqrt(  1.0 - 1.0 / gdp / gdp);
				gdm = gamma_cm * gds * ( 1.0 - bs * bds);
				bdm = sqrt( 1.0 -  1.0 / gdm / gdm);
				/* Define step function limits in eqn. (7) Dermer (1986) */
				g1 = gdp * gpip * ( 1.0 - bdp * bpip);
				g2 = gdp * gpip * ( 1.0 + bdp * bpip);
				g3 = gdm * gpip * ( 1.0 - bdm * bpip);
				g4 = gdm * gpip * ( 1.0 + bdm * bpip);
				if ((gamma_pion >= g1) && (gamma_pion <= g2))
				{
					fac1 =  1.0;
				}
				if ((gamma_pion >= g3) && (gamma_pion <= g4))
				{
					fac2 =  1.0;
				}
				/* Integrate eqn. (8)	*/
				fnc1 =   (fac1 * bw(md))/ ( 4.0 * bdp * gdp * bpip * gpip);
				fnc2 = (fac2 * bw(md)) / ( 4.0 * bdm * gdm * bpip * gpip);
				tot1 += ((fnca + fnc1) * delta /  2.0);
				tot2 += ((fncb + fnc2) * delta /  2.0);
				fnca = fnc1;
				fncb = fnc2;
			}
			double ePrimary = sqrt(TP*GeV*TP*GeV + mProton*mProton*c_squared*c_squared);
			sig = DermerXSection(ePrimary);
			fp1 = (eta * tot1 * sig) / mpi;
			fp2 = (eta * tot2 * sig) / mpi;
				/*  fpi = (dN/dT_pi) <zeta sigma_pi> according to Stecker */
			fpi = fp1 + fp2;
		}
		/*
		Piece the regimes together
		*/
		if (TP > 7.0)
		{
				dsig=sum;
		}
		else if (TP < 3.0)
		{
				dsig=1.0*fpi;                                  /************/
		}
		else if (TP >= 3.0 && TP <= 7.0)
		{
				dsig=(((TP -  3.0) * sum) + (( 7.0 - TP) * 1.0*fpi))/  4.0;     /************/
		}
	}
		
	qpion[ipion] = dsig; 
	}//xyclose
	/********** End pion source function
	 Do integral over epion to get gamma source, eqn. (2) Dermer (1986) *********/
	fac = pow( 10.0,  0.01);
		epi0 = Esec + mpi*mpi /  4.0 / Esec;
		qgam =  0.0;

		for (i=0; i<150; i++)        //xx
		{
			total_pion_energy = epi0 * (fac +  1.0) /  2.0;
			depi = epi0 * (fac -  1.0);
			momentum_pi = sqrt(  total_pion_energy*total_pion_energy
					- mpi*mpi);
			kinetic_energy_pion = total_pion_energy - mpi;
			ipion  = (int) floor(   30.0001 +  10.0 * log10(kinetic_energy_pion));

			if (ipion < 0)
				{
						qpi =  qpion[0]
						+  (qpion[1] - qpion[0])
					* (kinetic_energy_pion - pion_energy[0])
					/ (pion_energy[1] - pion_energy[0]);
				}
					else if (ipion >= NUMBER_OF_PION_BINS)
				{
						qpi =   qpion[NUMBER_OF_PION_BINS-2]
						+   (  qpion[NUMBER_OF_PION_BINS-1]
						- qpion[NUMBER_OF_PION_BINS-2])
					* (  kinetic_energy_pion
						- pion_energy[NUMBER_OF_PION_BINS-2])
					/ (  pion_energy[NUMBER_OF_PION_BINS-1]
						- pion_energy[NUMBER_OF_PION_BINS-2]);
				}
					else
				{
						qpi = ( qpion[ipion+1] * (kinetic_energy_pion - pion_energy[ipion])
					+ qpion[ipion] * (pion_energy[ipion+1] - kinetic_energy_pion) ) /
						(pion_energy[ipion+1] - pion_energy[ipion]);
				}

				qgam +=  2.0 * depi / momentum_pi * qpi;    // 2 photons per pi0
				epi0 *= fac;
		}//xxclose

	res = Esec * qgam;

 	return(res);	
}

double HadronicInteraction::pionb_pi0(double proton_gamma, double total_pion_energy, double cos_theta_max) const
{
  /*
    Invariant x-section E D3sig/d3p Stephens and Badhwar 1981, eqn(11)
    Dermer 1986
  */
  double a=140.0;
  double b=5.43;
  double c1=6.1;
  double c2=-3.3;
  double c3=0.6;
  double s;
  double momentum_pion;
  double gamma_cm; /* Lorentz factor of the center of mass wrt the lab frame */
  double beta_cm;
  double center_of_momentum_energy; /* energy of particle in the CM frame */
  double ps, total_proton_energy, sns, tth, cth, sth;
  long double p_perpendicular, q, esmax, psmax, xs1, xs2, x;
  long double fep, ex, u;
  double mpi = mPion/GeV*c_squared;
  double mp = mProton/GeV*c_squared;
  
  s = 2.0 * mp*mp * (  proton_gamma  +  1.0);
  momentum_pion = sqrt(  total_pion_energy*total_pion_energy
		       - mpi*mpi);
  gamma_cm = sqrt((proton_gamma +  1.0) /  2.0);
  beta_cm = sqrt( 1.0 -  1.0 / gamma_cm / gamma_cm);
  center_of_momentum_energy = gamma_cm * (   total_pion_energy
					   - (
					      beta_cm
					      * momentum_pion
					      * cos_theta_max
					     )
				         );
  ps = sqrt(  center_of_momentum_energy*center_of_momentum_energy
	    - mpi*mpi);
  total_proton_energy = mp * proton_gamma;

  sns = sqrt( 1.0 - cos_theta_max*cos_theta_max);
  tth = sns / (  gamma_cm
	       * (
		    cos_theta_max
		  - beta_cm * total_pion_energy / momentum_pion
		 )
	      );

  cth = sqrt( 1.0 / (tth*tth +  1.0));
  sth = sqrt( 1.0 - cth*cth);

  p_perpendicular = ps * sth;
  u =  1.0 +  4.0 * mp*mp / s;
  q = c1 + p_perpendicular * (c2 + c3 * p_perpendicular);
  q /= sqrt(u);

  esmax =   (s + mpi*mpi -  4.0 * mp*mp)
          /  2.0 / sqrt(s);
  psmax = sqrt(esmax*esmax - mpi*mpi);

  xs1 = ps * cth / psmax;
  xs2 =    4.0
        * (p_perpendicular*p_perpendicular + mpi*mpi) / s;
  x = sqrt(xs1*xs1 + xs2);

  if (x >= 1.0)
  {
    return 0.0;
  }

  fep = pow(( 1.0 -  4.0 * mp*mp / s),
	     2.0);
  fep *= (   1.0
	  +  23.0 / pow(total_proton_energy,  2.6));

  ex = exp(-b * p_perpendicular / u);

  return a * fep * ex * pow(( 1.0 - x), q);
}

double HadronicInteraction::bw(double md) const
{
  /* Calculates Breit-Wigner Distribution, eqn(9) Dermer 1986 */
	double mio = 1.232e0;
	double gamma = 0.0575e0;
  	double bw1 ;
  	bw1 = (md - mio);
  	bw1 = (bw1 * bw1) + gamma*gamma;
  	return (gamma /  M_PI / bw1);
}      

// ------------------------------------ general Function ----------------------------------------------------------------
// general Function are called in process and performInteraction. Here is seted which model is used.

double HadronicInteraction::xSection(double ePrimary) const {
	switch (flag_Function)
	{
	case 0:
		return KellnerXSection(ePrimary);
		break;
	case 1:
		return KamaeXSection(ePrimary);
		break;
	case 2:
		return DermerXSection(ePrimary);
		break;
	case 3:
		if(ePrimary < EnergySplit) {
			return KamaeXSection(ePrimary);
		}
		else{
			return KellnerXSection(ePrimary);
		}
		break;
	case 4:
		if(ePrimary < EnergySplit) {
			return DermerXSection(ePrimary);
		}
		else{
			return KellnerXSection(ePrimary);
		}
		break;
	case 5:
		if(ePrimary < 20*GeV) {
			return DermerXSection(ePrimary);
		}
		if(ePrimary < EnergySplit) {
			return KamaeXSection(ePrimary);
		}
		else{
			return KellnerXSection(ePrimary);
		}
		break;
	default:
		KISS_LOG_WARNING << "flag for Funiction not understood \nreturn 0";
		return 0;
		break;
	}
}

double HadronicInteraction::spectrumPion(double x, double ePrimary) const {
	return 0;
}

double HadronicInteraction::spectrumPhoton(double x, double ePrimary) const {
	switch (flag_Function)
	{
	case 0:
		return KellnerSpectrumPhoton(x, ePrimary);
		break;
	case 1:
		return KamaeSpectrumPhoton(x*ePrimary, ePrimary);
		break;
	case 2:
		return DermerSpectrumPhoton(x*ePrimary, ePrimary);
		break;
	case 3:
		if(ePrimary < EnergySplit) {
			return KamaeSpectrumPhoton(x*ePrimary, ePrimary);
		}
		else{
			return KellnerSpectrumPhoton(x,ePrimary);
		}
		break;
	case 4:
		if(ePrimary < EnergySplit) {
			return DermerSpectrumPhoton(x*ePrimary, ePrimary);
		}
		else{
			return KellnerSpectrumPhoton(x,ePrimary);
		}
		break;
	case 5:
		if(ePrimary < 20*GeV) {
			return DermerSpectrumPhoton(x*ePrimary, ePrimary);
		}
		if(ePrimary < EnergySplit) {
			return KamaeSpectrumPhoton(x*ePrimary, ePrimary);
		}
		else {
			return KellnerSpectrumPhoton(x,ePrimary);
		}
		break;
	default:
		KISS_LOG_WARNING << "flag for Funiction not understood \nreturn 0";
		return 0;
		break;
	}
	
}

// ----------------------------------- process and performInteraction ----------------------------------------------------
void HadronicInteraction::process(Candidate *candidate) const {

	const int id = candidate->current.getId();
	if (!isNucleus(id))
		return;

	const double ePrimary = candidate->current.getEnergy();
	if (ePrimary < 1. * GeV)
		return;

	const double _xSection = xSection(ePrimary);
	const double meanFreePath = 1. / (_xSection * massDensity->getNucleonDensity(candidate->current.getPosition()));
	const double stepLength = candidate->getCurrentStep();
	const double interactionProbability = stepLength / meanFreePath;

	Random &random = Random::instance();
	if (random.rand() > interactionProbability) {
		candidate->limitNextStep(this->limit * meanFreePath);
		return;
	}

	performInteraction(candidate);
}

void HadronicInteraction::performInteraction(Candidate *candidate) const {
	std::vector<int> outPartID;
	std::vector<double> outPartE;
	double eAvailable = candidate->current.getEnergy();

	Random &random = Random::instance();
	const double startPiont = random.rand();
	bool doPhoton = (startPiont >= 0. && startPiont < 0.5);
	bool doPiCharged = (startPiont >= 0.5 && startPiont <= 1.);
	bool donePhoton = false;
	bool donePiCharged = false;

	do {
		if (doPhoton && !donePhoton) {
			const double xMin = std::min(1 * GeV / eAvailable, 1e-3);
			const double xMax = 1.;
			int nPhoton = std::round(gaussInt([this, eAvailable](double x) { return this->spectrumPhoton(x, eAvailable); }, xMin, xMax));
			for (int i = 0; i < nPhoton; ++i) {
				const double ePhoton = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPhoton(x, eAvailable); }, eAvailable);
				if (eAvailable >= ePhoton) {
					eAvailable -= ePhoton;
					// create photons if allowed
					if (havePhotons) {
						outPartID.push_back(22);
						outPartE.push_back(ePhoton);
					}
				} else {
					break;
				}
			}
			doPiCharged = true;
			donePhoton = true;
		}

		if (doPiCharged && !donePiCharged) {
			const double xMin = 0.;
			const double xMax = 1.;
			int nPiCharged = 2 * std::round(gaussInt([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, xMin, xMax));
			for (int i = 0; i < nPiCharged / 2; ++i) {
				const double ePiPlus = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, eAvailable);
				const double ePiMinus = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, eAvailable);
				if (eAvailable >= ePiPlus + ePiMinus) {
					eAvailable -= ePiPlus;
					eAvailable -= ePiMinus;
					// create pions if allowed
					if (havePions) {
						outPartID.push_back(211);
						outPartE.push_back(ePiPlus);
						outPartID.push_back(-211);
						outPartE.push_back(ePiMinus);
					}
				} else {
					break;
				}
			}
			doPhoton = true;
			donePiCharged = true;
		}
	} while (!donePhoton || !donePiCharged);

	candidate->current.setEnergy(eAvailable);

	const Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	for (int i = 0; i < outPartID.size(); ++i) {
		candidate->addSecondary(outPartID[i], outPartE[i], pos);
	}
}

} // namespace CRPropa
