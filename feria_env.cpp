#include "feria.h"

/*//////////////////////////////////////
The values are expressed in units of 
    cm, g, sec, radian, Kelvin
/*//////////////////////////////////////


/*/////parameters
public:
	double mass_grav, rCB, inc, rot, PA, Rout, Rin, Height, Flare, tanFlare;
	double densCB, densProfIRE, densProfKEP, TempCB, tempProfIRE, tempProfKEP;
	double velrot_CB_xrCB;
/*/////

ireModel::ireModel(double mass_msun, double rCB_au, double inc_deg, double PA_deg, double rot_sign, double Rout_au, double Rin_au, double HeightIRE_au, double FlareIRE_deg, double HeightKEP_au, double FlareKEP_deg, double DensCB_in, double densityProfileIRE, double densityProfileKEP, double TempCB_in, double tempProfileIRE, double tempProfileKEP) {
    mass_grav = mass_msun * Msun * Grav;
    rCB = rCB_au * CMperAU;
    inc = inc_deg / DEGperRAD;
    PA = PA_deg / DEGperRAD;
    rot = rot_sign;
    Rout = Rout_au * CMperAU;
    Rin = Rin_au * CMperAU;
	
    HeightIRE = HeightIRE_au * CMperAU;
    tanFlareIRE = tan(FlareIRE_deg / DEGperRAD / 2.);
	HeightKEP = HeightKEP_au * CMperAU;
	tanFlareKEP = tan(FlareKEP_deg / DEGperRAD / 2.);

	DensCB = DensCB_in;
    densProfIRE = densityProfileIRE;
	densProfKEP = densityProfileKEP;
	
    TempCB = TempCB_in;
	tempProfIRE = tempProfileIRE;
	tempProfKEP = tempProfileKEP;

	velrot_CB_xrCB = sqrt(2. * mass_grav * rCB) * rot;
}


int ireModel::which_velStr(double dist_protostar, double r, double theta, double z) {
	if (dist_protostar > Rout || dist_protostar < Rin || dist_protostar <= EPS) return noGAS;
	//if (r > Rout || r < Rin || r <= EPS) return noGAS; // for vertical edge

	if (dist_protostar >= rCB) {
		if (fabs(z) < r * tanFlareIRE + HeightIRE / 2.) return isIRE;
	} else {
		if (fabs(z) < r * tanFlareKEP + HeightKEP / 2.) return isKEP;
	}
	
	return noGAS;
}


void ireModel::env_vel(double dist_protostar, double r, double theta, double z, double *vel_polar, int velStrID) {
    if (velStrID == noGAS) {
        for (int i = 0; i < naxis; ++i) vel_polar[i] = 0.;
        return ;
	}
	
	if (velStrID == isKEP) {
		vel_polar[tAXIS] = sqrt(mass_grav / dist_protostar) * rot;
		vel_polar[rAXIS] = 0.;
		vel_polar[ZAXIS] = 0.;
		return ;
	}
	
    vel_polar[tAXIS] = velrot_CB_xrCB / dist_protostar;
	double vel_inf2 = 2. * mass_grav / dist_protostar - pow(vel_polar[tAXIS], 2), vel_inf;
	if (vel_inf2 > 0) {
		vel_inf = sqrt(vel_inf2);
	} else {
		vel_inf = 0.;
	}
    vel_polar[rAXIS] = -vel_inf * r / dist_protostar;
    vel_polar[ZAXIS] = -vel_inf * z / dist_protostar;
	
    return ;
}

double ireModel::env_dens(double dist_protostar, double r, double theta, double z, int velStrID) {
    double ans = 0.;
    //*Edit below if you want to calculate the density. //
	
    if (velStrID == isIRE) ans = DensCB * pow(dist_protostar / rCB, densProfIRE);
	if (velStrID == isKEP) ans = DensCB * pow(dist_protostar / rCB, densProfKEP);

    //Edit above if you want to calculate the density. */
    return ans;
}

double ireModel::env_temp(double dist_protostar, double r, double theta, double z, int velStrID) {
    double ans = 0.;
    //*Edit below if you want to calculate the temperature. //
    
	if (velStrID == isIRE) ans = TempCB * pow(dist_protostar / rCB, tempProfIRE);
	if (velStrID == isKEP) ans = TempCB * pow(dist_protostar / rCB, tempProfKEP);

    //Edit above if you want to calculate the temperature. */
    return ans;
}


