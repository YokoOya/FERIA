#include "feria.h"

/*//////////////////////////////////////
 The values are expressed in units of
 cm, g, sec, radian, Kelvin
 /*//////////////////////////////////////

/*/////parameters
public:
	char outputfilename[lstr], overwrite, name_line[lstr], name_object[lstr], radesys[lstr], centRa_hms[lstr], centDec_dms[lstr];
	double Fldres_as, Velres_kmps, distance_pc, mass_msun, rCB_au, inc_deg, PA_deg, rot_sign, Rout_au, Rin_au;
	double HeightIRE_au, FlareIRE_deg, HeightKEP_au, FlareKEP_deg;
	double H2densityCB, densityProfileIRE, densityProfileKEP, fractionalDensity;
	double TempCB, tempProfileIRE, tempProfileKEP;
	double linewidth_cmps, beam_maj_as, beam_min_as, beam_pa_deg;
	double restfreq_Hz, centRa[3], centDec[3], centRa_deg, centDec_deg, Vsys_cmps;
/*/////

sourceParams::sourceParams(char *outputfilename_in, double Fldres_in, double Velres_in, double distance_pc_in, double mass_msun_in, double rCB_au_in, double inc_deg_in, double PA_deg_in, double rot_sign_in, double Rout_au_in, double Rin_au_in, double HeightIRE_au_in, double FlareIRE_deg_in, double HeightKEP_au_in, double FlareKEP_deg_in, double H2densityCB_in, double densityProfileIRE_in, double densityProfileKEP_in, double fractionalDensity_in, double TempCB_in, double tempProfileIRE_in, double tempProfileKEP_in, double linewidth_kmps_in, double beam_maj_as_in, double beam_min_as_in, double beam_pa_deg_in, char* name_line_in, double restfreq_GHz_in, char* name_object_in, char* radesys_in, char* center_ra_in, char* center_dec_in, double Vsys_kmps_in) {
    
    strcpy(outputfilename, outputfilename_in);
    Fldres_as = Fldres_in;
    Velres_kmps = Velres_in;
    distance_pc = distance_pc_in;
    mass_msun = mass_msun_in;
    rCB_au = rCB_au_in;
    inc_deg = inc_deg_in;
    PA_deg = PA_deg_in;
    rot_sign = rot_sign_in;
    Rout_au = Rout_au_in;
    Rin_au = Rin_au_in;
    HeightIRE_au = HeightIRE_au_in;
    FlareIRE_deg = FlareIRE_deg_in;
	HeightKEP_au = HeightKEP_au_in;
	FlareKEP_deg = FlareKEP_deg_in;
    H2densityCB = H2densityCB_in;
    densityProfileIRE = densityProfileIRE_in;
	densityProfileKEP = densityProfileKEP_in;
    fractionalDensity = fractionalDensity_in;
    TempCB = TempCB_in;
	tempProfileIRE = tempProfileIRE_in;
	tempProfileKEP = tempProfileKEP_in;
    linewidth_cmps = linewidth_kmps_in * CMperKM;
    beam_maj_as = beam_maj_as_in;
    beam_min_as = beam_min_as_in;
    beam_pa_deg = beam_pa_deg_in;
    strcpy(name_line, name_line_in);
    restfreq_Hz = restfreq_GHz_in * 1e9;
    strcpy(name_object, name_object_in);
    Vsys_cmps = Vsys_kmps_in * CMperKM;
    
	strcpy(radesys, radesys_in);
    strcpy(centRa_hms, center_ra_in);
    strcpy(centDec_dms, center_dec_in);
    sscanf(centRa_hms, "%lfh%lfm%lfs", &centRa[0], &centRa[1], &centRa[2]);
    sscanf(centDec_dms, "%lfd%lfm%lfs", &centDec[0], &centDec[1], &centDec[2]);
    if (centRa[0] < 0) for (int i = 1; i < 3; ++i) centRa[i] = -centRa[i];
    if (centDec[0] < 0) for (int i = 1; i < 3; ++i) centDec[i] = -centDec[i];
    centRa_deg = ((centRa[2] / 60. + centRa[1]) / 60. + centRa[0]) / 24. * 360.;
    centDec_deg = (centDec[2] / 60. + centDec[1]) / 60. + centDec[0];
	
	centRa_scaled_as = centRa_deg * 3600.;
	centDec_as = centDec_deg * 3600.;
	
}


sourceParams::~sourceParams() {
}

