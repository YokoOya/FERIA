#include "feria.h"

/*//////////////////////////////////////
 The values are expressed in units of
 cm, g, sec, radian, Kelvin
 /*//////////////////////////////////////

/*/////parameters
public:
	double **emission, *emission_fits;
 	int npix_PV;
	double crpix[naxis_PV], crval[naxis_PV], cdelt[naxis_PV], maxval, minval;
	double beam_maj_as, beam_min_as, beam_maj_rad, beam_min_rad, beam_maj_deg, beam_min_deg, beam_maj_cm, beam_min_cm, beam_pa_rad, beam_pa_deg, beam_cpa, beam_spa;
	double cent_offset_ra_au, cent_offset_dec_au,  PV_PA_deg, PV_PA_rad;
/*/////


PVdiagram::PVdiagram(int npix_in, double pvCrpix[naxis_PV], double pvCrval[naxis_PV], double pvCdelt[naxis_PV], sourceParams &source, double cent_offset_ra_auin, double cent_offset_dec_auin, double PV_PA_degin) {
	for (int i = 0; i < naxis_PV; ++i) {
		crpix[i] = pvCrpix[i];
		crval[i] = pvCrval[i] * (i < naxis_PV - 1 ? CMperAU : CMperKM);
		cdelt[i] = pvCdelt[i] * (i < naxis_PV - 1 ? CMperAU : CMperKM);
	}
	
	npix_PV = npix_in;
	
	maxval = -INF;
	minval = INF;
	
	beam_maj_as = source.beam_maj_as;
	beam_min_as = source.beam_min_as;
	beam_maj_rad = beam_maj_as / ASperRAD;
	beam_maj_deg = beam_maj_rad * DEGperRAD;
	beam_min_rad = beam_min_as / ASperRAD;
	beam_min_deg = beam_min_rad * DEGperRAD;
	beam_pa_deg = source.beam_pa_deg;
	beam_pa_rad = beam_pa_deg / DEGperRAD;
	
	emission = (double**) malloc(sizeof(double*) * npix_PV);
	emission[0] = (double*) malloc(sizeof(double) * npix_PV * nvel);
	emission_fits = (double*) malloc(sizeof(double) * npix_PV * nvel);
	
	for (int ip = 0; ip < npix_PV; ++ip) {
		emission[ip] = emission[0] + nvel * ip;
		for (int iv = 0; iv < nvel; ++iv) {
			emission[ip][iv] = 0.;
			emission_fits[ip * nvel + iv] = 0.;
		}
	}
	
	cent_offset_ra_au = cent_offset_ra_auin;
	cent_offset_dec_au = cent_offset_dec_auin;
	PV_PA_deg = PV_PA_degin;
	PV_PA_rad = PV_PA_degin / DEGperRAD;
	
	
	double theta_rad = PV_PA_rad - beam_pa_rad;
	beam_slice_as = 1. / sqrt(pow(sin(theta_rad) / beam_min_as, 2) + pow(cos(theta_rad) / beam_maj_as, 2));
	beam_slice_rad = beam_slice_as / ASperRAD;
	beam_slice_deg = beam_slice_rad * DEGperRAD;
	
	
	centcube_ra_scaled_as = source.centRa_scaled_as;
	centcube_dec_as = source.centDec_as;
	centpv_dec_as = centcube_dec_as + double(cent_offset_dec_au) / source.distance_pc;
	centpv_ra_scaled_as = centcube_ra_scaled_as + double(cent_offset_ra_au) / source.distance_pc / cos(pi / 180. * (centcube_dec_as + centpv_dec_as) / 2. / 3600.);
	
}

PVdiagram::~PVdiagram() {
	free(emission[0]);
	free(emission);
	free(emission_fits);
}


double PVdiagram::getIndex(int axis, double pos) {
	return crpix[axis] + (pos * (axis == VAXIS_PV ? CMperKM : CMperAU) - crval[axis]) / cdelt[axis];
}

double PVdiagram::getPos(int axis, double index) {
	return (crval[axis] + cdelt[axis] * (index - crpix[axis])) / (axis == VAXIS_PV ? CMperKM : CMperAU);
}

void PVdiagram::getIndexList(double *indList, double *posList) {
	for (int i = 0; i < naxis_PV; ++i) indList[i] = getIndex(i, posList[i]);
	return ;
}

void PVdiagram::getPosList(double *indList, double *posList) {
	for (int i = 0; i < naxis_PV; ++i) posList[i] = getPos(i, indList[i]);
	return ;
}



void PVdiagram::setMaxMin() {
	maxval = -INF;
	minval = INF;
	double val;
	for (int ip = 0; ip < npix_PV; ++ip) for (int iv = 0; iv < nvel; ++iv) {
		val = emission[ip][iv];
		if (maxval < val) maxval = val;
		if (minval > val) minval = val;
	}
	return ;
}


void PVdiagram::getMaxMin(double *max, double *min) {
	setMaxMin(); 
	if (max) *max = maxval;
	if (min) *min = minval;
	return ;
}

double PVdiagram::getSum() {
	double ans = 0;
	for (int ip = 0; ip < npix_PV; ++ip) for (int iv = 0; iv < nvel; ++iv) ans += emission[ip][iv];
	return ans;
}



void PVdiagram::setVal(int ip, int iv, double val) {
	emission[ip][iv] = val;
	return ;
}

double PVdiagram::getVal(int ip, int iv) {
	return emission[ip][iv];
}




void PVdiagram::setEmissionFits() {
	for (int iv = 0; iv < nvel; ++iv) for (int ip = 0; ip < npix_PV; ++ip) emission_fits[iv * npix_PV + ip] = emission[npix_PV - 1 - ip][iv];
	
	return ;
}



