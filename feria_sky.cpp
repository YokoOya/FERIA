#include "feria.h"

/*//////////////////////////////////////
 The values are expressed in units of
 cm, g, sec, radian, Kelvin
 /*//////////////////////////////////////

/*/////parameters
public:
	double ***emission, *emission_fits;
	double crpix[naxis], crval[naxis], cdelt[naxis], maxval, minval;
	double linewidth_FWHM_cmps, beam_maj_as, beam_min_as, beam_maj_rad, beam_min_rad, beam_maj_deg, beam_min_deg, beam_maj_cm, beam_min_cm, beam_pa_rad, beam_pa_deg, beam_cpa, beam_spa;
	double du, dv, dw;
	dcomp beam_uv[npix<<1][npix<<1], lw_uv[nvel<<1];
/*/////


skyPlane::skyPlane(double distance_pc, double skyCrpix[naxis], double skyCrval[naxis], double skyCdelt[naxis], double linewidth_kmps, double beam_maj_asin, double beam_min_asin, double beam_pa_degin) {
    for (int i = 0; i < naxis; ++i) {
        crpix[i] = skyCrpix[i];
        crval[i] = skyCrval[i] * (i == VAXIS ? CMperKM : CMperAU);
        cdelt[i] = skyCdelt[i] * (i == VAXIS ? CMperKM : CMperAU);
    }
    maxval = -INF;
    minval = INF;
    
    linewidth_FWHM_cmps = linewidth_kmps * CMperKM;
	beam_maj_as = beam_maj_asin;
	beam_min_as = beam_min_asin; 
    beam_maj_rad = beam_maj_as / ASperRAD;
	beam_maj_deg = beam_maj_rad * DEGperRAD;
    beam_min_rad = beam_min_as / ASperRAD;
	beam_min_deg = beam_min_rad * DEGperRAD;
    beam_maj_cm = beam_maj_as * distance_pc * CMperAU;
    beam_min_cm = beam_min_as * distance_pc * CMperAU;
	beam_pa_deg = beam_pa_degin;
    beam_pa_rad = beam_pa_deg / DEGperRAD;
	beam_cpa = cos(beam_pa_rad);
	beam_spa = sin(beam_pa_rad);
    
    du = 1. / cdelt[XAXIS] / npix / 2.;
    dv = 1. / cdelt[YAXIS] / npix / 2.;
    dw = 1. / cdelt[VAXIS] / nvel / 2.;
	
    make_beamuv();
    make_linewidthuv();
	
	
	emission = (double***) malloc(sizeof(double**) * npix);
	emission[0] = (double**) malloc(sizeof(double*) * npix * npix);
	emission[0][0] = (double*) malloc(sizeof(double) * npix * npix * nvel);
	
	emission_fits = (double*) malloc(sizeof(double) * npix * npix * nvel);

	for (int ix = 0; ix < npix; ++ix) {
		emission[ix] = emission[0] + npix * ix;
		for (int iy = 0; iy < npix; ++iy) {
			emission[ix][iy] = emission[0][0] + npix * nvel * ix + nvel * iy;
			for (int iv = 0; iv < nvel; ++iv) {
				emission[ix][iy][iv] = 0.;
				emission_fits[(ix * npix + iy) * nvel + iv] = 0.;
			}
		}
	}
	
}

skyPlane::~skyPlane() {
	
	free(emission[0][0]);
	free(emission[0]);
	free(emission);

	free(emission_fits);

}




void skyPlane::projection(int ix, int iy, double (*gasColumn)[ndata]) {
    double vel, iv, emit;
    int iv_low, iv_upp;
    for (int iz = 0; iz < npix; ++iz) {
		
        vel = gasColumn[iz][VzVAL];
        iv = (vel - crval[VAXIS]) / cdelt[VAXIS] + crpix[VAXIS];
		
		if ((iv + 1) * (iv - nvel) > 0) continue;
		iv_low = floor(iv); iv_upp = ceil(iv);
		if (iv_low == iv_upp) iv_upp += 1;
		
		
		//*Edit below if you want to calculate the emissivity. //
		
		emit = gasColumn[iz][NVAL]; // * gasColumn[iz][TVAL];
		
		//*Edit above if you want to calculate the emissivity. */

        if (iv_low >= 0) emission[ix][iy][iv_low] += emit * (iv_upp - iv);
        if (iv_upp <= nvel - 1) emission[ix][iy][iv_upp] += emit * (iv - iv_low);
    }
    return ;
}


void skyPlane::setCrval(int axis, double val) {
	crval[axis] = val * (axis == VAXIS ? CMperKM : CMperAU);
	return ;
}

double skyPlane::getIndex(int axis, double pos) {
	return crpix[axis] + (pos * (axis == VAXIS ? CMperKM : CMperAU) - crval[axis]) / cdelt[axis];
}

double skyPlane::getPos(int axis, double index) {
	return (crval[axis] + cdelt[axis] * (index - crpix[axis])) / (axis == VAXIS ? CMperKM : CMperAU);
}


void skyPlane::getIndexList(double *indList, double *posList) {
	for (int i = 0; i < naxis; ++i) indList[i] = getIndex(i, posList[i]);
	return ;
}

void skyPlane::getPosList(double *indList, double *posList) {
	for (int i = 0; i < naxis; ++i) posList[i] = getPos(i, indList[i]);
	return ;
}




void skyPlane::setMaxMin() {
    maxval = -INF;
    minval = INF;
	double val;
    for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) for (int iv = 0; iv < nvel; ++iv) {
        val = emission[ix][iy][iv];
        if (maxval < val) maxval = val;
        if (minval > val) minval = val;
    }
    return ;
}

void skyPlane::getMaxMin(double *max, double *min) {
	setMaxMin(); 
    if (max) *max = maxval;
    if (min) *min = minval;
    return ;
}

double skyPlane::getSum() {
    double ans = 0;
    for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) for (int iv = 0; iv < nvel; ++iv) ans += emission[ix][iy][iv];
    return ans; 
}

double skyPlane::getVal(double *indList) {
	double ans = 0, weight;
	
	for (int i = 0; i < naxis; ++i) if (indList[i] * (indList[i] - (i == VAXIS ? nvel : npix) + 1.) > 0.) return ans;
	
	for (int ix = 0; ix < 2; ++ix) for (int iy = 0; iy < 2; ++iy) for (int iv = 0; iv < 2; ++iv) {
		if (indList[XAXIS] + ix >= npix or indList[YAXIS] + iy >= npix or indList[VAXIS] >= nvel) continue;
		weight = (ix + 1) % 2 - (indList[XAXIS] - floor(indList[XAXIS] * (1. + db_erratio))) * pow(-1., ix);
		weight *= (iy + 1) % 2 - (indList[YAXIS] - floor(indList[YAXIS] * (1. + db_erratio))) * pow(-1., iy);
		weight *= (iv + 1) % 2 - (indList[VAXIS] - floor(indList[VAXIS] * (1. + db_erratio))) * pow(-1., iv);
		
		ans += emission[int(indList[XAXIS]) + ix][int(indList[YAXIS]) + iy][int(indList[VAXIS]) + iv] * weight;
	}
	
	
	
	return ans;
}

void skyPlane::setEmissionFits() {
	for (int ich = 0; ich < nvel; ++ich) for (int idec = 0; idec < npix; ++idec) for (int ira = 0; ira < npix; ++ira) emission_fits[(ich * npix + idec) * npix + ira] = emission[npix - 1 - ira][idec][ich];
	return ;
}

void skyPlane::readEmissionFits() {
	for (int ich = 0; ich < nvel; ++ich) for (int idec = 0; idec < npix; ++idec) for (int ira = 0; ira < npix; ++ira) emission[npix - 1 - ira][idec][ich] = emission_fits[(ich * npix + idec) * npix + ira];
	setMaxMin(); 
	return ;
}

//Convolution
void skyPlane::FFT_main(dcomp *a, int k, int parity) {
    if (k == 1) return ;
    
    k = k >> 1;
    FFT_main(a, k, parity);
    FFT_main(&a[k], k, parity);
    
    dcomp *q, *s;
    q = (dcomp *) malloc(sizeof(dcomp) * k);
    s = (dcomp *) malloc(sizeof(dcomp) * k);
    
    for (int i = 0; i < k; ++i) {
        q[i] = a[i];
        s[i] = a[i + k];
    }
    
    for (int i = 0; i < k; ++i) {
        dcomp w = polar(double(1.0), double(parity * pi * i / k));
        a[i] = q[i % k] + w * s[i % k];
        a[i + k] = q[i % k] - w * s[i % k];
    }
    
    return ;
}

void skyPlane::FFT(dcomp *a, int k, int parity) {
    dcomp *b;
    b = (dcomp *) malloc(sizeof(dcomp) * k);
    
    for (int i = 0; i < k; ++i) b[i] = a[i];
    for (int i = 0; i < k; ++i) {
        int num = 0;
        for (int j = 0; (1 << j) < k; ++j) {
            num *= 2;
            num += (i >> j) % 2;
        }
        a[i] = b[num];
    }
    
    FFT_main(a, k, parity);
    for (int i = 0; i < k; ++i) a[i] /= sqrt(k);
    
    return ;
}

void skyPlane::make_beamuv() {
    double u, v, val1, val2;
	double a_u2 = - pi * pi / 4. / log(2.) * (pow(beam_min_cm * beam_cpa, 2) + pow(beam_maj_cm * beam_spa, 2));
	double a_v2 = - pi * pi / 4. / log(2.) * (pow(beam_min_cm * beam_spa, 2) + pow(beam_maj_cm * beam_cpa, 2));
	double a_uv = - pi * pi / 2. / log(2.) * (beam_min_cm * beam_min_cm - beam_maj_cm * beam_maj_cm) * beam_cpa * beam_spa;
	
    int iiu, iiv;
	for (int iu = 0; iu < npix; ++iu) {
		u = du * iu; iiu = 2 * npix - 1 - iu;
		for (int iv = 0; iv < npix; ++iv) {
			v = dv * iv; iiv = 2 * npix - 1 - iv;
			val1 = exp(a_u2 * u * u + a_v2 * v * v + a_uv * u * v);
			val2 = exp(a_u2 * u * u + a_v2 * v * v - a_uv * u * v);
			beam_uv[iu][iv] = dcomp(val1);
			beam_uv[iu][iiv] = dcomp(val2);
			beam_uv[iiu][iv] = dcomp(val2);
			beam_uv[iiu][iiv] = dcomp(val1);
		}
    }
    return ;
}

void skyPlane::make_linewidthuv() {
	double w, val, a = - pow(pi * linewidth_FWHM_cmps, 2) / 4. / log(2.);
    for (int iw = 0; iw < nvel; ++iw) {
        w = dw * iw;
        val = 1. * exp(a * w * w);
        lw_uv[iw] = dcomp(val);
        lw_uv[2 * nvel - 1 - iw] = dcomp(val);
    }
    
    return ;
}

void skyPlane::convolution_beam(int iv) {
	dcomp **emission_uv, *emission_temp;
	
	emission_uv = (dcomp**) malloc(sizeof(dcomp*) * npix * 2);
	emission_uv[0] = (dcomp*) malloc(sizeof(dcomp) * npix * npix * 4);
	emission_temp = (dcomp*) malloc(sizeof(dcomp) * npix * 2);
	for (int iu = 0; iu < npix * 2; ++iu) {
		emission_uv[iu] = emission_uv[0] + npix * 2 * iu;
	}
	
	
	for (int iy = 0; iy < npix * 2; ++iy) {
        if (iy >= npix) {
            for (int iu = 0; iu < npix * 2; ++iu) emission_uv[iu][iy] = dcomp(0.);
            continue;
        }
		for (int ix = 0; ix < npix * 2; ++ix) emission_temp[ix] = (ix < npix ? dcomp(emission[ix][iy][iv]) : dcomp(0.));
        FFT(emission_temp, npix * 2, 1);
        for (int iu = 0; iu < npix * 2; ++iu) emission_uv[iu][iy] = emission_temp[iu];
	}
	for (int iu = 0; iu < npix * 2; ++iu) FFT(emission_uv[iu], npix * 2, 1);
    
	for (int iu = 0; iu < npix * 2; ++iu) for (int iv = 0; iv < npix * 2; ++iv) emission_uv[iu][iv] *=beam_uv[iu][iv];
    
    for (int iu = 0; iu < npix * 2; ++iu) FFT(emission_uv[iu], npix * 2, -1);   
    for (int iy = 0; iy < npix; ++iy) {
        for (int iu = 0; iu < npix * 2; ++iu) emission_temp[iu] = emission_uv[iu][iy];
        FFT(emission_temp, npix * 2, -1);
        for (int ix = 0; ix < npix; ++ix) emission[ix][iy][iv] = emission_temp[ix].real();
    }
	
	
	
	free(emission_uv[0]);
	free(emission_uv);
	free(emission_temp);
	
    return ;
}

void skyPlane::convolution_line(int ix, int iy) {
	
	dcomp *emission_uv;
	emission_uv = (dcomp*) malloc(sizeof(dcomp) * nvel * 2);
	
	
	for (int iv = 0; iv < nvel * 2; ++iv) emission_uv[iv] = (iv < nvel ? dcomp(emission[ix][iy][iv]) : dcomp(0.));
    FFT(emission_uv, nvel * 2, 1);
    for (int iw = 0; iw < nvel * 2; ++iw) emission_uv[iw] *= lw_uv[iw];
    
    FFT(emission_uv, nvel * 2, -1);
    
	for (int iv = 0; iv < nvel; ++iv) emission[ix][iy][iv] = emission_uv[iv].real();
	
	
	free(emission_uv);
	
    return;
}

void skyPlane::normalize() {
	setMaxMin();
	
	if (max(abs(maxval), abs(minval)) < EPS) {
		puts("\n!! WARN :: the intensity is 0 everywhere !!\n\n\t-> The values are not normalized.\n\t-> The value for the first pixel is set to be 1 as a dummy pixel.");
		emission[0][0][0] = 1.;
		return ; 
	}
	
	double norm = 1. / max(abs(maxval), abs(minval));
	
	for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) for (int iv = 0; iv < nvel; ++iv) emission[ix][iy][iv] *= norm;
	
	return; 
}

