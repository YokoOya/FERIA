#include "feria.h"

/*//////////////////////////////////////
 The values are expressed in units of
 cm, g, sec, radian, Kelvin
 /*//////////////////////////////////////

/*/////parameters
private:
	double ****meshdata;
	double crpix[naxis], crval[naxis], cdelt[naxis];
	double maxval[ndata], minval[ndata];
	double cinc, sinc, cpa, spa;
/*/////parameters

mesh::mesh(double meshCrpix[naxis], double meshCrval_au[naxis], double meshCdelt_au[naxis], double inc_deg, double pa_deg) {
	pa_deg *= -1;
	inc_deg -= 90.;
    for (int i = 0; i < naxis; ++i) {
        crpix[i] = meshCrpix[i];
        crval[i] = meshCrval_au[i] * CMperAU;
        cdelt[i] = meshCdelt_au[i] * CMperAU;
    }
    for (int idata = 0; idata < ndata; ++idata) {
        maxval[idata] = -INF;
        minval[idata] = INF;
    }
    
    cinc = cos(inc_deg / DEGperRAD);
    sinc = sin(inc_deg / DEGperRAD);
    cpa = cos(pa_deg / DEGperRAD);
    spa = sin(pa_deg / DEGperRAD);
	
	
	meshdata = (double****) malloc(sizeof(double***) * npix);
	meshdata[0] = (double***) malloc(sizeof(double**) * npix * npix);
	meshdata[0][0] = (double**) malloc(sizeof(double*) * npix * npix * npix);
	meshdata[0][0][0] = (double*) malloc(sizeof(double) * npix * npix * npix * ndata);

	for (int ix = 0; ix < npix; ++ix) {
		meshdata[ix] = meshdata[0] + npix * ix;
		for (int iy = 0; iy < npix; ++iy) {
			meshdata[ix][iy] = meshdata[0][0] + npix * npix * ix + npix * iy;
			for (int iz = 0; iz < npix; ++iz) {
				meshdata[ix][iy][iz] = meshdata[0][0][0] + npix * npix * ndata * ix + npix * ndata * iy + ndata * iz;
				for (int idata = 0; idata < ndata; ++idata) {
					meshdata[ix][iy][iz][idata] = 0.; 
				}
			}
		}
	}

}

mesh::~mesh() {	
	free(meshdata[0][0][0]);
	free(meshdata[0][0]);
	free(meshdata[0]);
	free(meshdata);
	
}

bool mesh::checkAxis(int ix, int iy, int iz, int idata) {
    if (ix * (npix - 1 - ix) < 0 || iy * (npix - 1 - iy) < 0 || iz * (npix - 1 - iz) < 0 || idata * (ndata - 1 - idata) < 0) return false;
    return true;
}

void mesh::setData(int ix, int iy, int iz, int idata, double val) {
    if (!checkAxis(ix, iy, iz, idata)) return ;
    meshdata[ix][iy][iz][idata] = val;
    return ;
}

double mesh::getData(int ix, int iy, int iz, int idata) {
    if (!checkAxis(ix, iy, iz, idata)) return 0.;
    return meshdata[ix][iy][iz][idata];
}

void mesh::getColumn(int ix, int iy, double (*gasColumn)[ndata]) {
	if (checkAxis(ix, iy, 0, 0)) for (int iz = 0; iz < npix; ++iz) for (int idata = 0; idata < ndata; ++idata) gasColumn[iz][idata] = meshdata[ix][iy][iz][idata];
	return ;
}

void mesh::setMaxMin(int idata) {
    if (!checkAxis(0, 0, 0, idata)) return ;
    maxval[idata] = -INF;
    minval[idata] = INF;
    for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) for (int iz = 0; iz < npix; ++iz) {
        double val = meshdata[ix][iy][iz][idata];
        if (maxval[idata] < val) maxval[idata] = val;
        if (minval[idata] > val) minval[idata] = val;
    }
    return ;
}

void mesh::getMaxMin(double *max, double *min, int idata) {
    bool check = checkAxis(0, 0, 0, idata);
    if (max) *max = (check ? maxval[idata] : -INF);
    if (min) *min = (check ? minval[idata] : INF);
    return ;
}

void mesh::getPos_cart(int ix, int iy, int iz, double *pos_cart) {
    pos_cart[XAXIS] = (ix - crpix[XAXIS]) * cdelt[XAXIS] + crval[XAXIS];
    pos_cart[YAXIS] = (iy - crpix[YAXIS]) * cdelt[YAXIS] + crval[YAXIS];
    pos_cart[ZAXIS] = (iz - crpix[ZAXIS]) * cdelt[ZAXIS] + crval[ZAXIS];
    return ;
}

void mesh::pos_cart2polar(double *pos_cart, double *pos_polar) {
    double x = -pos_cart[XAXIS] * spa + pos_cart[YAXIS] * cpa;
    double y = -(pos_cart[XAXIS] * cpa + pos_cart[YAXIS] * spa) * sinc + pos_cart[ZAXIS] * cinc;
    double z = -(pos_cart[XAXIS] * cpa + pos_cart[YAXIS] * spa) * cinc - pos_cart[ZAXIS] * sinc;
    
    pos_polar[rAXIS] = sqrt(pow(x, 2) + pow(y, 2));
	pos_polar[tAXIS] = atan2(y, x);
    pos_polar[ZAXIS] = z;
    return ;
}

void mesh::getPos_polar(int ix, int iy, int iz, double *pos_polar) {
    double pos_cart[naxis];
    getPos_cart(ix, iy, iz, pos_cart);
    pos_cart2polar(pos_cart, pos_polar);
    return ;
}

void mesh::setVel_cart(int ix, int iy, int iz, double *vel_cart) {
    setData(ix, iy, iz, VxVAL, vel_cart[VxVAL]);
    setData(ix, iy, iz, VyVAL, vel_cart[VyVAL]);
    setData(ix, iy, iz, VzVAL, vel_cart[VzVAL]);
    return ; 
}

void mesh::vel_polar2cart(double *pos_polar, double *vel_cart, double *vel_polar) {
    double theta = pos_polar[tAXIS];
    double cth = cos(theta);
    double sth = sin(theta);
    vel_cart[XAXIS] = vel_polar[rAXIS] * (-sth * sinc * cpa - cth * spa) + vel_polar[tAXIS] * (-cth * sinc * cpa + sth * spa) + vel_polar[ZAXIS] * (-cinc * cpa);
    vel_cart[YAXIS] = vel_polar[rAXIS] * (-sth * sinc * spa + cth * cpa) + vel_polar[tAXIS] * (-cth * sinc * spa - sth * cpa) + vel_polar[ZAXIS] * (-cinc * spa);
    vel_cart[ZAXIS] = vel_polar[rAXIS] * (sth * cinc) + vel_polar[tAXIS] * (cth * cinc) + vel_polar[ZAXIS] * (-sinc);
    return ;
}

void mesh::setVel_polar(int ix, int iy, int iz, double *pos_polar, double *vel_polar) {
    double vel_cart[naxis] = {0.};
    vel_polar2cart(pos_polar, vel_cart, vel_polar);
    setVel_cart(ix, iy, iz, vel_cart);
    return ;
}

void mesh::getVel_cart(int ix, int iy, int iz, double *vel_cart) {
    for (int i = 0; i < naxis; ++i) vel_cart[i] = meshdata[ix][iy][iz][i];
    return ;
}

void mesh::vel_cart2polar(double *pos_polar, double *vel_cart, double *vel_polar) {
    double theta = pos_polar[tAXIS];
    double cth = cos(theta);
    double sth = sin(theta);
    vel_polar[rAXIS] = vel_cart[XAXIS] * (-spa * cth - cpa * sinc * sth) + vel_cart[YAXIS] * (cpa * cth - spa * sinc * sth) + vel_cart[ZAXIS] * cinc * sth;
    vel_polar[tAXIS] = vel_cart[XAXIS] * (-spa * sth - cpa * sinc * cth) + vel_cart[YAXIS] * (cpa * sth - spa * sinc * cth) + vel_cart[ZAXIS] * cinc * cth;
    vel_polar[ZAXIS] = vel_cart[XAXIS] * (-cpa * cinc) + vel_cart[YAXIS] * (-spa * cinc) + vel_cart[ZAXIS] * (-sinc);
    return ; 
}

void mesh::getVel_polar(int ix, int iy, int iz, double *vel_polar) {
    double pos_polar[naxis], vel_cart[naxis];
    getPos_polar(ix, iy, iz, pos_polar);
    getVel_cart(ix, iy, iz, vel_cart);
    vel_cart2polar(pos_polar, vel_cart, vel_polar);
}
