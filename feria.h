#ifndef _HEAD_INCLUDED
#define _HEAD_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <ctime>
#include "fitsio.h"


using namespace std; 

typedef complex<double> dcomp; 

//---- Header ---//
const int lbNpix = 5;
const int lbNvel = 5;

//---- Constant ---//
const double EPS = 1e-308;
const double INF = 1e+308;
const double pi = acos(-1.);
const double db_erratio = 5e-15; 

const double c_light = 2.99792458e+10;	// cm s^-1
const double k_boltz = 1.3806488e-16;	// cm^2 g s^-2 K^-1
const double Grav = 6.67384e-8;			// cm^3 g^-1 s^-2
const double h_planck = 6.62606957e-27;	// cm^2 g s^-1
const double Msun = 1.989e+33;			// g
const double CMperM = 1e+2;				// cm
const double CMperKM = 1e+5;				// cm
const double CMperAU = 1.49597871e+13;	// cm
const double ASperRAD = 180. * 3600. / pi;
const double DEGperRAD = 180. / pi;
const double CMperPC = CMperAU * ASperRAD;	// cm

const char FITS[] = ".fits";
const char PARAMSNAME[] = "parameter";
const char EQUINOX[] = "J2000";

//---- Axis ---//
const int npix = 1 << lbNpix;
const int nvel = 1 << lbNvel;
const int XAXIS = 0;
const int rAXIS = 0;
const int YAXIS = 1;
const int tAXIS = 1; 
const int ZAXIS = 2;
const int VAXIS = 2;
const int naxis = 3;

const int PAXIS_PV = 0;
const int VAXIS_PV = 1;
const int naxis_PV = 2;

const int noGAS = 0;
const int isIRE = 1;
const int isKEP = 2;


const int VxVAL = 0;
const int VyVAL = 1;
const int VzVAL = 2;
const int NVAL = 3;
const int TVAL = 4;
const int ndata = 5;

const int lfilename = 256;
const int lstr = 256;
const int lcard = 80;

//---- Prototype ---//
class mesh {
private:
	double ****meshdata; 
    double crpix[naxis], crval[naxis], cdelt[naxis];
    double maxval[ndata], minval[ndata];
    double cinc, sinc, cpa, spa;
public:
    mesh(double meshCrpix[naxis], double meshCrval_au[naxis], double meshCdelt_au[naxis], double inc, double pa);
	~mesh(); 
    bool checkAxis(int ix, int iy, int iz, int idata);
    void setData(int ix, int iy, int iz, int idata, double val);
    double getData(int ix, int iy, int iz, int idata);
    void getColumn(int ix, int iy, double (*gasColumn)[ndata]);
    void setMaxMin(int idata);
    void getMaxMin(double *max, double *min, int idata); 
    void getPos_cart(int ix, int iy, int iz, double *pos_cart); 
    void pos_cart2polar(double *pos_cart, double *pos_polar);
    void getPos_polar(int ix, int iy, int iz, double *pos_polar);
    void setVel_cart(int ix, int iy, int iz, double *vel_cart);
    void vel_polar2cart(double *pos_polar, double *vel_cart, double *vel_polar);
    void setVel_polar(int ix, int iy, int iz, double *pos_polar, double *vel_polar); 
    void getVel_cart(int ix, int iy, int iz, double *vel_cart);
    void vel_cart2polar(double *pos_polar, double *vel_cart, double *vel_polar);
    void getVel_polar(int ix, int iy, int iz, double *vel_polar); 
};

class ireModel {
public:
	double mass_grav, rCB, inc, rot, PA, Rout, Rin, HeightIRE, tanFlareIRE, HeightKEP, tanFlareKEP;
	double DensCB, densProfIRE, densProfKEP, TempCB, tempProfIRE, tempProfKEP;
	double velrot_CB_xrCB;
	
    ireModel(double mass_msun, double rCB_au, double inc_deg, double PA_deg, double rot_sign, double Rout_au, double Rin_au, double HeightIRE_au, double FlareIRE_deg, double HeightKEP_au, double FlareKEP_deg, double DensCB_in, double densityProfileIRE, double densityProfileKEP, double TempCB_in, double tempProfileIRE, double tempProfileKEP);
	
	int which_velStr(double dist_protostar, double r, double theta, double z);
    void env_vel(double dist_protostar, double r, double theta, double z, double *vel_polar, int velStrID);
    double env_temp(double dist_protostar, double r, double theta, double z, int velStrID);
    double env_dens(double dist_protostar, double r, double theta, double z, int velStrID);
};


class sourceParams {
public:
	char outputfilename[lstr], name_line[lstr], name_object[lstr], radesys[lstr], centRa_hms[lstr], centDec_dms[lstr];
	double Fldres_as, Velres_kmps, distance_pc, mass_msun, rCB_au, inc_deg, PA_deg, rot_sign, Rout_au, Rin_au;
	double HeightIRE_au, FlareIRE_deg, HeightKEP_au, FlareKEP_deg;
	double DensCB, densityProfileIRE, densityProfileKEP;
	double TempCB, tempProfileIRE, tempProfileKEP;
	double linewidth_cmps, beam_maj_as, beam_min_as, beam_pa_deg;
	double restfreq_Hz, centRa[3], centDec[3], centRa_deg, centDec_deg, centRa_scaled_as, centDec_as, Vsys_cmps;
	
	sourceParams(char *outputfilename_in, double Fldres_in, double Velres_in, double distance_pc_in, double mass_msun_in, double rCB_au_in, double inc_deg_in, double PA_deg_in, double rot_sign_in, double Rout_au_in, double Rin_au_in, double HeightIRE_au_in, double FlareIRE_deg_in, double HeightKEP_au_in, double FlareKEP_deg_in, double DensCB_in, double densityProfileIRE_in, double densityProfileKEP_in, double TempCB_in, double tempProfileIRE_in, double tempProfileKEP_in, double linewidth_kmps_in, double beam_maj_as_in, double beam_min_as_in, double beam_pa_deg_in, char* name_line_in, double restfreq_GHz_in, char* name_object_in, char* radesys_in, char* center_ra_in, char* center_dec_in, double Vsys_kmps_in);
	~sourceParams();
};

class skyPlane {
public:
	double ***emission, *emission_fits;
	double crpix[naxis], crval[naxis], cdelt[naxis], maxval, minval;
    double linewidth_FWHM_cmps, beam_maj_as, beam_min_as, beam_maj_rad, beam_min_rad, beam_maj_deg, beam_min_deg, beam_maj_cm, beam_min_cm, beam_pa_rad, beam_pa_deg, beam_cpa, beam_spa;
    double du, dv, dw;
    dcomp beam_uv[npix<<1][npix<<1], lw_uv[nvel<<1];
    
    skyPlane(double distance_pc, double skyCrpix[naxis], double skyCrval[naxis], double skyCdelt[naxis], double linewidth_kmps, double beam_maj_as, double beam_min_as, double beam_pa_degin);
	~skyPlane();
	void setCrval(int axis, double val);
    void projection(int ix, int iy, double (*gasColumn)[ndata]);
	double getIndex(int axis, double pos);
	double getPos(int axis, double index);
	void getIndexList(double *indList, double *posList);
	void getPosList(double *indList, double *posList);
    void setMaxMin();
    void getMaxMin(double *max, double *min);
    double getSum();
    double getVal(double *indList);
	void setEmissionFits();
	void readEmissionFits(); 
    
    void FFT_main(dcomp *a, int k, int parity);
    void FFT(dcomp *a, int k, int parity);
    void make_beamuv();
    void make_linewidthuv();
    void convolution_beam(int iv);
    void convolution_line(int ix, int iy);
	void normalize();
};

class PVdiagram {
public:
	double **emission, *emission_fits;
	int npix_PV;
	double crpix[naxis_PV], crval[naxis_PV], cdelt[naxis_PV], maxval, minval;
	double beam_maj_as, beam_min_as, beam_maj_rad, beam_min_rad, beam_maj_deg, beam_min_deg, beam_maj_cm, beam_min_cm, beam_pa_rad, beam_pa_deg, beam_cpa, beam_spa;
	double cent_offset_ra_au, cent_offset_dec_au, PV_PA_deg, PV_PA_rad, beam_slice_as, beam_slice_rad, beam_slice_deg;
	double centcube_ra_scaled_as, centcube_dec_as, centpv_ra_scaled_as, centpv_dec_as;
	
	PVdiagram(int npix_in, double pvCrpix[naxis_PV], double pvCrval[naxis_PV], double pvCdelt[naxis_PV], sourceParams &source, double cent_offset_ra_auin, double cent_offset_dec_auin, double PV_PA_degin);
	~PVdiagram();
	double getIndex(int axis, double pos);
	double getPos(int axis, double index); 
	void getIndexList(double *indList, double *posList);
	void getPosList(double *indList, double *posList);
	void setMaxMin();
	void getMaxMin(double *max, double *min);
	double getSum(); 
	void setVal(int ip, int iv, double val);
	double getVal(int ip, int iv);
	void setEmissionFits();
};


void radec2arcsec(double *ra_hms, double *dec_dms, double *ra_scaled_as, double *dec_as);
void arcsec2radec(double *ra_hms, double *dec_dms, double ra_scaled_as, double dec_as);
void radecShift(double *ra_scaled_as, double *dec_as, double ra_shift_as, double dec_shift_as); 
void ire_fits_write_record(fitsfile* fptr, sourceParams &source, int *statusptr);
void pv_fits_write_record(fitsfile *fptr, sourceParams &source, PVdiagram &PVdiag, int *statusptr);
void ire_fitsout(fitsfile *fptr, long *fpixel, long nelements, double *valList, int *statusptr);
void ire_cube_fitsout(char* outfilename, skyPlane &sky, sourceParams &source);
void ire_cube_fitsin(char* infilename, skyPlane &sky, sourceParams &source, bool *f_makecube);
void ire_PV_fitsout(char* outfilename, PVdiagram &PVdiag, sourceParams &source, int npix_PV);

#endif
