#include "feria.h"

void radec2arcsec(double *ra_hms, double *dec_dms, double *ra_scaled_as, double *dec_as) {
	if (dec_dms[0] < 0) for (int i = 1; i < 3; ++i) dec_dms[i] *= -1.;
	
	*dec_as = (dec_dms[0] * 60. + dec_dms[1]) * 60. + dec_dms[2];
	*ra_scaled_as = ((ra_hms[0] * 60. + ra_hms[1]) * 60. + ra_hms[2]) * 360. / 24.;
	
	return ;
}

void arcsec2radec(double *ra_hms, double *dec_dms, double ra_scaled_as, double dec_as) {
	int sign_ra = (ra_scaled_as < 0 ? -1 : 1), sign_dec = (dec_as < 0 ? -1 : 1);
	
	ra_scaled_as = abs(ra_scaled_as) * 24. / 360.;
	dec_as = abs(dec_as);
	
	ra_hms[0] = floor(ra_scaled_as * (1. + db_erratio) / 3600.);
	ra_scaled_as -= ra_hms[0] * 3600.;
	ra_hms[1] = floor(ra_scaled_as * (1. + db_erratio) / 60.);
	ra_scaled_as -=  ra_hms[1] * 60.;
	ra_hms[2] = ra_scaled_as;
	
	
	dec_dms[0] = floor(dec_as * (1. + db_erratio) / 3600.);
	dec_as -= dec_dms[0] * 3600.;
	dec_dms[1] = floor(dec_as * (1. + db_erratio) / 60.);
	dec_as -=  dec_dms[1] * 60.;
	dec_dms[2] = dec_as;
	
	ra_hms[0] *= sign_ra;
	dec_dms[0] *= sign_dec;
		
	return ;
}


void radecShift(double *ra_scaled_as, double *dec_as, double ra_shift_as, double dec_shift_as) {
	double dec_shifted_as = *dec_as + dec_shift_as;
	
	*ra_scaled_as += ra_shift_as / cos(pi / 180. * (*dec_as + dec_shifted_as) / 2. / 3600.);
	*dec_as = dec_shifted_as;
	
	return ;
}




void ire_fits_write_record(fitsfile *fptr, sourceParams &source, int *statusptr) {
	char card_null[lcard] = "";
	
	char ireKey[][lstr] = {
		"OBJ",
		"LINE",
		"DIST",
		"MASS",
		"RCB",
		"INC",
		"PA",
		"ROT",
		"ROUT",
		"RIN",
		"TIRE",
		"FLAREIRE",
		"DPROIRE",
		"TPROIRE",
		"TKEP",
		"FLAREKEP",
		"DPROKEP",
		"TPROKEP",
		"H2",
		"DFRAC",
		"TEMP",
		"LW",
		"BMAJ",
		"BMIN",
		"BPA",
		"FREST",
		"VSYS"
	}, ireKeyname[lstr] = "";
	char ireComment[][lstr] = {
		"Target source",
		"Molecular line",
		"Distance (pc)",
		"Protostellar mass (Msun)",
		"Radius of the centrifugal barrier (au)",
		"Inclination angle (0 degree for face-on; deg)",
		"Position angle of the elongation (deg)",
		"Rotation ('1' for positive)",
		"Outer radius (au)",
		"Inner radius (au)",
		"Scale height (au) (IRE)",
		"Flared angle (deg) (IRE)",
		"Density profile of H2 (IRE)",
		"Temperature profile (IRE)",
		"Scale height (au) (Kep)",
		"Flared angle (deg) (Kep)",
		"Density profile of H2 (Kep)",
		"Temperature profile (Kep)",
		"n(H2) at the CB (cm-3)",
		"Fractional abundance of the target molecule",
		"Gas kinetic temperature at the CB (K)",
		"Intrinsic linewidth (km s-1)",
		"Beam (major) (arcsec)",
		"Beam (minor) (arcsec)",
		"Beam (PA) (deg)",
		"Rest frequency (Hz)",
		"Systemic velocity (km s-1)"
	};

	char* ireStr[] = {
		source.name_object,
		source.name_line
	};
	double ireVal[] = {
		source.distance_pc,
		source.mass_msun,
		source.rCB_au,
		source.inc_deg,
		source.PA_deg,
		source.rot_sign,
		source.Rout_au,
		source.Rin_au,
		source.HeightIRE_au,
		source.FlareIRE_deg,
		source.densityProfileIRE,
		source.tempProfileIRE,
		source.HeightKEP_au,
		source.FlareKEP_deg,
		source.densityProfileKEP,
		source.tempProfileKEP,
		source.H2densityCB,
		source.fractionalDensity,
		source.TempCB,
		source.linewidth_cmps / CMperKM,
		source.beam_maj_as,
		source.beam_min_as,
		source.beam_pa_deg,
		source.restfreq_Hz,
		source.Vsys_cmps / CMperKM
	};
	
	
	

	int nirenames = sizeof(ireStr)/sizeof(ireStr[0]);
	for (int i = 0; i < nirenames; ++i) {
		sprintf(ireKeyname, "IRE%s", ireKey[i]);
		fits_write_key(fptr, TSTRING, ireKeyname, ireStr[i], ireComment[i], statusptr);
	}
	int nireparams = sizeof(ireVal)/sizeof(ireVal[0]);
	for (int i = 0; i < nireparams; ++i) {
		sprintf(ireKeyname, "IRE%s", ireKey[i + nirenames]);
		fits_write_key(fptr, TDOUBLE, ireKeyname, &ireVal[i], ireComment[i + nirenames], statusptr);
	}


	fits_write_record(fptr, card_null, statusptr);
	fits_write_date(fptr, statusptr);
	
	
	return ;
}

void pv_fits_write_record(fitsfile *fptr, sourceParams &source, PVdiagram &PVdiag, int *statusptr) {
	char card_null[lcard] = "";
	
	char pvKey[][lstr] = {
		"CENTRA",
		"CENTDEC",
		"FLDRA",
		"FLDDEC",
		"OFFRA",
		"OFFDEC",
		"BEAM"
	}, pvKeyname[lstr] = "";
	char pvComment[][lstr] = {
		"Slice Center in RA",
		"Slice Center in DEC",
		"Center of Cube in RA",
		"Center of Cube in DEC",
		"Slice offset from the cube center (au)",
		"Sliece offset from the cube center (au)",
		"Beam size along the slice (deg)"
	};
	
	char centpv_str[][lstr] = {"", ""};
	double centpv_ra_hms[3], centpv_dec_dms[3];
	arcsec2radec(centpv_ra_hms, centpv_dec_dms, PVdiag.centpv_ra_scaled_as, PVdiag.centpv_dec_as);
	
	sprintf(centpv_str[0], "%02.0lfh%02.0lfm%09.6lfs", centpv_ra_hms[0], centpv_ra_hms[1], centpv_ra_hms[2]);
	sprintf(centpv_str[1], "%02.0lfd%02.0lfm%09.6lfs", centpv_dec_dms[0], centpv_dec_dms[1], centpv_dec_dms[2]);

	
	char* pvStr[] = {
		centpv_str[0],
		centpv_str[1],
		source.centRa_hms,
		source.centDec_dms,
	};
	double pvVal[] = {
		PVdiag.cent_offset_ra_au,
		PVdiag.cent_offset_dec_au,
		PVdiag.beam_slice_deg,
	};
	
	sprintf(pvComment[6], "%s = %lf as", pvComment[6], PVdiag.beam_slice_as); 
	
	int npvnames = sizeof(pvStr)/sizeof(pvStr[0]);
	for (int i = 0; i < npvnames; ++i) {
		sprintf(pvKeyname, "PV%s", pvKey[i]);
		fits_write_key(fptr, TSTRING, pvKeyname, pvStr[i], pvComment[i], statusptr);
	}
	int npvparams = sizeof(pvVal)/sizeof(pvVal[0]);
	for (int i = 0; i < npvparams; ++i) {
		sprintf(pvKeyname, "PV%s", pvKey[i + npvnames]);
		fits_write_key(fptr, TDOUBLE, pvKeyname, &pvVal[i], pvComment[i + npvnames], statusptr);
	}
	
	
	fits_write_record(fptr, card_null, statusptr);
	fits_write_date(fptr, statusptr);
	
	
	return ;
}




void printerror(int status) {
	if (status) {
		fits_report_error(stderr, status);
		exit(status);
	}
	return;
}


void ire_fitsout(char* outfilename, fitsfile *fptr, long *fpixel, long nelements, double *valList, int* statusptr) {
	if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, valList, statusptr)) printerror(*statusptr);
	
	char err_text[lcard];
	fits_get_errstatus(*statusptr, err_text);
	printf("\n\tStatus of fits_write_pix = %d -> %s\n\n", *statusptr, err_text);
	
	
	if (fits_close_file(fptr, statusptr)) printerror(*statusptr);
	printf("Fits file was created: %s\n", outfilename);
}


void ire_cube_fitsout(char* outfilename, skyPlane &sky, sourceParams &source) {
	fitsfile *fptr;
	int status = 0;
	
	if (fopen(outfilename, "r") != NULL) {
		printf("\nFits file will be removed: %s\n", outfilename);
		remove(outfilename);
	}
	if (fits_create_file(&fptr, outfilename, &status)) printerror(status);
	
	long npix_axis[naxis];
	for (int i = 0; i < naxis; ++i) npix_axis[i] = (i == VAXIS ? nvel : npix);
	if (fits_create_img(fptr, FLOAT_IMG, naxis, npix_axis, &status)) printerror(status);
	
	//---Move the COMMENT by 'fits_create_img' task to backward---//
	int lcomment = 2;
	char card_createimg[lcomment][lcard];
	for (int i = 0; i < lcomment; ++i) fits_read_card(fptr, "COMMENT", card_createimg[i], &status);
	for (int i = 0; i < lcomment; ++i) fits_delete_key(fptr, "COMMENT", &status);
	
	
	//---Create other records---//
	char card_null[lcard + 1] = "";
	for (int i = 0; i < lcard; ++i) strcat(card_null, " ");
	
	char dummy_str[2][lcard] = {"", ""};
	double dummy_flt = 0;
	dummy_flt = 1; //bscale
	fits_write_key(fptr, TDOUBLE, "BSCALE", &dummy_flt, "PHYSICAL = PIXEL*BSCALE + BZERO", &status);
	dummy_flt = 0.; //bzero
	fits_write_key(fptr, TDOUBLE, "BZERO", &dummy_flt, "", &status);
	sprintf(dummy_str[0], "%lf as", sky.beam_maj_as);
	fits_write_key(fptr, TDOUBLE, "BMAJ", &sky.beam_maj_deg, dummy_str[0], &status);
	sprintf(dummy_str[0], "%lf as", sky.beam_min_as);
	fits_write_key(fptr, TDOUBLE, "BMIN", &sky.beam_min_deg, dummy_str[0], &status);
	fits_write_key(fptr, TDOUBLE, "BPA", &sky.beam_pa_deg, "deg", &status);
	sprintf(dummy_str[0], "Intensity");
	fits_write_key(fptr, TSTRING, "BTYPE", &dummy_str[0], "", &status);
	sprintf(dummy_str[0], "IRE model");
	sprintf(dummy_str[1], "for %s", source.name_object);
	fits_write_key(fptr, TSTRING, "OBJECT", &dummy_str[0], dummy_str[1], &status);
	fits_write_record(fptr, card_null, &status);
	sprintf(dummy_str[0], "JY/BEAM");
	fits_write_key(fptr, TSTRING, "BUNIT", &dummy_str[0], "Brightness (pixel) unit", &status);
	
	bool fequinox = true;
	for (int i = 0; i < int(strlen(source.radesys)); ++i) {
		if (source.radesys[i] != EQUINOX[i]) fequinox = false;
	}
	if (fequinox) {
		dummy_flt = 2e3; //equinox
		fits_write_key(fptr, TDOUBLE, "EQUINOX", &dummy_flt, "", &status);
	} else {
		sprintf(dummy_str[0], "ICRS"); //radesys
		fits_write_key(fptr, TSTRING, "RADESYS", &dummy_str[0], "", &status);
	}
	/*
	 fits_write_key(fptr, TDOUBLE, "LONPOLE", &lonpole, "", &status);
	 fits_write_key(fptr, TDOUBLE, "LATPOLE", &latpole, "", &status);
	 //*/
	
	fits_write_record(fptr, card_null, &status);
	char ctype[naxis][9] = {"RA---SIN", "DEC--SIN", "VRAD    "}, cunit[naxis][9] = {"deg     ", "deg     ", "m/s    "};
	double crval[] = {source.centRa_deg, source.centDec_deg, source.Vsys_cmps / CMperM}; //sky.crval[VAXIS] / CMperM}; //
	char crvalComment[][lstr] = {"", "", ""};
	strcpy(crvalComment[XAXIS], source.centRa_hms);
	strcpy(crvalComment[YAXIS], source.centDec_dms);
	double cdelt[] = {sky.cdelt[XAXIS] / CMperAU / source.distance_pc / double(3600.), sky.cdelt[YAXIS] / CMperAU / source.distance_pc / double(3600.), sky.cdelt[VAXIS] / CMperM};
	//double crpix[] = {npix >> 1, npix >> 1, nvel >> 1};
	double crpix[naxis]; for (int i = 0; i < naxis; ++i) crpix[i] = sky.crpix[i] + 1;
	char cdeltComment[][lstr] = {"", "", ""};
	sprintf(cdeltComment[XAXIS], "%lf au", sky.cdelt[XAXIS] / CMperAU);
	sprintf(cdeltComment[YAXIS], "%lf au", sky.cdelt[YAXIS] / CMperAU);
	for (int i = 0; i < naxis; ++i) {
		char ckeyname[5][7] = {"CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT"};
		char id[] = "";
		sprintf(id, "%d", i+1);
		for (int j = 0; j < 5; ++j) strcat(ckeyname[j], id);
		fits_write_key(fptr, TSTRING, ckeyname[0], &ctype[i], "", &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[1], &crval[i], crvalComment[i], &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[2], &cdelt[i], cdeltComment[i], &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[3], &crpix[i], "", &status);
		fits_write_key(fptr, TSTRING, ckeyname[4], &cunit[i], "", &status);
	}
	
	fits_write_record(fptr, card_null, &status);
	char rfreqComment[] = "Rest Frequency (Hz) of ";
	strcat(rfreqComment, source.name_line);
	fits_write_key(fptr, TDOUBLE, "RESTFREQ", &source.restfreq_Hz, rfreqComment, &status);
	char specsys[] = "LSRK    ";
	fits_write_key(fptr, TSTRING, "SPECSYS", &specsys, "Spectral reference frame", &status);
	/*
	 fits_write_key(fptr, TSTRING, "ALTRVAL", &nakami_str, "Alternate frequency reference value", &status);
	 fits_write_key(fptr, TSTRING, "ALTRPIX", &nakami_str, "Alternate frequency reference pixel", &status);
	 fits_write_key(fptr, TSTRING, "VELREF", &nakami_str, "1 LSR, 2 HEL, 3 OBS, +256 Radio", &status);
	 //*/
	
	fits_write_record(fptr, card_null, &status);
	char comment_model[lstr] = {"This is a result of an infalling-rotating envelope model with the following physical parameters"};
	fits_write_comment(fptr, comment_model, &status);
	
	
	ire_fits_write_record(fptr, source, &status); // Parameters for IRE model
	
	
	//---Move the COMMENT by 'fits_create_img' task to backward---//
	fits_write_record(fptr, card_createimg[0], &status);
	fits_write_record(fptr, card_createimg[1], &status);
	
	
	//---Output the Results---//
	//use sky.emission[npix][npix][nvel]
	long fpixel[naxis];
	long nelements = 1;
	for (int i = 0; i < naxis; ++i) {
		fpixel[i] = 1;
		nelements *= npix_axis[i];
	}
	
	
	
	ire_fitsout(outfilename, fptr, fpixel, nelements, sky.emission_fits, &status);

	
	return ;
}




void ire_cube_fitsin(char* infilename, skyPlane &sky, sourceParams &source, bool *f_makecube) {
	fitsfile *fptr;
	int status = 0, nfound, anynull, firstelem = 1;
	double nulval = 0;
	long npix_axis[naxis], npixels = 1;
	
	if (fits_open_file(&fptr, infilename, READONLY, &status)) printerror(status);
	
	fits_read_keys_lng(fptr, "NAXIS", 1, naxis, npix_axis, &nfound, &status);
	for (int i = 0; i < naxis; ++i) {
		if (npix_axis[i] != (i == VAXIS ? nvel : npix)) {
			puts("\n!! WARN :: #Pixel mismatch !!");
			puts("\n\t-> The fits file will be replaced.\n");
			*f_makecube = true;
			return ;
		}
	}
	
	
	for (int i = 0; i < naxis; ++i) npixels *= npix_axis[i];
	if (fits_read_img(fptr, TDOUBLE, firstelem, npixels, &nulval, sky.emission_fits, &anynull, &status)) printerror(status);
	
	sky.readEmissionFits();
	
	return ;
}






void ire_PV_fitsout(char* outfilename, PVdiagram &PVdiag, sourceParams &source, int npix_PV) {
	fitsfile *fptr;
	int status = 0;
	
	if (fopen(outfilename, "r") != NULL) {
		printf("\nFits file will be removed: %s\n", outfilename);
		remove(outfilename);
	}
	
	
	if (fits_create_file(&fptr, outfilename, &status)) {
		printerror(status);
	}
	
	long npix_axis[naxis_PV];
	for (int i = 0; i < naxis_PV; ++i) npix_axis[i] = (i == VAXIS_PV ? nvel : npix_PV);
	if (fits_create_img(fptr, FLOAT_IMG, naxis_PV, npix_axis, &status)) printerror(status);

	//---Move the COMMENT by 'fits_create_img' task to backward---//
	int lcomment = 2;
	char card_createimg[lcomment][lcard];
	for (int i = 0; i < lcomment; ++i) fits_read_card(fptr, "COMMENT", card_createimg[i], &status);
	for (int i = 0; i < lcomment; ++i) fits_delete_key(fptr, "COMMENT", &status);
	
	
	//---Create other records---//
	char card_null[lcard + 1] = "";


	for (int i = 0; i < lcard; ++i) strcat(card_null, " ");
	

	char dummy_str[2][lcard] = {"", ""};
	double dummy_flt = 0;
	dummy_flt = 1; //bscale
	fits_write_key(fptr, TDOUBLE, "BSCALE", &dummy_flt, "PHYSICAL = PIXEL*BSCALE + BZERO", &status);
	dummy_flt = 0.; //bzero
	fits_write_key(fptr, TDOUBLE, "BZERO", &dummy_flt, "", &status);
	sprintf(dummy_str[0], "%lf as", PVdiag.beam_maj_as);
	fits_write_key(fptr, TDOUBLE, "BMAJ", &PVdiag.beam_maj_deg, dummy_str[0], &status);
	sprintf(dummy_str[0], "%lf as", PVdiag.beam_min_as);
	fits_write_key(fptr, TDOUBLE, "BMIN", &PVdiag.beam_min_deg, dummy_str[0], &status);
	fits_write_key(fptr, TDOUBLE, "BPA", &PVdiag.beam_pa_deg, "deg", &status);
	sprintf(dummy_str[0], "Intensity");
	fits_write_key(fptr, TSTRING, "BTYPE", &dummy_str[0], "", &status);
	sprintf(dummy_str[0], "PV of IRE model");
	sprintf(dummy_str[1], "for %s", source.name_object);
	fits_write_key(fptr, TSTRING, "OBJECT", &dummy_str[0], dummy_str[1], &status);
	fits_write_record(fptr, card_null, &status);
	sprintf(dummy_str[0], "JY/BEAM");
	fits_write_key(fptr, TSTRING, "BUNIT", &dummy_str[0], "Brightness (pixel) unit", &status);
	
	bool fequinox = true;
	for (int i = 0; i < int(strlen(source.radesys)); ++i) {
		if (source.radesys[i] != EQUINOX[i]) fequinox = false;
	}
	if (fequinox) {
		dummy_flt = 2e3; //equinox
		fits_write_key(fptr, TDOUBLE, "EQUINOX", &dummy_flt, "", &status);
	} else {
		sprintf(dummy_str[0], "ICRS"); //radesys
		fits_write_key(fptr, TSTRING, "RADESYS", &dummy_str[0], "", &status);
	}
	/*
	 fits_write_key(fptr, TDOUBLE, "LONPOLE", &lonpole, "", &status);
	 fits_write_key(fptr, TDOUBLE, "LATPOLE", &latpole, "", &status);
	 //*/
	
	fits_write_record(fptr, card_null, &status);
	char ctype[naxis_PV][9] = {"ANGLE   ", "VRAD    "}, cunit[naxis_PV][9] = {"deg     ", "m/s    "};
	
	double crpix[naxis_PV]; for (int i = 0; i < naxis_PV; ++i) crpix[i] = PVdiag.crpix[i] + 1;
	double crval[] = {PVdiag.crval[PAXIS_PV] / CMperAU / source.distance_pc / double(3600.), PVdiag.crval[VAXIS_PV] / CMperM};
	char crvalComment[][lstr] = {"", ""};
	double cdelt[] = {PVdiag.cdelt[PAXIS_PV] / CMperAU / source.distance_pc / double(3600.), PVdiag.cdelt[VAXIS_PV] / CMperM};
	char cdeltComment[][lstr] = {"", ""};
	sprintf(cdeltComment[0], "%lf au", PVdiag.cdelt[PAXIS_PV] / CMperAU);
	for (int i = 0; i < naxis_PV; ++i) {
		char ckeyname[5][7] = {"CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT"};
		char id[] = "";
		sprintf(id, "%d", i+1);
		for (int j = 0; j < 5; ++j) strcat(ckeyname[j], id);
		fits_write_key(fptr, TSTRING, ckeyname[0], &ctype[i], "", &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[1], &crval[i], crvalComment[i], &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[2], &cdelt[i], cdeltComment[i], &status);
		fits_write_key(fptr, TDOUBLE, ckeyname[3], &crpix[i], "", &status);
		fits_write_key(fptr, TSTRING, ckeyname[4], &cunit[i], "", &status);
	}
	
	fits_write_record(fptr, card_null, &status);
	char rfreqComment[] = "Rest Frequency (Hz) of ";
	strcat(rfreqComment, source.name_line);
	fits_write_key(fptr, TDOUBLE, "RESTFREQ", &source.restfreq_Hz, rfreqComment, &status);
	char specsys[] = "LSRK    ";
	fits_write_key(fptr, TSTRING, "SPECSYS", &specsys, "Spectral reference frame", &status);
	/*
	 fits_write_key(fptr, TSTRING, "ALTRVAL", &nakami_str, "Alternate frequency reference value", &status);
	 fits_write_key(fptr, TSTRING, "ALTRPIX", &nakami_str, "Alternate frequency reference pixel", &status);
	 fits_write_key(fptr, TSTRING, "VELREF", &nakami_str, "1 LSR, 2 HEL, 3 OBS, +256 Radio", &status);
	 //*/
	
	fits_write_record(fptr, card_null, &status);
	char comment_model[lstr] = {"This is a result of an infalling-rotating envelope model with the following physical parameters"};
	fits_write_comment(fptr, comment_model, &status);
	ire_fits_write_record(fptr, source, &status); // Parameters for IRE model

	fits_write_record(fptr, card_null, &status);
	sprintf(comment_model, "This is a result of a PV slice of the above IRE model with the following parameters");
	fits_write_comment(fptr, comment_model, &status);
	pv_fits_write_record(fptr, source, PVdiag, &status); // Parameters for PV diagram

	
	
	
	//---Move the COMMENT by 'fits_create_img' task to backward---//
	fits_write_record(fptr, card_createimg[0], &status);
	fits_write_record(fptr, card_createimg[1], &status);

	
	//---Output the Results---//
	long fpixel[naxis_PV];
	long nelements = 1;
	for (int i = 0; i < naxis_PV; ++i) {
		fpixel[i] = 1;
		nelements *= npix_axis[i];
	}
	

	ire_fitsout(outfilename, fptr, fpixel, nelements, PVdiag.emission_fits, &status);
	
	
	return ; 
}

