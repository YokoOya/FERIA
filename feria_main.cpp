#include "feria.h"


using namespace std;

int main() {
	
	puts("\n----------------------");
	clock_t start, end;
	start = clock();
	//

    //---Input---//
    char buff[lstr][lstr], filename_cubefits[lfilename], filename_PVfits[lfilename], name_line[lstr], name_object[lstr], radesys[lstr], center_ra[lstr], center_dec[lstr], dummy_char;
	double Fldres_as, Velres_kmps, distance_pc, mass_msun, rCB_au, inc_deg, PA_deg, rot_sign, Rout_au, Rin_au;
	double HeightIRE_au, FlareIRE_deg, HeightKEP_au, FlareKEP_deg;
	double DensCB, densityProfileIRE, densityProfileKEP;
    double TempCB, tempProfileIRE, tempProfileKEP;
    double linewidth_kmps, beam_maj_as, beam_min_as, beam_pa_deg;
    double restfreq_GHz, Vsys_kmps;
	int paramnameIDList[lstr] = {}, paramnameIDList_PV[lstr] = {}, velStrID = noGAS;
	bool f_makecube = true, f_overwrite = false, f_norm = false;
	
	double PV_PA_deg, PV_cent_offset_ra_au, PV_cent_offset_dec_au;
	double maxval = -INF, minval = INF;
	
    puts("\nSetting: ");
    for (int i = 0; ~scanf("%[^\n] ", buff[i]); ++i) {
        if (buff[i][0] == '#') {
			printf("\n\t%s\n", buff[i]);
            --i;
            continue;
        } 
        printf("\t%s\n", buff[i]);
    }
    
    int ind = -1, paramnameID = -1;
    sscanf(buff[++ind], "%s ", filename_cubefits);
    sscanf(buff[++ind], "%c ", &dummy_char);
	if (dummy_char == 'y') f_overwrite = true;
	
	sscanf(buff[++ind], "%s ", name_object);
	paramnameIDList[++paramnameID] = ind;	// parameter used for filename_cubefits
	sscanf(buff[++ind], "%s ", radesys);
	sscanf(buff[++ind], "%s ", center_ra);
	sscanf(buff[++ind], "%s ", center_dec);
	sscanf(buff[++ind], "%lf ", &Vsys_kmps);
	paramnameIDList[++paramnameID] = ind;
	
	sscanf(buff[++ind], "%s ", name_line);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &restfreq_GHz);
	
    sscanf(buff[++ind], "%lf ", &Fldres_as);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &Velres_kmps);
	paramnameIDList[++paramnameID] = ind;

    sscanf(buff[++ind], "%lf ", &distance_pc);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &mass_msun);
	paramnameIDList[++paramnameID] = ind;
	
    sscanf(buff[++ind], "%lf ", &rCB_au);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &inc_deg);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &PA_deg);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &rot_sign);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &Rout_au);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &Rin_au);
	paramnameIDList[++paramnameID] = ind;
	
    sscanf(buff[++ind], "%lf ", &HeightIRE_au);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &FlareIRE_deg);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &densityProfileIRE);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &tempProfileIRE);
	paramnameIDList[++paramnameID] = ind;

	sscanf(buff[++ind], "%lf ", &HeightKEP_au);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &FlareKEP_deg);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &densityProfileKEP);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &tempProfileKEP);
	paramnameIDList[++paramnameID] = ind;

    sscanf(buff[++ind], "%lf ", &DensCB);
    sscanf(buff[++ind], "%lf ", &TempCB);

    sscanf(buff[++ind], "%lf ", &linewidth_kmps);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &beam_maj_as);
	paramnameIDList[++paramnameID] = ind;
    sscanf(buff[++ind], "%lf ", &beam_min_as);
	paramnameIDList[++paramnameID] = ind;
	sscanf(buff[++ind], "%lf ", &beam_pa_deg);
	paramnameIDList[++paramnameID] = ind;
    if (beam_maj_as < beam_min_as) {
        double beam_temp = beam_maj_as;
        beam_maj_as = beam_min_as;
		beam_min_as = beam_temp;
		beam_pa_deg = 90 - beam_pa_deg;
    }
	
	sscanf(buff[++ind], "%c ", &dummy_char);
	if (dummy_char == 'y') f_norm = true;
	
	
	puts("\nSize of the Cube Data:");
	printf("\t%d x %d x %d \t\t# in pixel\n", npix, npix, nvel);
	printf("\t%lf arcsec x %lf arcsec x %lf km/s\n", Fldres_as * npix, Fldres_as * npix, Velres_kmps * nvel);
	
	puts("\n----------\n");
    
    //---Parameters---//
    sourceParams source(filename_cubefits, Fldres_as, Velres_kmps, distance_pc, mass_msun, rCB_au, inc_deg, PA_deg, rot_sign, Rout_au, Rin_au, HeightIRE_au, FlareIRE_deg, HeightKEP_au, FlareKEP_deg, DensCB, densityProfileIRE, densityProfileKEP, TempCB, tempProfileIRE, tempProfileKEP, linewidth_kmps, beam_maj_as, beam_min_as, beam_pa_deg, name_line, restfreq_GHz, name_object, radesys, center_ra, center_dec, Vsys_kmps);
    
    //---Output file (cube)---//
	char paramnameList[lstr][lstr] = {"", "-Vsys", "_Line", "_Pix", "as", "kmps_D", "M", "CB", "I", "PA", "Rot", "Rout", "Rin", "_IRE-T", "Flare", "Nprof", "Tprof", "_Kep-T", "Flare", "Nprof", "Tprof", "_LW", "_Beam", "x", "PA"};		// list of the parameters for the name of the filename_cubefits
	
	bool filename_cubefits_params = true;
	for (int i = 0; i < int(strlen(PARAMSNAME)); ++i) {		// is outputfile named after the parameters? or a given name?
		if (filename_cubefits[strlen(filename_cubefits) - strlen(PARAMSNAME) + i] != PARAMSNAME[i]) {
			filename_cubefits_params = false;
			break;
		}
	}
	if (filename_cubefits_params) {
		for (int i = 0; i <= paramnameID; ++i) {			// filename_cubefits = paramnameList[0] + buff[paramnameIDlist[0]] + ...
			int id = paramnameIDList[i];
			sprintf(filename_cubefits, "%s%s", filename_cubefits, paramnameList[i]);
			if (i == 0) sprintf(filename_cubefits, "%s", paramnameList[i]);
			for (int j = 0; j < int(strlen(buff[id])); ++j) {
				if (buff[id][j] == ' ' || buff[id][j] == '\t') break;
				if (buff[id][j] == '(' || buff[id][j] == ')' || buff[id][j] == '=') continue;	// avoid specific characters in the transition (e.g., CO(J=2-1))
				sprintf(filename_cubefits, "%s%c", filename_cubefits, buff[id][j]);
			}
		}
	}
	
	
	//---Output file (PV)---//
	char paramnameList_PV[lstr][lstr] = {"_PV-PA", "deg-CentRA", "Dec"};
	sscanf(buff[++ind], "%lf ", &PV_PA_deg);
	paramnameIDList_PV[0] = ind;
	sscanf(buff[++ind], "%lf ", &PV_cent_offset_ra_au);
	paramnameIDList_PV[1] = ind;
	sscanf(buff[++ind], "%lf ", &PV_cent_offset_dec_au);
	paramnameIDList_PV[2] = ind;

	sprintf(filename_PVfits, "%s", filename_cubefits);
	for (int i = 0; i < 3; ++i) {
		int id = paramnameIDList_PV[i];
		sprintf(filename_PVfits, "%s%s", filename_PVfits, paramnameList_PV[i]);
		for (int j = 0; j < int(strlen(buff[id])); ++j) {
			if (buff[id][j] == ' ' || buff[id][j] == '\t') break;
			if (buff[id][j] == '(' || buff[id][j] == ')' || buff[id][j] == '=') continue;	// avoid specific characters in the transition (e.g., CO(J=2-1))
			sprintf(filename_PVfits, "%s%c", filename_PVfits, buff[id][j]);
		}
	}
	
		 
		 
    for (int i = 0; i < int(strlen(FITS)); ++i) {			// filename_cubefits should end with ".fits"
        if (filename_cubefits[strlen(filename_cubefits) - strlen(FITS) + i] != FITS[i]) {
            strcat(filename_cubefits, FITS);
            break;
        }
    }
	for (int i = 0; i < int(strlen(FITS)); ++i) {			// filename_cubefits should end with ".fits"
		if (filename_PVfits[strlen(filename_PVfits) - strlen(FITS) + i] != FITS[i]) {
			strcat(filename_PVfits, FITS);
			break;
		}
	}

	
    if (fopen(filename_cubefits, "r") != NULL) {				// does the outputfile exist? if so, overwrite or not?
		printf("\n\nThere already exists the fits file:\n\t'%s'\n", filename_cubefits);
        if (f_overwrite) {
            puts("\n\t-> The fits file will be replaced.");
        } else {
            puts("\n\t-> The fits file will be read.");
			f_makecube = false;
        }
    }
	
	
	//---Make the plane of the sky---//
	double skyCrpix[naxis], skyCrval[naxis], skyCdelt[naxis];
	for (int i = 0; i < naxis; ++i) {		// make PPV data
		skyCrpix[i] = (i == VAXIS ? nvel / 2 - 1: npix / 2 - 1);
		skyCrval[i] = 0.;
		skyCdelt[i] = (i == VAXIS ? Velres_kmps : Fldres_as * distance_pc * pow(-1, i + 1));
	}
	
	skyPlane sky(distance_pc, skyCrpix, skyCrval, skyCdelt, linewidth_kmps, beam_maj_as, beam_min_as, beam_pa_deg);
	

	if (!f_makecube) {
		//---Read fits file---//
		printf("\n\nRead the fits file:\n\t'%s'\n", filename_cubefits);
		ire_cube_fitsin(filename_cubefits, sky, source, &f_makecube);
		
		sky.setMaxMin();
		sky.getMaxMin(&maxval, &minval);
		//printf("\n\t\tCube: max = %e, min = %e, sum = %e\n", maxval, minval, sky.getSum());

	}
	if (f_makecube) {
		//---Make 3D IRE---//
		
		printf("\n\nMake a fits file:\n\t'%s'\n", filename_cubefits);
		ireModel ire_model(mass_msun, rCB_au, inc_deg, PA_deg, rot_sign, Rout_au, Rin_au, HeightIRE_au, FlareIRE_deg, HeightKEP_au, FlareKEP_deg, DensCB, densityProfileIRE, densityProfileKEP, TempCB, tempProfileIRE, tempProfileKEP);
		
		double meshCrpix[naxis], meshCrval_au[naxis], meshCdelt_au[naxis];
		for (int i = 0; i < naxis; ++i) {
			meshCrpix[i] = npix / 2. - 1.;
			meshCrval_au[i] = 0.;
			meshCdelt_au[i] = Fldres_as * distance_pc; //Fldsize_as * distance_pc / npix;
		}
		mesh ire_mesh(meshCrpix, meshCrval_au, meshCdelt_au, inc_deg, PA_deg);
		
		double pos_polar[naxis] = {}, vel_polar[naxis] = {}, dist_protostar, density, temperature;
		
		for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) for (int iz = 0; iz < npix; ++iz) {		// make PPP data
			ire_mesh.getPos_polar(ix, iy, iz, pos_polar);
			dist_protostar = sqrt(pow(pos_polar[rAXIS], 2) + pow(pos_polar[ZAXIS], 2));
			
			velStrID = ire_model.which_velStr(dist_protostar, pos_polar[rAXIS], pos_polar[tAXIS], pos_polar[ZAXIS]);
			ire_model.env_vel(dist_protostar, pos_polar[rAXIS], pos_polar[tAXIS], pos_polar[ZAXIS], vel_polar, velStrID);
			density = ire_model.env_dens(dist_protostar, pos_polar[rAXIS], pos_polar[tAXIS], pos_polar[ZAXIS], velStrID);
			temperature = ire_model.env_temp(dist_protostar, pos_polar[rAXIS], pos_polar[tAXIS], pos_polar[ZAXIS], velStrID);
			
			ire_mesh.setVel_polar(ix, iy, iz, pos_polar, vel_polar);
			ire_mesh.setData(ix, iy, iz, NVAL, density);
			ire_mesh.setData(ix, iy, iz, TVAL, temperature);
		}
		
		for (int i = 0; i < ndata; ++i) ire_mesh.setMaxMin(i);
		
		/*Comment
		//printf("\ncdelt = %lf au\n", meshCdelt_au[0]);
		puts("\n\tValues in the PPP data (0:Vx (km/s), 1:Vy (km/s), 2:Vz (km/s), 3:N (cm-3), 4:T (K))");
		for (int i = 0; i < ndata; ++i) {
			ire_mesh.getMaxMin(&maxval, &minval, i);
			printf("\n\t\tvalue %d -> max = %e, min = %e", i, maxval / (i < 3 ? CMperKM : 1.), minval / (i < 3 ? CMperKM : 1.));
		}
		puts("");
		//*/
		
		
		//---Projection---//
		double gasColumn[npix][ndata];	// along z-axis
		for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) {		// make PPV data
			ire_mesh.getColumn(ix, iy, gasColumn);
			sky.projection(ix, iy, gasColumn);
		}

		/*Comment
		puts("\n\tIntensity in the PPV data:");
		double max1, min1;
		sky.setMaxMin();
		sky.getMaxMin(&max1, &min1);
		printf("\t\tCube (before convolution): max = %e, min = %e, sum = %e\n", max1, min1, sky.getSum());
		//*/
		
		//---Convolution---//
		for (int iv = 0; iv < nvel; ++iv) sky.convolution_beam(iv);
		for (int ix = 0; ix < npix; ++ix) for (int iy = 0; iy < npix; ++iy) sky.convolution_line(ix, iy);
		
		
		/*Comment
		double max2, min2;
		sky.setMaxMin();
		sky.getMaxMin(&max2, &min2);
		printf("\t\tCube (after convolution): max = %e, min = %e, sum = %e\n\n", max2, min2, sky.getSum());
		//*/
		
		sky.setMaxMin();
		//*Comment
		sky.getMaxMin(&maxval, &minval);
		printf("\n\t\tCube (before normalize): max = %e, min = %e, sum = %e\n", maxval, minval, sky.getSum());
		//*/
		
		
		if (f_norm) sky.normalize();		// normalize the intensities by the peak intensity
		
		
		sky.setMaxMin();
		/*Comment
		sky.getMaxMin(&maxval, &minval);
		printf("\n\t\tCube (after normalize): max = %e, min = %e, sum = %e\n\n", maxval, minval, sky.getSum());
		//*/
		

		
		//---Output fits file---//
		sky.setEmissionFits();	// match the data format (the order of the axes and pixels) to the fits format
		ire_cube_fitsout(filename_cubefits, sky, source);
		
		
	}
	sky.setCrval(VAXIS, Vsys_kmps);

	
	
	//---PV diagram---//
	printf("\n----------\n\nMake a PV diagram:\n\t'%s'\n", filename_PVfits);
	
	
	double PV_PA_rad = PV_PA_deg / 180. * pi;
	int npix_PV = ceil(npix / fmax(fabs(cos(PV_PA_rad)), fabs(sin(PV_PA_rad))) * (1. - db_erratio));
	if (!(npix_PV % 2)) npix_PV += 1;
	double posList[naxis] = {0., 0., 0.}, indList[naxis] = {0., 0., 0.};
	
	double PVline_posList[npix_PV][naxis];
	double dra_au = Fldres_as * distance_pc * sin(PV_PA_rad), ddec_au = Fldres_as * distance_pc * cos(PV_PA_rad) * -1.;
	
	
	double pvCrpix[naxis], pvCrval[naxis], pvCdelt[naxis];
	for (int i = 0; i < naxis_PV; ++i) {		// make PPV data
		pvCrpix[i] = (i == VAXIS_PV ? nvel / 2 - 1: npix_PV / 2 - 1);
		pvCrval[i] = (i == VAXIS_PV ? Vsys_kmps : 0.);
		pvCdelt[i] = (i == VAXIS_PV ? Velres_kmps : sqrt(dra_au * dra_au + ddec_au * ddec_au));
	}
	PVdiagram PVdiag(npix_PV, pvCrpix, pvCrval, pvCdelt, source, PV_cent_offset_ra_au, PV_cent_offset_dec_au, PV_PA_deg);
	
	
	sky.getPosList(indList, posList);
	PVline_posList[0][VAXIS] = posList[VAXIS];
	
	for (int ip = 0; ip < npix_PV; ++ip) {
		PVline_posList[ip][XAXIS] = PV_cent_offset_ra_au * -1. + dra_au * (ip - floor(npix_PV / 2));
		PVline_posList[ip][YAXIS] = PV_cent_offset_dec_au + ddec_au * (ip - floor(npix_PV / 2));
		PVline_posList[ip][VAXIS] = PVline_posList[0][VAXIS];
	}
	
	for (int ip = 0; ip < npix_PV; ++ip) {
		for (int i = 0; i < naxis; ++i) posList[i] = PVline_posList[ip][i];
		sky.getIndexList(indList, posList);
		for (int iv = 0; iv < nvel; ++iv) {
			indList[VAXIS] = sky.getIndex(VAXIS, PVdiag.getPos(VAXIS_PV, double(iv)));
			
			PVdiag.setVal(ip, iv, sky.getVal(indList));
		}
	}
	
	
	PVdiag.setMaxMin();
	/*Comment
	PVdiag.getMaxMin(&maxval, &minval);
	printf("\n\t\tPV diagram: max = %e, min = %e, sum = %e\n\n", maxval, minval, PVdiag.getSum());
	//*/
	
	
	//---Output fits file---//
	PVdiag.setEmissionFits();	// match the data format (the order of the axes and pixels) to the fits format
	ire_PV_fitsout(filename_PVfits, PVdiag, source, npix_PV);
	

	//
	end = clock();
	puts("\n----------------------");
	printf("Duration time: %lf [sec]\n", double(end - start) / CLOCKS_PER_SEC);
    puts("Done.\n\n");
	
	
	puts("\007"); 
	return 0;
}

