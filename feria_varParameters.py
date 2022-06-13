import os
import sys
import glob
import subprocess


fname_cubein = "feria.in"
fname_exec = "./feria"


Obj = 'IRAS16293-2422_A1'
Ra = '16h32m22.878s'
Dec = '-24d28m36.69s'
Vsys = '2.5'

line = 'H2CS707'
restfreq = '240.2668724'

pixSize = '0.004'
velRes = '0.1'


DList = ['140']
MList = ['0.1', '0.3', '1.0']
CBList = ['30.', '100.', '300.']
IList = ['0.', '30.', '60.', '90.']
PAList = ['50']
RotList = ['1']
RoutList = ['300.']
RinList = ['0']

ireHeightList = ['0']
ireFlareList = ['30']
ireNProfList = ['-1.5']
ireTProfList = ['0.0']

kepHeightList = ['0']
kepFlareList = ['30']
kepNProfList = ['-1.5']
kepTProfList = ['0.0']


CBDensList = ['1e-2']
CBTempList = ['10']

LWList = ['1.0']
BmajList = ['0.118']
BminList = ['0.077']
BpaList = ['-85.842']

PVpaList = ['50']	#, '80', '110', '140', '170', '200']
PVraList = ['0.']
PVdecList = ['0.']



####################

def makeCube(D, M, CB, I,
			PA, Rot,
			Rout, Rin,
			ireHeight, ireFlare, ireNProf, ireTProf,
			kepHeight, kepFlare, kepNProf, kepTProf,
			CBD, CBT,
			LW, Bmaj, Bmin, Bpa,
			PVpa, PVra, PVdec):
	f_cubein = open(fname_cubein, 'w')

	f_cubein.write("parameter\t# Name of the output fits file, 'parameter' to name after the given parameters\n")
	f_cubein.write("n\t\t# Overwrite? (y/n)\n")

	f_cubein.write("\n### Parameters for the source ###\n")
	f_cubein.write("%s\t# Name of the object\n" % Obj)
	f_cubein.write("ICRS\t# Equinox or Radesys ('J2000' or 'ICRS')\n")
	f_cubein.write("%s\t# Right Ascension of the field center\n" % Ra)
	f_cubein.write("%s\t# Declination of the field center\n" % Dec)
	f_cubein.write("%s\t# Systemic Velocity (km/s)\n" % Vsys)

	f_cubein.write("\n### Parameters for the line ###\n")
	f_cubein.write("%s\t# Name of the line\n" % line)
	f_cubein.write("%s\t# Restfrequency of the line (GHz)\n" % restfreq)

	f_cubein.write("\n### Parameters for the field ###\n")
	f_cubein.write("%s\t# Pixel size (arcsec)\n" % pixSize)
	f_cubein.write("%s\t# Velocity resolution (km/s)\n" % velRes)

	f_cubein.write("\n### Parameters for IRE ###\n")
	f_cubein.write("%s\t# Distance to the source (pc)\n" % D)
	f_cubein.write("%s\t# Protostellar mass (Msun)\n" % M)
	f_cubein.write("%s\t# Radius of the centrifugal barrier (au)\n" % CB)
	f_cubein.write("%s\t# Inclination angle of the disk/envelope system (degree; 0 deg for a face-on configuration)\n" % I)
	f_cubein.write("%s\t# Position angle of the line along which the disk/envelope system is extended\n" % PA)
	f_cubein.write("%s\t# Direction of the rotation: -1 for counterclockwise, 1 for clockwise for the inclination angle of 0 deg\n" % Rot)
	f_cubein.write("%s\t# Outer radius of the infalling-rotating envelope\n" % Rout)
	if Rin == 'CB':
		f_cubein.write("%s\t# Inner radius of the infalling-rotating envelope\n" % CB)
	else:
		f_cubein.write("%s\t# Inner radius of the infalling-rotating envelope\n" % Rin)

	f_cubein.write("%s\t# Scale height of the infalling-rotating envelope\n" % ireHeight)
	f_cubein.write("%s\t# Flared degree of the scale hight of the infalling-rotating envelope\n" % ireFlare)
	f_cubein.write("%s\t# Density profile (IRE): log(n(X)) / log(r)\n" % ireNProf)
	f_cubein.write("%s\t# Temperature profile (IRE): log(T) / log(r)\n" % ireTProf)
	f_cubein.write("%s\t# Scale height of the Keplerian disk\n" % kepHeight)
	f_cubein.write("%s\t# Flared degree of the scale hight of the Keplerian disk\n" % kepFlare)
	f_cubein.write("%s\t# Density profile (Kepler): log(n(X)) / log(r)\n" % kepNProf)
	f_cubein.write("%s\t# Temperature profile (Kepler): log(T) / log(r)\n" % kepTProf)
	f_cubein.write("%s\t# Molecular density at the centrifugal barrier (cm^-3)\n" % CBD)
	f_cubein.write("%s\t# Temperature at the centrifugal barrier (K)\n" % CBT)
	f_cubein.write("%s\t# Linewidth (km/s) (FWHM)\n" % LW)
	f_cubein.write("%s\t# Beam size: major axis (arcsec)\n" % Bmaj)
	f_cubein.write("%s\t# Beam size: minor axis (arcsec)\n" % Bmin)
	f_cubein.write("%s\t# Beam: position angle of the major axis (degree)\n" % Bpa)
	f_cubein.write("y\t# Normalize the intensity by the maximum value in the cube? (y/n)\n")

	f_cubein.write("\n### the position axis in the PV diagram ###\n")
	f_cubein.write("%s\t# position angle of the position axis (deg)\n" % PVpa)
	f_cubein.write("%s\t# offset in RA for the origin of the position axis (au)\n" % PVra)
	f_cubein.write("%s\t# offset in Dec for the origin of the position axis (au)\n" % PVdec)

	f_cubein.close()
	
	subprocess.call("%s < %s" % (fname_exec, fname_cubein), shell = True)





nCube = 1
nPV = 1
iCube = 0
iPV = 0
nCubeidList = [len(DList), len(MList), len(CBList), len(IList),
				len(PAList), len(RotList),
				len(RoutList), len(RinList),
				len(ireHeightList), len(ireFlareList), len(ireNProfList), len(ireTProfList),
				len(kepHeightList), len(kepFlareList), len(kepNProfList), len(kepTProfList),
				len(CBDensList), len(CBTempList),
				len(LWList), len(BmajList), len(BminList), len(BpaList)]
nPVidList = [len(PVpaList), len(PVraList), len(PVdecList)]

for item in nCubeidList:
	nCube = nCube * item
for item in nPVidList:
	nPV = nPV * item




Cubeid = [0 for i in range(len(nCubeidList))]

for iCube in range(nCube):
	print "\n<<<<< Make a Cube (%d / %d) >>>>>" % (iCube + 1, nCube)
	print Cubeid
	print ""
	
	PVid = [0 for i in range(len(nPVidList))]
	for iPV in range(nPV):
		
		print "<<< Make a PV (%d / %d) >>>\n" % (iPV + 1, nPV)
		print PVid
		
		makeCube(DList[Cubeid[0]], MList[Cubeid[1]], CBList[Cubeid[2]], IList[Cubeid[3]],
				PAList[Cubeid[4]], RotList[Cubeid[5]],
				RoutList[Cubeid[6]], RinList[Cubeid[7]],
				ireHeightList[Cubeid[8]], ireFlareList[Cubeid[9]], ireNProfList[Cubeid[10]], ireTProfList[Cubeid[11]],
				kepHeightList[Cubeid[12]], kepFlareList[Cubeid[13]], kepNProfList[Cubeid[14]], kepTProfList[Cubeid[15]],
				CBDensList[Cubeid[16]], CBTempList[Cubeid[17]],
				LWList[Cubeid[18]], BmajList[Cubeid[19]], BminList[Cubeid[20]], BpaList[Cubeid[21]],
				PVpaList[PVid[0]], PVraList[PVid[1]], PVdecList[PVid[2]])

		PVid[-1] = PVid[-1] + 1
		for i in range(1, len(nPVidList), 1):
			if PVid[-i] == nPVidList[-i]:
				PVid[-i] = 0
				PVid[-(i + 1)] = PVid[-(i + 1)] + 1


	Cubeid[-1] = Cubeid[-1] + 1
	for i in range(1, len(nCubeidList), 1):
		if Cubeid[-i] == nCubeidList[-i]:
			Cubeid[-i] = 0
			Cubeid[-(i + 1)] = Cubeid[-(i + 1)] + 1



print "\n\007\007"
print "\n\n\n----------\n'%s' has been done.\n\n\n" % sys.argv[0]


