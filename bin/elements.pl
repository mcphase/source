#!/usr/bin/perl
# 
# elements.pl
#
# Datatables for ions and elements from various sources. For use by cif2mcphas to generate McPhase inputfiles
#
# Data are stored as hashtables:
#
# element{} - Has basic information including: scattering lenghts, if it is magnetic, and possible valence states
# magions{} - Has information for each magnetic ions: spin S (High/Low spin) or total "J", radial integral <R^k> k=2,4,6
# magff{}   - Magnetic form factors for each magnetic ions
#
# This program is part of the McPhase package, licensed under the GNU GPL v2. Please see the COPYING file
#
# Sat Sep 27 17:15:19 KST 2014 - Duc Le - mducle@snu.ac.kr

# Scattering lengths, guess whether magnetic and possible valence states (from CHEMIX) of the elements
# Ionic radius and colours are taken from the file elements.ini from VESTA ( http://jp-minerals.org/vesta/ )
# Name, real_b, imag_b, ismag,  valence,	r_ion,	colour
%element = (
"H",	[-3.739,0,	0,	"1,-1",		0.46,	[255,204,204]],
"D",	[6.671,	0,	0,	"1,-1",		0.46,	[204,204,255]],
"He",	[3.26,	0,	0,	0,		1.22,	[252,232,206]],
"Li",	[-1.90,	0,	0,	1,		1.57,	[134,224,116]],
"Be",	[7.79,	0,	0,	2,		1.12,	[ 94,215,123]],
"B",	[5.30,	-0.213,	0,	3,		0.81,	[ 31,162, 15]],
"C",	[6.646,	0,	0,	"4,-4,2",	0.77,	[128, 73, 41]],
"N",	[9.36,	0,	0,	"3,-3,5,4,2",	0.74,	[176,185,230]],
"O",	[5.803,	0,	0,	-2,		0.74,	[254,  3,  0]],
"F",	[5.654,	0,	0,	-1,		0.72,	[176,185,230]],
"Ne",	[4.566,	0,	0,	0,		1.60,	[254, 55,181]],
"Na",	[3.63,	0,	0,	1,		1.91,	[249,220, 60]],
"Mg",	[5.375,	0,	0,	2,		1.60,	[251,123, 21]],
"Al",	[3.449,	0,	0,	3,		1.43,	[129,178,214]],
"Si",	[4.1491,0,	0,	4,		1.18,	[ 27, 59,250]],
"P",	[5.13,	0,	0,	"5,-3,3,4",	1.10,	[192,156,194]],
"S",	[2.847,	0,	0,	"6,-2,2,4",	1.04,	[255,250,  0]],
"Cl",	[9.577,	0,	0,	"1,-1,3,5,7",	0.99,	[ 49,252,  2]],
"Ar",	[1.909,	0,	0,	0,		1.92,	[207,254,196]],
"K",	[3.67,	0,	0,	1,		2.35,	[161, 33,246]],
"Ca",	[4.70,	0,	0,	2,		1.97,	[ 90,150,189]],
"Sc",	[12.29,	0,	0,	3,		1.64,	[181, 99,171]],
"Ti",	[-3.438,0,	0,	"4,3",		1.47,	[120,202,255]],
"V",	[-.3824,0,	1,	"5,4,3,2",	1.35,	[229, 25,  0]],
"Cr",	[3.635,	0,	1,	"3,2,6",	1.29,	[  0,  0,158]],
"Mn",	[-3.73,	0,	1,	"2,3,4,6,7",	1.37,	[168,  8,158]],
"Fe",	[9.45,	0,	1,	"2,3",		1.26,	[181,113,  0]],
"Co",	[2.49,	0,	1,	"2,3",		1.25,	[  0,  0,175]],
"Ni",	[10.3,	0,	1,	"2,3",		1.25,	[183,187,189]],
"Cu",	[7.718,	0,	1,	"2,1",		1.28,	[ 34, 71,220]],
"Zn",	[5.680,	0,	0,	2,		1.37,	[143,143,129]],
"Ga",	[7.288,	0,	0,	3,		1.53,	[158,227,115]],
"Ge",	[8.185,	0,	0,	4,		1.22,	[126,110,166]],
"As",	[6.58,	0,	0,	"3,-3,5",	1.21,	[116,208, 87]],
"Se",	[7.970,	0,	0,	"4,-2,6",	1.04,	[154,239, 15]],
"Br",	[6.795,	0,	0,	"1,-1,5",	1.14,	[126, 49,  2]],
"Kr",	[7.81,	0,	0,	0,		1.98,	[250,193,243]],
"Rb",	[7.09,	0,	0,	1,		2.50,	[255,  0,153]],
"Sr",	[7.02,	0,	0,	2,		2.15,	[  0,255, 38]],
"Y",	[7.75,	0,	0,	3,		1.82,	[102,152,142]],
"Zr",	[7.16,	0,	0,	4,		1.60,	[  0,255,  0]],
"Nb",	[7.054,	0,	0,	"5,3",		1.47,	[ 76,178,118]],
"Mo",	[6.715,	0,	0,	"6,5,4,3,2",	1.40,	[179,134,175]],
"Tc",	[6.8,	0,	0,	7,		1.35,	[205,175,202]],
"Ru",	[7.03,	0,	1,	"3,4,6,8,2",	1.34,	[207,183,173]],
"Rh",	[5.88,	0,	1,	"3,2,4",	1.34,	[205,209,171]],
"Pd",	[5.91,	0,	0,	"2,4",		1.37,	[193,195,184]],
"Ag",	[5.922,	0,	0,	1,		1.44,	[183,187,189]],
"Cd",	[4.87,	-0.70,	0,	2,		1.52,	[242, 30,220]],
"In",	[4.065,	-0.0539,0,	3,		1.67,	[215,128,187]],
"Sn",	[6.225,	0,	0,	"4,2",		1.58,	[154,142,185]],
"Sb",	[5.57,	0,	0,	"3,-3,5",	1.41,	[215,131, 79]],
"Te",	[5.80,	0,	0,	"4,-2,6",	1.37,	[173,162, 81]],
"I",	[5.28,	0,	0,	"1,-1,5,7",	1.33,	[142, 31,138]],
"Xe",	[4.92,	0,	0,	0,		2.18,	[154,161,248]],
"Cs",	[5.42,	0,	0,	1,		2.72,	[ 14,254,185]],
"Ba",	[5.07,	0,	0,	2,		2.24,	[ 30,239, 44]],
"La",	[8.24,	0,	0,	3,		1.88,	[ 90,196, 73]],
"Ce",	[4.84,	0,	1,	"3,4",		1.82,	[209,252,  6]],
"Pr",	[4.58,	0,	1,	"3,4",		1.82,	[252,225,  5]],
"Nd",	[7.69,	0,	1,	3,		1.82,	[251,141,  6]],
"Pm",	[12.6,	0,	1,	3,		1.81,	[  0,  0,244]],
"Sm",	[0.80,	-1.65,	1,	"3,2",		1.81,	[252,  6,125]],
"Eu",	[7.22,	-1.26,	1,	"3,2",		2.06,	[250,  7,213]],
"Gd",	[6.5,	-13.82,	1,	3,		1.79,	[192,  3,255]],
"Tb",	[7.38,	0,	1,	"3,4",		1.77,	[113,  4,254]],
"Dy",	[16.9,	-0.276,	1,	3,		1.77,	[ 49,  6,252]],
"Ho",	[8.01,	0,	1,	3,		1.76,	[  7, 65,251]],
"Er",	[7.79,	0,	1,	3,		1.75,	[ 73,114, 58]],
"Tm",	[7.07,	0,	1,	"3,2",		1.00,	[  0,  0,224]],
"Yb",	[12.43,	0,	1,	"3,2",		1.94,	[ 39,252,244]],
"Lu",	[7.21,	0,	0,	3,		1.72,	[ 38,253,181]],
"Hf",	[7.7,	0,	0,	4,		1.59,	[180,179, 89]],
"Ta",	[6.91,	0,	0,	5,		1.47,	[183,154, 86]],
"W",	[4.86,	0,	0,	"6,5,4,3,2",	1.41,	[141,138,127]],
"Re",	[9.2,	0,	0,	"7,-1,6,4,2",	1.37,	[179,176,142]],
"Os",	[10.7,	0,	1,	"4,2,3,6,8",	1.35,	[200,177,120]],
"Ir",	[10.6,	0,	1,	"4,2,3,6",	1.36,	[201,206,114]],
"Pt",	[9.60,	0,	0,	"4,2",		1.39,	[203,197,191]],
"Au",	[7.63,	0,	0,	"3,1",		1.44,	[254,178, 56]],
"Hg",	[12.692,0,	0,	"2,1",		1.55,	[211,183,203]],
"Tl",	[8.776,	0,	0,	"1,3",		1.71,	[149,137,108]],
"Pb",	[9.405,	0,	0,	"2,4",		1.75,	[ 82, 83, 91]],
"Bi",	[8.532,	0,	0,	"3,5",		1.82,	[210, 47,247]],
"Po",	["",	"",	0,	"4,2",		1.77,	[  0,  0,255]],
"At",	["",	"",	0,	"1,-1,3,5,7",	0.62,	[  0,  0,255]],
"Rn",	["",	"",	0,	0,		0.80,	[255,255,  0]],
"Fr",	["",	"",	0,	1,		1.00,	[  0,  0,  0]],
"Ra",	[10.0,	0,	0,	2,		2.35,	[109,169, 88]],
"Ac",	["",	"",	0,	3,		2.03,	[100,158,114]],
"Th",	[10.31,	0,	0,	4,		1.80,	[ 37,253,120]],
"Pa",	[9.1,	0,	0,	"5,4",		1.63,	[ 41,250, 53]],
"U",	[8.417,	0,	1,	"6,5,4,3",	1.56,	[121,161,170]],
"Np",	[10.55,	0,	1,	"5,4,3,6",	1.56,	[ 76, 76, 76]],
"Pu",	["",	"",	1,	"4,3,5,6",	1.64,	[ 76, 76, 76]],
"Am",	[8.3,	0,	1,	"3,4,5,6",	1.73,	[ 76, 76, 76]],
"Cm",	["",	"",	1,	3,		0.80,	[ 76, 76, 76]]
);

# atomic mass
#  Element Symbol AtomicMass 
%mass = (
"H", [1.007], 
"He", [4.002], 
"Li", [6.941], 
"Be", [9.012], 
"B", [10.811], 
"C", [12.011], 
"N", [14.007], 
"O", [15.999], 
"F", [18.998], 
"Ne", [20.18], 
"Na", [22.99], 
"Mg", [24.305], 
"Al", [26.982], 
"Si", [28.086], 
"P", [30.974], 
"S", [32.065], 
"Cl", [35.453], 
"Ar", [39.948], 
"K", [39.098], 
"Ca", [40.078], 
"Sc", [44.956], 
"Ti", [47.867], 
"V", [50.942], 
"Cr", [51.996], 
"Mn", [54.938], 
"Fe", [55.845], 
"Co", [58.933], 
"Ni", [58.693], 
"Cu", [63.546], 
"Zn", [65.38], 
"Ga", [69.723], 
"Ge", [72.64], 
"As", [74.922], 
"Se", [78.96], 
"Br", [79.904], 
"Kr", [83.798], 
"Rb", [85.468], 
"Sr", [87.62], 
"Y", [88.906], 
"Zr", [91.224], 
"Nb", [92.906], 
"Mo", [95.96], 
"Tc", [98], 
"Ru", [101.07], 
"Rh", [102.906], 
"Pd", [106.42], 
"Ag", [107.868], 
"Cd", [112.411], 
"In", [114.818], 
"Sn", [118.71], 
"Sb", [121.76], 
"Te", [127.6], 
"I", [126.904], 
"Xe", [131.293], 
"Cs", [132.905], 
"Ba", [137.327], 
"La", [138.905], 
"Ce", [140.116], 
"Pr", [140.908], 
"Nd", [144.242], 
"Pm", [145], 
"Sm", [150.36], 
"Eu", [151.964], 
"Gd", [157.25], 
"Tb", [158.925], 
"Ey", [162.5], 
"Ho", [164.93], 
"Er", [167.259], 
"Tm", [168.934], 
"Yb", [173.054], 
"Lu", [174.967], 
"Hf", [178.49], 
"Ta", [180.948], 
"W", [183.84], 
"Re", [186.207], 
"Os", [190.23], 
"Ir", [192.217], 
"Pt", [195.084], 
"Au", [196.967], 
"Hg", [200.59], 
"Tl", [204.383], 
"Pb", [207.2], 
"Bi", [208.98], 
"Po", [210], 
"At", [210], 
"Rn", [222], 
"Fr", [223], 
"Ra", [226], 
"Ac", [227], 
"Th", [232.038], 
"Pa", [231.036], 
"U", [238.029], 
"Np", [237], 
"Pu", [244], 
"Am", [243], 
"Cm", [247], 
"Bk", [247], 
"Cf", [251], 
"Es", [252], 
"Fm", [257], 
"Md", [258], 
"No", [259], 
"Lr", [262], 
"Rf", [261], 
"Eb", [262], 
"Sg", [266], 
"Bh", [264], 
"Hs", [267], 
"Mt", [268], 
"Es", [271], 
"Rg", [272], 
"Cn", [285], 
"Nh", [284], 
"Fl", [289], 
"Mc", [288], 
"Lv", [292], 
"Ts", [295], 
"Og", [294]
);


# Magnetic atoms table
# Name, conf,   S=HS,   S=LS,   J=?,    GJ,	R2, 	R4, 	R6, 	zeta,		F^2,		F^4,		F^6
# Radial expectation values are from CFIELD(theta.c) by Peter Fabi [Rare Earths] and Maurits Haverkort's thesis.
# Spin orbit zeta and Coulomb Fk parameters (in meV) taken from ic1ion (mostly from Carnall et al. J Chem Phys)
$a0=0.5292; $a2=$a0**2; $a4=$a0**4;
%magions = (
"V2+",	["3d3",	1.5,	"",	"",	2.,	0.595/$a2,0.793/$a4,0,	20.4574,	7131.45,	4460.08,	0],
"V3+",	["3d2",	1,	"",	"",	2.,	0.470/$a2,0.456/$a4,0,	26.2846,	8845.53,	6506.57,	0],
"V4+",	["3d1",	0.5,	"",	"",	2.,	0.392/$a2,0.299/$a4,0,	30.7481,	0,		0,		0],
"Cr2+",	["3d4",	2.,	1,	"",	2.,	0.515/$a2,0.605/$a4,0,	28.7643,	8019.3,		5358.35,	0],
"Cr3+",	["3d3",	1.5,	"",	"",	2.,	0.416/$a2,0.362/$a4,0,	33.8477,	9598.86,	6014.47,	0],
"Cr4+",	["3d2",	1,	"",	"",	2.,	0.351/$a2,0.243/$a4,0,	40.6668,	9990.27,	6620.61,	0],
"Mn2+",	["3d5",	2.5,	"",	"",	2.,	0.452/$a2,0.475/$a4,0,	42.508,		8717.95,	5194.32,	0],
"Mn3+",	["3d4",	2.,	1,	"",	2.,	0.371/$a2,0.293/$a4,0,	43.6424,	10115.2,	5741.09,	0],
"Mn4+",	["3d3",	1.5,	"",	"",	2.,	0.316/$a2,0.201/$a4,0,	49.8416,	10792.1,	6734.32,	0],
"Mn6+",	["3d1",	1,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Fe2+",	["3d6",	2.,	0,	"",	2.,	0.402/$a2,0.380/$a4,0,	51.0815,	9110.23,	4828.76,	0],
"Fe3+",	["3d5",	2.5,	"",	"",	2.,	0.334/$a2,0.240/$a4,0,	295.082,	12042.6,	7534.39,	0],
"Co2+",	["3d7",	1.5,	0.5,	"",	2.,	0.360/$a2,0.310/$a4,0,	66.2076,	10563.1,	6820.57,	0],
"Co3+",	["3d6",	2.,	0,	"",	2.,	0.302/$a2,0.200/$a4,0,	63.7775,	10913.7,	7998.47,	0],
"Co4+",	["3d5",	2.5,	"",	"",	2.,	0.262/$a2,0.141/$a4,0,	0,		0,		0,		0],
"Ni2+",	["3d8",	1,	0,	"",	2.,	0.325/$a2,0.256/$a4,0,	80.3418,	10778.3,	7546.99,	0],
"Ni3+",	["3d7",	1.5,	0.5,	"",	2.,	0.275/$a2,0.168/$a4,0,	101.171,	11623.6,	7975.03,	0],
"Cu2+",	["3d9",	0.5,	0,	"",	2.,	0.295/$a2,0.214/$a4,0,	102.907,	11564.6,	7278.29,	0],
"Ru+",	["4d7",	1.5,	0.5,	"",	2.,	0.890/$a2,1.605/$a4,0,	0,		0,		0,		0],
"Ru2+",	["4d6",	2.,	0,	"",	2.,	0.762/$a2,1.068/$a4,0,	133.531,	4158.61,	2717.67,	0],
"Ru3+",	["4d5",	2.5,	0,	"",	2.,	0.674/$a2,0.790/$a4,0,	148.409,	4445.81,	2928.44,	0],
"Ru4+",	["4d4",	2.,	1,	"",	2.,	0.609/$a2,0.619/$a4,0,	167.379,	3532.31,	1983.99,	0],
"Ru6+",	["4d2",	1,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Rh+",	["4d8",	1,	0,	"",	2.,	0.787/$a2,1.257/$a4,0,	0,		4050.53,	2621.54,	0],
"Rh2+",	["4d7",	1.5,	0.5,	"",	2.,	0.683/$a2,0.861/$a4,0,	206.31,		4364.48,	2847.64,	0],
"Rh3+",	["4d6",	2.,	0,	"",	2.,	0.610/$a2,0.649/$a4,0,	175.562,	4633.24,	3047.1,		0],
"Rh4+",	["4d5",	2.5,	0,	"",	2.,	0.554/$a2,0.515/$a4,0,	194.655,	0,		0,		0],
"Ce2+",	["4f2",	"",	"",	"Pr3+",	4./5,	0,	0,	0,	0,		0,		0,		0],
"Ce3+",	["4f1",	"",	"",	"Ce3+",	6./7,	1.309,	3.964,	23.31,	80.255,		0,		0,		0],
"Pr3+",	["4f2",	"",	"",	"Pr3+",	4./5,	1.1963,	3.3335,	18.353,	93.1989,	8539.78,	6242.23,	4079.2],
"Pr4+",	["4f1",	"",	"",	"Ce3+",	6./7,	0,	0,	0,	0,		0,		0,		0],
"Nd2+",	["4f4",	"",	"",	"Nd2+",	3./5,	1.392,	5.344,	45.450,	0,		0,		0,		0],
"Nd3+",	["4f3",	"",	"",	"Nd3+",	8./11,	1.114,	2.910,	15.03,	109.763,	9053.08,	6545,		4433.3],
"Pm3+",	["4f4",	"",	"",	"Pm3+",	3./5,	1.0353,	2.5390,	12.546,	127.084,	9472.39,	6806.73,	4674.2],
"Sm2+",	["4f6",	"",	"",	"Sm2+",	0.,	1.197,	3.861,	28.560,	0,		0,		0,		0],
"Sm3+",	["4f5",	"",	"",	"Sm3+",	2./7,	0.9743,	2.260,	10.55,	145.805,	9894.56,	7088.8,		4990.36],
"Eu2+",	["4f7",	"",	"",	"Eu2+",	2.,	1.098,	3.368,	23.580,	0,		0,		0,		0],
"Eu3+",	["4f6",	"",	"",	"Eu3+",	0.,	0.9175,	2.020,	9.039,	165.891,	10306.2,	7348.29,	5276.77],
"Gd2+",	["4f8",	"",	"",	"Gd2+",	3./2,	1.028,	2.975,	19.850,	0,		0,		0,		0],
"Gd3+",	["4f7",	"",	"",	"Gd3+",	2.,	0.8671,	1.820,	7.831,	186.968,	10621.6,	7541.34,	5551.52],
"Tb2+",	["4f9",	"",	"",	"Tb2+",	4./3,	0.968,	2.655,	16.980,	0,		0,		0,		0],
"Tb3+",	["4f8",	"",	"",	"Tb3+",	3./2,	0.8220,	1.651,	6.852,	211.641,	11034,		7800.96,	5858.5],
"Tb4+",	["4f7",	"",	"",	"Gd3+",	2.,	0,	0,	0,	0,		0,		0,		0],
"Dy2+",	["4f10","",	"",	"Dy2+",	5./4,	0.913,	2.391,	14.730,	0,		0,		0,		0],
"Dy3+",	["4f9",	"",	"",	"Dy3+",	4./3,	0.7814,	1.505,	6.048,	237.182,	11394.5,	7981.11,	6123.08],
"Ho2+",	["4f11","",	"",	"Ho2+",	6./5,	0.866,	2.169,	12.920,	0,		0,		0,		0],
"Ho3+",	["4f10","",	"",	"Ho3+",	5./4,	0.7446,	1.379,	5.379,	265.946,	11724.4,	8232.18,	6449.91],
"Er2+",	["4f12","",	"",	"Er2+",	7./6,	0.824,	1.979,	11.450,	0,		0,		0,		0],
"Er3+",	["4f11","",	"",	"Er3+",	6./5,	0.7111,	1.270,	4.816,	294.586,	12086.4,	8419.02,	6696.39],
"Tm2+",	["4f13","",	"",	"Tm2+",	8./7,	0.785,	1.819,	10.240,	0,		0,		0,		0],
"Tm3+",	["4f12","",	"",	"Tm3+",	7./6,	0.6804,	1.174,	4.340,	326.822,	12415,		8630.91,	6940.01],
"Yb3+",	["4f13","",	"",	"Yb3+",	8./7,	0.6522,	1.089,	3.932,	363.026,	0,		0,		0],
"Os2+",	["5d6",	2.,	0,	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Os3+",	["5d5",	2.5,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Os4+",	["5d4",	2.,	1,	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Os6+",	["5d2",	1,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Ir2+",	["5d7",	1.5,	0.5,	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Ir3+",	["5d6",	2.,	0,	"J=0.5",2.,	0,	0,	0,	0,		0,		0,		0],
"Ir4+",	["5d5",	2.5,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Ir6+",	["5d3",	1.5,	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"U2+",	["5f4",	"",	"",	"U2+",	3./5,	3.257,	26.82,	462.85,	0,		0,		0,		0],
"U3+",	["5f3",	"",	"",	"U3+",	8./11,	2.346,	10.906,	90.544,	201.598,	4911.14,	4086.52,	2862.05],
"U4+",	["5f2",	"",	"",	"U4+",	4./5,	2.042,	7.632,	47.774,	215.732,	5487.29,	4995.69,	3879.09],
"U5+",	["5f1",	"",	"",	"",	6./7,	0,	0,	0,	0,		0,		0,		0],
"Np3+",	["5f4",	"",	"",	"Np3+",	3./5,	2.297,	11.00,	98.63,	240.157,	5626.65,	4617.42,	3179.45],
"Np4+",	["5f3",	"",	"",	"Np4+",	8./11,	1.884,	6.504,	37.80,	258.879,	5944.92,	5139.76,	3289.05],
"Np5+",	["5f2",	"",	"",	"",	4./5,	0,	0,	0,	0,		0,		0,		0],
"Np6+",	["5f1",	"",	"",	"",	6./7,	0,	0,	0,	0,		0,		0,		0],
"Pu3+",	["5f5",	"",	"",	"Pu3+",	2./7,	2.1025,	9.1775,	75.30,	277.973,	6035.43,	4876.67,	3427.79],
"Pu4+",	["5f4",	"",	"",	"Pu4+",	3./5,	1.838,	6.401,	38.77,	293.347,	6124.07,	4896.76,	3804.33],
"Pu5+",	["5f3",	"",	"",	"",	8./11,	0,	0,	0,	0,		0,		0,		0],
"Pu6+",	["5f2",	"",	"",	"",	4./5,	0,	0,	0,	0,		0,		0,		0],
"Am2+",	["5f7",	"",	"",	"",	2.,	0,	0,	0,	0,		0,		0,		0],
"Am3+",	["5f6",	"",	"",	"",	0.,	0,	0,	0,	317.895,	6434.78,	5157.74,	3645.14],
"Am4+",	["5f5",	"",	"",	"",	2./7,	0,	0,	0,	349.759,	7869.59,	5415.47,	4146.46],
"Am5+",	["5f4",	"",	"",	"",	3./5,	0,	0,	0,	0,		0,		0,		0],
"Am6+",	["5f3",	"",	"",	"",	8./11,	0,	0,	0,	0,		0,		0,		0],
"Cm3+",	["5f7",	"",	"",	"",	2.,	0,	0,	0,	358.19,		6825.95,	5447.62,	4076.1],
);

# Magnetic formfactors tables
%magff = (
"V2+",	["3d3", "FFj0A= 0.4085 FFj0a=23.8526 FFj0B= 0.6091 FFj0b= 8.2456 FFj0C=-0.1676 FFj0c= 0.0415 FFj0D= 0.1496\nFFj2A= 3.4386 FFj2a=16.5303 FFj2B= 1.9638 FFj2b= 6.1415 FFj2C= 0.2997 FFj2c= 2.2669 FFj2D= 0.0009\nFFj4A=-1.1729 FFj4a=14.9732 FFj4B= 0.9092 FFj4b= 7.6131 FFj4C= 0.4105 FFj4c= 2.0391 FFj4D= 0.0067\n"],
"V3+",	["3d2", "FFj0A= 0.3598 FFj0a=19.3364 FFj0B= 0.6632 FFj0b= 7.6172 FFj0C=-0.3064 FFj0c= 0.0296 FFj0D= 0.2835\nFFj2A= 2.3005 FFj2a=14.6821 FFj2B= 2.0364 FFj2b= 6.1304 FFj2C= 0.4099 FFj2c= 2.3815 FFj2D= 0.0014\nFFj4A=-0.9417 FFj4a=14.2045 FFj4B= 0.5284 FFj4b= 6.6071 FFj4C= 0.4411 FFj4c= 1.9672 FFj4D= 0.0076\n"],
"V4+",	["3d1", "FFj0A= 0.3106 FFj0a=16.8160 FFj0B= 0.7198 FFj0b= 7.0487 FFj0C=-0.0521 FFj0c= 0.3020 FFj0D= 0.0221\nFFj2A= 1.8377 FFj2a=12.2668 FFj2B= 1.8247 FFj2b= 5.4578 FFj2C= 0.3979 FFj2c= 2.2483 FFj2D= 0.0012\nFFj4A=-0.7654 FFj4a=13.0970 FFj4B= 0.3071 FFj4b= 5.6739 FFj4C= 0.4476 FFj4c= 1.8707 FFj4D= 0.0081\n"],
"Cr2+",	["3d4", "FFj0A= 1.2024 FFj0a=-0.0055 FFj0B= 0.4158 FFj0b=20.5475 FFj0C= 0.6032 FFj0c= 6.9560 FFj0D=-1.2218\nFFj2A= 2.6422 FFj2a=16.0598 FFj2B= 1.9198 FFj2b= 6.2531 FFj2C= 0.4446 FFj2c= 2.3715 FFj2D= 0.0020\nFFj4A=-0.8930 FFj4a=15.6641 FFj4B= 0.5590 FFj4b= 7.0333 FFj4C= 0.4093 FFj4c= 1.9237 FFj4D= 0.0081\n"],
"Cr3+",	["3d3", "FFj0A=-0.3094 FFj0a= 0.0274 FFj0B= 0.3680 FFj0b=17.0355 FFj0C= 0.6559 FFj0c= 6.5236 FFj0D= 0.2856\nFFj2A= 1.6262 FFj2a=15.0656 FFj2B= 2.0618 FFj2b= 6.2842 FFj2C= 0.5281 FFj2c= 2.3680 FFj2D= 0.0023\nFFj4A=-0.7327 FFj4a=14.0727 FFj4B= 0.3268 FFj4b= 5.6741 FFj4C= 0.4114 FFj4c= 1.8101 FFj4D= 0.0085\n"],
"Cr4+",	["3d2", "FFj0A=-0.2320 FFj0a= 0.0433 FFj0B= 0.3101 FFj0b=14.9518 FFj0C= 0.7182 FFj0c= 6.1726 FFj0D= 0.2042\nFFj2A= 1.0293 FFj2a=13.9498 FFj2B= 1.9933 FFj2b= 6.0593 FFj2C= 0.5974 FFj2c= 2.3457 FFj2D= 0.0027\nFFj4A=-0.6748 FFj4a=12.9462 FFj4B= 0.1805 FFj4b= 6.7527 FFj4C= 0.4526 FFj4c= 1.7999 FFj4D= 0.0098\n"],
"Mn2+",	["3d5", "FFj0A= 0.4220 FFj0a=17.6840 FFj0B= 0.5948 FFj0b= 6.0050 FFj0C= 0.0043 FFj0c=-0.6090 FFj0D=-0.0219\nFFj2A= 2.0515 FFj2a=15.5561 FFj2B= 1.8841 FFj2b= 6.0625 FFj2C= 0.4787 FFj2c= 2.2323 FFj2D= 0.0027\nFFj4A=-0.7416 FFj4a=15.2555 FFj4B= 0.3831 FFj4b= 6.4693 FFj4C= 0.3935 FFj4c= 1.7997 FFj4D= 0.0093\n"],
"Mn3+",	["3d4", "FFj0A= 0.4198 FFj0a=14.2829 FFj0B= 0.6054 FFj0b= 5.4689 FFj0C= 0.9241 FFj0c=-0.0088 FFj0D=-0.9498\nFFj2A= 1.2427 FFj2a=14.9966 FFj2B= 1.9567 FFj2b= 6.1181 FFj2C= 0.5732 FFj2c= 2.2577 FFj2D= 0.0031\nFFj4A=-0.6603 FFj4a=13.6066 FFj4B= 0.2322 FFj4b= 6.2175 FFj4C= 0.4104 FFj4c= 1.7404 FFj4D= 0.0101\n"],
"Mn4+",	["3d3", "FFj0A= 0.3760 FFj0a=12.5661 FFj0B= 0.6602 FFj0b= 5.1329 FFj0C=-0.0372 FFj0c= 0.5630 FFj0D= 0.0011\nFFj2A= 0.7879 FFj2a=13.8857 FFj2B= 1.8717 FFj2b= 5.7433 FFj2C= 0.5981 FFj2c= 2.1818 FFj2D= 0.0034\nFFj4A=-0.5127 FFj4a=13.4613 FFj4B= 0.0313 FFj4b= 7.7631 FFj4C= 0.4282 FFj4c= 1.7006 FFj4D= 0.0113\n"],
"Mn5+",   ["3d2", "FFj0A= 0.2924 FFj0a=11.6655 FFj0B= 0.7405 FFj0b= 5.0741 FFj0C=-1.7883 FFj0c= 0.0059 FFj0D= 1.7557\nFFj2A=-0.2394 FFj2a=10.7309 FFj2B=-0.1190 FFj2b= 6.5989 FFj2C= 0.3505 FFj2c= 1.4912 FFj2D= 0.0078\nFFj4A= 1.6706 FFj4a= 6.6566 FFj4B=-1.8204 FFj4b= 6.1942 FFj4C= 0.1925 FFj4c= 0.3249 FFj4D=-0.0433\n"],
"Mn6+",	["3d1", ""],
"Fe2+",	["3d6", "FFj0A= 0.0263 FFj0a=34.9597 FFj0B= 0.3668 FFj0b=15.9435 FFj0C= 0.6188 FFj0c= 5.5935 FFj0D=-0.0119\nFFj2A= 1.6490 FFj2a=16.5593 FFj2B= 1.9064 FFj2b= 6.1325 FFj2C= 0.5206 FFj2c= 2.1370 FFj2D= 0.0035\nFFj4A=-0.5401 FFj4a=17.2268 FFj4B= 0.2865 FFj4b= 3.7422 FFj4C= 0.2658 FFj4c= 1.4238 FFj4D= 0.0076\n"],
"Fe3+",	["3d5", "FFj0A= 0.3972 FFj0a=13.2442 FFj0B= 0.6295 FFj0b= 4.9034 FFj0C=-0.0314 FFj0c= 0.3496 FFj0D= 0.0044\nFFj2A= 1.3602 FFj2a=11.9976 FFj2B= 1.5188 FFj2b= 5.0025 FFj2C= 0.4705 FFj2c= 1.9914 FFj2D= 0.0038\nFFj4A=-0.5507 FFj4a=11.4929 FFj4B= 0.2153 FFj4b= 4.9063 FFj4C= 0.3468 FFj4c= 1.5230 FFj4D= 0.0095\n"],
"Co2+",	["3d7", "FFj0A= 0.4332 FFj0a=14.3553 FFj0B= 0.5857 FFj0b= 4.6077 FFj0C=-0.0382 FFj0c= 0.1338 FFj0D= 0.0179\nFFj2A= 1.9049 FFj2a=11.6444 FFj2B= 1.3159 FFj2b= 4.3574 FFj2C= 0.3146 FFj2c= 1.6453 FFj2D= 0.0017\nFFj4A=-0.4759 FFj4a=14.0462 FFj4B= 0.2747 FFj4b= 3.7306 FFj4C= 0.2458 FFj4c= 1.2504 FFj4D= 0.0057\n"],
"Co3+",	["3d6", "FFj0A= 0.3902 FFj0a=12.5078 FFj0B= 0.6324 FFj0b= 4.4574 FFj0C=-0.1500 FFj0c= 0.0343 FFj0D= 0.1272\nFFj2A= 1.7058 FFj2a= 8.8595 FFj2B= 1.1409 FFj2b= 3.3086 FFj2C= 0.1474 FFj2c= 1.0899 FFj2D=-0.0025\nFFj4A=-0.4466 FFj4a=13.3912 FFj4B= 0.1419 FFj4b= 3.0110 FFj4C= 0.2773 FFj4c= 1.3351 FFj4D= 0.0093\n"],
"Co4+",   ["3d5", "FFj0A= 0.3515 FFj0a=10.7785 FFj0B= 0.6778 FFj0b= 4.2343 FFj0C=-0.0389 FFj0c= 0.2409 FFj0D= 0.0098\nFFj2A= 1.3110 FFj2a= 8.0252 FFj2B= 1.1551 FFj2b= 3.1792 FFj2C= 0.1608 FFj2c= 1.1301 FFj2D=-0.0011\nFFj4A=-0.4091 FFj4a=13.1937 FFj4B=-0.0194 FFj4b= 3.4169 FFj4C= 0.3534 FFj4c= 1.4214 FFj4D= 0.0112\n"],
"Ni2+",	["3d8", "FFj0A= 0.0163 FFj0a=35.8826 FFj0B= 0.3916 FFj0b=13.2233 FFj0C= 0.6052 FFj0c= 4.3388 FFj0D=-0.0133\nFFj2A= 1.7080 FFj2a=11.0160 FFj2B= 1.2147 FFj2b= 4.1031 FFj2C= 0.3150 FFj2c= 1.5334 FFj2D= 0.0018\nFFj4A=-0.3803 FFj4a=10.4033 FFj4B= 0.2838 FFj4b= 3.3780 FFj4C= 0.2108 FFj4c= 1.1036 FFj4D= 0.0050\n"],
"Ni3+",	["3d7", "FFj0A= 0.0012 FFj0a=34.9998 FFj0B= 0.3468 FFj0b=11.9874 FFj0C= 0.6667 FFj0c= 4.2518 FFj0D=-0.0148\nFFj2A= 1.4683 FFj2a= 8.6713 FFj2B= 1.1068 FFj2b= 3.2574 FFj2C= 0.1794 FFj2c= 1.1058 FFj2D=-0.0023\nFFj4A=-0.4014 FFj4a= 9.0462 FFj4B= 0.2314 FFj4b= 3.0753 FFj4C= 0.2192 FFj4c= 1.0838 FFj4D= 0.0060\n"],
"Cu2+",	["3d9", "FFj0A= 0.0232 FFj0a=34.9686 FFj0B= 0.4023 FFj0b=11.5640 FFj0C= 0.5882 FFj0c= 3.8428 FFj0D=-0.0137\nFFj2A= 1.5189 FFj2a=10.4779 FFj2B= 1.1512 FFj2b= 3.8132 FFj2C= 0.2918 FFj2c= 1.3979 FFj2D= 0.0017\nFFj4A=-0.3914 FFj4a=14.7400 FFj4B= 0.1275 FFj4b= 3.3840 FFj4C= 0.2548 FFj4c= 1.2552 FFj4D= 0.0103\n"],
"Ru+",    ["4d7", "FFj0A= 0.4410 FFj0a=33.3086 FFj0B= 1.4775 FFj0b= 9.5531 FFj0C=-0.9361 FFj0c= 6.7220 FFj0D= 0.0176\nFFj2A= 5.2826 FFj2a=23.6832 FFj2B= 3.5813 FFj2b= 8.1521 FFj2C=-0.0257 FFj2c= 0.4255 FFj2D= 0.0131\nFFj4A=-1.6278 FFj4a=18.5063 FFj4B= 1.1828 FFj4b=10.1886 FFj4C= 0.8138 FFj4c= 3.4180 FFj4D=-0.0009\n"],
"Ru2+",	["4d6", ""],
"Ru3+",	["4d5", ""],
"Ru4+",	["4d4", ""],
"Ru6+",	["4d2", ""],
"Rh+",	["4d8", "FFj0A= 0.3342 FFj0a=29.7564 FFj0B= 1.2209 FFj0b= 9.4384 FFj0C=-0.5755 FFj0c= 5.3320 FFj0D= 0.0210\nFFj2A= 4.0260 FFj2a=18.9497 FFj2B= 3.1663 FFj2b= 6.9998 FFj2C=-0.0296 FFj2c= 0.4862 FFj2D= 0.0127\nFFj4A=-1.4673 FFj4a=17.9572 FFj4B= 0.7381 FFj4b= 9.9444 FFj4C= 0.8485 FFj4c= 3.1263 FFj4D=-0.0012\n"],
"Rh2+",	["4d7", ""],
"Rh3+",	["4d6", ""],
"Rh4+",	["4d5", ""],
"Ce2+",   ["4f2", "FFj0A= 0.2953 FFj0a=17.6846 FFj0B= 0.2923 FFj0b= 6.7329 FFj0C= 0.4313 FFj0c= 5.3827 FFj0D=-0.0194\nFFj2A= 0.9809 FFj2a=18.0630 FFj2B= 1.8413 FFj2b= 7.7688 FFj2C= 0.9905 FFj2c= 2.8452 FFj2D= 0.0120\nFFj4A=-0.6468 FFj4a=10.5331 FFj4B= 0.4052 FFj4b= 5.6243 FFj4C= 0.3412 FFj4c= 1.5346 FFj4D= 0.0080\nFFj6A=-0.1212 FFj6a= 7.9940 FFj6B=-0.0639 FFj6b= 4.0244 FFj6C= 0.1519 FFj6c= 1.0957 FFj6D= 0.0078\n"],
"Ce3+",	["4f1", "FFj0A=+0.2291 FFj0a=+18.18 FFj0B=+0.7897 FFj0b=+5.807 FFj0C=-0.0191 FFj0c=0.0000 FFj0D=0.0000\nFFj2A=+2.1284 FFj2a=+8.9174 FFj2B=+1.1229 FFj2b=+2.8371 FFj2C=+0.01108 FFj2c=0.000 FFj2D=+0.0000\nFFj4A=+0.4221 FFj4a=+1.7572 FFj4B=-0.4087 FFj4b=+14.604 FFj4C=+0.01465 FFj4c=+0.000 FFj4D=+0.0000\nFFj6A=+0.13076 FFj6a=+0.8650 FFj6B=-0.15173 FFj6b=+5.6704 FFj6C=+0.00281 FFj6c=+0.000 FFj6D=+0.0000\n"],
"Pr3+",	["4f2", "FFj0A=+0.0504 FFj0a=+24.9989 FFj0B=+0.2572 FFj0b=+12.0377 FFj0C=+0.7142 FFj0c=+5.0039 FFj0D=-0.0219\nFFj2A=+0.8734 FFj2a=+18.9876 FFj2B=+1.5594 FFj2b=+6.0872 FFj2C=+0.8142 FFj2c=+2.4150 FFj2D=+0.0111\nFFj4A=-0.3970 FFj4a=+10.9919 FFj4B=+0.0818 FFj4b=+5.9897 FFj4C=+0.3656 FFj4c=+1.5021 FFj4D=+0.0110\nFFj6A=-0.0224 FFj6a=+7.9931 FFj6B=-0.1202 FFj6b=+3.9406 FFj6C=+0.1299 FFj6c=+0.8938 FFj6D=+0.0051\n"],
"Pr4+",	["4f1", ""],
"Nd2+",	["4f4", "FFj0A= 0.1645 FFj0a=25.0453 FFj0B= 0.2522 FFj0b=11.9782 FFj0C= 0.6012 FFj0c= 4.9461 FFj0D=-0.0180\nFFj2A= 1.4530 FFj2a=18.3398 FFj2B= 1.6196 FFj2b= 7.2854 FFj2C= 0.8752 FFj2c= 2.6224 FFj2D= 0.0126\nFFj4A=-0.5744 FFj4a=10.9304 FFj4B= 0.4210 FFj4b= 6.1052 FFj4C= 0.3124 FFj4c= 1.4654 FFj4D= 0.0081\nFFj6a= 8.0086 FFj6B= 0.0272 FFj6b= 4.0284 FFj6C= 0.1104 FFj6c= 1.0682 FFj6D= 0.0139\n"],
"Nd3+",	["4f3", "FFj0A=0.0540 FFj0a=25.0293 FFj0B=0.3101 FFj0b=12.1020 FFj0C=0.6575 FFj0c=4.7223 FFj0D=-0.0216\nFFj2A=0.6751 FFj2a=18.3421 FFj2B=1.6272 FFj2b=7.2600 FFj2C=0.9644 FFj2c=2.6016 FFj2D=0.0150\nFFj4A=-0.4053 FFj4a=+14.0141 FFj4B=+0.0329 FFj4b=+7.0046 FFj4C=+0.3759 FFj4c=+1.7074 FFj4D=+0.0209\nFFj6A=-0.0416 FFj6a=+8.0136 FFj6B=-0.1261 FFj6b=+4.0399 FFj6C=+0.1400 FFj6c=+1.0873 FFj6D=+0.0102\n"],
"Pm3+",	["4f4", ""],
"Sm2+",	["4f6", "FFj0A= 0.0909 FFj0a=25.2032 FFj0B= 0.3037 FFj0b=11.8562 FFj0C= 0.6250 FFj0c= 4.2366 FFj0D=-0.0200\nFFj2A= 1.0360 FFj2a=18.4249 FFj2B= 1.4769 FFj2b= 7.0321 FFj2C= 0.8810 FFj2c= 2.4367 FFj2D= 0.0152\nFFj4A=-0.4150 FFj4a=14.0570 FFj4B= 0.1368 FFj4b= 7.0317 FFj4C= 0.3272 FFj4c= 1.5825 FFj4D= 0.0192\nFFj6A=-0.1428 FFj6a= 6.0407 FFj6B= 0.0723 FFj6b= 2.0329 FFj6C= 0.0550 FFj6c= 0.5134 FFj6D= 0.0081\n"],
"Sm3+",	["4f5", "FFj0A= 0.0288 FFj0a=25.2068 FFj0B= 0.2973 FFj0b=11.8311 FFj0C= 0.6954 FFj0c= 4.2117 FFj0D=-0.0213\nFFj2A= 0.4707 FFj2a=18.4301 FFj2B= 1.4261 FFj2b= 7.0336 FFj2C= 0.9574 FFj2c= 2.4387 FFj2D= 0.0182\nFFj4A=-0.4288 FFj4a=10.0525 FFj4B= 0.1782 FFj4b= 5.0191 FFj4C= 0.2833 FFj4c= 1.2364 FFj4D= 0.0088\nFFj6A=-0.0944 FFj6a= 6.0299 FFj6B=-0.0498 FFj6b= 2.0743 FFj6C= 0.1372 FFj6c= 0.6451 FFj6D=-0.0132\n"],
"Eu2+",	["4f7", "FFj0A= 0.0755 FFj0a=25.2960 FFj0B= 0.3001 FFj0b=11.5993 FFj0C= 0.6438 FFj0c= 4.0252 FFj0D=-0.0196\nFFj2A= 0.8970 FFj2a=18.4429 FFj2B= 1.3769 FFj2b= 7.0054 FFj2C= 0.9060 FFj2c= 2.4213 FFj2D= 0.0190\nFFj4A=-0.4145 FFj4a=10.1930 FFj4B= 0.2447 FFj4b= 5.1644 FFj4C= 0.2661 FFj4c= 1.2054 FFj4D= 0.0065\nFFj6A=-0.1252 FFj6a= 6.0485 FFj6B= 0.0507 FFj6b= 2.0852 FFj6C= 0.0572 FFj6c= 0.6460 FFj6D= 0.0132\n"],
"Eu3+",	["4f6", "FFj0A= 0.0204 FFj0a=25.3078 FFj0B= 0.3010 FFj0b=11.4744 FFj0C= 0.7005 FFj0c= 3.9420 FFj0D=-0.0220\nFFj2A= 0.3985 FFj2a=18.4514 FFj2B= 1.3307 FFj2b= 6.9556 FFj2C= 0.9603 FFj2c= 2.3780 FFj2D= 0.0197\nFFj4A=-0.4095 FFj4a=10.2113 FFj4B= 0.1485 FFj4b= 5.1755 FFj4C= 0.2720 FFj4c= 1.2374 FFj4D= 0.0131\nFFj6A=-0.0817 FFj6a= 6.0389 FFj6B=-0.0596 FFj6b= 2.1198 FFj6C= 0.1243 FFj6c= 0.7639 FFj6D=-0.0001\n"],
"Gd2+",	["4f8", "FFj0A= 0.0636 FFj0a=25.3823 FFj0B= 0.3033 FFj0b=11.2125 FFj0C= 0.6528 FFj0c= 3.7877 FFj0D=-0.0199\nFFj2A= 0.7756 FFj2a=18.4695 FFj2B= 1.3124 FFj2b= 6.8990 FFj2C= 0.8956 FFj2c= 2.3383 FFj2D= 0.0199\nFFj4A=-0.3824 FFj4a=10.3436 FFj4B= 0.1955 FFj4b= 5.3057 FFj4C= 0.2622 FFj4c= 1.2032 FFj4D= 0.0097\nFFj6A=-0.1351 FFj6a= 5.0298 FFj6B= 0.0828 FFj6b= 2.0248 FFj6C= 0.0315 FFj6c= 0.5034 FFj6D= 0.0187\n"],
"Gd3+",	["4f7", "FFj0A=0.0186 FFj0a=25.3867 FFj0B=0.2895 FFj0b=11.1421 FFj0C=0.7135 FFj0c=3.7520 FFj0D=-0.0217\nFFj2A=0.3347 FFj2a=18.4758 FFj2B=1.2465 FFj2b=6.8767 FFj2C=0.9537 FFj2c=2.3184 FFj2D=0.0217\nFFj4A=-0.3621 FFj4a=10.3531 FFj4B=0.1016 FFj4b=5.3104 FFj4C=0.2649 FFj4c=1.2185 FFj4D=0.0147\nFFj6A=-0.0662 FFj6a=6.0308 FFj6B=-0.0850 FFj6b=2.1542 FFj6C=0.1323 FFj6c=0.8910 FFj6D=0.0048\n"],
"Tb2+",	["4f9", "FFj0A= 0.0547 FFj0a=25.5086 FFj0B= 0.3171 FFj0b=10.5911 FFj0C= 0.6490 FFj0c= 3.5171 FFj0D=-0.0212\nFFj2A= 0.6688 FFj2a=18.4909 FFj2B= 1.2487 FFj2b= 6.8219 FFj2C= 0.8888 FFj2c= 2.2751 FFj2D= 0.0215\nFFj4A=-0.3443 FFj4a=10.4686 FFj4B= 0.1481 FFj4b= 5.4156 FFj4C= 0.2575 FFj4c= 1.1824 FFj4D= 0.0104\nFFj6A=-0.0758 FFj6a= 6.0319 FFj6B=-0.0540 FFj6b= 2.1583 FFj6C= 0.1199 FFj6c= 0.8895 FFj6D= 0.0051\n"],
"Tb3+",	["4f8", "FFj0A= 0.0177 FFj0a=25.5095 FFj0B= 0.2921 FFj0b=10.5769 FFj0C= 0.7133 FFj0c= 3.5122 FFj0D=-0.0231\nFFj2A= 0.2892 FFj2a=18.4973 FFj2B= 1.1678 FFj2b= 6.7972 FFj2C= 0.9437 FFj2c= 2.2573 FFj2D= 0.0232\nFFj4A=-0.3228 FFj4a=10.4763 FFj4B= 0.0638 FFj4b= 5.4189 FFj4C= 0.2566 FFj4c= 1.1962 FFj4D= 0.0159\nFFj6A=-0.0559 FFj6a= 6.0311 FFj6B=-0.1020 FFj6b= 2.2365 FFj6C= 0.1264 FFj6c= 1.1066 FFj6D= 0.0167\n"],
"Tb4+",	["4f7", ""],
"Dy2+",	["4f10", "FFj0A= 0.1308 FFj0a=18.3155 FFj0B= 0.3118 FFj0b= 7.6645 FFj0C= 0.5795 FFj0c= 3.1469 FFj0D=-0.0226\nFFj2A= 0.5917 FFj2a=18.5114 FFj2B= 1.1828 FFj2b= 6.7465 FFj2C= 0.8801 FFj2c= 2.2141 FFj2D= 0.0229\nFFj4A=-0.3206 FFj4a=12.0714 FFj4B= 0.0904 FFj4b= 8.0264 FFj4C= 0.2616 FFj4c= 1.2296 FFj4D= 0.0143\nFFj6A=-0.0568 FFj6a= 6.0324 FFj6B=-0.1003 FFj6b= 2.2396 FFj6C= 0.1401 FFj6c= 1.1062 FFj6D= 0.0109\n"],
"Dy3+",	["4f9", "FFj0A= 0.1157 FFj0a=15.0732 FFj0B= 0.3270 FFj0b= 6.7991 FFj0C= 0.5821 FFj0c= 3.0202 FFj0D=-0.0249\nFFj2A= 0.2523 FFj2a=18.5172 FFj2B= 1.0914 FFj2b= 6.7362 FFj2C= 0.9345 FFj2c= 2.2082 FFj2D= 0.0250\nFFj4A=-0.2829 FFj4a= 9.5247 FFj4B= 0.0565 FFj4b= 4.4292 FFj4C= 0.2437 FFj4c= 1.0665 FFj4D= 0.0092\nFFj6A=-0.0423 FFj6a= 6.0376 FFj6B=-0.1248 FFj6b= 2.2437 FFj6C= 0.1359 FFj6c= 1.2002 FFj6D= 0.0188\n"],
"Ho2+",	["4f11", "FFj0A= 0.0995 FFj0a=18.1761 FFj0B= 0.3305 FFj0b= 7.8556 FFj0C= 0.5921 FFj0c= 2.9799 FFj0D=-0.0230\nFFj2A= 0.5094 FFj2a=18.5155 FFj2B= 1.1234 FFj2b= 6.7060 FFj2C= 0.8727 FFj2c= 2.1589 FFj2D= 0.0242\nFFj4A=-0.2976 FFj4a= 9.7190 FFj4B= 0.1224 FFj4b= 4.6345 FFj4C= 0.2279 FFj4c= 1.0052 FFj4D= 0.0063\nFFj6A=-0.0725 FFj6a= 6.0453 FFj6B=-0.0318 FFj6b= 2.2428 FFj6C= 0.0738 FFj6c= 1.2018 FFj6D= 0.0252\n"],
"Ho3+",	["4f10", "FFj0A= 0.0566 FFj0a=18.3176 FFj0B= 0.3365 FFj0b= 7.6880 FFj0C= 0.6317 FFj0c= 2.9427 FFj0D=-0.0248\nFFj2A= 0.2188 FFj2a=18.5157 FFj2B= 1.0240 FFj2b= 6.7070 FFj2C= 0.9251 FFj2c= 2.1614 FFj2D= 0.0268\nFFj4A=-0.2717 FFj4a= 9.7313 FFj4B= 0.0474 FFj4b= 4.6378 FFj4C= 0.2292 FFj4c= 1.0473 FFj4D= 0.0124\nFFj6A=-0.0289 FFj6a= 6.0504 FFj6B=-0.1545 FFj6b= 2.2305 FFj6C= 0.1550 FFj6c= 1.2605 FFj6D= 0.0177\n"],
"Er2+",	["4f12", "FFj0A= 0.1122 FFj0a=18.1223 FFj0B= 0.3462 FFj0b= 6.9106 FFj0C= 0.5649 FFj0c= 2.7614 FFj0D=-0.0235\nFFj2A= 0.4693 FFj2a=18.5278 FFj2B= 1.0545 FFj2b= 6.6493 FFj2C= 0.8679 FFj2c= 2.1201 FFj2D= 0.0261\nFFj4A=-0.2975 FFj4a= 9.8294 FFj4B= 0.1189 FFj4b= 4.7406 FFj4C= 0.2116 FFj4c= 1.0039 FFj4D= 0.0117\nFFj6A=-0.0648 FFj6a= 6.0559 FFj6B=-0.0515 FFj6b= 2.2303 FFj6C= 0.0825 FFj6c= 1.2638 FFj6D= 0.0250\n"],
"Er3+",	["4f11", "FFj0A= 0.0586 FFj0a=17.9802 FFj0B= 0.3540 FFj0b= 7.0964 FFj0C= 0.6126 FFj0c= 2.7482 FFj0D=-0.0251\nFFj2A= 0.1710 FFj2a=18.5337 FFj2B= 0.9879 FFj2b= 6.6246 FFj2C= 0.9044 FFj2c= 2.1004 FFj2D= 0.0278\nFFj4A=-0.2568 FFj4a= 9.8339 FFj4B= 0.0356 FFj4b= 4.7415 FFj4C= 0.2172 FFj4c= 1.0281 FFj4D= 0.0148\nFFj6A=-0.0110 FFj6a= 6.0609 FFj6B=-0.1954 FFj6b= 2.2242 FFj6C= 0.1818 FFj6c= 1.2958 FFj6D= 0.0149\n"],
"Tm2+",	["4f13", "FFj0A= 0.0983 FFj0a=18.3236 FFj0B= 0.3380 FFj0b= 6.9178 FFj0C= 0.5875 FFj0c= 2.6622 FFj0D=-0.0241\nFFj2A= 0.4198 FFj2a=18.5417 FFj2B= 0.9959 FFj2b= 6.6002 FFj2C= 0.8593 FFj2c= 2.0818 FFj2D= 0.0284\nFFj4A=-0.2677 FFj4a= 9.8883 FFj4B= 0.0925 FFj4b= 4.7838 FFj4C= 0.2056 FFj4c= 0.9896 FFj4D= 0.0124\nFFj6A=-0.0842 FFj6a= 4.0699 FFj6B= 0.0807 FFj6b= 0.8492 FFj6C=-0.2087 FFj6c= 0.0386 FFj6D= 0.2095\n"],
"Tm3+",	["4f12", "FFj0A= 0.0581 FFj0a=15.0922 FFj0B= 0.2787 FFj0b= 7.8015 FFj0C= 0.6854 FFj0c= 2.7931 FFj0D=-0.0224\nFFj2A= 0.1760 FFj2a=18.5417 FFj2B= 0.9105 FFj2b= 6.5787 FFj2C= 0.8970 FFj2c= 2.0622 FFj2D= 0.0294\nFFj4A=-0.2292 FFj4a= 9.8948 FFj4B= 0.0124 FFj4b= 4.7850 FFj4C= 0.2108 FFj4c= 1.0071 FFj4D= 0.0151\nFFj6A=-0.0727 FFj6a= 4.0730 FFj6B= 0.0243 FFj6b= 0.6888 FFj6C= 3.9459 FFj6c= 0.0023 FFj6D=-3.9076\n"],
"Yb3+",	["4f13", "FFj0A= 0.0416 FFj0a=16.0949 FFj0B= 0.2849 FFj0b= 7.8341 FFj0C= 0.6961 FFj0c= 2.6725 FFj0D=-0.0229\nFFj2A= 0.1570 FFj2a=18.5553 FFj2B= 0.8484 FFj2b= 6.5403 FFj2C= 0.8880 FFj2c= 2.0367 FFj2D= 0.0318\nFFj4A=-0.2121 FFj4a= 8.1967 FFj4B= 0.0325 FFj4b= 3.1533 FFj4C= 0.1975 FFj4c= 0.8842 FFj4D= 0.0093\nFFj6A=-0.0345 FFj6a= 5.0073 FFj6B=-0.0677 FFj6b= 2.0198 FFj6C= 0.0985 FFj6c= 0.5485 FFj6D=-0.0076\n"],
"Os2+",	["5d6", ""],
"Os3+",	["5d5", ""],
"Os4+",	["5d4", ""],
"Os6+",	["5d2", ""],
"Ir2+",	["5d7", ""],
"Ir3+",	["5d6", ""],
"Ir4+",	["5d5", ""],
"Ir6+",	["5d3", ""],
"U2+",	["5f4", ""],
"U3+",	["5f3", "FFj0A= 0.5058 FFj0a=23.2882 FFj0B= 1.3464 FFj0b= 7.0028 FFj0C=-0.8724 FFj0c= 4.8683 FFj0D= 0.0192\nFFj2A= 4.1582 FFj2a=16.5336 FFj2B= 2.4675 FFj2b= 5.9516 FFj2C=-0.0252 FFj2c= 0.7646 FFj2D= 0.0057\nFFj4A=-0.9859 FFj4a=16.6010 FFj4B= 0.6116 FFj4b= 6.5147 FFj4C= 0.6020 FFj4c= 2.5970 FFj4D=-0.0010\nFFj6A=-0.3797 FFj6a= 9.9525 FFj6B= 0.0459 FFj6b= 5.0379 FFj6C= 0.2748 FFj6c= 1.6072 FFj6D= 0.0016\n"],
"U4+",	["5f2", "FFj0A= 0.3291 FFj0a=23.5475 FFj0B= 1.0836 FFj0b= 8.4540 FFj0C=-0.4340 FFj0c= 4.1196 FFj0D= 0.0214\nFFj2A= 3.7449 FFj2a=13.8944 FFj2B= 2.6453 FFj2b= 4.8634 FFj2C=-0.5218 FFj2c= 3.1919 FFj2D= 0.0009\nFFj4A=-1.0540 FFj4a=16.6055 FFj4B= 0.4339 FFj4b= 6.5119 FFj4C= 0.6746 FFj4c= 2.5993 FFj4D=-0.0011\nFFj6A=-0.1793 FFj6a=11.8961 FFj6B=-0.2269 FFj6b= 5.4280 FFj6C= 0.3291 FFj6c= 1.7008 FFj6D= 0.0030\n"],
"U5+",	["5f1", "FFj0A= 0.3650 FFj0a=19.8038 FFj0B= 3.2199 FFj0b= 6.2818 FFj0C=-2.6077 FFj0c= 5.3010 FFj0D= 0.0233\nFFj2A= 3.0724 FFj2a=12.5460 FFj2B= 2.3076 FFj2b= 5.2314 FFj2C=-0.0644 FFj2c= 1.4738 FFj2D= 0.0035\nFFj4A=-0.9588 FFj4a=16.4851 FFj4B= 0.1576 FFj4b= 6.4397 FFj4C= 0.7785 FFj4c= 2.6402 FFj4D=-0.0010\nFFj6A=-0.0399 FFj6a=11.8909 FFj6B=-0.3458 FFj6b= 5.5803 FFj6C= 0.3340 FFj6c= 1.6448 FFj6D= 0.0029\n"],
"Np3+",	["5f4", "FFj0A= 0.5157 FFj0a=20.8654 FFj0B= 2.2784 FFj0b= 5.8930 FFj0C=-1.8163 FFj0c= 4.8457 FFj0D= 0.0211\nFFj2A= 3.7170 FFj2a=15.1333 FFj2B= 2.3216 FFj2b= 5.5025 FFj2C=-0.0275 FFj2c= 0.7996 FFj2D= 0.0052\nFFj4A=-0.9029 FFj4a=16.5858 FFj4B= 0.4006 FFj4b= 6.4699 FFj4C= 0.6545 FFj4c= 2.5631 FFj4D=-0.0004\nFFj6A=-0.2427 FFj6a=11.8444 FFj6B=-0.1129 FFj6b= 5.3774 FFj6C= 0.2848 FFj6c= 1.5676 FFj6D= 0.0022\n"],
"Np4+",	["5f3", "FFj0A= 0.4206 FFj0a=19.8046 FFj0B= 2.8004 FFj0b= 5.9783 FFj0C=-2.2436 FFj0c= 4.9848 FFj0D= 0.0228\nFFj2A= 2.9203 FFj2a=14.6463 FFj2B= 2.5979 FFj2b= 5.5592 FFj2C=-0.0301 FFj2c= 0.3669 FFj2D= 0.0141\nFFj4A=-0.9887 FFj4a=12.4415 FFj4B= 0.5918 FFj4b= 5.2941 FFj4C= 0.5306 FFj4c= 2.2625 FFj4D=-0.0021\nFFj6A=-0.2436 FFj6a= 9.5988 FFj6B=-0.1317 FFj6b= 4.1014 FFj6C= 0.3029 FFj6c= 1.5447 FFj6D= 0.0019\n"],
"Np5+",	["5f2", "FFj0A= 0.3692 FFj0a=18.1900 FFj0B= 3.1510 FFj0b= 5.8500 FFj0C=-2.5446 FFj0c= 4.9164 FFj0D= 0.0248\nFFj2A= 2.3308 FFj2a=13.6540 FFj2B= 2.7219 FFj2b= 5.4935 FFj2C=-0.1357 FFj2c= 0.0493 FFj2D= 0.1224\nFFj4A=-0.8146 FFj4a=16.5809 FFj4B=-0.0055 FFj4b= 6.4751 FFj4C= 0.7956 FFj4c= 2.5623 FFj4D=-0.0004\nFFj6A=-0.1157 FFj6a= 9.5649 FFj6B=-0.2654 FFj6b= 4.2599 FFj6C= 0.3298 FFj6c= 1.5494 FFj6D= 0.0025\n"],
"Np6+",	["5f1", "FFj0A= 0.2929 FFj0a=17.5611 FFj0B= 3.4866 FFj0b= 5.7847 FFj0C=-2.8066 FFj0c= 4.8707 FFj0D= 0.0267\nFFj2A= 1.8245 FFj2a=13.1803 FFj2B= 2.8508 FFj2b= 5.4068 FFj2C=-0.1579 FFj2c= 0.0444 FFj2D= 0.1438\nFFj4A=-0.6738 FFj4a=16.5531 FFj4B=-0.2297 FFj4b= 6.5055 FFj4C= 0.8513 FFj4c= 2.5528 FFj4D=-0.0003\nFFj6A=-0.0128 FFj6a= 9.5692 FFj6B=-0.3611 FFj6b= 4.3035 FFj6C= 0.3419 FFj6c= 1.5406 FFj6D= 0.0032\n"],
"Pu3+",	["5f5", "FFj0A= 0.3840 FFj0a=16.6793 FFj0B= 3.1049 FFj0b= 5.4210 FFj0C=-2.5148 FFj0c= 4.5512 FFj0D= 0.0263\nFFj2A= 2.0885 FFj2a=12.8712 FFj2B= 2.5961 FFj2b= 5.1896 FFj2C=-0.1465 FFj2c= 0.0393 FFj2D= 0.1343\nFFj4A=-0.7014 FFj4a=16.3687 FFj4B=-0.1162 FFj4b= 6.6971 FFj4C= 0.7778 FFj4c= 2.4502 FFj4D= 0.0000\nFFj6A=-0.0364 FFj6a= 9.5721 FFj6B=-0.3181 FFj6b= 4.3424 FFj6C= 0.3210 FFj6c= 1.5233 FFj6D= 0.0041\n"],
"Pu4+",	["5f4", "FFj0A= 0.4934 FFj0a=16.8355 FFj0B= 1.6394 FFj0b= 5.6384 FFj0C=-1.1581 FFj0c= 4.1399 FFj0D= 0.0248\nFFj2A= 2.7244 FFj2a=12.9262 FFj2B= 2.3387 FFj2b= 5.1633 FFj2C=-0.1300 FFj2c= 0.0457 FFj2D= 0.1177\nFFj4A=-0.9160 FFj4a=12.2027 FFj4B= 0.4891 FFj4b= 5.1274 FFj4C= 0.5290 FFj4c= 2.1487 FFj4D=-0.0022\nFFj6A=-0.2394 FFj6a= 7.8367 FFj6B=-0.0785 FFj6b= 4.0243 FFj6C= 0.2643 FFj6c= 1.3776 FFj6D= 0.0012\n"],
"Pu5+",	["5f3", "FFj0A= 0.3888 FFj0a=16.5592 FFj0B= 2.0362 FFj0b= 5.6567 FFj0C=-1.4515 FFj0c= 4.2552 FFj0D= 0.0267\nFFj2A= 2.1409 FFj2a=12.8319 FFj2B= 2.5664 FFj2b= 5.1522 FFj2C=-0.1338 FFj2c= 0.0457 FFj2D= 0.1210\nFFj4A=-0.7035 FFj4a=16.3601 FFj4B=-0.0979 FFj4b= 6.7057 FFj4C= 0.7726 FFj4c= 2.4475 FFj4D= 0.0000\nFFj6A=-0.1090 FFj6a= 7.8188 FFj6B=-0.2243 FFj6b= 4.1000 FFj6C= 0.2947 FFj6c= 1.4040 FFj6D= 0.0015\n"],
"Pu6+",	["5f2", "FFj0A= 0.3172 FFj0a=16.0507 FFj0B= 3.4654 FFj0b= 5.3507 FFj0C=-2.8102 FFj0c= 4.5133 FFj0D= 0.0281\nFFj2A= 1.7262 FFj2a=12.3240 FFj2B= 2.6652 FFj2b= 5.0662 FFj2C=-0.1695 FFj2c= 0.0406 FFj2D= 0.1550\nFFj4A=-0.5560 FFj4a=16.3215 FFj4B=-0.3046 FFj4b= 6.7685 FFj4C= 0.8146 FFj4c= 2.4259 FFj4D= 0.0001\nFFj6A=-0.0001 FFj6a= 7.8196 FFj6B=-0.3354 FFj6b= 4.1439 FFj6C= 0.3097 FFj6c= 1.4027 FFj6D= 0.0020\n"],
"Am2+",   ["5f7", "FFj0A= 0.4743 FFj0a=21.7761 FFj0B= 1.5800 FFj0b= 5.6902 FFj0C=-1.0779 FFj0c= 4.1451 FFj0D= 0.0218\nFFj2A= 3.5237 FFj2a=15.9545 FFj2B= 2.2855 FFj2b= 5.1946 FFj2C=-0.0142 FFj2c= 0.5853 FFj2D= 0.0033\nFFj4A=-0.7433 FFj4a=16.4163 FFj4B= 0.3481 FFj4b= 6.7884 FFj4C= 0.6014 FFj4c= 2.3465 FFj4D= 0.0000\nFFj6A=-0.3176 FFj6a= 7.8635 FFj6B= 0.0771 FFj6b= 4.1611 FFj6C= 0.2194 FFj6c= 1.3387 FFj6D= 0.0018\n"],
"Am3+",	["5f6", "FFj0A= 0.4239 FFj0a=19.5739 FFj0B= 1.4573 FFj0b= 5.8722 FFj0C=-0.9052 FFj0c= 3.9682 FFj0D= 0.0238\nFFj2A= 2.8622 FFj2a=14.7328 FFj2B= 2.4099 FFj2b= 5.1439 FFj2C=-0.1326 FFj2c= 0.0309 FFj2D= 0.1233\nFFj4A=-0.8092 FFj4a=12.8542 FFj4B= 0.4161 FFj4b= 5.4592 FFj4C= 0.5476 FFj4c= 2.1721 FFj4D=-0.0011\nFFj6A=-0.3159 FFj6a= 6.9821 FFj6B= 0.0682 FFj6b= 3.9948 FFj6C= 0.2141 FFj6c= 1.1875 FFj6D=-0.0015\n"],
"Am4+",	["5f5", "FFj0A= 0.3737 FFj0a=17.8625 FFj0B= 1.3521 FFj0b= 6.0426 FFj0C=-0.7514 FFj0c= 3.7199 FFj0D= 0.0258\nFFj2A= 2.4141 FFj2a=12.9478 FFj2B= 2.3687 FFj2b= 4.9447 FFj2C=-0.2490 FFj2c= 0.0215 FFj2D= 0.2371\nFFj4A=-0.8548 FFj4a=12.2257 FFj4B= 0.3037 FFj4b= 5.9087 FFj4C= 0.6173 FFj4c= 2.1881 FFj4D=-0.0016\nFFj6A=-0.1787 FFj6a= 7.8805 FFj6B=-0.1274 FFj6b= 4.0898 FFj6C= 0.2565 FFj6c= 1.3152 FFj6D= 0.0017\n"],
"Am5+",	["5f4", "FFj0A= 0.2956 FFj0a=17.3725 FFj0B= 1.4525 FFj0b= 6.0734 FFj0C=-0.7755 FFj0c= 3.6619 FFj0D= 0.0277\nFFj2A= 2.0109 FFj2a=12.0534 FFj2B= 2.4155 FFj2b= 4.8358 FFj2C=-0.2264 FFj2c= 0.0275 FFj2D= 0.2128\nFFj4A=-0.6538 FFj4a=15.4625 FFj4B=-0.0948 FFj4b= 5.9971 FFj4C= 0.7295 FFj4c= 2.2968 FFj4D= 0.0000\nFFj6A=-0.0927 FFj6a= 6.0727 FFj6B=-0.2227 FFj6b= 3.7840 FFj6C= 0.2916 FFj6c= 1.3723 FFj6D= 0.0026\n"],
"Am6+",	["5f3", "FFj0A= 0.2302 FFj0a=16.9533 FFj0B= 1.4864 FFj0b= 6.1159 FFj0C=-0.7457 FFj0c= 3.5426 FFj0D= 0.0294\nFFj2A= 1.6778 FFj2a=11.3372 FFj2B= 2.4531 FFj2b= 4.7247 FFj2C=-0.2043 FFj2c= 0.0337 FFj2D= 0.1892\nFFj4A=-0.5390 FFj4a=15.4491 FFj4B=-0.2689 FFj4b= 6.0169 FFj4C= 0.7711 FFj4c= 2.2970 FFj4D= 0.0002\nFFj6A= 0.0152 FFj6a= 6.0788 FFj6B=-0.3549 FFj6b= 3.8610 FFj6C= 0.3125 FFj6c= 1.4031 FFj6D= 0.0036\n"],
"Cm3+",	["5f7", ""]
);

# Table of Z(k) coefficients for beyond dipole calculations
%zktab = (
"Ce2+",	"",
"Ce3+",	"Z1c0=+1.07142857  Z1c2=+1.71428571\n                  Z3c2=+0.09583148  Z3c4=+0.31943828\n                                    Z5c4=+0.00360750  Z5c6=+0.04329004\n                                                      Z7c6=+0.00000000\n",
"Pr3+",	"Z1c0=+1.60000000  Z1c2=+2.63111111\n                  Z3c2=-0.09865803  Z3c4=-0.13453368\n                                    Z5c4=-0.01836547  Z5c6=-0.11678557\n                                                      Z7c6=+0.00223313\n",
"Pr4+",	"",
"Nd2+",	"",
"Nd3+",	"Z1c0=+1.63636364  Z1c2=+2.95041322\n                  Z3c2=-0.20896503  Z3c4=-0.25329095\n                                    Z5c4=+0.03820789  Z5c6=+0.14258681\n                                                      Z7c6=-0.00614959\n",
"Pm3+",	"Z1c0=+1.2         Z1c2=+2.71515152\n                  Z3c2=-0.10866182  Z3c4=-0.07024602\n                                    Z5c4=+0.00509439  Z5c6=-0.06449156\n                                                      Z7c6=+0.00232318\n",
"Sm2+",	"",
"Sm3+",	"Z1c0=+0.35714286  Z1c2=+1.93650794\n                  Z3c2=+0.04614109  Z3c4=+0.06291966\n                                    Z5c4=-0.00426341  Z5c6=+0.00335243\n                                                      Z7c6=0\n",
"Eu2+",	"",
"Eu3+",	"",
"Gd2+",	"",
"Gd3+",	"",
"Tb2+",	"",
"Tb3+",	"Z1c0=4.5          Z1c2=1.666666666\n                  Z3c2=0            Z3c4=0.28459047\n                                     Z5c4=-0.05050505  Z5c6=+0.03496503\n                                                      Z7c6=-0.01131382\n",
"Tb4+",	"",
"Dy2+",	"",
"Dy3+",	"Z1c0=5.0          Z1c2=2.666666667\n                  Z3c2=-0.22360680  Z3c4=-0.08131156\n                                    Z5c4=+0.02525253  Z5c6=-0.12432012\n                                                      Z7c6=+0.05656908\n",
"Ho2+",	"",
"Ho3+",	"Z1c0=5.0          Z1c2=3.066666667\n                  Z3c2=-0.31304952  Z3c4=-0.28459047\n                                    Z5c4=+0.12626263  Z5c6=+0.14763015\n                                                      Z7c6=-0.11313815\n",
"Er2+",	"",
"Er3+",	"Z1c0=4.5          Z1c2=2.933333333\n                  Z3c2=-0.13416408  Z3c4=-0.16262313\n                                    Z5c4=-0.02525253  Z5c6=-0.04662005\n                                                      Z7c6=+0.11313815\n",
"Tm2+",	"",
"Tm3+",	"Z1c0=3.5          Z1c2=2.333333333\n                  Z3c2=+0.22360680  Z3c4=+0.08131156\n                                    Z5c4=-0.17676768  Z5c6=-0.02719503\n                                                      Z7c6=-0.05656908\n",
"Yb3+",	"Z1c0=2.0          Z1c2=1.333333333\n                  Z3c2=+0.44721360  Z3c4=+0.16262313\n                                    Z5c4=+0.10101010  Z5c6=+0.01554002\n                                                      Z7c6=+0.01131382\n",
);
