# common parameters

from itertools import count
KEYFHOF = "FINAL HEAT OF FORMATION"
KEYCA = "COSMO AREA"
KEYDEFAULT = "PUT KEYWORDS HERE"
KEYCOSMO = "NR. ELEM.             COORDINATES     RADIUS  COSMO-CHARGE" #to position cursor in COSMO file
MOPACPATH = "/opt/mopac/MOPAC2016.exe"
DEFCORE = 3
DEFPARAMFILE = "list_param.in"
DEFREFFILE = "list_dgexp.in"
DEFREPORFILE = "report.log"
DEFSUMMARYFILE = "summaryfile.log"
DEFTEMPLATEDIR= "templates"
DEFINPUTEXT = ".mop"
DEFTEMPLATEEXT = ".mop"
DEFGASDIR = "gas"
DEFOPTIMIZEPARAM = "optimize.log"
RCONST=8.31 # gases constant in J K-1 mol-1
TEMP = 298 # temperature in K
#~ NS = 3.33679E+028 # numeral density in m-3
NS = 0.033330 # numeral density in m-3
PI = 3.14159
J2KCAL = 0.0002388458966275
step = ("%04i" % i for i in count(1)) #work with itertools to make 002, 003 format name
CONFIG_FILENAME = 'config.in'
PERIODICTABLE = {'H' :1, 'He' :2, 'Li' :3, 'Be' :4, 'B' :5, 'C' :6, 'N' :7, 'O' :8, 'F' :9, 'Ne' :10, 'Na' :11, 'Mg' :12, 'Al' :13, 'Si' :14, 'P' :15, 'S' :16, 'Cl' :17, 'Ar' :18, 'K' :19, 'Ca' :20, 'Sc' :21, 'Ti' :22, 'V' :23, 'Cr' :24, 'Mn' :25, 'Fe' :26, 'Co' :27, 'Ni' :28, 'Cu' :29, 'Zn' :30, 'Ga' :31, 'Ge' :32, 'As' :33, 'Se' :34, 'Br' :35, 'Kr' :36}
TMPDIR = "tmp/"
BESTCONFNUM = 10 #"automatic" #number of best conformer stored
ENERGYTHRESHOLD = 1.0 #difference in kcal/mol to calculate 
RMSDTHRESHOLD = 0.3
FINDROTATABLEPATH = "/home/mjl/Dropbox/CODIGOS/ConFind/findrotatable.pl"
