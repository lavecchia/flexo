# common parameters

KEYFHOF = "FINAL HEAT OF FORMATION"
KEYCA = "COSMO AREA"
MOPACPATH = "/opt/mopac/MOPAC2016.exe"
DEFCORE = 3
DEFREPORFILE = "report.log"
DEFSUMMARYFILE = "summaryfile.log"
DEFINPUTEXT = ".mop"
DEFTEMPLATEEXT = ".mop"
RCONST=8.31 # gases constant in J K-1 mol-1
TEMP = 298 # temperature in K
#~ NS = 3.33679E+028 # numeral density in m-3
NS = 0.033330 # numeral density in m-3
PI = 3.14159
J2KCAL = 0.0002388458966275
CONFIG_FILENAME = 'config.in'
PERIODICTABLE = {'H' :1, 'He' :2, 'Li' :3, 'Be' :4, 'B' :5, 'C' :6, 'N' :7, 'O' :8, 'F' :9, 'Ne' :10, 'Na' :11, 'Mg' :12, 'Al' :13, 'Si' :14, 'P' :15, 'S' :16, 'Cl' :17, 'Ar' :18, 'K' :19, 'Ca' :20, 'Sc' :21, 'Ti' :22, 'V' :23, 'Cr' :24, 'Mn' :25, 'Fe' :26, 'Co' :27, 'Ni' :28, 'Cu' :29, 'Zn' :30, 'Ga' :31, 'Ge' :32, 'As' :33, 'Se' :34, 'Br' :35, 'Kr' :36}
TMPDIR = "tmp/"
BESTCONFNUM = 10 #"automatic" #number of best conformer stored
ENERGYTHRESHOLD = 1.0 #difference in kcal/mol to calculate 
RMSDTHRESHOLD = 0.5

