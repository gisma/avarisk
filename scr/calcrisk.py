'''
Created on 01.09.2012

@author: creu
'''
import os
import sys
import csv 
import time
from osgeo import ogr

Startzeit = time.time()
# Directories
WorkHomePath='/home/creu/workspace/avalanche/riskcalc/'
ASTERDataPath = WorkHomePath+'ASTER/'
#SRTMDataPath = WorkHomePath+'SRTM/'
#ATKISDataPath = WorkHomePath+'ATKIS/'
WorkDataPath = WorkHomePath+'tmpdata/'
# Pfade SAGA GIS 
sgbin="/usr/local/bin/"
PROJ_ID='32632'
TIF='.tif'
SGRD='.sgrd'
mapset = PROJ_ID 
print mapset
print WorkHomePath

# globale variable (vorlaeufig)
# snow calculation
WINDDIR ='99'
WINDSPEED='99'
WIND=0
# woods parameter
s_tol= 0.9
c_tol= 0.0001
woodkernelsize=3
# obsolet
processall=0
# projizieren 0 =nein 1 =ja
projection=1
# delete tifs 1=yes
cleantif=0
# Max solar input angel
ALPHA_MAX=205
# choose filter of DEM default is 3
DEMfilter=2 # 1=mdenoise SUN et al, 2 =adaptive Lee et al. SAGA GIS 3=destripe PEREGO
# denoise factor for visualitation 0.76 
dnoi=0.98
#destripe params for ASTER
pereg_ang1='1.8'
pereg_rad1='1'
pereg_dist1='3'

remove=[]

# hol die datei vom server
os.system("wget http://giswerk.org/test/php/llb.csv | mv llb.csv "+WorkHomePath+'input/' )
#file = urllib2.urlopen('http://giswerk.org/test/php/llb.csv')

# import der LLB Daten
params=[]
reader = csv.reader(open(WorkHomePath+'input/llb.csv', "rb"), delimiter=";")
for row in reader:
    params.append(row)

LON=params [0][0]
LAT=params  [0][1]
Region=params [0][2]
Hoehengrenze=params [0][3]
LWStop=params [0][4]
LWSlow=params [0][5]
SchneeTypTop=params [0][6]
SchneeTypLow=params [0][7]
N=params [0][8]
NE=params [0][9]
E=params [0][10]
SE=params [0][11]
S=params [0][12]
SW=params [0][13]
W=params [0][14]
NW=params [0][15]    
WINDDIR =params [0][16] 
WINDSPEED=params [0][17] 

# generate expositon rules (aspect)
## safe exposition=0 dangerous expositions=1  
## due to integer restrictions the ranges of degree are multiplied by 100
if N=="1":
    bad_ex_rule_n1 = '0 thru 225 = 1'
    bad_ex_rule_n2 = '3375 thru 3600 = 1'
    good_ex_rule_n1 = '0 thru 225 = 0'
    good_ex_rule_n2 = '3375 thru 3600 = 0'
if N=="":
    bad_ex_rule_n1 = '0 thru 225 = 0'
    bad_ex_rule_n2 = '3375 thru 3600 = 0'
    good_ex_rule_n1 = '0 thru 225 = 1'
    good_ex_rule_n2 = '3375 thru 3600 = 1'    
if NE=="1":
    bad_ex_rule_ne = '226 thru 675  = 1'
    good_ex_rule_ne = '226 thru 675  = 0'
if NE=="":
    bad_ex_rule_ne = '226 thru 675 = 0'
    good_ex_rule_ne = '226 thru 675 = 1'
if E=="1":
    bad_ex_rule_e  = '676 thru 1125 = 1'
    good_ex_rule_e  = '676 thru 1125 = 0'
if E=="":
    bad_ex_rule_e  = '676 thru 1125 = 0'
    good_ex_rule_e  = '676 thru 1125 = 1'
if SE=="1":
    bad_ex_rule_se = '1126 thru 1575 = 1'
    good_ex_rule_se = '1126 thru 1575 = 0'
if SE=="":
    bad_ex_rule_se = '1126 thru 1575 = 0'
    good_ex_rule_se = '1126 thru 1575 = 1'
if S=="1":
    bad_ex_rule_s  = '1576 thru 2025 = 1'
    good_ex_rule_s  = '1576 thru 2025 = 0'
if S=="":
    bad_ex_rule_s   = '1576 thru 2025 = 0'
    good_ex_rule_s   = '1576 thru 2025 = 1'
if SW=="1":
    bad_ex_rule_sw  = '2026 thru 2475 = 1'
    good_ex_rule_sw  = '2026 thru 2475 = 0'
if SW=="":
    bad_ex_rule_sw  = '2026 thru 2475 = 0'
    good_ex_rule_sw  = '2026 thru 2475 = 1'
if W=="1":
    bad_ex_rule_w   = '2476 thru 2925 = 1'
    good_ex_rule_w   = '2476 thru 2925 = 0'
if W=="":
    bad_ex_rule_w   = '2476 thru 2925 = 0'
    good_ex_rule_w   = '2476 thru 2925 = 1'
if NW=="1":
    bad_ex_rule_nw  = '2926 thru 3375 = 1'
    good_ex_rule_nw  = '2926 thru 3375 = 0'
if NW=="":
    bad_ex_rule_nw  = '2926 thru 3375 = 0'
    good_ex_rule_nw  = '2926 thru 3375 = 1'
fobj = open(WorkHomePath+"input/rulesbadExpo.rules", "w")
fobj.write(bad_ex_rule_n1 + "\n")
fobj.write(bad_ex_rule_ne + "\n")
fobj.write(bad_ex_rule_e + "\n")
fobj.write(bad_ex_rule_se + "\n")
fobj.write(bad_ex_rule_s + "\n")
fobj.write(bad_ex_rule_sw + "\n")
fobj.write(bad_ex_rule_w + "\n")
fobj.write(bad_ex_rule_nw + "\n")
fobj.write(bad_ex_rule_n2 + "\n") 
fobj.write('end')
fobj.close()
fobj = open(WorkHomePath+"input/rulesgoodExpo.rules", "w")
fobj.write(good_ex_rule_n1 + "\n")
fobj.write(good_ex_rule_ne + "\n")
fobj.write(good_ex_rule_e + "\n")
fobj.write(good_ex_rule_se + "\n")
fobj.write(good_ex_rule_s + "\n")
fobj.write(good_ex_rule_sw + "\n")
fobj.write(good_ex_rule_w + "\n")
fobj.write(good_ex_rule_nw + "\n")
fobj.write(good_ex_rule_n2 + "\n") 
fobj.write('end')
fobj.close()


# generate LWS GRM TopZone Rules => GRM
## good exposition, potential spinning snow,potential radiation input last hours, north exposition, bad vibes
## generate reclass for slope
## dangerous slopes are true 
## due to integer restrictions the ranges of degree are multiplied by 10
if LWStop=="1":
    LWStop_rule_green   = "0 thru 40 = 1 "
    LWStop_rule_orange  = "41 thru 90 = 2 "
    LWStop_rule_red     = "91 thru 92 = 3 "
if LWStop=="2":
    LWStop_rule_green   = "0 thru 35 = 1 " 
    LWStop_rule_orange  = "36 thru 40 = 2 "
    LWStop_rule_red     = "41 thru 90 = 3 "
if LWStop=="3":
    LWStop_rule_green   = "0 thru 30 = 1 " 
    LWStop_rule_orange  = "31 thru 35 = 2 "
    LWStop_rule_red     = "36 thru 90 = 3 "
if LWStop=="4":
    LWStop_rule_green   = "0 thru 10 = 1 " 
    LWStop_rule_orange  = "11 thru 30 = 2 "
    LWStop_rule_red     = "31 thru 90 = 3 "
fobj = open(WorkHomePath+"input/GRMbadExpotop.rules", "w")
fobj.write(LWStop_rule_green + "\n")
fobj.write(LWStop_rule_orange + "\n")
fobj.write(LWStop_rule_red + "\n")
fobj.write('end')
fobj.close()

# generate LWS GRM lowZone Rules
if LWSlow=="1":
    LWSlow_rule_green   = "0 thru 40 = 1 "
    LWSlow_rule_orange  = "41 thru 90 = 2 "
    
if LWSlow=="2":
    LWSlow_rule_green   = "0 thru 35 = 1 " 
    LWSlow_rule_orange  = "36 thru 40 = 2 "
    LWSlow_rule_red     = "41 thru 90 = 3 "
if LWSlow=="3":
    LWSlow_rule_green   = "0 thru 30 = 1 " 
    LWSlow_rule_orange  = "31 thru 35 = 2 "
    LWSlow_rule_red     = "36 thru 90 = 3 "
if LWSlow=="4":
    LWSlow_rule_green   = "0 thru 10 = 1 " 
    LWSlow_rule_orange  = "11 thru 30 = 2 "
    LWSlow_rule_red     = "31 thru 90 = 3 "
fobj = open(WorkHomePath+"input/GRMbadExpolow.rules", "w")
fobj.write(LWSlow_rule_green + "\n")
fobj.write(LWSlow_rule_orange + "\n")
fobj.write(LWSlow_rule_red + "\n")
fobj.write('end')
fobj.close()


## good exposition = LWS -1
## dangerous slopes are true 
## due to integer restrictions the ranges of degree are multiplied by 10
if LWStop=="1":
    LWStop_rule_green   = "0 thru 40 = 1 "
    LWStop_rule_orange  = "41 thru 90 = 2 "
    
if LWStop=="2":
    LWStop_rule_green   = "0 thru 40 = 1 "
    LWStop_rule_orange  = "41 thru 90 = 2 "
    
if LWStop=="3":
    LWStop_rule_green   = "0 thru 35 = 1 " 
    LWStop_rule_orange  = "36 thru 40 = 2 "
    LWStop_rule_red     = "41 thru 90 = 3 "
if LWStop=="4":
    LWStop_rule_green   = "0 thru 30 = 1 " 
    LWStop_rule_orange  = "31 thru 35 = 2 "
    LWStop_rule_red     = "36 thru 90 = 3 "
fobj = open(WorkHomePath+"input/GRMgoodExpotop.rules", "w")
fobj.write(LWStop_rule_green + "\n")
fobj.write(LWStop_rule_orange + "\n")
fobj.write(LWStop_rule_red + "\n")
fobj.write('end')
fobj.close()

# generate LWS GRM lowZone Rules
if LWStop=="1":
    LWStop_rule_green   = "0 thru 40 = 1 "
    LWStop_rule_orange  = "41 thru 90 = 2 "
    
if LWStop=="2":
    LWStop_rule_green   = "0 thru 40 = 1 "
    LWStop_rule_orange  = "41 thru 90 = 2 "
    
if LWStop=="3":
    LWStop_rule_green   = "0 thru 35 = 1 " 
    LWStop_rule_orange  = "36 thru 40 = 2 "
    LWStop_rule_red     = "41 thru 90 = 3 "
if LWStop=="4":
    LWStop_rule_green   = "0 thru 30 = 1 " 
    LWStop_rule_orange  = "31 thru 35 = 2 "
    LWStop_rule_red     = "36 thru 90 = 3 "
fobj = open(WorkHomePath+"input/GRMgoodExpolow.rules", "w")
fobj.write(LWSlow_rule_green + "\n")
fobj.write(LWSlow_rule_orange + "\n")
fobj.write(LWSlow_rule_red + "\n")
fobj.write('end')
fobj.close()

# generate TopZoneMask
## 0=lowzone 1=topzone
low_rule1    = '0 thru '+ Hoehengrenze+' = 1 '
low_rule2    = Hoehengrenze+' thru 99999 = 0 '
top_rule1   = '0 thru '+ Hoehengrenze+' = 0 '
top_rule2   = Hoehengrenze+' thru 99999 = 1'

fobj = open(WorkHomePath+"input/rulesLowMask.rules", "w")
fobj.write(low_rule1  + "\n")
fobj.write(low_rule2  + "\n")
fobj.write('end')
fobj.close()
fobj = open(WorkHomePath+"input/rulesTopMask.rules", "w")
fobj.write(top_rule1  + "\n")
fobj.write(top_rule2  + "\n")
fobj.write('end')
fobj.close()

print " LLB Import successful"
print ' generated static reclass rules -> AltitudeMask -> GRM Scheme, Exposition' 

print " Start now DEM data preprocessing this takes a while - be calm :)"

#generate ASTERID according to the cutting window
ASTER_ID = mapset #str((float(LON)-0.3))+' '+str((float(LAT)+0.2))+' '+str((float(LON)+0.3))+' '+str((float(LAT)-0.2))
print  ASTER_ID
os.system('rm '+WorkDataPath+'*')
# cut a window 0.6x0,4 deg out of alps.dem
gdal_trans_cmd='gdal_translate -of GTiff -projwin '+str((float(LON)-0.3))+' '+str((float(LAT)+0.2))+' '+str((float(LON)+0.3))+' '+str((float(LAT)-0.2))+' '+ASTERDataPath+'/latlon/alpen.tiff '+WorkDataPath+'DEMlatlon.tif'

#gwarp='gdalwarp -s_srs EPSG:4326 -t_srs EPSG:32632 -r near -multi -of GTiff -te '' '+' '+ASTERDataPath+'/latlon/alpen.tiff '+WorkDataPath+'tmp.tif'
print gdal_trans_cmd
os.system(gdal_trans_cmd)
if projection==1:
    # projection from geographic to utm32n
    os.system('rm '+WorkDataPath+'DEMtmp.tif')
    GdalwarpCmd_4326to32632_tif ='gdalwarp -s_srs EPSG:4326 -t_srs EPSG:32632 -r near -multi -of GTiff '+WorkDataPath+'DEMlatlon.tif '+WorkDataPath+'DEMproj.tif'
    print GdalwarpCmd_4326to32632_tif
    os.system(GdalwarpCmd_4326to32632_tif)  

##import of Actual Subset of  DEM to SAGA 
sagacmd = sgbin+'saga_cmd libio_gdal  0 -GRIDS='+WorkDataPath+'DEMproj'+SGRD+' -FILES='+WorkDataPath+'DEMproj.tif' 
print sagacmd
os.system(sagacmd)


# try to destripe and despeckle ASTER data denoise SUN et al.
if DEMfilter==1:
    gdal_trans_cmd='gdal_translate -of AAIGrid  '+WorkDataPath+'DEMproj.tif '+WorkDataPath+'DEMproj.asc'
    print gdal_trans_cmd
    os.system(gdal_trans_cmd)  
    
    denoise_cmd='mdenoise -i '+WorkDataPath+'DEMproj.asc -e -n 20 -t '+str(dnoi)+'  -o ' +WorkDataPath+'DEMdenoise.asc'
    print denoise_cmd
    os.system(denoise_cmd)  
                                                                                 
    gdal_trans_cmd='gdal_translate -of Gtiff  '+WorkDataPath+'DEMdenoise.asc '+WorkDataPath+ASTER_ID+TIF
    print gdal_trans_cmd
    os.system(gdal_trans_cmd)  
    # export denoise grid to SAGA
    sagacmd = sgbin+'saga_cmd libio_gdal  0 -GRIDS='+WorkDataPath+ASTER_ID+SGRD+' -FILES='+WorkDataPath+ASTER_ID+TIF
    print sagacmd
    os.system(sagacmd)


# try to destripe and despeckle ASTER data Anisotropic Lee local statistics filter
if DEMfilter==2:
    sagacmd = sgbin+'saga_cmd libgrid_filter 3 -INPUT='+WorkDataPath+'DEMproj'+SGRD+' -RESULT='+WorkDataPath+ASTER_ID+SGRD+' -NOISE_ABS=1 -NOISE_REL=1 -METHOD=1'
    print sagacmd
    os.system(sagacmd)
    os.system('rm '+ WorkDataPath+'DEMproj'+SGRD)
    # export  filtered DEM -> GRA
    sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+SGRD+' -FILE='+WorkDataPath+ASTER_ID+TIF+' -FORMAT=1'
    print sagacmd
    os.system(sagacmd)    
    
    
    
if DEMfilter==3:  #destriping using the module of perego
    sagacmd = sgbin+'saga_cmd libcontrib_a_perego 5 -INPUT='+WorkDataPath+'DEMproj'+SGRD+' -RESULT3='+WorkDataPath+'DEMproj'+SGRD+' -RESULT2='+WorkDataPath+'lowpass1'+SGRD+' -RESULT1='+WorkDataPath+'lowpass2'+SGRD+' -ANG='+pereg_ang1+' -R='+pereg_rad1+' -D='+pereg_dist1
    print sagacmd
    os.system(sagacmd)
    sagacmd = sgbin+'saga_cmd libio_gdal  0 -GRIDS='+WorkDataPath+ASTER_ID+SGRD+' -FILES='+WorkDataPath+'DEMproj'+SGRD
    print sagacmd
    os.system(sagacmd)
    #gdal_trans_cmd='gdal_translate -of AAIGrid  '+WorkDataPath+'DEMproj.tif '+WorkDataPath+'DEMproj.asc'
    #print gdal_trans_cmd
    #os.system(gdal_trans_cmd)  
    #denoise_cmd='mdenoise -i '+WorkDataPath+'DEMproj.asc -e -n 20 -t '+str(dnoi)+'  -o ' +WorkDataPath+'DEMdenoise.asc'
    #print denoise_cmd
    #os.system(denoise_cmd)  
    #gdal_trans_cmd='gdal_translate -of Gtiff  '+WorkDataPath+'DEMdenoise.asc '+WorkDataPath+ASTER_ID+TIF
    #print gdal_trans_cmd
    #os.system(gdal_trans_cmd)  
    # remove temporary files
    os.system('rm '+WorkDataPath+'DEMproj.asc')
    os.system('rm '+WorkDataPath+'DEMdenoise.asc')
    os.system('rm '+WorkDataPath+ASTER_ID+'DEMproj.sdat')
    os.system('rm '+WorkDataPath+ASTER_ID+'DEMproj.mgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'DEMproj.sgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'DEMproj.prj')
    os.system('rm '+WorkDataPath+'lowpass1.*')
    os.system('rm '+WorkDataPath+'lowpass2.*')
    # export  filtered DEM -> GRA
#    sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+SGRD+' -FILE='+WorkDataPath+ASTER_ID+TIF+' -FORMAT=1'
#    print sagacmd
#    os.system(sagacmd)    
    
    
#------------------------------------

### setup der Variablen zur Initialisierung GRASS Umgebung
gisbase = os.environ['GISBASE'] = "/usr/lib/grass64"
gisdbase = '/home/creu/workspace/avalanche/riskcalc'
location = 'Allgau'
mapset = 'N47E010' #ASTER_ID
sys.path.append(os.path.join(os.environ['GISBASE']))

## jetzt import der GRASS Library
import grass.script as grass
import grass.script.setup as gsetup
# initialisierung von GRASS aus Python
gsetup.init(gisbase, gisdbase, location, mapset)
# check ob es stimmt
print grass.gisenv()
#grass version message
print grass.read_command("g.version", flags="g")

# Pfade SAGA GIS 
sgbin="/usr/local/bin/"

# import preprocessed DEM to Grass
grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+TIF, output=ASTER_ID,overwrite='true')

#------------------------------------------------------------------------------------------------

if processall ==1:
    # rauhigkeit
    sagacmd =  sgbin+'saga_cmd libta_morphometry 16 -DEM='+WorkDataPath+ASTER_ID+SGRD+ ' -TRI='+WorkDataPath+ASTER_ID+'TRI'+SGRD+ ' -RADIUS=3'
    print sagacmd
    os.system(sagacmd)
    sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+'TRI'+SGRD+ ' -FILE='+WorkDataPath+ASTER_ID+'TRI'+TIF+' -FORMAT=1'
    print sagacmd
    os.system(sagacmd)
    # convergence
    sagacmd =  sgbin+'saga_cmd libta_morphometry 2 -ELEVATION='+WorkDataPath+ASTER_ID+SGRD+ ' -CONVERGENCE='+WorkDataPath+ASTER_ID+'mconv'+SGRD+ ' -RADIUS=5'
    print sagacmd
    os.system(sagacmd)
    sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+'mconv'+SGRD+ ' -FILE='+WorkDataPath+ASTER_ID+'mconv'+TIF+' -FORMAT=1'
    print sagacmd
    os.system(sagacmd)
    os.system('rm '+WorkDataPath+ASTER_ID+'DAH.sdat')
    os.system('rm '+WorkDataPath+ASTER_ID+'DAH.mgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'DAH.sgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'DAH.prj')
    os.system('rm '+WorkDataPath+ASTER_ID+'TRI.sdat')
    os.system('rm '+WorkDataPath+ASTER_ID+'TRI.mgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'TRI.sgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'TRI.prj')
    os.system('rm '+WorkDataPath+ASTER_ID+'mconv.sdat')
    os.system('rm '+WorkDataPath+ASTER_ID+'mconv.mgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'mconv.sgrd')
    os.system('rm '+WorkDataPath+ASTER_ID+'mconv.prj')
    
    # grass Wood
    grass.run_command('r.param.scale', input=ASTER_ID, output=ASTER_ID+'wood', s_tol=s_tol, c_tol=c_tol, size=woodkernelsize, param='feature',overwrite='true')
    grass.run_command("r.out.gdal", input=ASTER_ID+'wood', output=WorkDataPath+ASTER_ID+"wood"+'.tif',overwrite='true')
    # windeffect snow 
    if WINDDIR != "99" and WINDSPEED !="99" and WIND==1:
        sagacmd = sgbin+'saga_cmd libta_morphometry 15 -DEM='+WorkDataPath+ASTER_ID+SGRD+' -DIR_CONST='+WINDDIR+' -LEN='+WINDSPEED+' -EFFECT='+WorkDataPath+ASTER_ID+'snowall'+SGRD+' -LUV='+WorkDataPath+ASTER_ID+'snowluv'+SGRD+' -LEE='+WorkDataPath+ASTER_ID+'snowlee'+SGRD
        print sagacmd
        os.system(sagacmd)
        sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+'snowall'+SGRD+ ' -FILE='+WorkDataPath+ASTER_ID+'snowall'+TIF+' -FORMAT=1'
        print sagacmd
        os.system(sagacmd)
        sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+'snowluv'+SGRD+ ' -FILE='+WorkDataPath+ASTER_ID+'snowluv'+TIF+' -FORMAT=1'
        print sagacmd
        os.system(sagacmd)
        sagacmd = sgbin+'saga_cmd libio_gdal  1 -GRIDS='+WorkDataPath+ASTER_ID+'snowlee'+SGRD+ ' -FILE='+WorkDataPath+ASTER_ID+'snowlee'+TIF+' -FORMAT=1'
        print sagacmd
        os.system(sagacmd)
        grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'snowall'+TIF, output=ASTER_ID+'snowall',overwrite='true')
        grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'snowlee'+TIF, output=ASTER_ID+'snowlee',overwrite='true')
        grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'snowluv'+TIF, output=ASTER_ID+'snowluv',overwrite='true')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowall.sdat')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowall.mgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowall.sgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowall.prj')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowluv.sdat')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowluv.mgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowluv.sgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowluv.prj')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowlee.sdat')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowlee.mgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowlee.sgrd')
        os.system('rm '+WorkDataPath+ASTER_ID+'snowlee.prj')


grass.run_command("r.slope.aspect", elevation=ASTER_ID, slope=ASTER_ID+'slope', aspect=ASTER_ID+'aspect',prec="int",overwrite='true')  
## Import von TIFF nach GRASS 
#    grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'aspect'+TIF, output='tmpaspect',overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+'aspect', formula="int(A*10)", outfile=ASTER_ID+'aspectrec', overwrite='true')
grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'TRI'+TIF, output=ASTER_ID+'TRI',overwrite='true')
grass.run_command("r.in.gdal", input=WorkDataPath+ASTER_ID+'DAH'+TIF, output=ASTER_ID+'DAH',overwrite='true')


    # loeschen der temporaeren TIFFs
if cleantif==1:
    os.system('rm '+WorkDataPath+ASTER_ID+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'slope'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'aspect'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'curv'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'hcurv'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'vcurv'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'mconv'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'TRI'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'DAH'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'snowall'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'snowlee'+TIF)
    os.system('rm '+WorkDataPath+ASTER_ID+'snowluv'+TIF)



# Die wichtigsten Daten und Analysen sind in GRASS jetzt kanns losgehen
# classify MASK Files for the Altitude Discrimination 0/1
grass.run_command('r.reclass', input=ASTER_ID, output=ASTER_ID+"UpperZoneMask", rules=WorkHomePath+"input/rulesTopMask.rules",overwrite='true')
grass.run_command('r.reclass', input=ASTER_ID, output=ASTER_ID+"LowerZoneMask", rules=WorkHomePath+"input/rulesLowMask.rules",overwrite='true')

# clasiffy the exposition 0 nonbad 1=bad
grass.run_command('r.reclass', input=ASTER_ID+'aspectrec', output=ASTER_ID+'badexposition', rules=WorkHomePath+"input/rulesbadExpo.rules",overwrite='true')
grass.run_command('r.reclass', input=ASTER_ID+'aspectrec', output=ASTER_ID+'goodexposition', rules=WorkHomePath+"input/rulesgoodExpo.rules",overwrite='true')


# GRM static reclassify
## reclassify GRM green,yellow red slopes BAD Exposition
grass.run_command('r.reclass', input=ASTER_ID+'slope', output=ASTER_ID+"GRM_badlow_SlopeRisk", rules=WorkHomePath+"input/GRMbadExpolow.rules",overwrite='true')
grass.run_command('r.reclass', input=ASTER_ID+'slope', output=ASTER_ID+"GRM_badtop_SlopeRisk", rules=WorkHomePath+"input/GRMbadExpotop.rules",overwrite='true')
# combine results bad exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRM_badtop_SlopeRisk",bmap=ASTER_ID+"UpperZoneMask", cmap=ASTER_ID+"GRM_badlow_SlopeRisk",dmap=ASTER_ID+"LowerZoneMask", emap=ASTER_ID+'badexposition', formula="(A*B+C*D)*E", outfile=ASTER_ID+"GRMgoodRisk", overwrite='true')
## reclassify GRM green,yellow red slopes good Exposition
grass.run_command('r.reclass', input=ASTER_ID+'slope', output=ASTER_ID+"GRM_goodlow_SlopeRisk", rules=WorkHomePath+"input/GRMgoodExpolow.rules",overwrite='true')
grass.run_command('r.reclass', input=ASTER_ID+'slope', output=ASTER_ID+"GRM_goodtop_SlopeRisk", rules=WorkHomePath+"input/GRMgoodExpotop.rules",overwrite='true')
# combine results good exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRM_goodtop_SlopeRisk",bmap=ASTER_ID+"UpperZoneMask", cmap=ASTER_ID+"GRM_goodlow_SlopeRisk",dmap=ASTER_ID+"LowerZoneMask", emap=ASTER_ID+'badexposition', formula="(A*B+C*D)*E", outfile=ASTER_ID+"GRMbadRisk", overwrite='true')
# GRM static Result combine good and bad Expositions 
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMgoodRisk",bmap=ASTER_ID+'GRMbadRisk', formula="A+B", outfile=ASTER_ID+"GRM_StaticRisk", overwrite='true')
grass.run_command("r.neighbors", input=ASTER_ID+"GRM_StaticRisk", output=ASTER_ID+"GRM_mStaticRisk" ,method="mode" ,size=3, overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'GRM_StaticRisk', output=WorkDataPath+ASTER_ID+"GRM_StaticRisk"+TIF,overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'GRM_mStaticRisk', output=WorkDataPath+ASTER_ID+"GRM_mStaticRisk"+TIF,overwrite='true')

# appends temp files for later removal 
remove.append(ASTER_ID+"GRM_goodlow_SlopeRisk")
remove.append(ASTER_ID+"GRM_goodtop_SlopeRisk")
remove.append(ASTER_ID+"GRM_badtop_SlopeRisk")
remove.append(ASTER_ID+"GRM_badlow_SlopeRisk")
remove.append(ASTER_ID+'GRMbadRisk')
remove.append(ASTER_ID+"GRMgoodRisk")


# GRM discrimant functions method
# high altitude bad exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWStop+"   <= -0.2 * A + 8.5   && B > "+Hoehengrenze+", 1,0)", outfile=ASTER_ID+"GRMtopgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWStop+"   >  -0.2 * A + 8.5   &&     "+LWStop+" <=-0.2 * A + 10.5 && B > "+Hoehengrenze+" , 2,0)", outfile=ASTER_ID+"GRMtopyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWStop+"   >   -0.2 * A + 10.5  && B > "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"GRMtopred", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMtopgreen",bmap=ASTER_ID+"GRMtopyellow",cmap=ASTER_ID+"GRMtopred",dmap=ASTER_ID+'badexposition', formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"GRMbadExpotop", overwrite='true')
# low altitude bad exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWSlow+"   <= -0.2 * A + 8.5   && B <= "+Hoehengrenze+", 1,0)", outfile=ASTER_ID+"GRMlowgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWSlow+"   >  -0.2 * A + 8.5   &&     "+LWSlow+" <=-0.2 * A + 10.5 && B <="+Hoehengrenze+" , 2,0)", outfile=ASTER_ID+"GRMlowyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+LWSlow+"   >   -0.2 * A + 10.5  && B <= "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"GRMlowred", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMlowgreen",bmap=ASTER_ID+"GRMlowyellow",cmap=ASTER_ID+"GRMlowred",dmap=ASTER_ID+'badexposition', formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"GRMbadExpolow", overwrite='true')
# high altitude good exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWStop)-1)+"   <= -0.2 * A + 8.5   && B > "+Hoehengrenze+", 1,0)", outfile=ASTER_ID+"GRMtopgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWStop)-1)+"   >  -0.2 * A + 8.5   &&     "+str(int(LWStop)-1)+" <=-0.2 * A + 10.5 && B > "+Hoehengrenze+" , 2,0)", outfile=ASTER_ID+"GRMtopyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWStop)-1)+"   >   -0.2 * A + 10.5  && B > "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"GRMtopred", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMtopgreen",bmap=ASTER_ID+"GRMtopyellow",cmap=ASTER_ID+"GRMtopred",dmap=ASTER_ID+'goodexposition', formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"GRMgoodExpotop", overwrite='true')
# low altitude good exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWSlow)-1)+"   <= -0.2 * A + 8.5   && B <= "+Hoehengrenze+", 1,0)", outfile=ASTER_ID+"GRMlowgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWSlow)-1)+"   >  -0.2 * A + 8.5   &&     "+str(int(LWSlow)-1)+" <=-0.2 * A + 10.5 && B <="+Hoehengrenze+" , 2,0)", outfile=ASTER_ID+"GRMlowyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if("+str(int(LWSlow)-1)+"   >   -0.2 * A + 10.5  && B <= "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"GRMlowred", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMlowgreen",bmap=ASTER_ID+"GRMlowyellow",cmap=ASTER_ID+"GRMlowred",dmap=ASTER_ID+'goodexposition', formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"GRMgoodExpolow", overwrite='true')
# GRM Result
grass.run_command("r.mapcalculator", amap=ASTER_ID+"GRMbadExpotop",bmap=ASTER_ID+'GRMbadExpolow', cmap=ASTER_ID+'GRMgoodExpolow',dmap=ASTER_ID+'GRMgoodExpotop',formula="A+B+C+D", outfile=ASTER_ID+"GRMrisk", overwrite='true')
grass.run_command("r.neighbors", input=ASTER_ID+"GRMrisk", output=ASTER_ID+"GRMmrisk" ,method="mode" ,size=5, overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'GRMrisk', output=WorkDataPath+ASTER_ID+"GRMrisk"+TIF,overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'GRMmrisk', output=WorkDataPath+ASTER_ID+"GRMmrisk"+TIF,overwrite='true')


# appends temp files for later removal 
remove.append(ASTER_ID+"GRMtopgreen")
remove.append(ASTER_ID+"GRMtopyellow")
remove.append(ASTER_ID+"GRMtopred")
remove.append(ASTER_ID+"GRMlowgreen")
remove.append(ASTER_ID+"GRMlowyellow")
remove.append(ASTER_ID+"GRMlowred")
remove.append(ASTER_ID+"GRMbadExpolow")
remove.append(ASTER_ID+"GRMgoodExpotop")
remove.append(ASTER_ID+"GRMgoodExpolow")
remove.append(ASTER_ID+"GRMbadExpotop")
# SC non linear discrimant funtions method 

# high altitude bad exposition 
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",formula="if( A <= 42.0  &&  A > 27.0 , 1,0)", outfile=ASTER_ID+"slopemask", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWStop+"   <= 262.31153- 31.561649*A+ 1.2864667* A^2- 0.01745925*A ^3   && B > "+Hoehengrenze+" , 1,0)", outfile=ASTER_ID+"badSCtopgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWStop+"   >  262.31153- 31.561649*A+ 1.2864667* A^2- 0.01745925*A ^3   &&     "+LWStop+" <= 17.67 - 0.605 * A + 0.0055 * A ^ 2 && B > "+Hoehengrenze+" , 2,0)", outfile=ASTER_ID+"badSCtopyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask", formula="if( C==1 && "+LWStop+"   >   17.67 - 0.605 * A + 0.0055 * A ^ 2  && B > "+Hoehengrenze+" , 3,0)", outfile=ASTER_ID+"badSCtopred", overwrite='true')
#grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if( A > 42.0  && B > "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"badSCtopredplus", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"badSCtopgreen",bmap=ASTER_ID+"badSCtopyellow",cmap=ASTER_ID+"badSCtopred", dmap=ASTER_ID+'badexposition',formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"badSCtopslope", overwrite='true')
# low altitude bad exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWSlow+"   <= 262.31153- 31.561649*A+ 1.2864667* A^2- 0.01745925*A ^3   && B <= "+Hoehengrenze+" , 1,0)", outfile=ASTER_ID+"badSClowgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask", formula="if( C==1 && "+LWSlow+"   >  262.31153- 31.561649*A+ 1.2864667* A^2- 0.01745925*A ^3  &&     "+LWSlow+" <= 17.67 - 0.605 * A + 0.0055 * A ^ 2 && B <= "+Hoehengrenze+"  , 2,0)", outfile=ASTER_ID+"badSClowyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWSlow+"   >   17.67 - 0.605 * A + 0.0055 * A ^ 2  && B <= "+Hoehengrenze+"   , 3,0)", outfile=ASTER_ID+"badSClowred", overwrite='true')
#grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if( A > 42.0  && B <= "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"badSClowredplus", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"badSClowgreen",bmap=ASTER_ID+"badSClowyellow",cmap=ASTER_ID+"badSClowred",dmap=ASTER_ID+'badexposition',formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"badSClowslope", overwrite='true')
# high altitude good exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWStop+"   <= 12.287 - 0.3909 * A + 0.0030 * A ^ 2   && B > "+Hoehengrenze+"  , 1,0)", outfile=ASTER_ID+"goodSCtopgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask", formula="if( C==1 && "+LWStop+"   >  12.287 - 0.3909 * A + 0.0030 * A ^ 2   &&     "+LWStop+" <=  10.166558 * 0.97367151 ^ A && B > "+Hoehengrenze+"  , 2,0)", outfile=ASTER_ID+"goodSCtopyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWStop+"   >    10.166558 * 0.97367151 ^ A  && B > "+Hoehengrenze+"  , 3,0)", outfile=ASTER_ID+"goodSCtopred", overwrite='true')
#grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, formula="if( A > 42.0  && B > "+Hoehengrenze+", 3,0)", outfile=ASTER_ID+"goodSCtopredplus", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"goodSCtopgreen",bmap=ASTER_ID+"goodSCtopyellow",cmap=ASTER_ID+"goodSCtopred", dmap=ASTER_ID+'goodexposition', formula="(A*D+B*D+C*D)", outfile=ASTER_ID+"goodSCtopslope", overwrite='true')
# low altitude good exposition
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID, cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWSlow+"   <= 12.287 - 0.3909 * A + 0.0030 * A ^ 2   && B <= "+Hoehengrenze+"  , 1,0)", outfile=ASTER_ID+"goodSClowgreen", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWSlow+"   >  12.287 - 0.3909 * A + 0.0030 * A ^ 2   &&     "+LWSlow+" <=  10.166558 * 0.97367151 ^ A && B <= "+Hoehengrenze+"  , 2,0)", outfile=ASTER_ID+"goodSClowyellow", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope",bmap=ASTER_ID,cmap=ASTER_ID+"slopemask",formula="if( C==1 && "+LWSlow+"   >    10.166558 * 0.97367151 ^ A  && B <= "+Hoehengrenze+"  , 3,0)", outfile=ASTER_ID+"goodSClowred", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"slope", formula="if(A > 42.0 , 3,0)", outfile=ASTER_ID+"redplus42", overwrite='true')
grass.run_command("r.mapcalculator", amap=ASTER_ID+"goodSClowgreen",bmap=ASTER_ID+"goodSClowyellow",cmap=ASTER_ID+"goodSClowred", dmap=ASTER_ID+'goodexposition',emap=ASTER_ID+"redplus42", formula="(A+B+C+E)*D", outfile=ASTER_ID+"goodSClowslope", overwrite='true')
# SC Result
grass.run_command("r.mapcalculator", amap=ASTER_ID+"badSCtopslope",bmap=ASTER_ID+"badSClowslope",cmap=ASTER_ID+"goodSCtopslope",dmap=ASTER_ID+"goodSClowslope", formula="A+B+C+D", outfile=ASTER_ID+"SCrisk", overwrite='true')
grass.run_command("r.neighbors", input=ASTER_ID+"SCrisk", output=ASTER_ID+"SCmrisk" ,method="mode" ,size=5, overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'SCrisk', output=WorkDataPath+ASTER_ID+"SCrisk"+TIF,overwrite='true')
grass.run_command("r.out.gdal", input=ASTER_ID+'SCmrisk', output=WorkDataPath+ASTER_ID+"SCmrisk"+TIF,overwrite='true')


remove.append(ASTER_ID+"badSCtopgreen")
remove.append(ASTER_ID+"badSCtopyellow")
remove.append(ASTER_ID+"badSCtopred")
remove.append(ASTER_ID+"badSClowgreen")
remove.append(ASTER_ID+"badSClowyellow")
remove.append(ASTER_ID+"badSClowred")
remove.append(ASTER_ID+"badSCtopslope")
remove.append(ASTER_ID+"badSClowslope")
remove.append(ASTER_ID+"goodSCtopgreen")
remove.append(ASTER_ID+"goodSCtopyellow")
remove.append(ASTER_ID+"goodSCtopred")
remove.append(ASTER_ID+"goodSClowgreen")
remove.append(ASTER_ID+"goodSClowyellow")
remove.append(ASTER_ID+"goodSClowred")
remove.append(ASTER_ID+"goodSCtopslope")
remove.append(ASTER_ID+"goodSClowslope")

#grass.run_command("g.remove", rast=remove) 

Endzeit = time.time()
print str((Endzeit-Startzeit))+' Sekunden '
print "Ok SnowCard + GRM  berechnet"
 



  

 
