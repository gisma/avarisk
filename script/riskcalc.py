import os
import sys
import rpy2
import csv
import datetime





gisbase = os.environ['GISBASE']
gisdbase = '/home/creu/workspace/avalanche/riskcalc'
location = 'Allgau'
mapset = 'DGM'
override = 'TRUE'
gsetup.init(gisbase, gisdbase, location, mapset)
print grass.gisenv()

grass.message('Raster maps:')
for rast in grass.list_strings(type = 'rast'):
    print rast
print grass.read_command("g.version", flags="c") 
grass.run_command('r.param.scale', input='DGM_aktuell', output='test', s_tol=1.0, size=3, param='feature')

