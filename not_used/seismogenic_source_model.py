###########################################################
# Convert shapefile to nrml
###########################################################
#convert csv file returned by smoothed seismicity to shapefile
import csv

filename = 'smoothed_seismicity.csv'

srcm = {
        'src_id':[],
        'src_name':[],
        'lon':[],
        'lat':[],
        'tect_reg':[],
        'usd':[],
        'lsd':[],
        'msr':[],
        'rar':[],
        'mfd_type':[],
        'min_mag':[],
        'max_mag':[],
        'a_val':[],
        'b_val':[],
        'num_npd':[]
    }

with open(filename,'r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    #skip header
    for i in range(1):
        next(reader,None)
    for row in reader:
    #Store data
        #remove whitespace
        row = [x.strip() for x in row]

        srcm['lon'].append(row[0])
        srcm['lat'].append(row[1])
        srcm['a_val'].append(row[4])
        srcm['b_val'].append(row[5])

#create depending on 'num_npd'
#weight_i,strike_i,rake_i,dip_i
#create depending on 'num_hdd'
#hdd_d_i,hdd_w_i

#store to shapefile
import osgeo.ogr
spatialReference = osgeo.osr.SpatialReference()
spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
shapeData = driver.CreateDataSource('points-shifted.shp')

layer = shapeData.CreateLayer('point sources', spatialReference, osgeo.ogr.wkbPoint)
layerDefinition = layer.GetLayerDefn()

point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
point.SetPoint(0, 474595, 4429281)

featureIndex = 0
feature = osgeo.ogr.Feature(layerDefinition)
feature.SetGeometry(point)
feature.SetFID(featureIndex)

layer.CreateFeature(feature)

shapeData.Destroy()

#import osgeo.ogr
##required fields
# src_id,src_name,tect_reg,usd,lsd,msr,rar,mfd_type,min_mag,max_mag,a_val,   b_val,num_npd,weight_1,strike_1,rake_1,dip_1,weight_2,strike_2,rake_2,dip_2, num_hdd,hdd_d_1,hdd_w_1,hdd_d_2,hdd_w_2
#
#filename='t.shp'
## Open a Shapefile, and get field names
#source = osgeo.ogr.Open(filename, 1)
#layer = source.GetLayer()
#layer_defn = layer.GetLayerDefn()
#field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
#
## Add a new field
#new_field = osgeo.ogr.FieldDefn('MYFLD', osgeo.ogr.OFTInteger)
#layer.CreateField(new_field)
##
### Close the Shapefile
##source = None
##
##import oq_input.source_model_converter as smc
##smc.shp2nrml('test_point.shp','test')
#
