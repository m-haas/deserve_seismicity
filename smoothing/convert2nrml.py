###########################################################
# Convert model to nrml
###########################################################

import csv
import smoothed_seismicity2nrml as ss

input_file = 'grid_rate_0_50km_50_0.1_grid.xyz'
output_file = 'smoothed_seismicity_f.csv'
fieldnames = list(['Longitude','Latitude','Depth','Observed Count','Smoothed Rate','b-value'])
bval=0.84
mmax=7.8
mmin=4.0
ud = 1.
ld = 50
hdd = [(15,1)]

#xmin,xmax,ymin,ymax
bbox = [31.9,38.8,26.5,35.8]

cat = {'Longitude':[] ,
       'Latitude': [] ,
       'Depth': [] ,
       'Observed Count': []   ,
       'Smoothed Rate': [] ,
       'b-value': [] ,
}

with open(input_file,'r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    #has no header
    #for i in range(1):
    #    next(reader,None)
    for row in reader:
    #Store data
        #remove whitespace
        row = [x.strip() for x in row]
        #only rates within source region (boundaries included)
        lon=float(row[0])
        lat=float(row[1])
        rate=float(row[2])
        if lon >= bbox[0] and lon <= bbox[1] and lat >= bbox[2] and lat <= bbox[3]:
            cat['Longitude'].append(row[0])
            cat['Latitude'].append(row[1])
            cat['Depth'].append(0)
            cat['Observed Count'].append(0)
            cat['Smoothed Rate'].append(row[2])
            cat['b-value'].append(bval)


with open(output_file, 'w') as f:
    writer = csv.DictWriter(f,fieldnames=fieldnames)
    writer.writeheader()
    for i in range(len(cat[fieldnames[0]])):
        row = {}
        for key in cat:
            row[key] = cat[key][i]
        writer.writerow(row)

ss2n = ss.SmoothedSeismicity2Nrml(output_file, bvalue = bval, mmax = mmax)
ss2n.write_to_nrml('source_model_f.xml',mmin=mmin,upper_depth=ud,lower_depth=ld,hdd=hdd)


