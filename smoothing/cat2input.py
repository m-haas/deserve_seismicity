#converts an input hmtk catalogue csv to a ascii year,lon,lat,z,mag
import csv

input_file = 'DESERVE_catalogue.csv'
output_file = 'input_cat.csv'

mmin=4.0

cat = {'year':[] ,
       'lat': [] ,
       'lon': [] ,
       'z': []   ,
       'mag': []
}

with open(input_file,'r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    #skip header
    for i in range(1):
        next(reader,None)
    for row in reader:
    #Store data
        #remove whitespace
        row = [x.strip() for x in row]
        #store only if larger equal mmin
        if round(float(row[16])*10)/10 >= mmin:
            cat['year'].append(int(row[2]))
            cat['lon'].append(float(row[9]))
            cat['lat'].append(float(row[10]))
            # add default depth for unknown depth
            if int(float(row[14]))!=999:
                cat['z'].append(int(float(row[14])))
            else:
                cat['z'].append(15)
            cat['mag'].append(round(float(row[16])*10)/10)


with open(output_file, 'w') as f:
    fieldnames = list(['year','lat','lon','z','mag'])
    writer = csv.DictWriter(f,fieldnames=fieldnames)
    #writer.writeheader()
    for i in range(len(cat[fieldnames[0]])):
        row = {}
        for key in cat:
            row[key] = cat[key][i]
        writer.writerow(row)
