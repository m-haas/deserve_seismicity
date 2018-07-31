#######################################################
# Seismicity using EMEC catalogue
#######################################################

import numpy
import hmtk.parsers.catalogue as cp
import hmtk.plotting.seismicity.catalogue_plots as cplt
import hmtk.seismicity.declusterer as dc
import os

####################################################################
# Read in harmonized catalogue created with own catalogue.py routine
####################################################################

#catalogue csv with hmtk-style columns
filename = 'harmonized_catalogue.csv'
parser = cp.CsvCatalogueParser(filename)
catalogue = parser.read_file()

##plot depth distribution
depth_bin=5.0
##depth distribution
dd = catalogue.get_depth_distribution(numpy.arange(0.,max(catalogue.data['depth'])+depth_bin,depth_bin))
##largest depth are values set for unknown depth:999
#unknown = dd[-1]
#cplt.plot_depth_histogram(catalogue,depth_bin,filename="/home/ubuntu/deserve/regional/depth_distribution.eps",filetype="eps")

########################
#decluster the catalogue
########################
dc_config = {'time_distance_window': dc.distance_time_windows.GruenthalWindow(),#UhrhammerWindow(),#dc.distance_time_windows.GruenthalWindow(),
	     'fs_time_prop': .5}

declustering = dc.dec_gardner_knopoff.GardnerKnopoffType1()

cluster_index,cluster_flag = declustering.decluster(catalogue,dc_config)
#store if a event belongs to a cluster (cluster_flag) and if yes which cluster (cluster_index)
catalogue.data['cluster_index']=cluster_index
catalogue.data['cluster_flag']=cluster_flag

#write the clustered catalogue to a file with the cluster indices (fur further analysis)
import hmtk.parsers.catalogue.csv_catalogue_parser as cp2

filename = 'cluster_catalogue.csv'
if os.path.exists(filename):
    os.remove(filename)
cw = cp2.CsvCatalogueWriter(filename)
cw.OUTPUT_LIST.append('cluster_index')
cw.write_file(catalogue)

#write the declustered catalogue to a file
mainshock_flag = cluster_flag == 0
catalogue.purge_catalogue(mainshock_flag)
filename = 'declustered_catalogue.csv'
if os.path.exists(filename):
    os.remove(filename)
cw = cp2.CsvCatalogueWriter('declustered_catalogue.csv')
cw.write_file(catalogue)

################
#Completeness
#################
comp_config={'magnitude_bin': 0.5 , 'time_bin': 10. ,'increment_lock': True }
from hmtk.seismicity.completeness.comp_stepp_1971 import Stepp1971
completeness_algorithm = Stepp1971()
completeness_table = completeness_algorithm.completeness(catalogue, comp_config )
print completeness_table

def compl_plot(catalogue,completeness_table,zone_id,show=False):
    '''
    Creates a completeness plot for the catalogue/subcatalogue
    '''
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    from hmtk.seismicity.utils import decimal_time

    dyear = decimal_time(
        catalogue.data['year'],
        catalogue.data['month'],
        catalogue.data['day'],
        catalogue.data['hour'],
        catalogue.data['minute'],
        catalogue.data['second'])
    mag = catalogue.data['magnitude']
    # add completeness line
    yc = [i[0] for i in completeness_table]
    mc = [i[1] for i in completeness_table]
    yc.sort()
    mc.sort(reverse=True)
    comp = list()
    year = numpy.linspace(yc[0],dyear[-1],10000)
    for y in year:
        tmp = 0
        for i in range(len(yc)):
            if y > yc[i]:
                tmp = mc[i]
        comp.append(tmp)

    #uncomment to plot from completeness year
    #ll = (yc[0],mc[-1])
    ll = (1900,mc[-1])
    x1 = [y for y in dyear if y >= ll[0]]
    less = len(mag)-len(x1)
    x2 = year
    y1 = mag[less:]
    y2 = comp
    event = plt.plot(x1,y1,'ob',label='Events')
    compl = plt.plot(x2,y2)
    plt.setp(compl,color = 'r',linewidth=2.0)
    #plt.title(str(zone_id))
    plt.xlim(min(x1),max(x1))
    plt.xlabel('Year')
    plt.ylabel('Mw')
    rline = mlines.Line2D([],[],color='red',label='Completeness')
    plt.legend(handles=[rline],loc='lower left')
    if show:
        plt.show()
    else:
        plt.savefig('completeness_plot.eps')
    plt.clf()


########################
# Smoothed seismicity
########################

########################
#1) estimate b-values
########################

#estimate b-value for whole area using different approaches
mle_config={'magnitude_interval': 0.3 ,
            'Average Type': 'Weighted',
            'reference_magnitude': None }


weichert_config={'magnitude_interval': 0.3 ,
                 'reference_magnitude': None ,
                 # The remaining parameters are
                 'bvalue': 0.8 ,
                 'itstab': 1E-5 ,
                 'maxiter': 1000}

kijko_smit_config = mle_config

import hmtk.seismicity.occurrence as occ

class parameter_estimate():
    '''
    Object to keep an estimate
    '''
    def __init__(self,result):
        self.bval = result[0]
        self.sigmab = result[1]
        self.aval = result[2]
        self.sigmaa = result[3]

# function to estimate b-value for a given shapefile
import osgeo.ogr
import openquake.hazardlib.geo.point
import openquake.hazardlib.geo.polygon
import hmtk.seismicity.selector


#def b_for_subset(layer,catalogue):
#    '''
#    Given shapefile layer and a catalogue computes GR a,b for subset following 3 approaches
#    Returns dictionary whith values for each polygon in layer and each method
#    '''
#    #Initialize dictionary to hold estimates for subcatalogues
#    subcatalogues = {
#            'zone_id':list(),
#            'cat':list(),
#            'comp':list(),
#            'bmlh_aki':list(),
#            'weichert':list(),
#            'kijko_smit':list()
#            }
#
#    #get points for all polygon features
#    for feature in layer:
#    	#subcatalogue selector
#    	selector = hmtk.seismicity.selector.CatalogueSelector(catalogue,create_copy=True)
#    	#get the geometries from the shapefile
#    	geom = feature.GetGeometryRef()
#    	for ring in geom:
#    		points = []
#    		for point in range(ring.GetPointCount()):
#    			#lon/lat of the point (node of zone)
#    			lon = ring.GetPoint(point)[0]
#    			lat = ring.GetPoint(point)[1]
#    			points.append(openquake.hazardlib.geo.point.Point(lon,lat))
#    		#add subcatalogue for the zone
#    		#id/oq-polygon of the zone
#    		zone_id = feature.GetField('ID')
#    		zone_poly = openquake.hazardlib.geo.polygon.Polygon(points)
#    		subcatalogue=selector.within_polygon(zone_poly,distance=None)
#    		subcatalogues['zone_id'].append(zone_id)
#    		subcatalogues['cat'].append(selector.within_polygon(zone_poly,distance=None))
#                #get completeness table for subcatalogue
#                completeness_table = completeness_algorithm.completeness(subcatalogue, comp_config )
#                subcatalogues['comp'].append(completeness_table)
#                compl_plot(subcatalogue,completeness_table,zone_id)
#                #estimate the a,b values for all three methods
#                #Modified Aki
#                recurrence = occ.b_maximum_likelihood.BMaxLikelihood()
#                subcatalogues['bmlh_aki'].append(parameter_estimate(recurrence.calculate(catalogue,mle_config,completeness = completeness_table)))
#                #Weichert
#                recurrence = occ.weichert.Weichert()
#                subcatalogues['weichert'].append(parameter_estimate(recurrence.calculate(catalogue,weichert_config,completeness = completeness_table)))
#                #Kijko Smit
#                #recurrence = occ.kijko_smit.KijkoSmit()
#                #subcatalogues['kijko_smit'].append(parameter_estimate(recurrence.calculate(catalogue,kijko_smit_config,completeness = completeness_table)))
#    return subcatalogues

#estimate for whole region
#Get zones from a shapefile
#filename='DESERVE_catalogue_region.shp'
#shapefile = osgeo.ogr.Open(filename)
#layer = shapefile.GetLayer()
#whole = b_for_subset(layer,catalogue)

#estimate b_value for catalogue
selector = hmtk.seismicity.selector.CatalogueSelector(catalogue,create_copy=True)
import datetime
st = datetime.datetime(1900,1,1,0,0,0,0)
#catalogue = selector.within_time_period(start_time=st)
completeness_table = completeness_algorithm.completeness(catalogue, comp_config)

compl_plot(catalogue,completeness_table,1,show=False)

#reduce catalogue for b-value estimate
catalogue_reduced = selector.within_time_period(start_time=st)
selector = hmtk.seismicity.selector.CatalogueSelector(catalogue_reduced,create_copy=False)
catalogue_reduced = selector.within_magnitude_range(4,8)
recurrence = occ.b_maximum_likelihood.BMaxLikelihood()
aki = parameter_estimate(recurrence.calculate(catalogue_reduced,mle_config,completeness = completeness_table))
recurrence = occ.weichert.Weichert()
weichert = parameter_estimate(recurrence.calculate(catalogue_reduced,weichert_config,completeness = completeness_table))
#recurrence = occ.kijko_smit.KijkoSmit()
KS = parameter_estimate(recurrence.calculate(catalogue_reduced,kijko_smit_config,completeness = completeness_table))


#test = b_for_subset(layer,test_cat)
#
##write final catalogue to csv
#filename = 'DESERVE_catalogue.csv'
#if os.path.exists(filename):
#    os.remove(filename)
#cw = cp2.CsvCatalogueWriter(filename)
#cw.write_file(whole['cat'][0])
#
##completeness plot (GEM style)
#from hmtk.plotting.seismicity.catalogue_plots import plot_magnitude_time_density
#plot_magnitude_time_density(whole['cat'][0],0.1,1.0)
#
###estimate for quadrants
###Get zones from a shapefile
##filename='psha_quadrants.shp'
##shapefile = osgeo.ogr.Open(filename)
##layer = shapefile.GetLayer()
##quadrants = b_for_subset(layer,catalogue)
#
#
##pretty print the results
#def pretty_print(dictionary):
#    print 'zone_id','aki_b','aki_sigma','weichert_b','weichert_sigma'#,'kijko_smit_b','kijko_smit_sigma'
#    for q in range(len(dictionary['zone_id'])):
#        print dictionary['zone_id'][q],dictionary['bmlh_aki'][q].bval,dictionary['bmlh_aki'][q].sigmab,dictionary['weichert'][q].bval,dictionary['weichert'][q].sigmab#,dictionary['kijko_smit'][q].bval,dictionary['kijko_smit'][q].sigmab
#
#create a GR plot
#def GR_plot(catalogue,row_comp_table,binwidth,estimate):
def GR_plot(catalogue,binwidth,estimate,fname='GR_plot',show=False):
    '''
    Plots actual cumulative magnitude histogram against GR frequency law
    '''
    import scipy.stats
    import math
    import matplotlib.pyplot as plt
    mags = catalogue.data['magnitude']
    #filter for completeness
    #yc = row_comp_table[0]
    #years  = [y for y in catalogue.data['year'] if y >= yc]
    #less = len(mags)-len(years)
    #mags = mags[less:]
    #round magnitudes to binwidth (lower bounds of the bins)
    mags = [int(m*1/binwidth)/(1./binwidth) for m in mags]
    nrbins = int(((max(mags)-min(mags))/binwidth)+1)
    #observation period
    period = max(catalogue.data['year'])-min(catalogue.data['year'])+1
    observed = {'bins':[],'counts':[]}
    estimates = {'bins':[],'counts':[]}
    for i in range(nrbins):
        #binstart
        val = min(mags) + i*binwidth
        #store binlefts and observed counts
        observed['bins'].append(val)
        observed['counts'].append(len([m for m in mags if m >= val]))
        #store estimated values
        estimates['bins'].append(val+binwidth/2)
        sim = estimate.aval-estimate.bval*(val+binwidth/2)
        sim = period*10**sim
        estimates['counts'].append(sim)

    #cf = scipy.stats.cumfreq(mags,defaultreallimits=(min(mags),max(mags)+binwidth),numbins=nrbins)
    #get result in GR fashion (nr(m) >= m_bin)
    #GR_hist = []
    #for i in range(len(cf[0])):
    #    if i==0:
    #        GR_hist.append(cf[0][-1])
    #    else:
    #        GR_hist.append(cf[0][-1] - cf[0][i-1])
    #x = numpy.arange(min(mags),max(mags)+binwidth,binwidth)
    #Plot cumulative count
    GR_hist = plt.bar(left=observed['bins'],height=observed['counts'],width=binwidth,label='Observed')
    GR_est = plt.plot(estimates['bins'],estimates['counts'],'-r',label='Estimated')
    plt.setp(GR_est,color = 'r',linewidth=2.0)
    #plt.title(str(zone_id))
    plt.xlabel('Mw')
    plt.ylabel('N')
    plt.yscale('log')
    plt.legend(handles=[GR_hist,GR_est],labels=['Observed','Estimated'],loc='best')
    text = 'a={:3.2f}, b={:3.2f}, T={}a'.format(estimate.aval,estimate.bval,period)
    plt.annotate(text,xy=(0.05,0.05),xycoords='axes fraction',backgroundcolor='white')
    if show:
        plt.show()
    else:
        plt.savefig(fname+'.eps')
    plt.clf()

#Create GR plot
GR_plot(catalogue_reduced,.5,weichert,'GR_weichert',show=False)

## Imports the smoothed seismicity algorithm
#from hmtk.seismicity.smoothing.smoothed_seismicity import SmoothedSeismicity
## Imports the Kernel function
#from hmtk.seismicity.smoothing.kernels.isotropic_gaussian import IsotropicGaussian
## Assuming a homogeneous b - value
#bval = whole['weichert'][0].bval
#bsigma = whole['weichert'][0].sigmab
#import scipy.stats
#b = scipy.stats.norm(bval,bsigma)
#grid_spacing = 0.2
#
##create seismicity models for 10 samples from the bvalue distribution percentiles
#smooth_seis = SmoothedSeismicity ( grid_spacing ,use_3d = False ,bvalue = whole['weichert'][0].bval)
## Set up config ( e.g. 50 km band width , up to 3 bandwidths )
#config = {'Length_Limit': 3. ,'BandWidth': 50. ,'increment': False}
## Run the analysis !
#output_data = smooth_seis.run_analysis (whole['cat'][0] ,config ,whole['comp'][0] ,smoothing_kernel = IsotropicGaussian())
## To write the resulting data to a csv file
#smooth_seis.write_to_csv('smoothed_seismicity.csv')
#
#####################################
##depth distribution
#####################################
##binwidth=5.0
#depths = whole['cat'][0].data['depth']
###round depths to binwidth
#depths = [d for d in depths if d < 999]
##depths = [int(d*1/binwidth)/(1./binwidth) for d in depths]
##nrbins = int(((max(depths)-min(depths))/binwidth))
#nrbins=20
#dd = scipy.stats.histogram(depths,defaultlimits=(min(depths),max(depths)),numbins=nrbins)
#dpdf = dd[0]/sum(dd[0])
#depthbins = [i*dd[2]+dd[1] for i in range(len(dpdf))]
#hdd = []
##round the probabilities and include only > 0.01
#psum=[round(p*100)/100 for p in dpdf if p >= 0.01]
#psum=sum(psum)
#for i in range(len(dpdf)):
#    # include only depthbin where probability >= 0.01
#    if dpdf[i] >= 0.01:
#        #midpoints of the depthbins and probability
#        hdd.append((depthbins[i]+dd[2]/2,round(dpdf[i]*100)/100/psum))
#
##round lower seismogenic depth to next largest 10km step
#lsd = round(max(hdd)[0]/10)*10+10
#
###########################################
## Maximum magnitude
###########################################
#mmax_config={ 'input_mmax': None ,
#        'input_mmax_uncertainty': 0.2 ,
#        'b-value': b.mean(),
#        'sigma-b': b.std(),
#        'input_mmin': 4.0 ,
#        'tolerance': 1.0E-5 , # Default
#        'maximum_iterations': 1000}  # Defaults
#
#from hmtk.seismicity.max_magnitude.kijko_sellevol_bayes import KijkoSellevolBayes
#
#mmax_estimator = KijkoSellevolBayes()
#mmax,mmax_uncertainty = mmax_estimator.get_mmax(whole['cat'][0],mmax_config)
#mmax = round(mmax*10)/10
#
############################################################
## Convert model to nrml
############################################################
#
#import smoothed_seismicity2nrml as ss
#
#ss2n = ss.SmoothedSeismicity2Nrml('smoothed_seismicity.csv', bvalue = whole['weichert'][0].bval, mmax = mmax)
#ss2n.write_to_nrml('source_model.xml',mmin=4.0,upper_depth=1.,lower_depth=lsd,hdd=hdd)
#

