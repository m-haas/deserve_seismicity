#!/usr/bin/env/python

"""
Smoothed seismicity output to NRML
"""
import argparse
import numpy as np
from decimal import Decimal
from openquake.nrmllib import models
from openquake.nrmllib.hazard.writers import SourceModelXMLWriter

DEFAULT_MSR = "WC1994"
DEFAULT_TRT = "Active Shallow Crust"
ASPECT_RATIO = 1
MIN_MAGNITUDE = 4
DEFAULT_RAKE = 0.

#Probability, Strike, Dip, Rake
NODAL_PLANES = [models.NodalPlane(Decimal("0.25"), 0., 90., DEFAULT_RAKE),
                models.NodalPlane(Decimal("0.25"), 45., 90., DEFAULT_RAKE),
                models.NodalPlane(Decimal("0.25"), 90., 90., DEFAULT_RAKE),
                models.NodalPlane(Decimal("0.25"), 135., 90., DEFAULT_RAKE)]


class SmoothedSeismicity2Nrml(object):
    """
    Converts hmtk smoothed seismicity output to OQ Point Source Model
    """
    def __init__(self, input_filename, bvalue=None, mmax=None):
        """
        """
        self.data = np.genfromtxt(input_filename, delimiter=',', skip_header=1,dtype=None)
        self.npts = np.shape(self.data)[0]
        self.bvalue = bvalue
        self.mmax = mmax

    def write_to_nrml(self, output_filename, mmin, upper_depth=0.,
            lower_depth=50., source_model_name="POINT SOURCE MODEL",
            trt=DEFAULT_TRT, msr=DEFAULT_MSR, aspect=ASPECT_RATIO, hdd=None):
        """
        Converts the smoothed seismicity data to a set of oq-nrml point
        sources and writes to the output file
        """
        writer = SourceModelXMLWriter(output_filename)
        source_model = models.SourceModel(source_model_name)
        source_model.sources = []

        print 'Building source model ...'
        for iloc, row in enumerate(self.data):
            # Geometry
            #trt =(row[18])

            geom = models.PointGeometry("POINT (%9.4f %9.4f)" % (row[0], row[1]),
                upper_depth,
                lower_depth)
            if hdd:
                src_hdd=[]
                # each tuple with probality and depth
                for d in hdd:
                    src_hdd.append(models.HypocentralDepth(d[1], d[0]))
            else:
                src_hdd = [models.HypocentralDepth(Decimal("1.0"), row[2])]

            npd = [models.NodalPlane(1,0,90,0)]
            #if row[5]==1:
            #    npd = [models.NodalPlane(row[6], row[7], row[9], row[8])]
            #elif row[5] ==2:
            #    npd = [models.NodalPlane(row[6], row[7], row[9], row[8]), models.NodalPlane(row[10], row[11], row[13], row[12])]
            #else:
            #    npd =  [models.NodalPlane(row[6], row[7], row[9], row[8]), models.NodalPlane(row[10], row[11], row[13], row[12]), models.NodalPlane(row[14], row[15], row[17], row[16])]


            source_model.sources.append(models.PointSource(
                str(iloc),
                "PNT_%s" % str(iloc),
                trt,
                geom,
                msr,
                aspect,
                self.get_mfd(iloc, row, mmin),
                npd,
                src_hdd))
        print 'done!'
        print 'Writing to file ...'
        writer.serialize(source_model)
        print 'done!'

#       def get_npd(self, row):
#               '''
#               Creates nodal plane arrays from read line
#               '''
#               if row[5]==1:
#                       return [models.NodalPlane(Decimal(row[6]), row[7], row[9], row[8])]
#               elif row[5] ==2:
#                       return [models.NodalPlane(Decimal(row[6]), row[7], row[9], row[8]),
#                               models.NodalPlane(Decimal(row[10]), row[11], row[13], row[12])]
#               else:
#                       return [models.NodalPlane(Decimal(row[6]), row[7], row[9], row[8]),
#                               models.NodalPlane(Decimal(row[10]), row[11], row[13], row[12]),
#                                       models.NodalPlane(Decimal(row[14]), row[15], row[17], row[16])]
#

    def get_mfd(self, iloc, row, mmin):
        """
        Builds the truncated Gutenberg Richter distribution for the point
        from the smoothed seismicity data
        """
        if self.bvalue is None:
            bval = row[-1]
        else:
            bval = self.bvalue
        if row[4]!=0:
            aval = row[4]
            #aval = np.log10(row[4])
        else:
            aval = -2

        if self.mmax is None:
            raise ValueError('Maximum magnitude must be specified!')
        elif isinstance(self.mmax, np.ndarray):
            return models.TGRMFD(aval, bval, mmin, self.mmax[iloc])
        else:
            return models.TGRMFD(aval, bval, mmin, self.mmax)

    def get_bvalue_mmax_from_file(self, filename):
        """
        If opting for a spatially variable b-value and/or Mmax, these
        can be specified in another csv file of format
        [Longitude, Latitude, b-value, Mmax]

        """
        new_data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        if np.shape(new_data)[0] != self.npts:
            raise ValueError('Additional data not aligned with original!')
        self.bvalue = new_data[:, 2]
        self.mmax = new_data[:, 3]


#def set_up_arg_parser():
#    """
#    Can run as executable. To do so, set up the command line parser
#    """
#    parser = argparse.ArgumentParser(
#        description='Convert Hazard Map file from Nrml to Something Readable'
#            'To run just type: python hazard_map_converter.py '
#            '--input-file=/PATH/INPUT_FILE_NAME '
#            '--output-file=/PATH/OUTPUT_FILE_NAME.xml')
#    parser.add_argument('--input-file',
#                        help='path to NRML hazard map file',
#                        default=None)
#    parser.add_argument('--output-file',
#                        help='path to output csv file',
#    parser.add_argument('--scale-rel',
#                        help='magnitude scaling relation',
#                        default="WC1994")
#    parser.add_argument('--additional-data-file',
#                        help='path to additional data file (b-value, mmax)',
#                        default=None)
#    parser.add_argument('--b-value',
#                        help='Override with fixed b-value',
#                        default=None)
#    parser.add_argument('--mmax',
#                        help='Default maximum magnitude',
#                        default=None)
#    parser.add_argument('--mmin',
#                        help='Default minimum magnitude',
#                        default=None)
#    parser.add_argument('--upper-depth',
#                        help='Upper seismogenic depth',
#                        default=0.)
#    parser.add_argument('--lower-depth',
#                        help='Lower Seismogenic depth',
#                        default=50.)
#    parser.add_argument('--aspect-ratio',
#                        help='Default aspect ratio',
#                        default=1.0)
#    parser.add_argument('--tectonic-region',
#                        help='Default tectonic region type',
#                        default="Active Shallow Crust")
#    parser.add_argument('--hypocentral-depth',
#                        help='Default tectonic region type',
#                        default=None)
#    return parser
#
#
#if __name__=="__main__":
#    # Parse inputs
#    parser = set_up_arg_parser()
#    args = parser.parse_args()
#    # Instantiate model
#    builder = SmoothedSeismicity2Nrml(
#        args.input_file,
#        args.b_value,
#        args.mmax)
#    # If additional data is specified then update bvalue and mmax
#    if args.additional_data_file:
#        builder.get_bvalue_mmax_from_file(args.additional_data_file)
#    # Build model and serialise
#    builder.write_to_nrml(
#        args.output_file,
#        args.mmin,
#        args.upper_depth,
#        args.lower_depth,
#        trt=args.tectonic_region,
#        msr=args.scale_rel,
#        aspect=args.aspect_ratio,
#        hdd=args.hypocentral_depth)
#
