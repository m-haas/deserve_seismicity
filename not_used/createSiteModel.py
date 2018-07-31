import sys

from sitemodelwriter_0pt4 import WriteSiteModel0pt4

def main(argv):

    siteModelWriter = WriteSiteModel0pt4(vs30Type='measured',z1pt0=500,z2pt5=2.0)
    siteModelWriter.loader_csv(argv[1])
    siteModelWriter.write('siteModel.xml')

if __name__ == '__main__':
    """
    usage: python siteModel_creator.py <xyz_filename>
    """
    main(sys.argv)

