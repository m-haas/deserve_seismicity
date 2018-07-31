mport re
from lxml import etree

GML_NS = 'http://www.opengis.net/gml'
NRML04_NS = 'http://openquake.org/xmlns/nrml/0.4'
NRML04 = "{%s}" % NRML04_NS
NRML04_ROOT_TAG = '%snrml' % NRML04
NRML04_SITE_MODEL = etree.QName(NRML04_NS,'siteModel')
NRML04_SITE = etree.QName(NRML04_NS,'site')

NSMAP = {None: NRML04_NS, "gml": GML_NS}

class WriteSiteModel0pt4(object):

    def __init__(self,vs30Type='measured',z1pt0=500,z2pt5=1.5):
        """
        :param vs30_list:
            A list of (lon,lat,vs30,vs30Type,z1pt0,z2pt5) tuples. The last three
            parameters are optional.
        """
        self.site_list = None
        self.vs30Type = vs30Type
        self.z1pt0 = z1pt0
        self.z2pt5 = z2pt5
        self._nrml_init()

    def loader_csv(self,filename):
        self.site_list = []
        fin = open(filename,'r')
        for line in fin:
            aa = re.split(',',line)
            site = ()
            for i, elem in enumerate(aa):
                if i != 3:
                    site += (float(elem),)
                else:
                    site += (elem,)
            self.site_list.append(site)

    def write(self, outfile):
        """
        Write sources to the output file
        """
        # Check this
        for point in self.site_list:
            self._add_site(point)

        fou = open(outfile,'w')
        str = etree.tostring(self.nrml, pretty_print=True,
                xml_declaration=True, encoding="UTF-8")
        fou.write(str)
        fou.close()

    def _nrml_init(self):
        """
        Create the nrml header
        """
        self.nrml = etree.Element(NRML04_ROOT_TAG, nsmap=NSMAP)
        elem0 = etree.SubElement(self.nrml, NRML04_SITE_MODEL)

    def _add_site(self, point):
        print point
        site = etree.SubElement(self.nrml.find(NRML04_SITE_MODEL), NRML04_SITE)
        site.set('lon','%.2f' % point[0])
        site.set('lat','%.2f' % point[1])
        site.set('vs30','%.2f' % point[2])
        if len(point) > 3:
            site.set('vs30Type',point[3])
        else:
            site.set('vs30Type',self.vs30Type)
        if len(point) > 4:
            site.set('z1pt0','%.2f' % point[4])
        else:
            site.set('z1pt0','%.2f' % self.z1pt0)
        if len(point) > 5:
            site.set('z2pt5','%.2f' % point[5])
        else:
            site.set('z2pt5','%.2f' % self.z2pt5)

