
import xml.etree.ElementTree as ET
from filter import SingleTiltWedge

class BaseClass(object):
    """docstring for BaseClass"""
    def __init__(self):
        super(BaseClass, self).__init__()


    def __str__(self):
        return self.toStr()


    def fromStr(self, text):
        tree = ET.fromstring(text)
        self.fromXMLObj(tree)


    def toStr(self):
        return ET.tostring(self.toXMLObj(), encoding="UTF-8")


    def fromXMLFile(self, filename):
        with open(filename, 'rb') as f:
            text = f.read()

        self.fromStr(text)


    def toXMLFile(self, filename):
        with open(filename, 'wb') as f:
            f.write(self.toStr())


    def fromXMLObj(self, xml_obj):
        raise NotImplementedError('This method needs to be implemented by the child class!')


    def toXMLObj(self):
        raise NotImplementedError('This method needs to be implemented by the child class!')

    def copy(self):
        import copy
        return copy.deepcopy(self)
        

class Particle(BaseClass):
    def __init__(self, filename=None):
        super(Particle, self).__init__()

        self.filename = filename # particle filename
        
        self.rotation_phi = 0 # particle orientation
        self.rotation_psi = 0
        self.rotation_the = 0

        self.shift_x = 0 # particle shift
        self.shift_y = 0
        self.shift_z = 0

        self.pick_pos_x = 0 # particle pick position
        self.pick_pos_y = 0
        self.pick_pos_z = 0

        self.score = 0 # particle score

        self.class_label = 0 # particle class label

        self.wedge = SingleTiltWedge() # wedge


    def getTransformedVolume(self):
        from pytom.agnostic.io import read
        from pytom.agnostic.transform import translate3d, rotate3d

        v = read(self.filename)
        v2 = translate3d(v, -self.shift_x, -self.shift_y, -self.shift_z)
        v3 = rotate3d(v2, -self.rotation_psi, -self.rotation_phi, -self.rotation_the)

        return v3


    def fromXMLObj(self, xml_obj):
        self.filename = xml_obj.attrib['Filename']
        for child in xml_obj:
            if child.tag == 'Rotation':
                self.rotation_phi = float(child.get('Z1'))
                self.rotation_psi = float(child.get('Z2'))
                self.rotation_the = float(child.get('X'))
            elif child.tag == 'Shift':
                self.shift_x = float(child.get('X'))
                self.shift_y = float(child.get('Y'))
                self.shift_z = float(child.get('Z'))
            elif child.tag == 'PickPosition':
                self.pick_pos_x = float(child.get('X'))
                self.pick_pos_y = float(child.get('Y'))
                self.pick_pos_z = float(child.get('Z'))
            elif child.tag == 'Score':
                self.score = float(child.get('Value'))
            elif child.tag == 'Class':
                self.class_label = child.get('Class')
            elif child.tag == 'SingleTiltWedge':
                self.wedge = SingleTiltWedge(-90+float(child.get('Angle1')), 90-float(child.get('Angle2')))
            else:
                pass


    def toXMLObj(self):
        xml_obj = ET.Element('Particle', attrib={'Filename': self.filename})
        ET.SubElement(xml_obj, 'Rotation', attrib={'Z1': str(self.rotation_phi), 'Z2': str(self.rotation_psi), 'X': str(self.rotation_the)})
        ET.SubElement(xml_obj, 'Shift', attrib={'X': str(self.shift_x), 'Y': str(self.shift_y), 'Z': str(self.shift_z)})
        ET.SubElement(xml_obj, 'PickPosition', attrib={'X': str(self.pick_pos_x), 'Y': str(self.pick_pos_y), 'Z': str(self.pick_pos_z)})
        ET.SubElement(xml_obj, 'Score', attrib={'Value': str(self.score)})
        ET.SubElement(xml_obj, 'Class', attrib={'Class': str(self.class_label)})

        return xml_obj


class ParticleList(BaseClass):
    def __init__(self, filename=None):
        super(ParticleList, self).__init__()
        self._list = []

        if filename:
            self.fromXMLFile(filename)


    def average(self):
        v_sum = None
        wedge_sum = None

        for p in self._list:
            v = p.getTransformedVolume()
            # apply wedge to the volume itself
            v = p.wedge.apply(v, [-p.rotation_psi, -p.rotation_phi, -p.rotation_the])
            w = p.wedge.returnWedgeVolume(v.shape, [-p.rotation_psi, -p.rotation_phi, -p.rotation_the])
            if v_sum is None:
                v_sum = v
                wedge_sum = w
            else:
                v_sum += v
                wedge_sum += w

        # weighting
        wedge_sum[wedge_sum < 1] = 0 # prevent boosting by division
        from pytom.agnostic.transform import rfft, irfft, ifftshift, fourier_full2reduced
        fv_sum = rfft(v_sum)
        wedge_sum = fourier_full2reduced(ifftshift(wedge_sum))
        res = fv_sum/wedge_sum
        res[wedge_sum==0] = 0
        avg = irfft(res, v_sum.shape)

        return avg


    def fromXMLObj(self, xml_obj):
        for child in xml_obj:
            p = Particle()
            p.fromXMLObj(child)
            self._list.append(p)


    def toXMLObj(self):
        xml_obj = ET.Element('ParticleList')
        for p in self._list:
            xml_obj.append(p.toXMLObj())

        return xml_obj

    def __len__(self):
        return len(self._list)

    def __getitem__(self, key):
        if key.__class__ == int:
            return self._list[key]
        elif key.__class__ == slice:
            res = ParticleList()
            for i in xrange(key.start, key.stop, key.step):
                res.append(self._list[i])
            return res
        else:
            raise TypeError()

    def __setitem__(self, key, value):
        if key.__class__ == int:
            if value.__class__ != Particle:
                raise TypeError()
            else:
                self._list[key] = value
        else:
            raise TypeError()

    def __add__(self, pl):
        if pl.__class__ != ParticleList:
            raise TypeError()
        else:
            self._list.extend(pl._list)


