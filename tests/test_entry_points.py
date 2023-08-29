import unittest
import tempfile
import pathlib
import pytom.convert.coords2PL as coords2PL
from pytom.tools.parse_script_options import RequiredError
from pytom.basic.structures import ParticleList
from lxml import etree
import sys

class TestCoords2PL(unittest.TestCase):
    def setUp(self):
        self.old_argv = sys.argv.copy()
        self.own_dir = pathlib.Path(__file__).parent
        self.coord_file = self.own_dir / 'testData' / 'coords_tomogram_000_WBP.txt'
        self.pl = ParticleList()
        self.pl.loadCoordinateFile(self.coord_file)

    def tearDown(self):
        sys.argv = self.old_argv

    def test_missing_particle_list(self):
        sys.argv = ['test','-c', self.coord_file]
        with self.assertRaisesRegex(
            RequiredError, 
            r"Required flag not passed.*-p --particleList.*"
            ):
            coords2PL.entry_point()

    def test_missing_coordinate_file(self):
        sys.argv = ['test','-p', 'test.xml']
        with self.assertRaisesRegex(
            RequiredError, 
            r"Required flag not passed.*-c --coords.*"
            ):
            coords2PL.entry_point()

    def test_minimal(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = pathlib.Path(tmpdirname)
            out_file = tmpdirname / 'test.xml'
            sys.argv = ['test','-p', str(out_file.resolve()),'-c', str(self.coord_file.resolve())]
            coords2PL.entry_point()
            assert out_file.exists()
            out_pl = ParticleList()
            xml = etree.parse(out_file)
            out_pl.fromXML(xml.getroot())
            assert out_pl == self.pl

    def test_tomogram_name(self):
        tomo_name = 'TestWasHere'
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = pathlib.Path(tmpdirname)
            out_file = tmpdirname / 'test.xml'
            sys.argv = ['test','-p', str(out_file.resolve()),'-c', str(self.coord_file.resolve()),
                    '-t', tomo_name]
            coords2PL.entry_point()
            assert out_file.exists()
            out_pl = ParticleList()
            xml = etree.parse(out_file)
            out_pl.fromXML(xml.getroot())
            assert out_pl != self.pl
            for p in out_pl:
                assert p.getPickPosition().getOriginFilename() == tomo_name

