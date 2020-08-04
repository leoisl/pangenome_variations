from unittest import TestCase
from src.VarifierDataframe import VarifierDataframe

info="LENGTH_QRY=4725092;LENGTH_REF=5155066;QNAME=1_pilon_pilon_pilon_pilon_pilon_pilon;QSTART=273255;QSTRAND=+;"
class TestVarifierDataframe(TestCase):
    def test___parse_field_from_info___LENGTH_QRY(self):
        expected = 4725092
        actual = VarifierDataframe.parse_field_from_info(info, "LENGTH_QRY", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___LENGTH_REF(self):
        expected = 5155066
        actual = VarifierDataframe.parse_field_from_info(info, "LENGTH_REF", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QNAME(self):
        expected = "1_pilon_pilon_pilon_pilon_pilon_pilon"
        actual = VarifierDataframe.parse_field_from_info(info, "QNAME")
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QSTART(self):
        expected = 273255
        actual = VarifierDataframe.parse_field_from_info(info, "QSTART", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QSTRAND(self):
        expected = "+"
        actual = VarifierDataframe.parse_field_from_info(info, "QSTRAND")
        self.assertEqual(expected, actual)