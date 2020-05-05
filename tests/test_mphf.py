from unittest import TestCase

from src.MPHF import MPHF

class TestMPHF(TestCase):
    def test___add_three_objects_with_copies(self):
        mphf = MPHF()

        mphf.add_object(1)
        mphf.add_object(1)
        mphf.add_object("str")
        mphf.add_object(1)
        mphf.add_object(1.5)
        mphf.add_object("str")
        mphf.add_object(1.5)

        self.assertEqual(3, mphf.get_number_of_objects())

        self.assertEqual(0, mphf.get_id(1))
        self.assertEqual(1, mphf.get_id("str"))
        self.assertEqual(2, mphf.get_id(1.5))

        self.assertEqual(1, mphf.get_object(0))
        self.assertEqual("str", mphf.get_object(1))
        self.assertEqual(1.5, mphf.get_object(2))
