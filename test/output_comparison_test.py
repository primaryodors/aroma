import unittest
from compare_output import equivalent


class OutputComparisonTest(unittest.TestCase):
    def test_printed_rounding_is_accepted(self):
        self.assertTrue(equivalent("Initial eclipses: 2.84655\n", "Initial eclipses: 2.84656\n"))
        self.assertTrue(equivalent("energy: -0.0666979", "energy: -0.0666978"))

    def test_material_energy_change_is_rejected(self):
        self.assertFalse(equivalent("energy: 1.0000", "energy: 1.0100"))

    def test_counts_and_labels_are_exact(self):
        self.assertFalse(equivalent("2 atoms", "3 atoms"))
        self.assertFalse(equivalent("His155: 1.2345", "His156: 1.2345"))
        self.assertFalse(equivalent("energy: 1.2345", "distance: 1.2345"))


if __name__ == "__main__":
    unittest.main()
