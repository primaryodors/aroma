import unittest
from data.configutils import state_atom_constraints


class StateConfigTest(unittest.TestCase):
    def test_inactive_only(self):
        self.assertEqual(state_atom_constraints({"atomtoi": "1:CA 2:CA"}, "inactive"),
                         ["1:CA 2:CA"])

    def test_both_states(self):
        pocket = {"atomtoa": ["active target"], "atomtoi": ["inactive target"]}
        self.assertEqual(state_atom_constraints(pocket, "inactive"), ["inactive target"])
        self.assertEqual(state_atom_constraints(pocket, "active"), ["active target"])
        self.assertEqual(state_atom_constraints(pocket, "L"), ["active target"])

    def test_no_constraints(self):
        self.assertEqual(state_atom_constraints({}, "inactive"), [])


if __name__ == "__main__":
    unittest.main()
