"""Analytical checks for posterior integration and group ranking."""
import unittest

import numpy as np
import pandas as pd

from src.compare_shet import posterior_mean, rank_and_compare


class ShetTests(unittest.TestCase):
    def test_mean_on_linear_scale_with_unnormalized_log_density(self):
        x = np.linspace(-4, 0, 10001)
        mean, mass = posterior_mean(x, np.full_like(x, 2.0))
        expected = (1 - 1e-4) / (4 * np.log(10))
        self.assertAlmostEqual(mean, expected, places=7)
        self.assertAlmostEqual(mass, 8.0)
        self.assertNotAlmostEqual(mean, 10**x.mean(), places=2)

    def test_invalid_density_rejected(self):
        for x, density in [([0, -1], [1, 1]), ([-1, 0], [0, 0]),
                           ([-1, 0], [1, -1]), ([-1, 0], [1, np.nan])]:
            with self.assertRaises(ValueError):
                posterior_mean(x, density)

    def test_ties_missing_and_u_orientation(self):
        frame = pd.DataFrame({"st6_symbol": list("abcde"),
                              "group": ["Dog", "Dog", "Dog", "Human", "Human"],
                              "shet_posterior_mean": [0.2, 0.3, np.nan, 0.1, 0.2]})
        ranked, stats = rank_and_compare(frame)
        self.assertEqual(stats["n_dog"], 2)
        self.assertEqual(stats["U_dog"], 3.5)
        self.assertEqual(stats["probability_dog_greater"], 0.875)
        self.assertEqual(ranked.set_index("st6_symbol").loc["a", "shet_rank"], 2.5)
        self.assertTrue(np.isnan(ranked.set_index("st6_symbol").loc["c", "shet_rank"]))


if __name__ == "__main__":
    unittest.main()
