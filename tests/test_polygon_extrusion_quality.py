import unittest

from scripts.evaluate_polygon_extrusion_quality import summarize_model


def square_extrusion(height=1, top_shift=0):
    controls = []
    for u in range(2):
        for v in range(2):
            for w in range(2):
                controls.append([u + w * top_shift, v, w * height, 1])
    return {"blocks": [{"id": "square", "degrees": [1, 1, 1], "controls": controls}]}


class PolygonExtrusionQualityTests(unittest.TestCase):
    def test_thin_square_separates_planar_and_thickness_quality(self):
        result = summarize_model(square_extrusion(height=0.1))
        self.assertAlmostEqual(result["planar_spectral_kappa"], 1.0)
        self.assertAlmostEqual(result["minimum_corner_sine"], 1.0)
        self.assertAlmostEqual(result["thickness_direction_aspect"], 10.0)
        self.assertAlmostEqual(result["max_3d_spectral_kappa"], 10.0)

    def test_skew_translation_does_not_look_like_low_normal_thickness(self):
        result = summarize_model(square_extrusion(top_shift=4))
        self.assertAlmostEqual(result["thickness_direction_aspect"], 1.0)
        self.assertGreater(result["max_3d_spectral_kappa"], 1.0)

    def test_rejects_nonuniform_weights(self):
        invalid = square_extrusion()
        invalid["blocks"][0]["controls"][0][3] = 2
        with self.assertRaisesRegex(ValueError, "common control weight"):
            summarize_model(invalid)


if __name__ == "__main__":
    unittest.main()
