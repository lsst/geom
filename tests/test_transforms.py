#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import unittest

import numpy as np

import lsst.utils.tests
from lsst.geom import LinearTransform, AffineTransform, Point2D, Extent2D


class TransformTests:

    def testVectorization(self):
        x1 = np.random.randn(4, 3)
        y1 = np.random.randn(4, 3)
        x2, y2 = self.transform(x1, y1)
        self.assertEqual(x1.shape, x2.shape)
        self.assertEqual(x2.shape, y2.shape)
        for i in range(4):
            for j in range(3):
                x3, y3 = self.transform(Point2D(x1[i, j], y1[i, j]))
                self.assertAlmostEqual(x3, x2[i, j], 15)
                self.assertAlmostEqual(y3, y2[i, j], 15)


class LinearTransformTestCase(lsst.utils.tests.TestCase, TransformTests):

    def setUp(self):
        np.random.seed(5)
        matrix = np.random.randn(2, 2)
        self.transform = LinearTransform(matrix)


class AffineTransformTestCase(lsst.utils.tests.TestCase, TransformTests):

    def setUp(self):
        np.random.seed(6)
        matrix = np.random.randn(2, 2)
        vector = np.random.randn(2)
        self.transform = AffineTransform(LinearTransform(matrix), Extent2D(vector))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
