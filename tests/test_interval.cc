/*
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE IntervalCpp

#include "boost/test/unit_test.hpp"

#include "lsst/cpputils/tests.h"

#include "lsst/geom/Interval.h"
#include "ndarray/arange.h"

/*
 * Unit tests for C++-only functionality in IntervalI and IntervalD.
 *
 * See test_interval.py for remaining unit tests.
 */
namespace lsst {
namespace geom {

BOOST_AUTO_TEST_CASE(IntervalISlice) {
    auto array = ndarray::copy(ndarray::arange(5));
    auto interval = IntervalI::fromMinMax(2, 4);
    auto subarray = array[interval.getSlice()];
    BOOST_TEST(subarray.getShape() == ndarray::makeVector(3));
    BOOST_TEST(subarray[0] == 2);
    BOOST_TEST(subarray[1] == 3);
    BOOST_TEST(subarray[2] == 4);
}

BOOST_AUTO_TEST_CASE(IntervalIHash) {
    cpputils::assertValidHash<IntervalI>();
    cpputils::assertHashesEqual(IntervalI::fromMinMax(2, 5), IntervalI::fromMinSize(2, 4));
    cpputils::assertHashesEqual(IntervalI(), IntervalI::fromMinMax(1, -1));
}

BOOST_AUTO_TEST_CASE(IntervalDHash) {
    cpputils::assertValidHash<IntervalD>();
    cpputils::assertHashesEqual(IntervalD::fromMinMax(-2.0, 2.0), IntervalD::fromMinSize(-2.0, 4.0));
    cpputils::assertHashesEqual(IntervalD(), IntervalD::fromMinMax(1.0, -1.0));
}

}  // namespace geom
}  // namespace lsst
