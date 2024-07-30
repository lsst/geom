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
#define BOOST_TEST_MODULE AngleCpp
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop

#include "lsst/cpputils/tests.h"

#include "lsst/geom/Angle.h"

namespace lsst {
namespace geom {

BOOST_AUTO_TEST_CASE(Hash) {
    cpputils::assertValidHash<Angle>();
    cpputils::assertHashesEqual(0.0 * radians, 0.0 * degrees);

    cpputils::assertValidHash<AngleUnit>();
    cpputils::assertHashesEqual(degrees, AngleUnit(PI / 180.0));
}

}  // namespace geom
}  // namespace lsst
