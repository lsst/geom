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

#include "nanobind/nanobind.h"
#include "nanobind/stl/tuple.h"
#include "nanobind/stl/vector.h"
#include "nanobind/ndarray.h"
#include "lsst/sphgeom/python.h"

#include <cmath>
#include <memory>

#include "lsst/cpputils/python.h"
#include "lsst/sphgeom/UnitVector3d.h"
#include "lsst/sphgeom/Vector3d.h"
#include "lsst/sphgeom/LonLat.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/SpherePoint.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace lsst {
namespace geom {

namespace {

	
double toUnitX(double longitude, double latitude) {
    return std::cos(longitude) * std::cos(latitude);
}

double toUnitY(double longitude, double latitude) {
    return std::sin(longitude) * std::cos(latitude);
}

double toUnitZ(double longitude, double latitude) {
    return std::sin(latitude);
}

} // anonymous


using PySpherePoint = nb::class_<SpherePoint>;

void wrapSpherePoint(cpputils::python::WrapperCollection & wrappers) {
    wrappers.addSignatureDependency("lsst.sphgeom");
    wrappers.wrapType(
        PySpherePoint(wrappers.module, "SpherePoint"),
        [](auto & mod, auto & cls) mutable {
            /* Constructors */
            cls.def(nb::init<>());
            cls.def(nb::init<Angle const &, Angle const &>(), "longitude"_a, "latitude"_a);
            cls.def(nb::init<double, double, AngleUnit>(), "longitude"_a, "latitude"_a, "units"_a);
            cls.def(nb::init<sphgeom::Vector3d const &>(), "vector"_a);
            cls.def(nb::init<sphgeom::UnitVector3d const &>(), "unitVector"_a);
            cls.def(nb::init<sphgeom::LonLat const &>(), "lonLat"_a);
            cls.def(nb::init<SpherePoint const &>(), "other"_a);
            nb::implicitly_convertible<SpherePoint, sphgeom::LonLat>();
            nb::implicitly_convertible<sphgeom::LonLat, SpherePoint>();

            /* Operators */
            cls.def("__getitem__",
                    [](SpherePoint const &self, std::ptrdiff_t i) {
                        return self[cpputils::python::cppIndex(2, i)];
                    });
            cls.def("__eq__", &SpherePoint::operator==, nb::is_operator());
            cls.def("__ne__", &SpherePoint::operator!=, nb::is_operator());

            /* Members */
            cls.def("getLongitude", &SpherePoint::getLongitude);
            cls.def("getLatitude", &SpherePoint::getLatitude);
            cls.def("getRa", &SpherePoint::getRa);
            cls.def("getDec", &SpherePoint::getDec);
            cls.def("getVector", &SpherePoint::getVector);
            cls.def("getPosition", &SpherePoint::getPosition, "units"_a);
            cls.def("atPole", &SpherePoint::atPole);
            cls.def("isFinite", &SpherePoint::isFinite);
            cls.def("bearingTo", &SpherePoint::bearingTo, "other"_a);
            cls.def("separation", &SpherePoint::separation, "other"_a);
            cls.def("rotated", &SpherePoint::rotated, "axis"_a, "amount"_a);
            cls.def("offset", &SpherePoint::offset, "bearing"_a, "amount"_a);
            cls.def("getTangentPlaneOffset", &SpherePoint::getTangentPlaneOffset, "other"_a);
            cpputils::python::addOutputOp(cls, "__str__");
            cls.def("__len__", [](SpherePoint const &) { return 2; });
            cls.def("__reduce__", [cls](SpherePoint const &self) {
                return nb::make_tuple(cls, nb::make_tuple(nb::cast(self.getLongitude()),
                                                          nb::cast(self.getLatitude())));
            });


            /* Module level */
            mod.def("averageSpherePoint", averageSpherePoint);

            // Used only by pure-Python extension toUnitXYZ in _SpherePoint.py.
            mod.def("_toUnitX", nb::vectorize(&toUnitX), "longitude"_a, "latitude"_a);
            mod.def("_toUnitY", nb::vectorize(&toUnitY), "longitude"_a, "latitude"_a);
            mod.def("_toUnitZ", nb::vectorize(&toUnitZ), "longitude"_a, "latitude"_a);
        }
    );
}

}  // namespace geom
}  // namespace lsst
