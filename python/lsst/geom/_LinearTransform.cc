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
#include "nanobind/eigen/dense.h"
#include "nanobind/stl/pair.h"
#include "nanobind/stl/tuple.h"
#include "nanobind/stl/string.h"
#include "nanobind/ndarray.h"
#include "lsst/sphgeom/python.h"

#include "ndarray/nanobind.h"

#include "lsst/cpputils/python.h"

#include "lsst/geom/Extent.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/LinearTransform.h"

namespace nb = nanobind;
using namespace nanobind::literals;

namespace lsst {
namespace geom {

void wrapLinearTransform(cpputils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        nb::class_<LinearTransform>(wrappers.module, "LinearTransform"),
        [](auto & mod, auto & cls) {

            // Parameters enum is really only used as integer constants.
            cls.attr("XX") = nb::cast(int(LinearTransform::Parameters::XX));
            cls.attr("YX") = nb::cast(int(LinearTransform::Parameters::YX));
            cls.attr("XY") = nb::cast(int(LinearTransform::Parameters::XY));
            cls.attr("YY") = nb::cast(int(LinearTransform::Parameters::YY));

            /* Constructors */
            cls.def(nb::init<>());
            cls.def(nb::init<LinearTransform::Matrix const &>(), "matrix"_a);

            /* Operators */
            cls.def("__call__",
                    nb::overload_cast<Point2D const &>(&LinearTransform::operator(), nb::const_));
            cls.def("__call__",
                    nb::overload_cast<Extent2D const &>(&LinearTransform::operator(), nb::const_));
            cls.def("__call__",
                    // We use nanobind's wrappers for the Python C API to
                    // delegate to other wrapped methods because:
                    //  - defining this in pure Python is tricky because it's
                    //    an overload, not a standalone method;
                    //  - we'd rather not add a new pure-Python file just for
                    //    this;
                    //  - using nb::vectorize internal to the method would
                    //    involve defining a new internal callable every time
                    //    this method is called.
                    // The other viable alternative would be to define
                    // applyX and applyY as Python callables with nb::vectorize
                    // outside the lambda as C++ local variables, and then
                    // capture them by value in the lambda.  This just seems
                    // slightly cleaner, as it's closer to how one would
                    // implement this in pure Python, if it wasn't an overload.
                    [](nb::object self, nb::object x, nb::object y) {
                        return nb::make_tuple(self.attr("applyX")(x, y),
                                              self.attr("applyY")(x, y));
                    },
                    "x"_a, "y"_a);
            cls.def("__getitem__",
                    [](LinearTransform const &self, int i) { return self[cpputils::python::cppIndex(4, i)]; });
            cls.def("__getitem__", [](LinearTransform const &self, std::pair<int, int> i) {
                auto row = cpputils::python::cppIndex(2, i.first);
                auto col = cpputils::python::cppIndex(2, i.second);
                return self.getMatrix()(row, col);
            });
            cls.def("__mul__", &LinearTransform::operator*, nb::is_operator());
            cls.def("__add__", &LinearTransform::operator+, nb::is_operator());
            cls.def("__sub__", &LinearTransform::operator-, nb::is_operator());
            cls.def("__iadd__", &LinearTransform::operator+=);
            cls.def("__isub__", &LinearTransform::operator-=);

            /* Members */
            cls.def_static("makeScaling",
                           nb::overload_cast<double>(&LinearTransform::makeScaling),
                           "scale"_a);
            cls.def_static("makeScaling",
                           nb::overload_cast<double, double>(&LinearTransform::makeScaling));
            cls.def_static("makeRotation",
                           nb::overload_cast<Angle>(LinearTransform::makeRotation),
                           "angle"_a);
            cls.def("getParameterVector", &LinearTransform::getParameterVector);
            cls.def("getMatrix",
                    nb::overload_cast<>(& LinearTransform::getMatrix, nb::const_));
            cls.def("inverted", &LinearTransform::inverted);
            cls.def("computeDeterminant", &LinearTransform::computeDeterminant);
            cls.def("isIdentity", &LinearTransform::isIdentity);
            cls.def("applyX", nb::vectorize(&LinearTransform::applyX), "x"_a, "y"_a);
            cls.def("applyY", nb::vectorize(&LinearTransform::applyY), "x"_a, "y"_a);

            cls.def("set",
                    [](LinearTransform &self, double xx, double yx, double xy, double yy) {
                        self[LinearTransform::XX] = xx;
                        self[LinearTransform::XY] = xy;
                        self[LinearTransform::YX] = yx;
                        self[LinearTransform::YY] = yy;
                    },
                    "xx"_a, "yx"_a, "xy"_a, "yy"_a);

            cls.def("__str__", [](LinearTransform const &self) {
                return nb::str(nb::cast(self.getMatrix()));
            });
            cls.def("__repr__", [](LinearTransform const &self) {
                return nb::str("LinearTransform(\n{}\n)").format(nb::cast(self.getMatrix()));
            });
            cls.def("__reduce__", [cls](LinearTransform const &self) {
                return nb::make_tuple(cls, nb::make_tuple(nb::cast(self.getMatrix())));
            });
        }
    );
}

}  // namespace geom
}  // namespace lsst
