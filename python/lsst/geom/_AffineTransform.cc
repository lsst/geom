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
#include "lsst/sphgeom/python.h"
//#include "nanobind/stl.h"

#include "ndarray/nanobind.h"

#include "lsst/geom/AffineTransform.h"
#include "lsst/cpputils/python.h"

namespace nb = nanobind;
using namespace nanobind::literals;

namespace lsst {
namespace geom {

void wrapAffineTransform(cpputils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        nb::class_<AffineTransform>(wrappers.module, "AffineTransform"),
        [](auto & mod, auto & cls) mutable {

            // Parameters enum is really only used as integer constants.
            cls.attr("XX") = nb::cast(int(AffineTransform::Parameters::XX));
            cls.attr("YX") = nb::cast(int(AffineTransform::Parameters::YX));
            cls.attr("XY") = nb::cast(int(AffineTransform::Parameters::XY));
            cls.attr("YY") = nb::cast(int(AffineTransform::Parameters::YY));
            cls.attr("X") = nb::cast(int(AffineTransform::Parameters::X));
            cls.attr("Y") = nb::cast(int(AffineTransform::Parameters::Y));

            /* Constructors */
            cls.def(nb::init<>());
            cls.def(nb::init<Eigen::Matrix3d const &>(), "matrix"_a);
            cls.def(nb::init<Eigen::Matrix2d const &>(), "linear"_a);
            cls.def(nb::init<Eigen::Vector2d const &>(), "translation"_a);
            cls.def(nb::init<Eigen::Matrix2d const &, Eigen::Vector2d const &>(),
                    "linear"_a, "translation"_a);
            cls.def(nb::init<LinearTransform const &>(), "linear"_a);
            cls.def(nb::init<Extent2D const &>(), "translation"_a);
            cls.def(nb::init<LinearTransform const &, Extent2D const &>(), "linear"_a, "translation"_a);

            /* Operators and special methods */
            cls.def("__mul__", &AffineTransform::operator*, nb::is_operator());
            cls.def("__call__",
                    nb::overload_cast<Point2D const &>(&AffineTransform::operator(), nb::const_));
            cls.def("__call__",
                    nb::overload_cast<Extent2D const &>(&AffineTransform::operator(), nb::const_));
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
            cls.def("__setitem__", [](AffineTransform &self, int i, double value) {
                if (i < 0 || i > 5) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
                    throw nb::python_error();
                }
                self[i] = value;
            });
            cls.def("__getitem__", [](AffineTransform const &self, int row, int col) {
                if (row < 0 || row > 2 || col < 0 || col > 2) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d, %d", row, col);
                    throw nb::python_error();
                }
                return (self.getMatrix())(row, col);
            });
            cls.def("__getitem__", [](AffineTransform const &self, int i) {
                if (i < 0 || i > 5) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
                    throw nb::python_error();
                }
                return self[i];
            });
            cls.def("__str__", [](AffineTransform const &self) {
                return nb::str(nb::cast(self.getMatrix())); }
            );
            cls.def("__repr__", [](AffineTransform const &self) {
                return nb::str("AffineTransform(\n{}\n)").format(nb::cast(self.getMatrix()));
            });
            cls.def("__reduce__", [cls](AffineTransform const &self) {
                return nb::make_tuple(cls, nb::make_tuple(nb::cast(self.getMatrix())));
            });

            /* Members */
            cls.def("inverted", &AffineTransform::inverted);
            cls.def("isIdentity", &AffineTransform::isIdentity);
            cls.def("getTranslation", (Extent2D & (AffineTransform::*)()) & AffineTransform::getTranslation);
            cls.def("getLinear", (LinearTransform & (AffineTransform::*)()) & AffineTransform::getLinear);
            cls.def("getMatrix", &AffineTransform::getMatrix, nb::rv_policy::reference);
            cls.def("getParameterVector", &AffineTransform::getParameterVector);
            cls.def("setParameterVector", &AffineTransform::setParameterVector);
            cls.def("applyX", nb::vectorize(&AffineTransform::applyX), "x"_a, "y"_a);
            cls.def("applyY", nb::vectorize(&AffineTransform::applyY), "x"_a, "y"_a);
            cls.def_static("makeScaling", nb::overload_cast<double>(&AffineTransform::makeScaling));
            cls.def_static("makeScaling", nb::overload_cast<double, double>(&AffineTransform::makeScaling));
            cls.def_static("makeRotation", &AffineTransform::makeRotation, "angle"_a);
            cls.def_static("makeTranslation", &AffineTransform::makeTranslation, "translation"_a);

            /* Non-members */
            mod.def("makeAffineTransformFromTriple", makeAffineTransformFromTriple);
        }
    );
}

}  // namespace geom
}  // namespace lsst
