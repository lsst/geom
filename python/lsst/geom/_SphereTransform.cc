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

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"

#include "ndarray/pybind11.h"

#include "lsst/cpputils/python.h"

#include "lsst/geom/SphereTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace geom {

void wrapSphereTransform(cpputils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        py::class_<SphereTransform, std::shared_ptr<SphereTransform>>(wrappers.module, "SphereTransform"),
        [](auto & mod, auto & cls) {

            cls.def(py::init<>());
            cls.def(py::init<SphereTransform::Matrix const &>(), "matrix"_a);
            cls.def_static(
                "fit_unit_vectors", &SphereTransform::fit_unit_vectors,
                "from"_a, "to"_a, "weights"_a = py::none()
            );
            cls.def("__call__",
                    py::overload_cast<SpherePoint const &>(&SphereTransform::operator(), py::const_));
            cls.def("__call__",
                    py::overload_cast<sphgeom::UnitVector3d const &>(&SphereTransform::operator(), py::const_));
            cls.def("__call__",
                    [](py::object self, py::object x, py::object y, py::object z) {
                        return py::make_tuple(self.attr("applyX")(x, y, z),
                                              self.attr("applyY")(x, y, z),
                                              self.attr("applyZ")(x, y, z));
                    },
                    "x"_a, "y"_a, "z"_a);
            cls.def("__mul__", &SphereTransform::operator*, py::is_operator());
            cls.def("getMatrix", &SphereTransform::getMatrix, py::return_value_policy::reference_internal);
            cls.def_property_readonly("matrix", &SphereTransform::getMatrix, py::return_value_policy::reference_internal);
            cls.def("inverted", &SphereTransform::inverted);
            cls.def("applyX", py::vectorize(&SphereTransform::applyX), "x"_a, "y"_a, "z"_a);
            cls.def("applyY", py::vectorize(&SphereTransform::applyY), "x"_a, "y"_a, "z"_a);
            cls.def("applyZ", py::vectorize(&SphereTransform::applyZ), "x"_a, "y"_a, "z"_a);

            cls.def("__str__", [](SphereTransform const &self) {
                return py::str(py::cast(self.getMatrix()));
            });
            cls.def("__repr__", [](SphereTransform const &self) {
                return py::str("SphereTransform(\n{}\n)").format(py::cast(self.getMatrix()));
            });
            cls.def("__reduce__", [cls](SphereTransform const &self) {
                return py::make_tuple(cls, py::make_tuple(py::cast(self.getMatrix())));
            });
        }
    );
}

}  // namespace geom
}  // namespace lsst
