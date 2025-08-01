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
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "ndarray/pybind11.h"

#include "lsst/geom/Interval.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace geom {

namespace {

template <typename PyClass>
void declareCommonIntervalInterface(PyClass &cls) {
    using T = typename PyClass::type;
    using Element = typename T::Element;
    cls.def(py::init<>());
    // It's not clear why py::overload cast doesn't work for the next two
    // declarations - maybe the templated context matters?
    cls.def_static("fromSpannedPoints", [](ndarray::Array<Element const, 1> const &elements) {
        return T::fromSpannedPoints(elements);
    });
    cls.def_static("fromSpannedPoints",
                   [](std::vector<Element> const &elements) { return T::fromSpannedPoints(elements); });
    cls.def(py::init([](py::kwargs kw) -> T {
        if (kw.size() == 2) {
            if (kw.contains("min")) {
                if (kw.contains("max")) {
                    return T::fromMinMax(py::cast<Element>(kw["min"]), py::cast<Element>(kw["max"]));
                }
                if (kw.contains("size")) {
                    return T::fromMinSize(py::cast<Element>(kw["min"]), py::cast<Element>(kw["size"]));
                }
            }
            if (kw.contains("max") && kw.contains("size")) {
                return T::fromMaxSize(py::cast<Element>(kw["max"]), py::cast<Element>(kw["size"]));
            }
            if (kw.contains("center") && kw.contains("size")) {
                return T::fromCenterSize(py::cast<Element>(kw["center"]), py::cast<Element>(kw["size"]));
            }
        }
        PyErr_SetString(PyExc_TypeError,
                        "General constructor requires one of the following keyword-only "
                        "argument pairs: (min, max), (center, size), (min, size), (max, size).");
        throw py::error_already_set();
    }));
    cls.def(py::init<T const &>());
    cls.def("__eq__", [](T const &self, T const &other) { return self == other; }, py::is_operator());
    cls.def("__ne__", [](T const &self, T const &other) { return self != other; }, py::is_operator());
    cls.def("getMin", &T::getMin);
    cls.def_property_readonly("min", &T::getMin);
    cls.def("getMax", &T::getMax);
    cls.def_property_readonly("max", &T::getMax);
    cls.def("getSize", &T::getSize);
    cls.def_property_readonly("size", &T::getSize);
    cls.def("isEmpty", &T::isEmpty);
    cls.def("contains", py::overload_cast<T const &>(&T::contains, py::const_));
    cls.def("contains", py::vectorize(static_cast<bool (T::*)(Element) const>(&T::contains)));
    cls.def("__contains__", py::overload_cast<Element>(&T::contains, py::const_));
    cls.def("__contains__", py::overload_cast<T const &>(&T::contains, py::const_));
    cls.def("overlaps", &T::overlaps);
    cls.def("intersects", &T::intersects);
    cls.def("isDisjointFrom", &T::isDisjointFrom);
    cls.def("dilatedBy", &T::dilatedBy);
    cls.def("erodedBy", &T::erodedBy);
    cls.def("shiftedBy", &T::shiftedBy);
    cls.def("reflectedAbout", &T::reflectedAbout);
    cls.def("expandedTo", py::overload_cast<Element>(&T::expandedTo, py::const_));
    cls.def("expandedTo", py::overload_cast<T const &>(&T::expandedTo, py::const_));
    cls.def("clippedTo", &T::clippedTo);
    cls.def("__str__", &T::toString);
    cpputils::python::addOutputOp(cls, "__repr__");
    cls.def("__reduce__", [cls](IntervalD const &self) {
        return py::make_tuple(cls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
    });
}

}  // namespace

void wrapInterval(cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(py::classh<IntervalI>(wrappers.module, "IntervalI"),
                      [](auto &mod, auto &cls) {
                          py::enum_<IntervalI::EdgeHandlingEnum>(cls, "EdgeHandlingEnum")
                                  .value("EXPAND", IntervalI::EdgeHandlingEnum::EXPAND)
                                  .value("SHRINK", IntervalI::EdgeHandlingEnum::SHRINK);
                          cls.def(py::init<IntervalD const &, IntervalI::EdgeHandlingEnum>(), "other"_a,
                                  "edgeHandling"_a = IntervalI::EdgeHandlingEnum::EXPAND);
                          cls.def("getBegin", &IntervalI::getBegin);
                          cls.def_property_readonly("begin", &IntervalI::getBegin);
                          cls.def("getEnd", &IntervalI::getEnd);
                          cls.def_property_readonly("end", &IntervalI::getEnd);
                          declareCommonIntervalInterface(cls);
                      });

    wrappers.wrapType(py::classh<IntervalD>(wrappers.module, "IntervalD"),
                      [](auto &mod, auto &cls) {
                          cls.def(py::init<IntervalI const &>());
                          cls.def("getCenter", &IntervalD::getCenter);
                          cls.def_property_readonly("center", &IntervalD::getCenter);
                          cls.def("isFinite", &IntervalD::isFinite);
                          declareCommonIntervalInterface(cls);
                      });
}

}  // namespace geom
}  // namespace lsst
