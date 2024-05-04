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
#include "nanobind/stl/vector.h"
#include "nanobind/stl/pair.h"
#include "nanobind/ndarray.h"
#include "lsst/sphgeom/python.h"
#include "ndarray/nanobind.h"

#include "lsst/geom/Interval.h"
#include "lsst/cpputils/python.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace lsst {
namespace geom {

namespace {

template <typename PyClass>
void declareCommonIntervalInterface(PyClass &cls) {
    using T = typename PyClass::Type;
    using Element = typename T::Element;
    cls.def(nb::init<>());
    // It's not clear why nb::overload cast doesn't work for the next two
    // declarations - maybe the templated context matters?
    cls.def_static("fromSpannedPoints", [](ndarray::Array<Element const, 1> const &elements) {
        return T::fromSpannedPoints(elements);
    });
    cls.def_static("fromSpannedPoints",
                   [](std::vector<Element> const &elements) { return T::fromSpannedPoints(elements); });
    cls.def("__init__", [](T &self, nb::kwargs kw) {
        if (kw.size() == 2) {
            if (kw.contains("min")) {
                if (kw.contains("max")) {
                     self = T::fromMinMax(nb::cast<Element>(kw["min"]), nb::cast<Element>(kw["max"]));
                   return;
		}
                if (kw.contains("size")) {
                    self =  T::fromMinSize(nb::cast<Element>(kw["min"]), nb::cast<Element>(kw["size"]));
                   return;
		}
            }
            if (kw.contains("max") && kw.contains("size")) {
                self = T::fromMaxSize(nb::cast<Element>(kw["max"]), nb::cast<Element>(kw["size"]));
               return;
	    }
            if (kw.contains("center") && kw.contains("size")) {
                self = T::fromCenterSize(nb::cast<Element>(kw["center"]), nb::cast<Element>(kw["size"]));
               return;
	    }
        }
        PyErr_SetString(PyExc_TypeError,
                        "General constructor requires one of the following keyword-only "
                        "argument pairs: (min, max), (center, size), (min, size), (max, size).");
        throw nb::python_error();
    });
    cls.def(nb::init<T const &>());
    cls.def("__eq__", [](T const &self, T const &other) { return self == other; }, nb::is_operator());
    cls.def("__ne__", [](T const &self, T const &other) { return self != other; }, nb::is_operator());
    cls.def("getMin", &T::getMin);
    cls.def_prop_ro("min", &T::getMin);
    cls.def("getMax", &T::getMax);
    cls.def_prop_ro("max", &T::getMax);
    cls.def("getSize", &T::getSize);
    cls.def_prop_ro("size", &T::getSize);
    cls.def("isEmpty", &T::isEmpty);
    cls.def("contains", nb::overload_cast<T const &>(&T::contains, nb::const_));
    cls.def("contains", nb::vectorize((bool (T::*)(Element) const) &T::contains));
    cls.def("contains", static_cast<bool (T::*)(Element) const>(&T::contains));
    cls.def("__contains__", nb::overload_cast<Element>(&T::contains, nb::const_));
    cls.def("__contains__", nb::overload_cast<T const &>(&T::contains, nb::const_));
    cls.def("overlaps", &T::overlaps);
    cls.def("intersects", &T::intersects);
    cls.def("isDisjointFrom", &T::isDisjointFrom);
    cls.def("dilatedBy", &T::dilatedBy);
    cls.def("erodedBy", &T::erodedBy);
    cls.def("shiftedBy", &T::shiftedBy);
    cls.def("reflectedAbout", &T::reflectedAbout);
    cls.def("expandedTo", nb::overload_cast<Element>(&T::expandedTo, nb::const_));
    cls.def("expandedTo", nb::overload_cast<T const &>(&T::expandedTo, nb::const_));
    cls.def("clippedTo", &T::clippedTo);
    cls.def("__str__", &T::toString);
    cpputils::python::addOutputOp(cls, "__repr__");
    cls.def("__reduce__", [cls](IntervalD const &self) {
        return nb::make_tuple(cls, make_tuple(nb::cast(self.getMin()), nb::cast(self.getMax())));
    });
}

}  // namespace

void wrapInterval(cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(nb::class_<IntervalI>(wrappers.module, "IntervalI"),
                      [](auto &mod, auto &cls) {
                          nb::enum_<IntervalI::EdgeHandlingEnum>(cls, "EdgeHandlingEnum")
                                  .value("EXPAND", IntervalI::EdgeHandlingEnum::EXPAND)
                                  .value("SHRINK", IntervalI::EdgeHandlingEnum::SHRINK);
                          cls.def(nb::init<IntervalD const &, IntervalI::EdgeHandlingEnum>(), "other"_a,
                                  "edgeHandling"_a = IntervalI::EdgeHandlingEnum::EXPAND);
                          cls.def("getBegin", &IntervalI::getBegin);
                          cls.def_prop_ro("begin", &IntervalI::getBegin);
                          cls.def("getEnd", &IntervalI::getEnd);
                          cls.def_prop_ro("end", &IntervalI::getEnd);
                          declareCommonIntervalInterface(cls);
                      });

    wrappers.wrapType(nb::class_<IntervalD>(wrappers.module, "IntervalD"),
                      [](auto &mod, auto &cls) {
                          cls.def(nb::init<IntervalI const &>());
                          cls.def("getCenter", &IntervalD::getCenter);
                          cls.def_prop_ro("center", &IntervalD::getCenter);
                          cls.def("isFinite", &IntervalD::isFinite);
                          declareCommonIntervalInterface(cls);
                      });
}

}  // namespace geom
}  // namespace lsst
