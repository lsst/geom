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
#include "nanobind/ndarray.h"
#include "nanobind/stl/vector.h"
#include "nanobind/stl/list.h"
#include "nanobind/stl/array.h"
#include "nanobind/stl/tuple.h"

#include "ndarray/nanobind.h"

#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"
#include "lsst/geom/CoordinateBase.h"
#include "lsst/geom/CoordinateExpr.h"
#include "lsst/cpputils/python.h"

namespace nb = nanobind;
using namespace nanobind::literals;

namespace lsst {
namespace geom {
namespace {
template<std::size_t N, typename T>
struct TupleGenerator {
    using type = decltype(std::tuple_cat(std::declval<std::tuple<T>>(), typename TupleGenerator<N-1, T>::type{}));
};

template<typename T>
struct TupleGenerator<1, T> {
    using type = std::tuple<T>;
};

template <typename Derived, typename T, int N>
using PyCoordinateBase = nb::class_<CoordinateBase<Derived, T, N>>;

template <int N>
using PyCoordinateExpr = nb::class_<CoordinateExpr<N>, CoordinateBase<CoordinateExpr<N>, bool, N>>;

template <typename T, int N>
using PyExtentBase = nb::class_<ExtentBase<T, N>, CoordinateBase<Extent<T, N>, T, N>>;

template <typename T, int N>
using PyExtent = nb::class_<Extent<T, N>, ExtentBase<T, N>>;

template <typename T, int N>
using PyPointBase = nb::class_<PointBase<T, N>, CoordinateBase<Point<T, N>, T, N>>;

template <typename T, int N>
using PyPoint = nb::class_<Point<T, N>, PointBase<T, N>>;

template <typename Derived, typename T, int N>
void declareCoordinateBase(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "CoordinateBase" + suffix;
    wrappers.wrapType(
        PyCoordinateBase<Derived, T, N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            cls.def("__getitem__", [](CoordinateBase<Derived, T, N> &self, int i) -> T {
                return self[cpputils::python::cppIndex(N, i)];
            });
            cls.def("__setitem__", [](CoordinateBase<Derived, T, N> &self, int i, T value) {
                self[cpputils::python::cppIndex(N, i)] = value;
            });
            cls.def("__len__", [](CoordinateBase<Derived, T, N> &c) -> int { return N; });
        }
    );
}

template <int N>
void declareCoordinateExpr(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "CoordinateExpr" + suffix;
    declareCoordinateBase<CoordinateExpr<N>, bool, N>(wrappers, name);
    wrappers.wrapType(
        PyCoordinateExpr<N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            cls.def(nb::init<bool>(), "val"_a = false);
            cls.def("and_", &CoordinateExpr<N>::and_);
            cls.def("or_", &CoordinateExpr<N>::or_);
            cls.def("not_", &CoordinateExpr<N>::not_);
            mod.def("all", all<N>);
            mod.def("any", any<N>);
        }
    );
}

template <typename T, int N>
void declareExtentBase(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "ExtentBase" + suffix;
    declareCoordinateBase<Extent<T, N>, T, N>(wrappers, "Extent" + suffix);
    wrappers.wrapType(
        PyExtentBase<T, N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            // These are not the usual Python double-underscore operators - they do elementwise comparisons,
            // returning a CoordinateExpr object with boolean x and y values.  NumPy examples to the contrary
            // notwithstanding, true Python comparison operators are expected to return scalar bools.
            cls.def("eq", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.eq(other); });
            cls.def("ne", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.ne(other); });
            cls.def("lt", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.lt(other); });
            cls.def("le", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.le(other); });
            cls.def("gt", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.gt(other); });
            cls.def("ge", [](ExtentBase<T, N> const &self, Extent<T, N> other) { return self.ge(other); });
            cls.def("eq", [](ExtentBase<T, N> const &self, T other) { return self.eq(other); });
            cls.def("ne", [](ExtentBase<T, N> const &self, T other) { return self.ne(other); });
            cls.def("lt", [](ExtentBase<T, N> const &self, T other) { return self.lt(other); });
            cls.def("le", [](ExtentBase<T, N> const &self, T other) { return self.le(other); });
            cls.def("gt", [](ExtentBase<T, N> const &self, T other) { return self.gt(other); });
            cls.def("ge", [](ExtentBase<T, N> const &self, T other) { return self.ge(other); });

            cls.def("asPoint", &ExtentBase<T, N>::asPoint);
            cls.def("computeNorm", &ExtentBase<T, N>::computeNorm);
            cls.def("computeSquaredNorm", &ExtentBase<T, N>::computeSquaredNorm);
        }
    );
}

// Common functionality for all Extents, including declaring base classes.
template <typename T, int N>
PyExtent<T, N> declareExtent(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "Extent" + suffix;
    declareExtentBase<T, N>(wrappers, suffix);
    return wrappers.wrapType(
        PyExtent<T, N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            cls.def(nb::init<T>(), "value"_a = static_cast<T>(0));
            cls.def(nb::init<Point<int, N> const &>());
            cls.def(nb::init<Point<T, N> const &>());
            cls.def(nb::init<Extent<int, N> const &>());
            cls.def(nb::init<Extent<T, N> const &>());
            cls.def(nb::init<typename Extent<T, N>::EigenVector>());
            cls.def("__neg__", [](Extent<T, N> const &self) { return -self; });
            cls.def("__pos__", [](Extent<T, N> const &self) { return self; });
            cls.def("__mul__", [](Extent<T, N> const &self, int other) { return self * other; },
                    nb::is_operator());
            cls.def("__mul__", [](Extent<T, N> const &self, double other) { return self * other; },
                    nb::is_operator());
            cls.def("__rmul__", [](Extent<T, N> const &self, int other) { return self * other; },
                    nb::is_operator());
            cls.def("__rmul__", [](Extent<T, N> const &self, double other) { return self * other; },
                    nb::is_operator());
            cls.def("__add__",
                    [](Extent<T, N> const &self, Extent<int, N> const &other) { return self + other; },
                    nb::is_operator());
            cls.def("__add__",
                    [](Extent<T, N> const &self, Extent<double, N> const &other) { return self + other; },
                    nb::is_operator());
            cls.def("__add__",
                    [](Extent<T, N> const &self, Point<int, N> const &other) {
                        return self + Point<T, N>(other);
                    },
                    nb::is_operator());
            cls.def("__add__",
                    [](Extent<T, N> const &self, Point<double, N> const &other) { return self + other; },
                    nb::is_operator());
            cls.def("__sub__",
                    [](Extent<T, N> const &self, Extent<int, N> const &other) {
                        return self - Extent<T, N>(other);
                    },
                    nb::is_operator());
            cls.def("__sub__",
                    [](Extent<T, N> const &self, Extent<double, N> const &other) { return self - other; },
                    nb::is_operator());
            cls.def("__eq__",
                    [](Extent<T, N> const &self, Extent<T, N> const &other) { return self == other; },
                    nb::is_operator());
            cls.def("__ne__",
                    [](Extent<T, N> const &self, Extent<T, N> const &other) { return self != other; },
                    nb::is_operator());
            cls.def("clone", [](Extent<T, N> const &self) { return Extent<T, N>{self}; });
        }
    );
}

// Add functionality only found in N=2 Extents (and delgate to declareExtent for the rest)
template <typename T>
PyExtent<T, 2> declareExtent2(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    return wrappers.wrapType(
        declareExtent<T, 2>(wrappers, std::string("2") + suffix),
        [](auto & mod, auto & cls) {
            /* Members types and enums */
            cls.def_prop_ro_static("dimensions", [](nb::object /* cls */) { return 2; });
            /* Constructors */
            cls.def(nb::init<int, int>(), "x"_a, "y"_a);
            cls.def(nb::init<double, double>(), "x"_a, "y"_a);
            /* Members */
            auto getX = nb::overload_cast<>(&Extent<T, 2>::getX, nb::const_);
            auto getY = nb::overload_cast<>(&Extent<T, 2>::getY, nb::const_);
            cls.def("getX", getX);
            cls.def("getY", getY);
            cls.def("setX", &Extent<T, 2>::setX);
            cls.def("setY", &Extent<T, 2>::setY);
            cls.def_prop_rw("x", getX, &Extent<T, 2>::setX);
            cls.def_prop_rw("y", getY, &Extent<T, 2>::setY);
        }
    );
}

// Add functionality only found in N=3 Extents (and delgate to declareExtent for the rest)
template <typename T>
PyExtent<T, 3> declareExtent3(cpputils::python::WrapperCollection & wrappers, const std::string &suffix) {
    return wrappers.wrapType(
        declareExtent<T, 3>(wrappers, std::string("3") + suffix),
        [](auto & mod, auto & cls) mutable {
            /* Member types and enums */
            cls.def_prop_ro_static("dimensions", [](nb::object /* cls */) { return 3; });
            /* Constructors */
            cls.def(nb::init<int, int, int>(), "x"_a, "y"_a, "z"_a);
            cls.def(nb::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
            /* Members */
            auto getX = nb::overload_cast<>(&Extent<T, 3>::getX, nb::const_);
            auto getY = nb::overload_cast<>(&Extent<T, 3>::getY, nb::const_);
            auto getZ = nb::overload_cast<>(&Extent<T, 3>::getZ, nb::const_);
            cls.def("getX", getX);
            cls.def("getY", getY);
            cls.def("getZ", getZ);
            cls.def("setX", &Extent<T, 3>::setX);
            cls.def("setY", &Extent<T, 3>::setY);
            cls.def("setZ", &Extent<T, 3>::setZ);
            cls.def_prop_rw("x", getX, &Extent<T, 3>::setX);
            cls.def_prop_rw("y", getY, &Extent<T, 3>::setY);
            cls.def_prop_rw("z", getZ, &Extent<T, 3>::setZ);
        }
    );
}

// Declare mixed-type and type-overloaded operators for Extent with dimension
// N for both int and double. Because nanobind tries operators (like any
// overload) `in order', int has to come before double in any overloaded
// operators that dispatch on a scalar, and hence they have to be defined here
// instead of declareExtent.
template <int N>
void declareExtentOperators(cpputils::python::WrapperCollection & wrapper,
                         PyExtent<int, N> &clsI, PyExtent<double, N> &clsD) {
    wrapper.wrap(
        [clsI, clsD](auto & mod) mutable {
            // Python's integer division works differently than C++'s for negative numbers - Python
            // uses floor (rounds towards more negative), while C++ truncates (rounds towards zero).
            // Therefore one needs to be careful in the definition of division operators.
            clsI.def("__floordiv__",
                     [](Extent<int, N> const &self, int other) -> Extent<int, N> {
                         return floor(self / static_cast<double>(other));
                     },
                     nb::is_operator());

            clsI.def("__truediv__", [](Extent<int, N> const &self, double other) { return self / other; },
                     nb::is_operator());
            clsD.def("__truediv__", [](Extent<double, N> const &self, double other) { return self / other; },
                     nb::is_operator());

            clsI.def("__ifloordiv__", [](Extent<int, N> &self, int other) -> Extent<int, N> & {
                self = floor(self / static_cast<double>(other));
                return self;
            }, nb::rv_policy::none);

            clsI.def("__itruediv__", [](Extent<int, N> &self, double other) {
                PyErr_SetString(PyExc_TypeError, "In-place true division not supported for Extent<int,N>.");
                throw nb::python_error();
            }, nb::is_operator());
            clsD.def("__itruediv__", [](Extent<double, N> &self, double other) -> Extent<double, N> & {
                self /= other;
                return self;
            }, nb::rv_policy::none);

            clsI.def("__iadd__", [](Extent<int, N> &self, Extent<int, N> const &other) -> Extent<int, N> & {
                self += other;
                return self;
            }, nb::rv_policy::none);
            clsD.def(
               "__iadd__",
               [](Extent<double, N> &self, Extent<double, N> const &other) -> Extent<double, N> & {
                   self += other;
                    return self;
                } ,nb::rv_policy::none);
            clsD.def(
                "__iadd__",
                [](Extent<double, N> &self, Extent<int, N> const &other) -> Extent<double, N> & {
                    self += other;
                    return self;
                }, nb::rv_policy::none);
            clsI.def("__isub__", [](Extent<int, N> &self, Extent<int, N> const &other) -> Extent<int, N> & {
                self -= other;
                return self;
            }, nb::rv_policy::none);
            clsD.def("__isub__", [](Extent<double, N> &self, Extent<double, N> const &other) -> Extent<double, N> & {
                self -= other;
                return self;
            }, nb::rv_policy::none);
            clsD.def("__isub__", [](Extent<double, N> &self, Extent<int, N> const &other) -> Extent<double, N> & {
                self -= other;
                return self;
            }, nb::rv_policy::none);

            clsI.def("__imul__", [](Extent<int, N> &self, int other) -> Extent<int, N> & {
                self *= other;
                return self;
            }, nb::rv_policy::none);
            clsD.def("__imul__", [](Extent<double, N> &self, int other) -> Extent<double, N> & {
                self *= other;
                return self;
            }, nb::rv_policy::none);
            clsD.def("__imul__", [](Extent<double, N> &self, double other) -> Extent<double, N> & {
                self *= other;
                return self;
            },  nb::rv_policy::none);

            // Operator-like free functions
            mod.def("truncate", truncate<N>);
            mod.def("floor", floor<N>);
            mod.def("ceil", ceil<N>);
            // And method versions, since Python doens't have ADL
            clsD.def("truncate", truncate<N>);
            clsD.def("floor", floor<N>);
            clsD.def("ceil", ceil<N>);
        }
    );
}

template <typename T, int N>
void declarePointBase(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "PointBase" + suffix;
    declareCoordinateBase<Point<T, N>, T, N>(wrappers, "Point" + suffix);
    wrappers.wrapType(
        PyPointBase<T, N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            // These are not the usual Python double-underscore operators - they do elementwise comparisons,
            // returning a CoordinateExpr object with boolean x and y values.  NumPy examples to the contrary
            // notwithstanding, true Python comparison operators are expected to return scalar bools.
            cls.def("eq", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.eq(rhs); });
            cls.def("ne", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.ne(rhs); });
            cls.def("lt", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.lt(rhs); });
            cls.def("le", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.le(rhs); });
            cls.def("gt", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.gt(rhs); });
            cls.def("ge", [](PointBase<T, N> const &self, Point<T, N> const &rhs) { return self.ge(rhs); });
            cls.def("eq", [](PointBase<T, N> const &self, T rhs) { return self.eq(rhs); });
            cls.def("ne", [](PointBase<T, N> const &self, T rhs) { return self.ne(rhs); });
            cls.def("lt", [](PointBase<T, N> const &self, T rhs) { return self.lt(rhs); });
            cls.def("le", [](PointBase<T, N> const &self, T rhs) { return self.le(rhs); });
            cls.def("gt", [](PointBase<T, N> const &self, T rhs) { return self.gt(rhs); });
            cls.def("ge", [](PointBase<T, N> const &self, T rhs) { return self.ge(rhs); });
            /* Members */
            cls.def("asExtent", &PointBase<T, N>::asExtent);
            cls.def("shift", &PointBase<T, N>::shift);
            cls.def("scale", &PointBase<T, N>::scale);
            cls.def("distanceSquared", &PointBase<T, N>::distanceSquared);
            cls.def("toString", &PointBase<T, N>::toString);
        }
    );
}

// Common functionality
template <typename T, int N>
PyPoint<T, N> declarePoint(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    static std::string const name = "Point" + suffix;
    declarePointBase<T, N>(wrappers, suffix);
    using tuple_t = typename TupleGenerator<N, T>::type;
    return wrappers.wrapType(
        PyPoint<T, N>(wrappers.module, name.c_str()),
        [](auto & mod, auto & cls) {
            /* Constructors */
            cls.def(nb::init<T>(), "value"_a = static_cast<T>(0));
            // Note that we can't use T here because both types are needed
            cls.def(nb::init<Point<double, N> const &>());
            cls.def(nb::init<Point<int, N> const &>());
            cls.def(nb::init<Extent<int, N> const &>());
            cls.def(nb::init<Extent<double, N> const &>());
            cls.def(nb::init<typename Point<T, N>::EigenVector>());
            /* Operators */
            cls.def("__add__", [](Point<T, N> const &self, Extent<double, N> &other) { return self + other; },
                    nb::is_operator());
            cls.def("__add__", [](Point<T, N> const &self, Extent<int, N> &other) { return self + other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Point<T, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Extent<T, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Point<double, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Point<int, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Extent<double, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__sub__", [](Point<T, N> const &self, Extent<int, N> &other) { return self - other; },
                    nb::is_operator());
            cls.def("__eq__", [](Point<T, N> const &self, Point<T, N> const &other) { return self == other; },
                    nb::is_operator());
            cls.def("__ne__", [](Point<T, N> const &self, Point<T, N> const &other) { return self != other; },
                    nb::is_operator());
            /* Members */
            cls.def("clone", [](Point<T, N> const &self) { return Point<T, N>{self}; });
        }
    );
}

// Add functionality only found in N=2 Points
template <typename T>
PyPoint<T, 2> declarePoint2(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    return wrappers.wrapType(
        declarePoint<T, 2>(wrappers, std::string("2") + suffix),
        [](auto & mod, auto & cls) {
            /* Member types and enums */
            cls.def_prop_ro_static("dimensions", [](nb::object /* cls */) { return 2; });
            /* Constructors */
            cls.def(nb::init<int, int>(), "x"_a, "y"_a);
            cls.def(nb::init<double, double>(), "x"_a, "y"_a);

            cls.def("__init__", [](Point<T, 2> *point, nb::tuple arg) {
                    new (point) Point<T, 2>(nb::cast<T>(arg[0]), nb::cast<T>(arg[1]));
            });
            /* Members */
            auto getX = nb::overload_cast<>(&Point<T, 2>::getX, nb::const_);
            auto getY = nb::overload_cast<>(&Point<T, 2>::getY, nb::const_);
            cls.def("getX", getX);
            cls.def("getY", getY);
            cls.def("setX", &Point<T, 2>::setX);
            cls.def("setY", &Point<T, 2>::setY);
            cls.def_prop_rw("x", getX, &Point<T, 2>::setX);
            cls.def_prop_rw("y", getY, &Point<T, 2>::setY);
        }
    );
}

// Add functionality only found in N=3 Points
template <typename T>
PyPoint<T, 3> declarePoint3(cpputils::python::WrapperCollection & wrappers, std::string const &suffix) {
    return wrappers.wrapType(
        declarePoint<T, 3>(wrappers, std::string("3") + suffix),
        [](auto & mod, auto & cls) {
            /* Member types and enums */
            cls.def_prop_ro_static("dimensions", [](nb::object /* cls */) { return 3; });
            /* Constructors */
            cls.def(nb::init<int, int, int>(), "x"_a, "y"_a, "z"_a);
            cls.def(nb::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
            cls.def("__init__", [](Point<T, 3> *point, nb::tuple arg) {
                new (point) Point<T, 3>(nb::cast<T>(arg[0]), nb::cast<T>(arg[1]), nb::cast<T>(arg[2]));
            });
            /* Members */
            auto getX = nb::overload_cast<>(&Point<T, 3>::getX, nb::const_);
            auto getY = nb::overload_cast<>(&Point<T, 3>::getY, nb::const_);
            auto getZ = nb::overload_cast<>(&Point<T, 3>::getZ, nb::const_);
            cls.def("getX", getX);
            cls.def("getY", getY);
            cls.def("getZ", getZ);
            cls.def("setX", &Point<T, 3>::setX);
            cls.def("setY", &Point<T, 3>::setY);
            cls.def("setZ", &Point<T, 3>::setZ);
            cls.def_prop_rw("x", getX, &Point<T, 3>::setX);
            cls.def_prop_rw("y", getY, &Point<T, 3>::setY);
            cls.def_prop_rw("z", getZ, &Point<T, 3>::setZ);
        }
    );
}

// Declare mixed-type and type-overloaded operators for Point with dimension
// N for both int and double. Because nanobind tries operators (like any
// overload) `in order', int has to come before double in any overloaded
// operators that dispatch on a scalar, and hence they have to be defined here
// instead of declareExtent.
template <int N>
void declarePointOperators(cpputils::python::WrapperCollection & wrappers,
                        PyPoint<int, N> &clsI, PyPoint<double, N> &clsD) {
    wrappers.wrap(
        [clsI, clsD](auto & mod) mutable {
            clsI.def("__iadd__", [](Point<int, N> &self, Extent<int, N> const &other) {
                self += other;
                return &self;
            }, nb::rv_policy::none);
            clsD.def("__iadd__", [](Point<double, N> &self, Extent<int, N> const &other) {
                self += other;
                return &self;
            }, nb::rv_policy::none);
            clsD.def("__iadd__", [](Point<double, N> &self, Extent<double, N> const &other) {
                self += other;
                return &self;
            }, nb::rv_policy::none);
            clsI.def("__isub__", [](Point<int, N> &self, Extent<int, N> const &other) {
                self -= other;
                return &self;
            }, nb::rv_policy::none);
            clsD.def("__isub__", [](Point<double, N> &self, Extent<int, N> const &other) {
                self -= other;
                return &self;
            }, nb::rv_policy::none);
            clsD.def("__isub__", [](Point<double, N> &self, Extent<double, N> const &other) {
                self -= other;
                return &self;
            }, nb::rv_policy::none);
        }
    );
}

}  // anonymous

void wrapCoordinates(cpputils::python::WrapperCollection & wrappers) {
    // Only the interface-level classes are defined here, and these functions
    // call others to define their base classes, since those are never shared.

    declareCoordinateExpr<2>(wrappers, "2");
    declareCoordinateExpr<3>(wrappers, "3");

    auto clsExtent2I = declareExtent2<int>(wrappers, "I");
    auto clsExtent2D = declareExtent2<double>(wrappers, "D");

    auto clsExtent3I = declareExtent3<int>(wrappers, "I");
    auto clsExtent3D = declareExtent3<double>(wrappers, "D");

    auto clsPoint2I = declarePoint2<int>(wrappers, "I");
    auto clsPoint2D = declarePoint2<double>(wrappers, "D");

    auto clsPoint3I = declarePoint3<int>(wrappers, "I");
    auto clsPoint3D = declarePoint3<double>(wrappers, "D");

    declareExtentOperators(wrappers, clsExtent2I, clsExtent2D);
    declareExtentOperators(wrappers, clsExtent3I, clsExtent3D);
    declarePointOperators(wrappers, clsPoint2I, clsPoint2D);
    declarePointOperators(wrappers, clsPoint3I, clsPoint3D);

}

}  // namespace geom
}  // namespace lsst
