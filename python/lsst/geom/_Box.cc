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
#include "nanobind/stl/tuple.h"
#include "nanobind/ndarray.h"
#include "lsst/sphgeom/python.h"

#include "lsst/geom/Box.h"
#include "lsst/cpputils/python.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace lsst {
namespace geom {

void wrapBox(cpputils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        nb::class_<Box2I>(wrappers.module, "Box2I"),
        [](auto & mod, auto & cls) mutable {

            cls.attr("Point") = mod.attr("Point2I");
            cls.attr("Extent") = mod.attr("Extent2I");

            nb::enum_<Box2I::EdgeHandlingEnum>(cls, "EdgeHandlingEnum")
                    .value("EXPAND", Box2I::EdgeHandlingEnum::EXPAND)
                    .value("SHRINK", Box2I::EdgeHandlingEnum::SHRINK)
                    .export_values();

            cls.def(nb::init<>());
            cls.def(nb::init<Point2I const &, Point2I const &, bool>(), "minimum"_a, "maximum"_a,
                    "invert"_a = true);
            cls.def(nb::init<Point2I const &, Extent2I const &, bool>(), "corner"_a, "dimensions"_a,
                    "invert"_a = true);
            cls.def(nb::init<IntervalI const &, IntervalI const &>(), "x"_a, "y"_a);
            cls.def(nb::init<Box2D const &, Box2I::EdgeHandlingEnum>(), "other"_a,
                    "edgeHandling"_a = Box2I::EXPAND);
            cls.def(nb::init<Box2I const &>(), "other"_a);

            cls.def("__eq__", [](Box2I const &self, Box2I const &other) { return self == other; },
                    nb::is_operator());
            cls.def("__ne__", [](Box2I const &self, Box2I const &other) { return self != other; },
                    nb::is_operator());

            cls.def_static("makeCenteredBox", &Box2I::makeCenteredBox, "center"_a, "size"_a);
            cls.def("swap", &Box2I::swap);
            cls.def("getMin", &Box2I::getMin);
            cls.def("getMinX", &Box2I::getMinX);
            cls.def("getMinY", &Box2I::getMinY);
            cls.def("getMax", &Box2I::getMax);
            cls.def("getMaxX", &Box2I::getMaxX);
            cls.def("getMaxY", &Box2I::getMaxY);
            cls.def_prop_ro("minX", &Box2I::getMinX);
            cls.def_prop_ro("minY", &Box2I::getMinY);
            cls.def_prop_ro("maxX", &Box2I::getMaxX);
            cls.def_prop_ro("maxY", &Box2I::getMaxY);
            cls.def("getBegin", &Box2I::getBegin);
            cls.def("getBeginX", &Box2I::getBeginX);
            cls.def("getBeginY", &Box2I::getBeginY);
            cls.def("getEnd", &Box2I::getEnd);
            cls.def("getEndX", &Box2I::getEndX);
            cls.def("getEndY", &Box2I::getEndY);
            cls.def_prop_ro("beginX", &Box2I::getBeginX);
            cls.def_prop_ro("beginY", &Box2I::getBeginY);
            cls.def_prop_ro("endX", &Box2I::getEndX);
            cls.def_prop_ro("endY", &Box2I::getEndY);
            cls.def("getDimensions", &Box2I::getDimensions);
            cls.def("getWidth", &Box2I::getWidth);
            cls.def("getHeight", &Box2I::getHeight);
            cls.def("getArea", &Box2I::getArea);
            cls.def_prop_ro("width", &Box2I::getWidth);
            cls.def_prop_ro("height", &Box2I::getHeight);
            cls.def_prop_ro("area", &Box2I::getArea);
            cls.def("getCenter", &Box2I::getCenter);
            cls.def("getCenterX", &Box2I::getCenterX);
            cls.def("getCenterY", &Box2I::getCenterY);
            cls.def_prop_ro("centerX", &Box2I::getCenterX);
            cls.def_prop_ro("centerY", &Box2I::getCenterY);
            cls.def("getX", &Box2I::getX);
            cls.def("getY", &Box2I::getY);
            cls.def_prop_ro("x", &Box2I::getX);
            cls.def_prop_ro("y", &Box2I::getY);
            cls.def("isEmpty", &Box2I::isEmpty);
            cls.def("contains", nb::overload_cast<Point2I const &>(&Box2I::contains, nb::const_));
            cls.def("contains", nb::overload_cast<Box2I const &>(&Box2I::contains, nb::const_));
            cls.def("contains",
                    nb::vectorize(static_cast<bool (Box2I::*)(int x, int y) const>(&Box2I::contains)),
                    "x"_a, "y"_a);
	    cls.def("contains",
                    static_cast<bool (Box2I::*)(int x, int y) const>(&Box2I::contains),
                    "x"_a, "y"_a);
            cls.def("__contains__", nb::overload_cast<Point2I const &>(&Box2I::contains, nb::const_));
            cls.def("__contains__", nb::overload_cast<Box2I const &>(&Box2I::contains, nb::const_));
            cls.def("overlaps", &Box2I::overlaps);
            cls.def("intersects", &Box2I::intersects);
            cls.def("isDisjointFrom", &Box2I::isDisjointFrom);
            cls.def("grow", nb::overload_cast<int>(&Box2I::grow));
            cls.def("grow", nb::overload_cast<Extent2I const &>(&Box2I::grow));
            cls.def("shift", &Box2I::shift);
            cls.def("flipLR", &Box2I::flipLR);
            cls.def("flipTB", &Box2I::flipTB);
            cls.def("include", nb::overload_cast<Point2I const &>(&Box2I::include));
            cls.def("include", nb::overload_cast<Box2I const &>(&Box2I::include));
            cls.def("clip", &Box2I::clip);
            cls.def("dilatedBy", nb::overload_cast<int>(&Box2I::dilatedBy, nb::const_));
            cls.def("dilatedBy", nb::overload_cast<Extent2I const &>(&Box2I::dilatedBy, nb::const_));
            cls.def("erodedBy", nb::overload_cast<int>(&Box2I::erodedBy, nb::const_));
            cls.def("erodedBy", nb::overload_cast<Extent2I const &>(&Box2I::erodedBy, nb::const_));
            cls.def("shiftedBy", &Box2I::shiftedBy);
            cls.def("reflectedAboutX", &Box2I::reflectedAboutX);
            cls.def("reflectedAboutY", &Box2I::reflectedAboutY);
            cls.def("expandedTo", nb::overload_cast<Point2I const &>(&Box2I::expandedTo, nb::const_));
            cls.def("expandedTo", nb::overload_cast<Box2I const &>(&Box2I::expandedTo, nb::const_));
            cls.def("clippedTo", &Box2I::clippedTo);
            cls.def("getCorners", &Box2I::getCorners);
            cls.def("toString", &Box2I::toString);
            cls.def("__repr__", [](Box2I const &self) {
                return nb::str("Box2I(corner={}, dimensions={})")
                    .format(nb::repr(nb::cast(self.getMin())), nb::repr(nb::cast(self.getDimensions())));
            });
            cls.def("__str__", [](Box2I const &self) {
                return nb::str("(minimum={}, maximum={})")
                    .format(nb::str(nb::cast(self.getMin())), nb::str(nb::cast(self.getMax())));
            });
            cls.def("__reduce__", [cls](Box2I const &self) {
                return nb::make_tuple(cls, make_tuple(nb::cast(self.getMin()), nb::cast(self.getMax())));
            });
            auto getSlices = [](Box2I const &self) {
                return nb::make_tuple(nb::slice(self.getBeginY(), self.getEndY(), 1),
                                      nb::slice(self.getBeginX(), self.getEndX(), 1));
            };
            cls.def("getSlices", getSlices);
            cls.def_prop_ro("slices", getSlices);

            mod.attr("BoxI") = cls;
        }
    );

    wrappers.wrapType(
        nb::class_<Box2D>(wrappers.module, "Box2D"),
        [](auto & mod, auto & cls) mutable {

            cls.attr("Point") = mod.attr("Point2D");
            cls.attr("Extent") = mod.attr("Extent2D");

            cls.attr("EPSILON") = nb::float_(Box2D::EPSILON);
            cls.attr("INVALID") = nb::float_(Box2D::INVALID);

            cls.def(nb::init<>());
            cls.def(nb::init<Point2D const &, Point2D const &, bool>(), "minimum"_a, "maximum"_a,
                    "invert"_a = true);
            cls.def(nb::init<Point2D const &, Extent2D const &, bool>(), "corner"_a, "dimensions"_a,
                    "invert"_a = true);
            cls.def(nb::init<IntervalD const &, IntervalD const &>(), "x"_a, "y"_a);
            cls.def(nb::init<Box2I const &>());
            cls.def(nb::init<Box2D const &>());

            cls.def("__eq__", [](Box2D const &self, Box2D const &other) { return self == other; },
                    nb::is_operator());
            cls.def("__ne__", [](Box2D const &self, Box2D const &other) { return self != other; },
                    nb::is_operator());

            cls.def_static("makeCenteredBox", &Box2D::makeCenteredBox, "center"_a, "size"_a);
            cls.def("swap", &Box2D::swap);
            cls.def("getMin", &Box2D::getMin);
            cls.def("getMinX", &Box2D::getMinX);
            cls.def("getMinY", &Box2D::getMinY);
            cls.def("getMax", &Box2D::getMax);
            cls.def("getMaxX", &Box2D::getMaxX);
            cls.def("getMaxY", &Box2D::getMaxY);
            cls.def_prop_ro("minX", &Box2D::getMinX);
            cls.def_prop_ro("minY", &Box2D::getMinY);
            cls.def_prop_ro("maxX", &Box2D::getMaxX);
            cls.def_prop_ro("maxY", &Box2D::getMaxY);
            cls.def("getDimensions", &Box2D::getDimensions);
            cls.def("getWidth", &Box2D::getWidth);
            cls.def("getHeight", &Box2D::getHeight);
            cls.def("getArea", &Box2D::getArea);
            cls.def_prop_ro("width", &Box2D::getWidth);
            cls.def_prop_ro("height", &Box2D::getHeight);
            cls.def_prop_ro("area", &Box2D::getArea);
            cls.def("getX", &Box2D::getX);
            cls.def("getY", &Box2D::getY);
            cls.def_prop_ro("x", &Box2D::getX);
            cls.def_prop_ro("y", &Box2D::getY);
            cls.def("getCenter", &Box2D::getCenter);
            cls.def("getCenterX", &Box2D::getCenterX);
            cls.def("getCenterY", &Box2D::getCenterY);
            cls.def_prop_ro("centerX", &Box2D::getCenterX);
            cls.def_prop_ro("centerY", &Box2D::getCenterY);
            cls.def("isEmpty", &Box2D::isEmpty);
            cls.def("contains", nb::overload_cast<Point2D const &>(&Box2D::contains, nb::const_));
            cls.def("contains", nb::overload_cast<Box2D const &>(&Box2D::contains, nb::const_));
            cls.def("contains",
                    nb::vectorize(static_cast<bool (Box2D::*)(double x, double y) const>(&Box2D::contains)),
                    "x"_a, "y"_a);
            cls.def("__contains__", nb::overload_cast<Point2D const &>(&Box2D::contains, nb::const_));
            cls.def("__contains__", nb::overload_cast<Box2D const &>(&Box2D::contains, nb::const_));
            cls.def("intersects", &Box2D::intersects);
            cls.def("isDisjointFrom", &Box2D::isDisjointFrom);
            cls.def("overlaps", &Box2D::overlaps);
            cls.def("grow", nb::overload_cast<double>(&Box2D::grow));
            cls.def("grow", nb::overload_cast<Extent2D const &>(&Box2D::grow));
            cls.def("shift", &Box2D::shift);
            cls.def("flipLR", &Box2D::flipLR);
            cls.def("flipTB", &Box2D::flipTB);
            cls.def("include", nb::overload_cast<Point2D const &>(&Box2D::include));
            cls.def("include", nb::overload_cast<Box2D const &>(&Box2D::include));
            cls.def("clip", &Box2D::clip);
            cls.def("dilatedBy", nb::overload_cast<double>(&Box2D::dilatedBy, nb::const_));
            cls.def("dilatedBy", nb::overload_cast<Extent2D const &>(&Box2D::dilatedBy, nb::const_));
            cls.def("erodedBy", nb::overload_cast<double>(&Box2D::erodedBy, nb::const_));
            cls.def("erodedBy", nb::overload_cast<Extent2D const &>(&Box2D::erodedBy, nb::const_));
            cls.def("shiftedBy", &Box2D::shiftedBy);
            cls.def("reflectedAboutX", &Box2D::reflectedAboutX);
            cls.def("reflectedAboutY", &Box2D::reflectedAboutY);
            cls.def("expandedTo", nb::overload_cast<Point2D const &>(&Box2D::expandedTo, nb::const_));
            cls.def("expandedTo", nb::overload_cast<Box2D const &>(&Box2D::expandedTo, nb::const_));
            cls.def("clippedTo", &Box2D::clippedTo);
            cls.def("getCorners", &Box2D::getCorners);
            cls.def("toString", &Box2D::toString);
            cls.def("__repr__", [](Box2D const &self) {
                return nb::str("Box2D(corner={}, dimensions={})")
                    .format(nb::repr(nb::cast(self.getMin())), nb::repr(nb::cast(self.getDimensions())));
            });
            cls.def("__str__", [](Box2D const &self) {
                return nb::str("(minimum={}, maximum={})")
                    .format(nb::str(nb::cast(self.getMin())), nb::str(nb::cast(self.getMax())));
            });
            cls.def("__reduce__", [cls](Box2D const &self) {
                return nb::make_tuple(cls, make_tuple(nb::cast(self.getMin()), nb::cast(self.getMax())));
            });

            mod.attr("BoxD") = cls;
        }
    );
}

}  // namespace geom
}  // namespace lsst
