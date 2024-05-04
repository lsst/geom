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

#include "lsst/cpputils/python.h"

#include "lsst/geom/Angle.h"

namespace nb = nanobind;

namespace lsst {
namespace geom {

using PyAngle = nb::class_<Angle>;
using PyAngleUnit = nb::class_<AngleUnit>;

namespace {

template <typename OtherT>
void declareAngleComparisonOperators(PyAngle& cls) {
    cls.def("__eq__", [](Angle const& self, OtherT const& other) { return self == other; },
            nb::is_operator());
    cls.def("__ne__", [](Angle const& self, OtherT const& other) { return self != other; },
            nb::is_operator());
    cls.def("__le__", [](Angle const& self, OtherT const& other) { return self <= other; },
            nb::is_operator());
    cls.def("__ge__", [](Angle const& self, OtherT const& other) { return self >= other; },
            nb::is_operator());
    cls.def("__lt__", [](Angle const& self, OtherT const& other) { return self < other; }, nb::is_operator());
    cls.def("__gt__", [](Angle const& self, OtherT const& other) { return self > other; }, nb::is_operator());
}

} // anonymous

void wrapAngle(cpputils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        PyAngleUnit(wrappers.module, "AngleUnit"),
        [](auto & mod, auto & cls) mutable {
            cls.def("__eq__", [](AngleUnit const& self, AngleUnit const& other) { return self == other; },
                    nb::is_operator());
            cls.def("__ne__", [](AngleUnit const& self, AngleUnit const& other) { return !(self == other); },
                    nb::is_operator());
            cls.def("_mul", [](AngleUnit const& self, double other) { return other * self; },
                    nb::is_operator());
            cls.def("_rmul", [](AngleUnit const& self, double other) { return other * self; },
                    nb::is_operator());
            mod.attr("radians") = nb::cast(radians);
            mod.attr("degrees") = nb::cast(degrees);
            mod.attr("hours") = nb::cast(hours);
            mod.attr("arcminutes") = nb::cast(arcminutes);
            mod.attr("arcseconds") = nb::cast(arcseconds);
            mod.attr("milliarcseconds") = nb::cast(milliarcseconds);
        }
    );

    wrappers.wrapType(
        PyAngle(wrappers.module, "Angle"),
        [](auto & mod, auto & cls) mutable {
            cls.def(nb::init<double, AngleUnit>(), nb::arg("val"), nb::arg("units") = radians);
            cls.def(nb::init<>());
            declareAngleComparisonOperators<Angle>(cls);
            declareAngleComparisonOperators<double>(cls);
            declareAngleComparisonOperators<int>(cls);
            cls.def("__mul__", [](Angle const& self, double other) { return self * other; },
                    nb::is_operator());
            cls.def("__mul__", [](Angle const& self, int other) { return self * other; },
                    nb::is_operator());
            cls.def("__rmul__", [](Angle const& self, double other) { return self * other; },
                    nb::is_operator());
            cls.def("__rmul__", [](Angle const& self, int other) { return self * other; },
                    nb::is_operator());
            cls.def("__imul__", [](Angle& self, double other) { return self *= other; });
            cls.def("__imul__", [](Angle& self, int other) { return self *= other; });
            cls.def("__add__", [](Angle const& self, Angle const& other) { return self + other; },
                    nb::is_operator());
            cls.def("__sub__", [](Angle const& self, Angle const& other) { return self - other; },
                    nb::is_operator());
            cls.def("__neg__", [](Angle const& self) { return -self; }, nb::is_operator());
            cls.def("__iadd__", [](Angle& self, Angle const& other) { return self += other; });
            cls.def("__isub__", [](Angle& self, Angle const& other) { return self -= other; });
            cls.def("__truediv__", [](Angle const& self, double other) { return self / other; },
                    nb::is_operator());
            // Without an explicit wrapper, Python lets Angle / Angle -> Angle
            cls.def("__truediv__", [](Angle const& self, Angle const& other) {
                throw nb::type_error("unsupported operand type(s) for /: 'Angle' and 'Angle'");
            });
            cls.def("__float__", &Angle::operator double);
            cls.def("__abs__", [](Angle const& self) { return std::abs(self.asRadians()) * radians; });

            cls.def("__reduce__", [cls](Angle const& self) {
                return nb::make_tuple(cls, nb::make_tuple(nb::cast(self.asRadians())));
            });
            cpputils::python::addOutputOp(cls, "__str__");
            cls.def("__repr__", [](Angle const & self) {
                if (std::isfinite(self.asDegrees())) {
                    return nb::str("Angle({:0.17g}, degrees)").format(self.asDegrees());
                } else {
                    return nb::str("Angle(float('{}'), degrees)").format(self.asDegrees());
                }
            });
            cls.def("asAngularUnits", &Angle::asAngularUnits);
            cls.def("asRadians", &Angle::asRadians);
            cls.def("asDegrees", &Angle::asDegrees);
            cls.def("asHours", &Angle::asHours);
            cls.def("asArcminutes", &Angle::asArcminutes);
            cls.def("asArcseconds", &Angle::asArcseconds);
            cls.def("asMilliarcseconds", &Angle::asMilliarcseconds);
            cls.def("wrap", &Angle::wrap);
            cls.def("wrapCtr", &Angle::wrapCtr);
            cls.def("wrapNear", &Angle::wrapNear);
            cls.def("separation", &Angle::separation);
            mod.def("isAngle", isAngle<Angle>);
            mod.def("isAngle", isAngle<double>);
            nb::implicitly_convertible<Angle, sphgeom::Angle>();
            nb::implicitly_convertible<sphgeom::Angle, Angle>();
        }
    );

    wrappers.wrap(
        [](auto & mod) mutable {
            mod.attr("PI") = nb::float_(PI);
            mod.attr("TWOPI") = nb::float_(TWOPI);
            mod.attr("HALFPI") = nb::float_(HALFPI);
            mod.attr("ONE_OVER_PI") = nb::float_(ONE_OVER_PI);
            mod.attr("SQRTPI") = nb::float_(SQRTPI);
            mod.attr("INVSQRTPI") = nb::float_(INVSQRTPI);
            mod.attr("ROOT2") = nb::float_(ROOT2);
            mod.def("degToRad", degToRad);
            mod.def("radToDeg", radToDeg);
            mod.def("radToArcsec", radToArcsec);
            mod.def("radToMas", radToMas);
            mod.def("arcsecToRad", arcsecToRad);
            mod.def("masToRad", masToRad);
        }
    );
}

}  // namespace geom
}  // namespace lsst
