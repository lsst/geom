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

/*
 * A coordinate class intended to represent absolute positions.
 */
#ifndef LSST_GEOM_POINT_H
#define LSST_GEOM_POINT_H

#include <tuple>

#include "lsst/geom/CoordinateBase.h"
#include "lsst/geom/CoordinateExpr.h"
#include "lsst/geom/Extent.h"

namespace lsst {
namespace geom {

template <typename T, int N>
class PointBase : public CoordinateBase<Point<T, N>, T, N> {
    typedef CoordinateBase<Point<T, N>, T, N> Super;

public:
    PointBase(PointBase const &) = default;
    PointBase(PointBase &&) = default;
    PointBase &operator=(PointBase const &) = default;
    PointBase &operator=(PointBase &&) = default;
    ~PointBase() = default;

    /**
     *  Standard equality comparison.
     *
     *  Returns true iff all(this->eq(other));
     */
    bool operator==(Point<T, N> const &other) const noexcept { return all(this->eq(other)); }

    /**
     *  Standard inequality comparison.
     *
     *  Returns true iff any(this->ne(other));
     */
    bool operator!=(Point<T, N> const &other) const noexcept { return any(this->ne(other)); }

    /**
     *  @name Named vectorized comparison functions
     *
     *  Note that these return CoordinateExpr, not bool.
     *
     *  Unlike most arithmetic and assignment operators, scalar interoperability is provided
     *  for comparisons; expressions like
     *
     *      if (all(point.gt(0))) ...
     *
     *  are both ubiquitous and easy to interpret.
     */
    //@{
    CoordinateExpr<N> eq(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> ne(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> lt(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> le(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> gt(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> ge(Point<T, N> const &other) const noexcept;
    CoordinateExpr<N> eq(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->eq(Point<T, N>(scalar));
    }
    CoordinateExpr<N> ne(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->ne(Point<T, N>(scalar));
    }
    CoordinateExpr<N> lt(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->lt(Point<T, N>(scalar));
    }
    CoordinateExpr<N> le(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->le(Point<T, N>(scalar));
    }
    CoordinateExpr<N> gt(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->gt(Point<T, N>(scalar));
    }
    CoordinateExpr<N> ge(T scalar) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return this->ge(Point<T, N>(scalar));
    }
    //@}

    /**
     *  @name Arithmetic operators
     *
     *  No scalar interoperability is provided for Point arithmetic operations.
     */
    //@{
    Extent<T, N> operator-(Point<T, N> const &other) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return Extent<T, N>(this->_vector - other._vector);
    }
    Point<T, N> operator+(Extent<T, N> const &other) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return Point<T, N>(this->_vector + other.asEigen());
    }
    Point<T, N> operator-(Extent<T, N> const &other) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return Point<T, N>(this->_vector - other.asEigen());
    }
    Point<T, N> &operator+=(Extent<T, N> const &other) noexcept(Super::IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        this->_vector += other.asEigen();
        return static_cast<Point<T, N> &>(*this);
    }
    Point<T, N> &operator-=(Extent<T, N> const &other) noexcept(Super::IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        this->_vector -= other.asEigen();
        return static_cast<Point<T, N> &>(*this);
    }
    //@}

    /// Cast this object to an Extent of the same numeric type and dimensionality.
    Extent<T, N> asExtent() const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        return Extent<T, N>(static_cast<Point<T, N> const &>(*this));
    }

    /// Shift the point by the given offset.
    void shift(Extent<T, N> const &offset) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        this->_vector += offset.asEigen();
    }

    void scale(double factor) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) { this->_vector *= factor; }

    double distanceSquared(PointBase<T, N> const &other) const noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) {
        // the cast to double is lame but Eigen seems to require they be the same type
        return (this->asEigen() - other.asEigen()).squaredNorm();
    }

    std::string toString() const {
        std::stringstream out;
        out << "Point(";
        for (size_t i = 0; i < N; ++i) {
            if (i != 0) {
                out << ",";
            }
            out << (*this)[i];
        }
        out << ")";
        return out.str();
    }

protected:
    explicit PointBase(T val = static_cast<T>(0)) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(val) {}

    template <typename Vector>
    explicit PointBase(Eigen::MatrixBase<Vector> const &vector) : Super(vector) {}
};

/**
 *  A coordinate class intended to represent absolute positions.
 *
 *  See @ref geomOps for mathematical operators on Point.
 */
template <typename T, int N>
class Point : public PointBase<T, N> {
    typedef PointBase<T, N> Super;

public:
    typedef typename Super::EigenVector EigenVector;

    /// Construct a Point with all elements set to the same scalar value.
    explicit Point(T val = static_cast<T>(0)) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(val) {}

    Point(Point const &) = default;
    Point(Point &&) = default;
    ~Point() = default;

    Point &operator=(Point const &) = default;
    Point &operator=(Point &&) = default;

    /**
     *  Explicit converting constructor.
     *
     *  Converting from floating point to integer rounds to the nearest integer instead of truncating.
     *  This ensures that a floating-point pixel coordinate converts to the coordinate of the pixel
     *  it lies on (assuming the floating point origin is the center of the first pixel).
     */
    template <typename U>
    explicit Point(Point<U, N> const &other) noexcept(IS_NOTHROW_CONVERTIBLE<T, U>);

    /// Construct a Point from an Eigen vector.
    explicit Point(EigenVector const &vector) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(vector) {}

    /// Explicit constructor from Extent.
    explicit Point(Extent<T, N> const &other) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(other.asEigen()) {}

    void swap(Point &other) noexcept { this->_swap(other); }
};

/**
 *  A coordinate class intended to represent absolute positions (2-d specialization).
 *
 *  See @ref geomOps for mathematical operators on Point.
 */
template <typename T>
class Point<T, 2> : public PointBase<T, 2> {
    typedef PointBase<T, 2> Super;

public:
    typedef typename Super::EigenVector EigenVector;

    /// Construct a Point with all elements set to the same scalar value.
    explicit Point(T val = static_cast<T>(0)) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(val) {}

    Point(Point const &) = default;
    Point(Point &&) = default;
    ~Point() = default;

    Point &operator=(Point const &) = default;
    Point &operator=(Point &&) = default;

    /**
     *  Explicit converting constructor.
     *
     *  Converting from floating point to integer rounds to the nearest integer instead of truncating.
     *  This ensures that a floating-point pixel coordinate converts to the coordinate of the pixel
     *  it lies on (assuming the floating point origin is the center of the first pixel).
     */
    template <typename U>
    explicit Point(Point<U, 2> const &other) noexcept(IS_NOTHROW_CONVERTIBLE<T, U>);

    template <typename U>
    explicit Point(Extent<U, 2> const &other) noexcept(IS_NOTHROW_CONVERTIBLE<T, U>)
        : Super(Extent<T, 2>(other).asEigen()) {}

    /// Construct a Point from an Eigen vector.
    explicit Point(EigenVector const &vector) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(vector) {}

    /// Explicit constructor from Extent.
    explicit Point(Extent<T, 2> const &other) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(other.asEigen()) {}

    /// Explicit constructor from a pair of doubles.
    explicit Point(T x, T y) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(EigenVector(x, y)) {}

    /// Construct from a two-element array.
    explicit Point(T const xy[2]) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(xy[0], xy[1])) {}

    /// Construct from a std::pair.
    explicit Point(std::pair<T, T> const &xy) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(xy.first, xy.second)) {}

    /// Construct from std::tuple.
    explicit Point(std::tuple<T, T> const &xy) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(std::get<0>(xy), std::get<1>(xy))) {}

    void swap(Point &other) noexcept { this->_swap(other); }
};

/**
 *  A coordinate class intended to represent absolute positions (3-d specialization).
 *
 *  See @ref geomOps for mathematical operators on Point.
 */
template <typename T>
class Point<T, 3> : public PointBase<T, 3> {
    typedef PointBase<T, 3> Super;

public:
    typedef typename Super::EigenVector EigenVector;

    /// Construct a Point with all elements set to the same scalar value.
    explicit Point(T val = static_cast<T>(0)) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(val) {}

    Point(Point const &) = default;
    Point(Point &&) = default;
    ~Point() = default;

    Point &operator=(Point const &) = default;
    Point &operator=(Point &&) = default;

    /**
     *  Explicit converting constructor.
     *
     *  Converting from floating point to integer rounds to the nearest integer instead of truncating.
     *  This ensures that a floating-point pixel coordinate converts to the coordinate of the pixel
     *  it lies on (assuming the floating point origin is the center of the first pixel).
     */
    template <typename U>
    explicit Point(Point<U, 3> const &other) noexcept(IS_NOTHROW_CONVERTIBLE<T, U>);

    /// Construct a Point from an Eigen vector.
    explicit Point(EigenVector const &vector) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE) : Super(vector) {}

    /// Explicit constructor from Extent.
    explicit Point(Extent<T, 3> const &other) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(other.asEigen()) {}

    template <typename U>
    explicit Point(Extent<U, 3> const &other) noexcept(IS_NOTHROW_CONVERTIBLE<T, U>)
        : Super(Extent<T, 3>(other).asEigen()) {}

    /// Explicit constructor from a sequence of doubles.
    explicit Point(T x, T y, T z) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(x, y, z)) {}

    /// Construct from a two-element array.
    explicit Point(T const xyz[3]) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(xyz[0], xyz[1], xyz[2])) {}

    /// Construct from std::tuple.
    explicit Point(std::tuple<T, T, T> const &xyz) noexcept(Super::IS_ELEMENT_NOTHROW_COPYABLE)
            : Super(EigenVector(std::get<0>(xyz), std::get<1>(xyz), std::get<2>(xyz))) {}

    void swap(Point &other) noexcept { this->_swap(other); }
};

// Hash functions
template <typename T, int N>
std::size_t hash_value(Point<T, N> const &point) noexcept;

typedef Point<int, 2> PointI;
typedef Point<int, 2> Point2I;
typedef Point<int, 3> Point3I;
typedef Point<double, 2> PointD;
typedef Point<double, 2> Point2D;
typedef Point<double, 3> Point3D;

template <int N>
Point<double, N> operator+(Point<double, N> const &lhs, Extent<int, N> const &rhs) noexcept {
    return lhs + Extent<double, N>(rhs);
}

template <int N>
Point<double, N> operator+(Extent<int, N> const &rhs, Point<double, N> const &lhs) noexcept {
    return Point<double, N>(lhs) + rhs;
}

template <int N>
Point<double, N> &operator+=(Point<double, N> &lhs, Extent<int, N> const &rhs) noexcept {
    return lhs += Extent<double, N>(rhs);
}

template <int N>
Point<double, N> operator+(Point<int, N> const &lhs, Extent<double, N> const &rhs) noexcept {
    return Point<double, N>(lhs) + rhs;
}

template <int N>
Point<double, N> operator-(Point<double, N> const &lhs, Extent<int, N> const &rhs) noexcept {
    return lhs - Extent<double, N>(rhs);
}

template <int N>
Point<double, N> &operator-=(Point<double, N> &lhs, Extent<int, N> const &rhs) noexcept {
    return lhs -= Extent<double, N>(rhs);
}

template <int N>
Point<double, N> operator-(Point<int, N> const &lhs, Extent<double, N> const &rhs) noexcept {
    return Point<double, N>(lhs) - rhs;
}

template <int N>
Extent<double, N> operator-(Point<double, N> const &lhs, Point<int, N> const &rhs) noexcept {
    return lhs - Point<double, N>(rhs);
}

template <int N>
Extent<double, N> operator-(Point<int, N> const &lhs, Point<double, N> const &rhs) noexcept {
    return Point<double, N>(lhs) - rhs;
}

}  // namespace geom
}  // namespace lsst

namespace std {
template <typename T, int N>
struct hash<lsst::geom::Point<T, N>> {
    using argument_type = lsst::geom::Point<T, N>;
    using result_type = std::size_t;
    result_type operator()(argument_type const &x) const noexcept { return lsst::geom::hash_value(x); }
};
}  // namespace std

#endif
