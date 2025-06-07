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

#ifndef LSST_GEOM_SPHERETRANSFORM_H_
#define LSST_GEOM_SPHERETRANSFORM_H_

#include "ndarray.h"
#include "Eigen/Core"

#include "lsst/sphgeom/UnitVector3d.h"
#include "lsst/geom/SpherePoint.h"

namespace lsst {
namespace geom {


/**
 *  A rigid 3-d transform on the sphere.
 */
class SphereTransform {
public:
    using Matrix = Eigen::Matrix<double, 3, 3, Eigen::DontAlign>;

    /** Construct an empty (identity) SphereTransform. */
    SphereTransform() noexcept;

    /**
     *  Construct an SphereTransform from a matrix.
     *
     *  The caller is responsible for ensuring that this is an orthogonal
     *  matrix with positive determinant.
    */
    explicit SphereTransform(Matrix const& matrix) noexcept;

    SphereTransform(SphereTransform const& other) noexcept = default;
    SphereTransform(SphereTransform&& other) noexcept = default;
    ~SphereTransform() noexcept = default;

    /**
     *  Fit for the transform that best maps (in a least-squares sense) one set
     *  of unit vectors onto another.
     *
     *  @param[in] from    Unit vectors the returned transform would take as
     *                     inputs, with shape `(N, 3)`.
     *  @param[in] to      Unit vectors the returned transform would return as
     *                     output, with shape `(N, 3)`.
     *  @param[in] weights Weights for the vectors, with shape `(N,)`.
     */
    static SphereTransform fit_unit_vectors(
            ndarray::Array<double const, 2> const& from, ndarray::Array<double const, 2> const& to,
            ndarray::Array<double const, 1> const& weights = ndarray::Array<double const, 1>());

    /**
     *  Return the matrix representation of the transform.
     *
     *  This is an orthogonal matrix with positive determinant that acts on
     *  3-d unit vectors.
     */
    Matrix const& getMatrix() const noexcept { return _matrix; }

    /**
     *  Compose two transforms.
     *
     *  The returned transform first applies `other` and then `this` (i.e.
     *  transform multiply just like their matrices do).
     */
    SphereTransform operator*(SphereTransform const& other) const noexcept;

    /**
     *  Apply the transform to a `SpherePoint`.
     */
    SpherePoint operator()(SpherePoint const& point) const noexcept;

    /**
     *  Apply the transform to a unit vector.
     */
    sphgeom::UnitVector3d operator()(sphgeom::UnitVector3d const& vector) const noexcept;

    //@{
    /**
     *  Transform a 3-d unit vector given and returned as separate double values.
     *
     *  This interface is intended primarily for use in Python (where it is
     *  vectorized to support NumPy array arguments).
     */
    double applyX(double x, double y, double z) const noexcept;
    double applyY(double x, double y, double z) const noexcept;
    double applyZ(double x, double y, double z) const noexcept;
    //@}

    /** Return the inverse transform. */
    SphereTransform const inverted() const noexcept;

private:
    Matrix _matrix;
};

}  // namespace geom
}  // namespace lsst

#endif /* LSST_GEOM_SPHERETRANSFORM_H_ */
