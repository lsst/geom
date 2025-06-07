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

class SphereTransform {
public:
    using Matrix = Eigen::Matrix<double, 3, 3, Eigen::DontAlign>;

    /** Construct an empty (identity) SphereTransform. */
    SphereTransform() noexcept;

    /** Construct an SphereTransform from an Eigen::Matrix. */
    explicit SphereTransform(Matrix const& matrix) noexcept;

    SphereTransform(SphereTransform const& other) noexcept = default;
    SphereTransform(SphereTransform&& other) noexcept = default;
    ~SphereTransform() noexcept = default;

    static SphereTransform fit_unit_vectors(
            ndarray::Array<double const, 2> const& from, ndarray::Array<double const, 2> const& to,
            ndarray::Array<double const, 1> const& weights = ndarray::Array<double const, 1>());

    Matrix const& getMatrix() const noexcept { return _matrix; }

    SphereTransform operator*(SphereTransform const& other) const noexcept;

    SpherePoint operator()(SpherePoint const& point) const noexcept;

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
