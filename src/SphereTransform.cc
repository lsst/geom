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

#include <iostream>

#include "boost/format.hpp"
#include "ndarray/eigen.h"
#include "Eigen/SVD"
#include "Eigen/LU"

#include "lsst/pex/exceptions.h"
#include "lsst/geom/SphereTransform.h"
#include "lsst/geom/sphgeomUtils.h"

namespace lsst {
namespace geom {

SphereTransform::SphereTransform() noexcept : _matrix(Matrix::Identity()) {}

SphereTransform::SphereTransform(Matrix const& matrix) noexcept : _matrix(matrix) {}

SphereTransform SphereTransform::fit_unit_vectors(ndarray::Array<double const, 2> const& from,
                                                  ndarray::Array<double const, 2> const& to,
                                                  ndarray::Array<double const, 1> const& weights) {
    if (from.getSize<0>() != to.getSize<0>()) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("'from' array first dimension (%d) does not match 'to' array first "
                                         "dimension (%d).") %
                           from.getSize<0>() % to.getSize<0>())
                                  .str());
    }
    if (from.getSize<1>() != 3) {
        throw LSST_EXCEPT(
                pex::exceptions::LengthError,
                (boost::format("'from' array second dimension size must be 3; got %d.") % from.getSize<1>())
                        .str());
    }
    if (to.getSize<1>() != 3) {
        throw LSST_EXCEPT(
                pex::exceptions::LengthError,
                (boost::format("'to' array second dimension size must be 3; got %d.") % to.getSize<1>())
                        .str());
    }
    Eigen::Matrix<double, Eigen::Dynamic, 3> y = ndarray::asEigenMatrix(to);
    Eigen::Matrix<double, Eigen::Dynamic, 3> x = ndarray::asEigenMatrix(from);
    if (!weights.isEmpty()) {
        auto w = ndarray::asEigenMatrix(weights);
        y.array().colwise() *= w.array();
    }
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(x.transpose() * y, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d u_t = svd.matrixU().transpose();
    Matrix result = svd.matrixV() * u_t;
    if (result.determinant() < 0) {
        u_t.row(2) *= -1.0;
        result = svd.matrixV() * u_t;
    }
    return SphereTransform(result);
}

SphereTransform SphereTransform::operator*(SphereTransform const& other) const noexcept {
    return SphereTransform(_matrix * other._matrix);
}

SpherePoint SphereTransform::operator()(SpherePoint const& point) const noexcept {
    return SpherePoint(this->operator()(point.getVector()));
}

sphgeom::UnitVector3d SphereTransform::operator()(sphgeom::UnitVector3d const& vector) const noexcept {
    Eigen::Vector3d v = _matrix * asEigen(vector);
    return sphgeom::UnitVector3d(v.x(), v.y(), v.z());
}

double SphereTransform::applyX(double x, double y, double z) const noexcept {
    return _matrix.row(0).dot(Eigen::Vector3d(x, y, z));
}
double SphereTransform::applyY(double x, double y, double z) const noexcept {
    return _matrix.row(1).dot(Eigen::Vector3d(x, y, z));
}
double SphereTransform::applyZ(double x, double y, double z) const noexcept {
    return _matrix.row(2).dot(Eigen::Vector3d(x, y, z));
}

SphereTransform const SphereTransform::inverted() const noexcept {
    return SphereTransform(_matrix.transpose());
}

}  // namespace geom
}  // namespace lsst
