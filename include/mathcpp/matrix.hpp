#pragma once
#include <cstring>

#include "vector.hpp"

namespace mathcpp {

/*
3d rotation standart: X:F,T; Y:R,F; Z:R,T
X,Y,Z means the rotation axis, and 2 second letter represent local X and Y
vectors for rotation in 2d space R-RightVec,T-TopVec,F-FrontVec
*/
constexpr size_t RotationStandart3D[6] = {2, 1, 0, 2, 0, 1};

template <size_t SizeX, size_t SizeY, typename T>
class Matrix {
    T Nums[SizeX * SizeY];

public:
    template <size_t, size_t, typename>
    friend class Matrix;

    constexpr Matrix(
        std::initializer_list<Vector<SizeY, T>> const& vecs) noexcept
        : Matrix(vecs, std::make_integer_sequence<size_t, SizeX>{}) {}

private:
    constexpr T _MatrixConstrFromVectorsHelper(
        std::initializer_list<Vector<SizeY, T>> const& vecs,
        size_t ind) noexcept {
        return vecs[ind / SizeY][ind % SizeY];
    }
    template <size_t... Inds>
    constexpr Matrix(std::initializer_list<Vector<SizeY, T>> const& vecs,
                     std::integer_sequence<size_t, Inds...>) noexcept
        : Nums{_MatrixConstrFromVectorsHelper(vecs, Inds)...} {}

public:
    constexpr Matrix(std::initializer_list<T> const& nums) noexcept
        : Matrix(nums, std::make_integer_sequence<size_t, SizeX * SizeY>{}) {}

private:
    template <size_t... Inds>
    constexpr Matrix(std::initializer_list<T> const& nums,
                     std::integer_sequence<size_t, Inds...>) noexcept
        : Nums{nums.begin()[Inds]...} {}

public:
    constexpr Matrix(const T& num) noexcept {
        for (size_t i = 0; i < SizeX * SizeY; i++) Nums[i] = num;
    }
    // be aware that this constructor will not call move contructors for T, it
    // will use copy constructors
    constexpr Matrix(T const* const numsArrPtr) noexcept {
        for (size_t i = 0; i < SizeX * SizeY; i++) Nums[i] = numsArrPtr[i];
    }
    constexpr Matrix() noexcept {}

    ~Matrix() = default;

    constexpr Matrix(const Matrix<SizeX, SizeY, T>& matToCopy) noexcept {
        for (size_t i = 0; i < SizeX * SizeY; i++) Nums[i] = matToCopy.Nums[i];
    }
    constexpr Matrix(Matrix<SizeX, SizeY, T>&& matToCopy) noexcept {
        for (size_t i = 0; i < SizeX * SizeY; i++)
            Nums[i] = std::move(matToCopy.Nums[i]);
    }
    constexpr Matrix<SizeX, SizeY, T>& operator=(
        const Matrix<SizeX, SizeY, T>& matToCopy) noexcept {
        this->~Matrix();
        new (this) Matrix(matToCopy);
        return *this;
    }
    constexpr Matrix<SizeX, SizeY, T>& operator=(
        Matrix<SizeX, SizeY, T>&& matToCopy) noexcept {
        this->~Matrix();
        new (this) Matrix(std::move(matToCopy));
        return *this;
    }

    constexpr T& operator[](const size_t ind) noexcept { return Nums[ind]; }
    constexpr const T& operator[](const size_t ind) const noexcept {
        return Nums[ind];
    }

    constexpr T GetDeterminant() const noexcept {
        static_assert(SizeX == SizeY,
                      "Determinant is only defined for square matrix");

        if constexpr (SizeX == 1) {
            return Nums[0];
        } else {
            T tempNums[(SizeX - 1) * (SizeY - 1)] = {T()};
            int multiplier = 1;
            T det = 0;
            for (size_t ai = 0; ai < SizeX; ai++) {
                for (size_t x = 0; x < SizeX; x++) {
                    if (x == ai) continue;
                    for (size_t y = 1; y < SizeY; y++) {
                        size_t insertInd =
                            y - 1 + (SizeY - 1) * (x > ai ? (x - 1) : x);
                        tempNums[insertInd] = Nums[y + x * SizeY];
                    }
                }
                det += Nums[ai * SizeY] * multiplier *
                       Matrix<SizeX - 1, SizeX - 1, T>(&tempNums[0])
                           .GetDeterminant();
                multiplier = -multiplier;
            }
            return det;
        }
    }

    constexpr Matrix<SizeX, SizeY, T> GetInversedMatrix(
        const T& det) const noexcept {
        static_assert(SizeX == SizeY,
                      "Inversed matrix is only defined for square matrix "
                      "becouse of determinant");

        Matrix<SizeX, SizeY, T> retMat(0.f);

        T tempNums[(SizeX - 1) * (SizeY - 1)] = {0};
        for (size_t x = 0; x < SizeX; x++) {
            for (size_t y = 0; y < SizeY; y++) {
                for (size_t lx = 0; lx < SizeX; lx++) {
                    if (lx == y) continue;
                    for (size_t ly = 0; ly < SizeY; ly++) {
                        if (ly == x) continue;
                        tempNums[((lx > y) ? (lx - 1) : lx) * (SizeX - 1) +
                                 ((ly > x) ? (ly - 1) : ly)] =
                            Nums[ly + lx * SizeY];
                    }
                }
                int multipler = ((x + y) % 2 == 0) ? 1 : -1;
                retMat.Nums[y + x * SizeY] =
                    multipler *
                    Matrix<SizeX - 1, SizeY - 1, T>(tempNums).GetDeterminant() /
                    det;
                multipler = -multipler;
            }
        }

        return retMat;
    }

    // used to normalize whole matrix, since after some rotations matrix can get
    // a little broken becouse of Ting point precision error
    void Normalize() noexcept {
        for (size_t x = 0; x < SizeX; x++) {
            Vector<SizeY, T> v(&Nums[x * SizeY]);
            v = v.Normalize();
            for (size_t y = 0; y < SizeY; y++) Nums[x * SizeY + y] = v[y];
        }
    }

    // used to fix the problem when matrix becomes not Cartesian after some
    // rotations no operations regarding length of vectors will be done! axis
    // are defined as folows X=0,Y=1,Z=2 the way cross fix works is defined by
    // order of axis, first axis WONT change, the cross will be taken between
    // axis1 and axis2 so this will be axis3, and after that axis2 will be a
    // cross product between axis1 and axis3
    template <size_t axis1, size_t axis2>
    constexpr void CrossFix3D() noexcept {
        static_assert(SizeX == 3 and SizeY == 3,
                      "Matrix is not 3x3, cross fix is not defined");
        static_assert(axis1 != axis2,
                      "You cant have repeating axis rotation order");
        static_assert(axis1 >= 0 and axis1 <= 3 and axis2 >= 0 and axis2 <= 3,
                      "Incorrect axis set up");

        Vector<3, T> vecs[3] = {{&Nums[axis1 * 3]}, {&Nums[axis2 * 3]}, {0.f}};

        vecs[2] = vecs[0].CrossFL(vecs[1]);
        vecs[1] = vecs[2].CrossFL(vecs[0]);

        size_t axis3 = RotationStandart3D[axis1 * 2];
        if (axis3 == axis2) axis3 = RotationStandart3D[axis1 * 2 + 1];

        for (size_t i = 0; i < SizeX * SizeY; i++)
            Nums[i] = vecs[i / 3u][i % 3];
    }

    // all vectors are supposed to be length of 1
    // works only for Cartesian coordinates system
    template <size_t axisInd>
    constexpr T GetLocalCordForAxisC(
        const Vector<SizeY, T>& vec) const noexcept {
        static_assert(SizeX == SizeY, "Matrix should be square");
        static_assert(axisInd >= 0 and axisInd < SizeX, "Invalid axis");

        return Vector<SizeY, T>(&Nums[SizeY * axisInd]).DotFL(vec);
    }

    // all vectors are supposed to be length of 1
    // works only for Cartesian coordinates system
    template <size_t axisInd>
    constexpr Vector<SizeX, T> GetLocalCordsForAxisC(
        const Vector<SizeY, T>& vec) const noexcept {
        static_assert(SizeX == SizeY, "Matrix should be square");
        static_assert(axisInd >= 0 and axisInd < SizeX, "Invalid axis");

        T nums[SizeX] = {T()};
        for (size_t i = 0; i < SizeX; i++) nums[i] = GetLocalCordForAxisC(vec);

        return Vector<SizeX, T>(nums);
    }

    // cant be used for just 1 axis at a time, have to deal with whole matrix,
    // this will work with any gived matrix space but it comes at a cost... U
    // stands for Universal, meaning that it can work in any coordinates system,
    // not only Cartesian
    constexpr Vector<SizeX, T> GetLocalCordsForAxisU(
        const Vector<SizeY, T>& vec) const noexcept {
        static_assert(
            SizeX == SizeY,
            "Getting local cords when matrix is not a square is not defined");

        return (Vector<SizeX, T>)(GetInversedMatrix() * vec);
    }

    template <size_t axisForX, size_t axisForY>
    Vector<SizeY, T> RotateVectorByTwoVectorsU(
        const Vector<SizeY, T>& vec, const float angle) const noexcept {
        return RotateVectorByTwoVectorsU<axisForX, axisForY>(
            vec, angle, this->GetInversedMatrix(this->GetDeterminant()));
    }
    template <size_t axisForX, size_t axisForY>
    Vector<SizeY, T> RotateVectorByTwoVectorsU(const Vector<SizeY, T>& vec,
                                               const float angle,
                                               const T det) const noexcept {
        return RotateVectorByTwoVectorsU<axisForX, axisForY>(
            vec, angle, this->GetInversedMatrix(det));
    }
    template <size_t axisForX, size_t axisForY>
    Vector<SizeY, T> RotateVectorByTwoVectorsU(
        const Vector<SizeY, T>& vec, const float angle,
        const Matrix<SizeX, SizeY, T>& invMat) const noexcept {
        static_assert(
            SizeX == SizeY,
            "Matrix should be square, otherwise it dosent make sence");
        static_assert(axisForX < SizeX and axisForX >= 0 and
                          axisForY < SizeX and axisForY >= 0,
                      "Invalid axis indexes");
        static_assert(axisForX != axisForY, "Axes cant be same");

        Matrix<SizeX, SizeY, T> RotMat;
        for (size_t x = 0; x < SizeX; x++)
            if (x != axisForX && x != axisForY) RotMat[x * SizeY + x] = 1.f;

        RotMat[axisForX * SizeY + axisForX] = cosf(angle);
        RotMat[axisForX * SizeY + axisForY] = sinf(angle);
        RotMat[axisForY * SizeY + axisForX] = -sinf(angle);
        RotMat[axisForY * SizeY + axisForY] = cosf(angle);

        return (*this) * RotMat * invMat * vec;
    }

    //"C" stands for Cartesian coordinates system
    template <size_t axisX, size_t axisY>
    Vector<SizeY, T> RotateVectorC(const Vector<SizeY, T>& vec,
                                   const T angle) const noexcept {
        static_assert(
            axisX < SizeX and axisX >= 0 and axisY < SizeX and axisY >= 0,
            "Invalid axis indexes");
        static_assert(axisX != axisY, "Axes cant be same");

        T vecLen = vec.Length();

        T localCordsX = GetLocalCordForAxisC<axisX>(vec) / vecLen;
        if (localCordsX < -1)
            localCordsX = -1;
        else if (localCordsX > 1)
            localCordsX = 1;
        T localCordsY = GetLocalCordForAxisC<axisY>(vec) / vecLen;

        T angleToX = acosf(localCordsX);
        // check if angle is negative
        if (localCordsY < 0) angleToX = -angleToX;

        angleToX += angle;

        Vector<SizeY, T> xv(&Nums[SizeY * axisX]);
        Vector<SizeY, T> yv(&Nums[SizeY * axisY]);

        return (xv * cosf(angleToX) + yv * sinf(angleToX)) * vecLen;
    }

    template <size_t axis1, size_t axis2, size_t axis3>
    Matrix<SizeX, SizeY, T> RotateIn3DByAnglesC(
        const Vector3F& rot) const noexcept {
        return RotateIn3DByAnglesC<axis1, axis2, axis3>(rot[0], rot[1], rot[2]);
    }
    // axis1,axis2,axis3 are representing the order of rotation, X=0,Y=1,Z=2, so
    // for example rotation by XYZ will be 0,1,2 will work only for Cartesian
    // coordinate system(the "C" at end means Cartesian), it means that angle
    // between all axes is pi/2 axis vectors are supposed to be length of 1
    template <size_t axis1, size_t axis2, size_t axis3>
    Matrix<SizeX, SizeY, T> RotateIn3DByAnglesC(const float xr, const float yr,
                                                const float zr) const noexcept {
        static_assert(SizeX == 3 and SizeY == 3,
                      "Matrix is not 3x3, rotation is not defined");
        static_assert(axis1 != axis2 and axis1 != axis3 and axis2 != axis3,
                      "You cant have repeating axis rotation order");
        static_assert(axis1 >= 0 and axis1 <= 3 and axis2 >= 0 and
                          axis2 <= 3 and axis3 >= 0 and axis3 <= 3,
                      "Incorrect axis set up");

        Matrix<SizeX, SizeY, T> retMat(*this);

        const float rots[3] = {xr, yr, zr};
        constexpr size_t order[3] = {axis1, axis2, axis3};

#define rotVecMacr(i)                                                   \
    {                                                                   \
        constexpr size_t xvi = RotationStandart3D[order[i] * 2];        \
        constexpr size_t yvi = RotationStandart3D[order[i] * 2 + 1];    \
        Vector<3, T> xv(&retMat.Nums[SizeY * xvi]);                     \
        Vector<3, T> yv(&retMat.Nums[SizeY * yvi]);                     \
        Vector<3, T> nxv =                                              \
            xv * cosf(rots[order[i]]) + yv * sinf(rots[order[i]]);      \
        Vector<3, T> nyv =                                              \
            xv * -sinf(rots[order[i]]) + yv * cosf(rots[order[i]]);     \
        std::memcpy(&retMat.Nums[xvi * SizeX], &nxv[0], sizeof(T) * 3); \
        std::memcpy(&retMat.Nums[yvi * SizeX], &nyv[0], sizeof(T) * 3); \
    }

        rotVecMacr(0);
        rotVecMacr(1);
        rotVecMacr(2);

#undef rotVecMacr

        return retMat;
    }

    constexpr operator Vector<SizeX * SizeY, T>() noexcept {
        return Vector<SizeX * SizeY, T>(&Nums[0]);
    }

    constexpr Matrix<SizeX, SizeY, T> operator+(const T& num) const noexcept {
        Matrix<SizeX, SizeY, T> retMat;
        for (size_t i = 0; i < SizeX * SizeY; i++)
            retMat.Nums[i] = Nums[i] + num;
        return retMat;
    }
    constexpr Matrix<SizeX, SizeY, T> operator-(const T& num) const noexcept {
        Matrix<SizeX, SizeY, T> retMat;
        for (size_t i = 0; i < SizeX * SizeY; i++)
            retMat.Nums[i] = Nums[i] - num;
        return retMat;
    }
    constexpr Matrix<SizeX, SizeY, T> operator*(const T& num) const noexcept {
        Matrix<SizeX, SizeY, T> retMat;
        for (size_t i = 0; i < SizeX * SizeY; i++)
            retMat.Nums[i] = Nums[i] * num;
        return retMat;
    }
    constexpr Matrix<SizeX, SizeY, T> operator/(const T& num) const noexcept {
        Matrix<SizeX, SizeY, T> retMat;
        for (size_t i = 0; i < SizeX * SizeY; i++)
            retMat.Nums[i] = Nums[i] / num;
        return retMat;
    }

    template <size_t SizeX2>
    constexpr Matrix<SizeX2, SizeY, T> operator*(
        const Matrix<SizeX2, SizeX, T>& mat) const noexcept {
        Matrix<SizeX2, SizeY, T> retMat;
        for (size_t y = 0; y < SizeY; y++) {
            for (size_t x = 0; x < SizeX2; x++) {
                T sum = 0.f;
                for (size_t lo = 0; lo < SizeX; lo++) {
                    sum += Nums[SizeY * lo + y] * mat.Nums[SizeX * x + lo];
                }
                retMat[x * SizeY + y] = sum;
            }
        }
        return retMat;
    }

    constexpr Vector<SizeY, T> operator*(
        const Vector<SizeY, T>& vec) const noexcept {
        return (*this) * Matrix<1, SizeY, T>(&vec[0]);
    }
};

typedef Matrix<1, 1, float> Matrix11F;
typedef Matrix<1, 2, float> Matrix12F;
typedef Matrix<1, 3, float> Matrix13F;
typedef Matrix<1, 4, float> Matrix14F;
typedef Matrix<2, 1, float> Matrix21F;
typedef Matrix<2, 2, float> Matrix22F;
typedef Matrix<2, 3, float> Matrix23F;
typedef Matrix<2, 4, float> Matrix24F;
typedef Matrix<3, 1, float> Matrix31F;
typedef Matrix<3, 2, float> Matrix32F;
typedef Matrix<3, 3, float> Matrix33F;
typedef Matrix<3, 4, float> Matrix34F;
typedef Matrix<4, 1, float> Matrix41F;
typedef Matrix<4, 2, float> Matrix42F;
typedef Matrix<4, 3, float> Matrix43F;
typedef Matrix<4, 4, float> Matrix44F;

typedef Matrix<1, 1, double> Matrix11D;
typedef Matrix<1, 2, double> Matrix12D;
typedef Matrix<1, 3, double> Matrix13D;
typedef Matrix<1, 4, double> Matrix14D;
typedef Matrix<2, 1, double> Matrix21D;
typedef Matrix<2, 2, double> Matrix22D;
typedef Matrix<2, 3, double> Matrix23D;
typedef Matrix<2, 4, double> Matrix24D;
typedef Matrix<3, 1, double> Matrix31D;
typedef Matrix<3, 2, double> Matrix32D;
typedef Matrix<3, 3, double> Matrix33D;
typedef Matrix<3, 4, double> Matrix34D;
typedef Matrix<4, 1, double> Matrix41D;
typedef Matrix<4, 2, double> Matrix42D;
typedef Matrix<4, 3, double> Matrix43D;
typedef Matrix<4, 4, double> Matrix44D;

constexpr inline Matrix33F Mat3D_RotBaseF =
    Matrix33F{1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
constexpr inline Matrix33D Mat3D_RotBaseD =
    Matrix33D{1., 0., 0., 0., 1., 0., 0., 0., 1.};
constexpr inline Matrix22F Mat2D_RotBaseF = Matrix22F{1.f, 0.f, 0.f, 1.f};
constexpr inline Matrix22D Mat2D_RotBaseD = Matrix22D{1., 0., 0., 1.};

}  // namespace mathcpp
