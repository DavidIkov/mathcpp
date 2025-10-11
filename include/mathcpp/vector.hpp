#pragma once
#include <algorithm>
#include <initializer_list>
#include <utility>

#include "math.h"

namespace mathcpp {

template <size_t AxesAmount, typename T>
class Vector {
    T Axes[AxesAmount];

public:
    template <size_t, typename>
    friend class Vector;

    constexpr Vector(const T &num) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] = num;
    }
    constexpr Vector(T const *numsArrPtr) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] = numsArrPtr[i];
    }
    constexpr Vector(std::initializer_list<T> const &nums) noexcept
        : Vector(nums, std::make_integer_sequence<size_t, AxesAmount>{}) {}

private:
    template <size_t... Inds>
    constexpr Vector(std::initializer_list<T> const &nums,
                     std::integer_sequence<size_t, Inds...>) noexcept
        : Axes{nums.begin()[Inds]...} {}

public:
    constexpr Vector() noexcept {}

    ~Vector() = default;

    constexpr Vector(const Vector<AxesAmount, T> &vecToCopy) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] = vecToCopy[i];
    }
    constexpr Vector(Vector<AxesAmount, T> &&vecToCopy) noexcept {
        for (size_t i = 0; i < AxesAmount; i++)
            Axes[i] = std::move(vecToCopy[i]);
    }
    constexpr Vector<AxesAmount, T> &operator=(
        const Vector<AxesAmount, T> &vecToCopy) noexcept {
        this->~Vector();
        new (this) Vector(vecToCopy);
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator=(
        Vector<AxesAmount, T> &&vecToCopy) noexcept {
        this->~Vector();
        new (this) Vector(std::move(vecToCopy));
        return *this;
    }

    constexpr T &operator[](const size_t ind) noexcept { return Axes[ind]; }
    constexpr const T &operator[](const size_t ind) const noexcept {
        return Axes[ind];
    }

    // length but squared, square root didnt happen
    constexpr T LengthSQ() const noexcept {
        T sum = T();
        for (size_t i = 0; i < AxesAmount; i++) sum += Axes[i] * Axes[i];
        return sum;
    }
    T Length() const noexcept { return sqrtf(LengthSQ()); }

    constexpr Vector<AxesAmount, T> Normalize() const noexcept {
        return (*this) / Length();
    }
    // dot product with full length, means that it wasnt divided by length of
    // both vectors
    constexpr T DotFL(const Vector<AxesAmount, T> &vec) const noexcept {
        T sum = T();
        for (size_t i = 0; i < AxesAmount; i++) sum += vec.Axes[i] * Axes[i];
        return sum;
    }
    T Dot(const Vector<AxesAmount, T> &vec) const noexcept {
        return DotFL(vec) / vec.Length() / Length();
    }

    // first value is dist from p1, and second one is from p2
    constexpr static Vector<2, T> GetDistsToShortestDistBetweenVecs(
        const Vector<AxesAmount, T> &p1, const Vector<AxesAmount, T> &v1,
        const Vector<AxesAmount, T> &p2,
        const Vector<AxesAmount, T> &v2) noexcept {
        Vector<AxesAmount, T> diff = p1 - p2;
        T v1l = v1.Length();
        T v2l = v2.Length();

        // variables for system like
        //  {ax+by+c=0
        //  {dx+ey+f=0
        T a = v1.DotFL(v2);
        T b = -v2l * v2l;
        T c = diff.DotFL(v2);
        T d = v1l * v1l;
        T e = -v2.DotFL(v1);
        T f = diff.DotFL(v1);

        return Vector<2, T>((b * f - e * c) / (a * e - d * b),
                            (d * c - f * a) / (a * e - d * b));
    };
    constexpr static T GetShortestDistBetweenVecs(
        const Vector<AxesAmount, T> &p1, const Vector<AxesAmount, T> &v1,
        const Vector<AxesAmount, T> &p2,
        const Vector<AxesAmount, T> &v2) noexcept {
        Vector<2, T> dists = GetDistsToStortestDistBetweenVecs(p1, v1, p2, v2);

        Vector<AxesAmount, T> mp1 = p1 + dists[0] * v1;
        Vector<AxesAmount, T> mp2 = p2 + dists[1] * v2;

        return (mp1 - mp2).Length();
    };

    // cross product of full length, means that it wasnt deivided by length of
    // both vectors
    constexpr Vector<AxesAmount, T> CrossFL(
        const Vector<AxesAmount, T> &vec) const noexcept {
        static_assert(AxesAmount == 3,
                      "Cross is defined only for 3-axes vectors!");
        return Vector<AxesAmount, T>(
            Axes[1] * vec.Axes[2] - Axes[2] * vec.Axes[1],
            Axes[2] * vec.Axes[0] - Axes[0] * vec.Axes[2],
            Axes[0] * vec.Axes[1] - Axes[1] * vec.Axes[0]);
    }
    constexpr Vector<AxesAmount, T> Cross(
        const Vector<AxesAmount, T> &vec) const noexcept {
        return CrossFL(vec) / Length() / vec.Length();
    }

    constexpr Vector<AxesAmount, T> Clamp(const Vector<AxesAmount, T> &low,
                                          const Vector<AxesAmount, T> &high) {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++)
            retVec.Axes[i] = std::clamp(Axes[i], low[i], high[i]);
        return retVec;
    }

    constexpr Vector<AxesAmount, T> operator+(const T &num) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++) retVec.Axes[i] = Axes[i] + num;
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator-(const T &num) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++) retVec.Axes[i] = Axes[i] - num;
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator*(const T &num) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++) retVec.Axes[i] = Axes[i] * num;
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator/(const T &num) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++) retVec.Axes[i] = Axes[i] / num;
        return retVec;
    }

    constexpr Vector<AxesAmount, T> operator+(
        const Vector<AxesAmount, T> &vec) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++)
            retVec.Axes[i] = Axes[i] + vec.Axes[i];
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator-(
        const Vector<AxesAmount, T> &vec) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++)
            retVec.Axes[i] = Axes[i] - vec.Axes[i];
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator*(
        const Vector<AxesAmount, T> &vec) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++)
            retVec.Axes[i] = Axes[i] * vec.Axes[i];
        return retVec;
    }
    constexpr Vector<AxesAmount, T> operator/(
        const Vector<AxesAmount, T> &vec) const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++)
            retVec.Axes[i] = Axes[i] / vec.Axes[i];
        return retVec;
    }

    constexpr Vector<AxesAmount, T> operator-() const noexcept {
        Vector<AxesAmount, T> retVec;
        for (size_t i = 0; i < AxesAmount; i++) retVec.Axes[i] = -Axes[i];
        return retVec;
    }

    constexpr Vector<AxesAmount, T> &operator+=(const T &num) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] += num;
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator-=(const T &num) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] -= num;
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator*=(const T &num) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] *= num;
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator/=(const T &num) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] /= num;
        return *this;
    }

    constexpr Vector<AxesAmount, T> &operator+=(
        const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] += vec.Axes[i];
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator-=(
        const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] -= vec.Axes[i];
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator*=(
        const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] *= vec.Axes[i];
        return *this;
    }
    constexpr Vector<AxesAmount, T> &operator/=(
        const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++) Axes[i] /= vec.Axes[i];
        return *this;
    }

    template <typename T2>
    constexpr operator Vector<AxesAmount, T2>() const noexcept {
        Vector<AxesAmount, T2> vec;
        for (size_t i = 0; i < AxesAmount; i++) vec.Axes[i] = (T2)Axes[i];
        return vec;
    }

    constexpr bool operator==(const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++)
            if (Axes[i] != vec.Axes[i]) return false;
        return true;
    }
    constexpr bool operator!=(const Vector<AxesAmount, T> &vec) noexcept {
        for (size_t i = 0; i < AxesAmount; i++)
            if (Axes[i] != vec.Axes[i]) return true;
        return false;
    }
};

template <size_t AxesAmount, typename T>
constexpr Vector<AxesAmount, T> operator+(
    const T &num, const Vector<AxesAmount, T> &vec) noexcept {
    Vector<AxesAmount, T> retVec;
    for (size_t i = 0; i < AxesAmount; i++) retVec[i] = num + vec[i];
    return retVec;
}
template <size_t AxesAmount, typename T>
constexpr Vector<AxesAmount, T> operator-(
    const T &num, const Vector<AxesAmount, T> &vec) noexcept {
    Vector<AxesAmount, T> retVec;
    for (size_t i = 0; i < AxesAmount; i++) retVec[i] = num - vec[i];
    return retVec;
}
template <size_t AxesAmount, typename T>
constexpr Vector<AxesAmount, T> operator*(
    const T &num, const Vector<AxesAmount, T> &vec) noexcept {
    Vector<AxesAmount, T> retVec;
    for (size_t i = 0; i < AxesAmount; i++) retVec[i] = num * vec[i];
    return retVec;
}
template <size_t AxesAmount, typename T>
constexpr Vector<AxesAmount, T> operator/(
    const T &num, const Vector<AxesAmount, T> &vec) noexcept {
    Vector<AxesAmount, T> retVec;
    for (size_t i = 0; i < AxesAmount; i++) retVec[i] = num / vec[i];
    return retVec;
}

typedef Vector<1, int> Vector1I;
typedef Vector<2, int> Vector2I;
typedef Vector<3, int> Vector3I;
typedef Vector<4, int> Vector4I;

typedef Vector<1, unsigned> Vector1U;
typedef Vector<2, unsigned> Vector2U;
typedef Vector<3, unsigned> Vector3U;
typedef Vector<4, unsigned> Vector4U;

typedef Vector<1, float> Vector1F;
typedef Vector<2, float> Vector2F;
typedef Vector<3, float> Vector3F;
typedef Vector<4, float> Vector4F;

typedef Vector<1, double> Vector1D;
typedef Vector<2, double> Vector2D;
typedef Vector<3, double> Vector3D;
typedef Vector<4, double> Vector4D;

typedef Vector<1, bool> Vector1B;
typedef Vector<2, bool> Vector2B;
typedef Vector<3, bool> Vector3B;
typedef Vector<4, bool> Vector4B;

constexpr inline Vector2F Vec2D_X_F = Vector2F{1.f, 0.f};
constexpr inline Vector2F Vec2D_Y_F = Vector2F{0.f, 1.f};
constexpr inline Vector2D Vec2D_X_D = Vector2D{1., 0.};
constexpr inline Vector2D Vec2D_Y_D = Vector2D{0., 1.};

constexpr inline Vector3F Vec3D_X_F = Vector3F{1.f, 0.f, 0.f};
constexpr inline Vector3F Vec3D_Y_F = Vector3F{0.f, 1.f, 0.f};
constexpr inline Vector3F Vec3D_Z_F = Vector3F{0.f, 0.f, 1.f};

constexpr inline Vector3D Vec3D_X_D = Vector3D{1., 0., 0.};
constexpr inline Vector3D Vec3D_Y_D = Vector3D{0., 1., 0.};
constexpr inline Vector3D Vec3D_Z_D = Vector3D{0., 0., 1.};

}  // namespace mathcpp
