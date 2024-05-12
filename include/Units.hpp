#pragma once

#include <cmath>
#include <ratio>

template <typename MassDim, typename LengthDim, typename TimeDim, typename AngleDim>
class RQuantity {
  public:
  explicit constexpr RQuantity() : value(0.0) {
  }

  explicit constexpr RQuantity(double val) : value(val) {
  }

  explicit constexpr RQuantity(long double val) : value(static_cast<double>(val)) {
  }

  // The intrinsic operations for a quantity with a unit is addition and subtraction
  constexpr RQuantity const &operator+=(const RQuantity &rhs) {
    value += rhs.value;
    return *this;
  }

  constexpr RQuantity const &operator-=(const RQuantity &rhs) {
    value -= rhs.value;
    return *this;
  }

  constexpr RQuantity operator-() {
    return RQuantity(value * -1);
  }

  constexpr RQuantity const &operator*=(const double rhs) {
    value *= rhs;
    return *this;
  }

  constexpr RQuantity const &operator/=(const double rhs) {
    value /= rhs;
    return *this;
  }

  // Returns the value of the quantity in multiples of the specified unit
  constexpr double convert(const RQuantity &rhs) const {
    return value / rhs.value;
  }

  // returns the raw value of the quantity (should not be used)
  constexpr double getValue() const {
    return value;
  }

  constexpr RQuantity<MassDim, LengthDim, TimeDim, AngleDim> abs() const {
    return RQuantity<MassDim, LengthDim, TimeDim, AngleDim>(std::fabs(value));
  }

  constexpr RQuantity<std::ratio_divide<MassDim, std::ratio<2>>,
                      std::ratio_divide<LengthDim, std::ratio<2>>,
                      std::ratio_divide<TimeDim, std::ratio<2>>,
                      std::ratio_divide<AngleDim, std::ratio<2>>>
  sqrt() const {
    return RQuantity<std::ratio_divide<MassDim, std::ratio<2>>,
                     std::ratio_divide<LengthDim, std::ratio<2>>,
                     std::ratio_divide<TimeDim, std::ratio<2>>,
                     std::ratio_divide<AngleDim, std::ratio<2>>>(std::sqrt(value));
  }

  private:
  double value;
};

// Predefined (physical unit) quantity types:
// ------------------------------------------
#define QUANTITY_TYPE(_Mdim, _Ldim, _Tdim, _Adim, name)                                            \
  typedef RQuantity<std::ratio<_Mdim>, std::ratio<_Ldim>, std::ratio<_Tdim>, std::ratio<_Adim>>    \
    name;

// Unitless
QUANTITY_TYPE(0, 0, 0, 0, Number)
constexpr Number number(1.0);

// Standard arithmetic operators:
// ------------------------------
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> operator+(const RQuantity<M, L, T, A> &lhs,
                                          const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(lhs.getValue() + rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> operator-(const RQuantity<M, L, T, A> &lhs,
                                          const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(lhs.getValue() - rhs.getValue());
}
template <typename M1,
          typename L1,
          typename T1,
          typename A1,
          typename M2,
          typename L2,
          typename T2,
          typename A2>
constexpr RQuantity<std::ratio_add<M1, M2>,
                    std::ratio_add<L1, L2>,
                    std::ratio_add<T1, T2>,
                    std::ratio_add<A1, A2>>
operator*(const RQuantity<M1, L1, T1, A1> &lhs, const RQuantity<M2, L2, T2, A2> &rhs) {
  return RQuantity<std::ratio_add<M1, M2>,
                   std::ratio_add<L1, L2>,
                   std::ratio_add<T1, T2>,
                   std::ratio_add<A1, A2>>(lhs.getValue() * rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> operator*(const double &lhs, const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(lhs * rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> operator*(const RQuantity<M, L, T, A> &lhs, const double &rhs) {
  return RQuantity<M, L, T, A>(lhs.getValue() * rhs);
}
template <typename M1,
          typename L1,
          typename T1,
          typename A1,
          typename M2,
          typename L2,
          typename T2,
          typename A2>
constexpr RQuantity<std::ratio_subtract<M1, M2>,
                    std::ratio_subtract<L1, L2>,
                    std::ratio_subtract<T1, T2>,
                    std::ratio_subtract<A1, A2>>
operator/(const RQuantity<M1, L1, T1, A1> &lhs, const RQuantity<M2, L2, T2, A2> &rhs) {
  return RQuantity<std::ratio_subtract<M1, M2>,
                   std::ratio_subtract<L1, L2>,
                   std::ratio_subtract<T1, T2>,
                   std::ratio_subtract<A1, A2>>(lhs.getValue() / rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_subtract<std::ratio<0>, M>,
                    std::ratio_subtract<std::ratio<0>, L>,
                    std::ratio_subtract<std::ratio<0>, T>,
                    std::ratio_subtract<std::ratio<0>, A>>
operator/(const double &x, const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio_subtract<std::ratio<0>, M>,
                   std::ratio_subtract<std::ratio<0>, L>,
                   std::ratio_subtract<std::ratio<0>, T>,
                   std::ratio_subtract<std::ratio<0>, A>>(x / rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> operator/(const RQuantity<M, L, T, A> &rhs, const double &x) {
  return RQuantity<M, L, T, A>(rhs.getValue() / x);
}

// Comparison operators for quantities:
// ------------------------------------
template <typename M, typename L, typename T, typename A>
constexpr bool operator==(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() == rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator!=(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() != rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator<=(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() <= rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator>=(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() >= rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator<(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() < rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator>(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return (lhs.getValue() > rhs.getValue());
}

// Common math functions:
// ------------------------------

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> abs(const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::abs(rhs.getValue()));
}

template <typename R, typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_multiply<M, R>,
                    std::ratio_multiply<L, R>,
                    std::ratio_multiply<T, R>,
                    std::ratio_multiply<A, R>>
pow(const RQuantity<M, L, T, A> &lhs) {
  return RQuantity<std::ratio_multiply<M, R>,
                   std::ratio_multiply<L, R>,
                   std::ratio_multiply<T, R>,
                   std::ratio_multiply<A, R>>(std::pow(lhs.getValue(), double(R::num) / R::den));
}

template <int R, typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_multiply<M, std::ratio<R>>,
                    std::ratio_multiply<L, std::ratio<R>>,
                    std::ratio_multiply<T, std::ratio<R>>,
                    std::ratio_multiply<A, std::ratio<R>>>
pow(const RQuantity<M, L, T, A> &lhs) {
  return RQuantity<std::ratio_multiply<M, std::ratio<R>>,
                   std::ratio_multiply<L, std::ratio<R>>,
                   std::ratio_multiply<T, std::ratio<R>>,
                   std::ratio_multiply<A, std::ratio<R>>>(std::pow(lhs.getValue(), R));
}

template <int R, typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_divide<M, std::ratio<R>>,
                    std::ratio_divide<L, std::ratio<R>>,
                    std::ratio_divide<T, std::ratio<R>>,
                    std::ratio_divide<A, std::ratio<R>>>
root(const RQuantity<M, L, T, A> &lhs) {
  return RQuantity<std::ratio_divide<M, std::ratio<R>>,
                   std::ratio_divide<L, std::ratio<R>>,
                   std::ratio_divide<T, std::ratio<R>>,
                   std::ratio_divide<A, std::ratio<R>>>(std::pow(lhs.getValue(), 1.0 / R));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_divide<M, std::ratio<2>>,
                    std::ratio_divide<L, std::ratio<2>>,
                    std::ratio_divide<T, std::ratio<2>>,
                    std::ratio_divide<A, std::ratio<2>>>
sqrt(const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio_divide<M, std::ratio<2>>,
                   std::ratio_divide<L, std::ratio<2>>,
                   std::ratio_divide<T, std::ratio<2>>,
                   std::ratio_divide<A, std::ratio<2>>>(std::sqrt(rhs.getValue()));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_divide<M, std::ratio<3>>,
                    std::ratio_divide<L, std::ratio<3>>,
                    std::ratio_divide<T, std::ratio<3>>,
                    std::ratio_divide<A, std::ratio<3>>>
cbrt(const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio_divide<M, std::ratio<3>>,
                   std::ratio_divide<L, std::ratio<3>>,
                   std::ratio_divide<T, std::ratio<3>>,
                   std::ratio_divide<A, std::ratio<3>>>(std::cbrt(rhs.getValue()));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_multiply<M, std::ratio<2>>,
                    std::ratio_multiply<L, std::ratio<2>>,
                    std::ratio_multiply<T, std::ratio<2>>,
                    std::ratio_multiply<A, std::ratio<2>>>
square(const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio_multiply<M, std::ratio<2>>,
                   std::ratio_multiply<L, std::ratio<2>>,
                   std::ratio_multiply<T, std::ratio<2>>,
                   std::ratio_multiply<A, std::ratio<2>>>(std::pow(rhs.getValue(), 2));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_multiply<M, std::ratio<3>>,
                    std::ratio_multiply<L, std::ratio<3>>,
                    std::ratio_multiply<T, std::ratio<3>>,
                    std::ratio_multiply<A, std::ratio<3>>>
cube(const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio_multiply<M, std::ratio<3>>,
                   std::ratio_multiply<L, std::ratio<3>>,
                   std::ratio_multiply<T, std::ratio<3>>,
                   std::ratio_multiply<A, std::ratio<3>>>(std::pow(rhs.getValue(), 3));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> hypot(const RQuantity<M, L, T, A> &lhs,
                                      const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::hypot(lhs.getValue(), rhs.getValue()));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> mod(const RQuantity<M, L, T, A> &lhs,
                                    const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::fmod(lhs.getValue(), rhs.getValue()));
}

template <typename M1,
          typename L1,
          typename T1,
          typename A1,
          typename M2,
          typename L2,
          typename T2,
          typename A2>
constexpr RQuantity<M1, L1, T1, A1> copysign(const RQuantity<M1, L1, T1, A1> &lhs,
                                             const RQuantity<M2, L2, T2, A2> &rhs) {
  return RQuantity<M1, L1, T1, A1>(std::copysign(lhs.getValue(), rhs.getValue()));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> ceil(const RQuantity<M, L, T, A> &lhs,
                                     const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::ceil(lhs.getValue() / rhs.getValue()) * rhs.getValue());
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> floor(const RQuantity<M, L, T, A> &lhs,
                                      const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::floor(lhs.getValue() / rhs.getValue()) * rhs.getValue());
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> trunc(const RQuantity<M, L, T, A> &lhs,
                                      const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::trunc(lhs.getValue() / rhs.getValue()) * rhs.getValue());
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A> round(const RQuantity<M, L, T, A> &lhs,
                                      const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<M, L, T, A>(std::round(lhs.getValue() / rhs.getValue()) * rhs.getValue());
}

// Common trig functions:
// ------------------------------

constexpr Number
sin(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::sin(rhs.getValue()));
}

constexpr Number
cos(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::cos(rhs.getValue()));
}

constexpr Number
tan(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::tan(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
asin(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::asin(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
acos(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::acos(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
atan(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::atan(rhs.getValue()));
}

constexpr Number
sinh(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::sinh(rhs.getValue()));
}

constexpr Number
cosh(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::cosh(rhs.getValue()));
}

constexpr Number
tanh(const RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>> &rhs) {
  return Number(std::tanh(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
asinh(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::asinh(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
acosh(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::acosh(rhs.getValue()));
}

constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
atanh(const Number &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::atanh(rhs.getValue()));
}

template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>
atan2(const RQuantity<M, L, T, A> &lhs, const RQuantity<M, L, T, A> &rhs) {
  return RQuantity<std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>>(
    std::atan2(lhs.getValue(), rhs.getValue()));
}

inline namespace literals {
constexpr long double operator"" _pi(long double x) {
  return static_cast<double>(x) * 3.1415926535897932384626433832795;
}
constexpr long double operator"" _pi(unsigned long long int x) {
  return static_cast<double>(x) * 3.1415926535897932384626433832795;
}
} // namespace literals

// Conversion macro, which utilizes the string literals
#define ConvertTo(_x, _y) (_x).convert(1.0_##_y)

QUANTITY_TYPE(0, 0, 1, 0, QTime)

constexpr QTime second(1.0); // SI base unit
constexpr QTime millisecond = second / 1000;
constexpr QTime minute = 60 * second;
constexpr QTime hour = 60 * minute;
constexpr QTime day = 24 * hour;

inline namespace literals {
constexpr QTime operator"" _s(long double x) {
  return QTime(x);
}
constexpr QTime operator"" _ms(long double x) {
  return static_cast<double>(x) * millisecond;
}
constexpr QTime operator"" _min(long double x) {
  return static_cast<double>(x) * minute;
}
constexpr QTime operator"" _h(long double x) {
  return static_cast<double>(x) * hour;
}
constexpr QTime operator"" _day(long double x) {
  return static_cast<double>(x) * day;
}
constexpr QTime operator"" _s(unsigned long long int x) {
  return QTime(static_cast<double>(x));
}
constexpr QTime operator"" _ms(unsigned long long int x) {
  return static_cast<double>(x) * millisecond;
}
constexpr QTime operator"" _min(unsigned long long int x) {
  return static_cast<double>(x) * minute;
}
constexpr QTime operator"" _h(unsigned long long int x) {
  return static_cast<double>(x) * hour;
}
constexpr QTime operator"" _day(unsigned long long int x) {
  return static_cast<double>(x) * day;
}
} // namespace literals

QUANTITY_TYPE(0, 1, 0, 0, QLength)

constexpr QLength meter(1.0); // SI base unit
constexpr QLength decimeter = meter / 10;
constexpr QLength centimeter = meter / 100;
constexpr QLength millimeter = meter / 1000;
constexpr QLength kilometer = 1000 * meter;
constexpr QLength inch = 2.54 * centimeter;
constexpr QLength foot = 12 * inch;
constexpr QLength yard = 3 * foot;
constexpr QLength mile = 5280 * foot;
constexpr QLength tile = 24 * inch;

inline namespace literals {
constexpr QLength operator"" _mm(long double x) {
  return static_cast<double>(x) * millimeter;
}
constexpr QLength operator"" _cm(long double x) {
  return static_cast<double>(x) * centimeter;
}
constexpr QLength operator"" _m(long double x) {
  return static_cast<double>(x) * meter;
}
constexpr QLength operator"" _km(long double x) {
  return static_cast<double>(x) * kilometer;
}
constexpr QLength operator"" _mi(long double x) {
  return static_cast<double>(x) * mile;
}
constexpr QLength operator"" _yd(long double x) {
  return static_cast<double>(x) * yard;
}
constexpr QLength operator"" _ft(long double x) {
  return static_cast<double>(x) * foot;
}
constexpr QLength operator"" _in(long double x) {
  return static_cast<double>(x) * inch;
}
constexpr QLength operator"" _tile(long double x) {
  return static_cast<double>(x) * tile;
}
constexpr QLength operator"" _mm(unsigned long long int x) {
  return static_cast<double>(x) * millimeter;
}
constexpr QLength operator"" _cm(unsigned long long int x) {
  return static_cast<double>(x) * centimeter;
}
constexpr QLength operator"" _m(unsigned long long int x) {
  return static_cast<double>(x) * meter;
}
constexpr QLength operator"" _km(unsigned long long int x) {
  return static_cast<double>(x) * kilometer;
}
constexpr QLength operator"" _mi(unsigned long long int x) {
  return static_cast<double>(x) * mile;
}
constexpr QLength operator"" _yd(unsigned long long int x) {
  return static_cast<double>(x) * yard;
}
constexpr QLength operator"" _ft(unsigned long long int x) {
  return static_cast<double>(x) * foot;
}
constexpr QLength operator"" _in(unsigned long long int x) {
  return static_cast<double>(x) * inch;
}
constexpr QLength operator"" _tile(unsigned long long int x) {
  return static_cast<double>(x) * tile;
}
} // namespace literals

QUANTITY_TYPE(0, 1, -1, 0, QSpeed)

constexpr QSpeed mps = meter / second;
constexpr QSpeed miph = mile / hour;
constexpr QSpeed kmph = kilometer / hour;

inline namespace literals {
constexpr QSpeed operator"" _mps(long double x) {
  return static_cast<double>(x) * mps;
}
constexpr QSpeed operator"" _miph(long double x) {
  return static_cast<double>(x) * mile / hour;
}
constexpr QSpeed operator"" _kmph(long double x) {
  return static_cast<double>(x) * kilometer / hour;
}
constexpr QSpeed operator"" _mps(unsigned long long int x) {
  return static_cast<double>(x) * mps;
}
constexpr QSpeed operator"" _miph(unsigned long long int x) {
  return static_cast<double>(x) * mile / hour;
}
constexpr QSpeed operator"" _kmph(unsigned long long int x) {
  return static_cast<double>(x) * kilometer / hour;
}
} // namespace literals


QUANTITY_TYPE(0, 1, -3, 0, QJerk)

QUANTITY_TYPE(0, 1, -2, 0, QAcceleration)

constexpr QAcceleration mps2 = meter / (second * second);
constexpr QAcceleration G = 9.80665 * mps2;

inline namespace literals {
constexpr QAcceleration operator"" _mps2(long double x) {
  return QAcceleration(x);
}
constexpr QAcceleration operator"" _mps2(unsigned long long int x) {
  return QAcceleration(static_cast<double>(x));
}
constexpr QAcceleration operator"" _G(long double x) {
  return static_cast<double>(x) * G;
}
constexpr QAcceleration operator"" _G(unsigned long long int x) {
  return static_cast<double>(x) * G;
}
} // namespace literals

// Physical quantity types
QUANTITY_TYPE(0, -1, 0, 1, QCurvature);

// Predegined Number unit
constexpr Number pct = number / 100;

// Predefined Length Units
constexpr QLength field = 12 * foot;

// Predefined Speed Unit
constexpr QSpeed ftps = foot / second;

// Predefined Acceleration Unit
constexpr QAcceleration ftps2 = ftps / second;

// Predefined Jerk Unit
constexpr QJerk mps3 = mps2 / second;
constexpr QJerk ftps3 = ftps2 / second;

inline namespace literals {
// number unit literals
constexpr Number operator"" _pct(long double x) {
    return static_cast<double>(x) * pct;
}
constexpr Number operator"" _pct(unsigned long long int x) {
    return static_cast<double>(x) * pct;
}

// Length Unit Literals
constexpr QLength operator"" _field(long double x) {
    return static_cast<double>(x) * field;
}
constexpr QLength operator"" _field(unsigned long long int x) {
    return static_cast<double>(x) * field;
}

// Predefined Speed Unit
constexpr QSpeed operator"" _ftps(long double x) {
    return static_cast<double>(x) * ftps;
}
constexpr QSpeed operator"" _ftps(unsigned long long int x) {
    return static_cast<double>(x) * ftps;
}


// Predefined Acceleration Unit
constexpr QAcceleration operator"" _ftps2(long double x) {
    return static_cast<double>(x) * ftps2;
}
constexpr QAcceleration operator"" _ftps2(unsigned long long int x) {
    return static_cast<double>(x) * ftps2;
}

// Jerk Literals
constexpr QJerk operator"" _mps3(long double x) {
    return static_cast<double>(x) * mps3;
}
constexpr QJerk operator"" _mps3(unsigned long long int x) {
    return static_cast<double>(x) * mps3;
}
constexpr QJerk operator"" _ftps3(long double x) {
    return static_cast<double>(x) * ftps3;
}
constexpr QJerk operator"" _ftps3(unsigned long long int x) {
    return static_cast<double>(x) * ftps3;
}

} // namespace literals
