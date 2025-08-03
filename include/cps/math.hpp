#ifndef CPS_MATH_HPP_INCLUDED
#   define CPS_MATH_HPP_INCLUDED

#   include <Eigen/Dense>
#   include <cmath>
#   include <type_traits>

namespace cps {


using Scalar = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;


inline bool valid(Scalar const s) { return !std::isnan(s) && std::isfinite(s); }
inline bool valid(Vector const& u) { return valid(u.dot(u)); }


template<typename R>
R cast(Scalar value)
{
    if constexpr (std::is_same<R, float>::value == false && std::is_same<R, double>::value == false)
    {
        if (std::isnan(value))
            value = std::numeric_limits<Scalar>::max();
        else if (value < std::numeric_limits<Scalar>::lowest())
            value = std::numeric_limits<Scalar>::lowest();
        else if (value > std::numeric_limits<Scalar>::max())
            value = std::numeric_limits<Scalar>::max();
    }

    if constexpr (std::is_integral<R>::value)
        value = std::round(value);

    if (value <= (Scalar)std::numeric_limits<R>::lowest())
        return std::numeric_limits<R>::lowest();
    if (value >= (Scalar)std::numeric_limits<R>::max())
        return std::numeric_limits<R>::max();

    return (R)value;
}


template<typename T>
Scalar epsilon_around(T const x)
{
    using Type = std::decay_t<T>;
    static_assert(std::is_same<Type, float>::value || std::is_same<Type, double>::value);
    int x_exponent;
    std::frexp(x, &x_exponent);
    return std::pow(2.0, x_exponent - (std::numeric_limits<Type>::digits >> 2));
}


template<typename T>
inline Scalar epsilon_step(T const x, Scalar const dx)
{
    using Type = std::decay_t<T>;
    if constexpr (std::numeric_limits<Type>::is_integer)
        return 1.0 / std::fabs(dx);
    else
        return epsilon_around(x);
}


}

#endif
