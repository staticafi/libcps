#ifndef CPS_MATH_HPP_INCLUDED
#   define CPS_MATH_HPP_INCLUDED

#   include <Eigen/Dense>
#   include <cmath>
#   include <type_traits>

namespace cps {


template<typename T> struct is_floating_point : std::is_floating_point<T> {};

// using namespace boost::multiprecision;
// template<> struct is_floating_point<boost::multiprecision::cpp_bin_float_quad> : std::true_type {};
// template<> struct is_floating_point<boost::multiprecision::cpp_dec_float_100> : std::true_type {};


using Scalar = long double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;


inline bool valid(Scalar const s) { return !std::isnan(s) && std::isfinite(s); }
inline bool valid(Vector const& u) { return valid(u.dot(u)); }


template<typename R>
R cast(Scalar value)
{
    if constexpr (is_floating_point<R>::value == false)
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
Scalar epsilon_around(Scalar const x)
{
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    int x_exponent;
    std::frexp(x, &x_exponent);
    return std::pow(2.0, x_exponent - (std::numeric_limits<T>::digits >> 1));
}


template<typename T>
Scalar real_epsilon_step_along_vector(Vector const& v)
{
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    Scalar best_step{ 0.0 };
    for (std::size_t i{ 0ULL }; i != v.size(); ++i)
    {
        Scalar const step{ epsilon_around<T>(v(i)) };
        if (step > best_step)
            best_step = step;
    }
    return best_step;
}


Scalar integral_epsilon_step_along_vector(
    Vector const& v,
    std::uint16_t num_steps_per_unit = 2U,
    std::uint16_t max_steps = 1000U,
    Scalar epsilon = 1e-6f
    );


bool subspace_orthogonal_to_vector(Matrix const& space, Vector const& vector, Matrix& sub_space);


}

#endif
