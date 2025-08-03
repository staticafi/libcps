#include <cps/math.hpp>

namespace cps {


Scalar real_epsilon_step_along_vector(Vector const& v)
{
    Scalar best_step{ 0.0 };
    for (std::size_t i{ 0ULL }; i != v.size(); ++i)
    {
        Scalar const step{ epsilon_around(v(i)) };
        if (step > best_step)
            best_step = step;
    }
    return best_step;
}


Scalar integral_epsilon_step_along_vector(
    Vector const& v,
    std::uint16_t const num_steps_per_unit,
    std::uint16_t const max_steps,
    Scalar const epsilon
    )
{
    Scalar best_step{ 0.0 };
        Scalar min_error = std::numeric_limits<Scalar>::max();
    if (v.size() > 0)
    {
        for (std::uint16_t i{ 1U }; i <= max_steps; ++i)
        {
            Vector const w{ (((Scalar)i / (Scalar)num_steps_per_unit) * v).array().round() };
            if (w.norm() > 0.9)
            {
                Scalar const t{ w.dot(v) };
                Scalar const error{ (w - t * v).norm() };
                if (error < min_error)
                {
                    min_error = error;
                    best_step = t;
                    if (min_error < epsilon)
                        break;
                }
            }
        }
    }
    return best_step;
}


}
