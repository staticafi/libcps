#include <cps/math.hpp>

namespace cps {


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


bool subspace_orthogonal_to_vector(Matrix const& space, Vector const& vector, Matrix& sub_space)
{
    if (!valid(vector) || vector.norm() < 1e-9)
        return false;
    Vector const g{ vector.normalized() };
    Matrix M(space.cols(), 0);
    for (std::size_t i{ 0ULL }; i < space.cols(); ++i)
    {
        Vector w{ Vector::Unit(space.cols(), i) };
        w -= w.dot(g) * g;
        for (std::size_t j{ 0ULL }; j != M.cols(); ++j)
            w -= w.dot(M.col(j)) * M.col(j);
        if (valid(w) && w.norm() >= 1e-9)
        {
            M.conservativeResize(Eigen::NoChange, M.cols() + 1);
            M.col(M.cols() - 1) = w.normalized();
        }
    }
    sub_space = space * M;
    return true;
}


}
