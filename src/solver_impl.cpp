#include <cps/solver_impl.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>
#include <set>
#include <limits>
#include <type_traits>

namespace cps {


SolverImpl::SolverImpl(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config_
        )
    : Component{}
    , config{ config_ }
    , constants{ parameter_indices, comparators, {}, {} }
    , round_constants{ seed_input, seed_output }
    , sample{}
    , best_io{}
    , state{ State::ROUND_BEGIN }
    , state_round_begin{ this }
    , state_local_space{ this }
    , state_constraints{ this }
    , state_gradient{ this }
    , state_fuzzing_gradient_descent{ this }
    , state_round_end{ this }
    , state_success{}
    , state_failure{}
    , state_processors{
        { State::ROUND_BEGIN, &state_round_begin },
        { State::LOCAL_SPACE, &state_local_space },
        { State::CONSTRAINTS, &state_constraints },
        { State::GRADIENT, &state_gradient },
        { State::FUZZING_GRADIENT_DESCENT, &state_fuzzing_gradient_descent },
        { State::ROUND_END, &state_round_end },
        { State::SUCCESS, &state_success },
        { State::FAILURE, &state_failure },
    }
    , origin{ Vector(0) }
    , matrix{ Matrix(0,0) }
    , constraints{}
    , gradient{ Vector(0) }
    , statistics{}
{
    if (config.build_local_space)
    {
        std::set<std::size_t> variable_indices;
        std::set<std::size_t> function_indices;
        variable_indices.insert(constants.parameter_indices.back().begin(), constants.parameter_indices.back().end());
        function_indices.insert(constants.parameter_indices.size() - 1ULL);
        std::unordered_set<std::size_t>  work_set{};
        for (std::size_t  i = 1UL; i < constants.parameter_indices.size(); ++i)
            work_set.insert(i - 1UL);
        while (true)
        {
            bool changed{ false };
            for (auto it = work_set.begin(); it != work_set.end(); )
            {
                auto const& params{ constants.parameter_indices.at(*it) };
                bool intersection{ false };
                for (auto i : params)
                    if (variable_indices.contains(i))
                    {
                        intersection = true;
                        break;
                    }
                if (intersection)
                {
                    variable_indices.insert(params.begin(), params.end());
                    function_indices.insert(*it);

                    changed = true;

                    it = work_set.erase(it);
                }
                else
                    ++it;
            }
            if (changed == false)
                break;
        }
        constants.active_variable_indices.assign(variable_indices.begin(), variable_indices.end());
        constants.active_function_indices.assign(function_indices.begin(), function_indices.end());
    }
    else
    {
        constants.active_variable_indices = constants.parameter_indices.back();
        std::sort(constants.active_variable_indices.begin(), constants.active_variable_indices.end());
        constants.active_function_indices.push_back(constants.parameter_indices.size() - 1ULL);
    }

    state_processors.at(state)->enter();
}


void SolverImpl::compute_next_input(std::vector<Variable>& input)
{
    input = round_constants.seed_input;

    sample.ready = false;
    sample.vector.resize(matrix.cols());
    sample.vector.setZero();
    while (!sample.ready)
    {
        if (is_finished())
            return;
        auto& processor{ state_processors.at(state) };
        processor->update();
        State const next_state{ processor->transition() };
        if (next_state != state)
        {
            state = next_state;
            state_processors.at(state)->enter();
        }
    }

    Vector const u{ origin + matrix * sample.vector };
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        std::visit(
            [i, &u](auto&& x) { x = cast<std::decay_t<decltype(x)> >(u(i)); },
            input.at(constants.active_variable_indices.at(i))
        );
    best_io.candidate = input;

    ++statistics.insert({ to_string(state), 0ULL }).first->second;
}


void SolverImpl::process_output(std::vector<Evaluation> const& output_)
{
    if (is_finished())
        return;

    std::vector<Evaluation> output;
    {
        output.reserve(round_constants.seed_output.size());
        for (std::size_t i{ 0ULL }; i < output_.size() && i < round_constants.seed_output.size(); ++i)
        {
            output.push_back(output_.at(i));
            if (output.at(i).predicate != round_constants.seed_output.at(i).predicate)
                break;
        }
    }

    if (output.size() == round_constants.seed_output.size())
    {
        struct local
        {
            static bool is_better_evaluation(Comparator const comparator, Scalar const current_valuation, Scalar const new_valuation)
            {
                switch (comparator)
                {
                    case Comparator::EQUAL:
                        return std::fabs(new_valuation) < std::fabs(current_valuation);
                    case Comparator::UNEQUAL:
                        return std::fabs(new_valuation) > std::fabs(current_valuation);
                    case Comparator::LESS:
                    case Comparator::LESS_EQUAL:
                        return new_valuation < current_valuation;
                    case Comparator::GREATER:
                    case Comparator::GREATER_EQUAL:
                        return new_valuation > current_valuation;
                    default: { UNREACHABLE(); } break;
                }
            }
        };

        std::size_t const last{ output.size() - 1ULL };
        if (output.at(last).predicate != round_constants.seed_output.at(last).predicate)
        {
            state = State::SUCCESS;
            best_io.input = best_io.candidate;
            best_io.output = output;
            return;
        }
        else if (local::is_better_evaluation(constants.comparators.back(), round_constants.seed_output.at(last).function, output.at(last).function)
                    && (best_io.input.empty() ||
                        local::is_better_evaluation(constants.comparators.back(), best_io.output.at(last).function, output.at(last).function)))
        {
            best_io.input = best_io.candidate;
            best_io.output = output;
        }
    }

    state_processors.at(state)->update(output);
}


void SolverImpl::StateRoundBegin::enter()
{
    std::size_t const n{ solver().constants.active_variable_indices.size() };
    solver().origin.resize(n);
    for (std::size_t i{ 0ULL }; i != n; ++i)
        std::visit(
            [i, this](auto x) { solver().origin(i) = (double)x; },
            solver().round_constants.seed_input.at(solver().constants.active_variable_indices.at(i))
        );
    solver().matrix.setIdentity(n,n);
    solver().best_io.clear();
}


Vector SolverImpl::GradientComputationBase::compute_partial_step_vector()
{
    Vector const u{ Vector::Unit(solver().matrix.cols(), column_index) };
    if (step_coeffs.size() == STEP_COEFFS.size())
    {
        Scalar const epsilon = solver().epsilon_step_along_vector(solver().matrix * u);
        for (Scalar& coeff : step_coeffs)
            coeff *= epsilon;
    }
    current_coeff = step_coeffs.back();
    step_coeffs.pop_back();
    return current_coeff * u;
}


void SolverImpl::StateLocalSpace::enter()
{
    active_function_index = 0ULL;
    reset_gradient_computation();
    gradient = Vector::Zero(solver().matrix.cols());
}


void SolverImpl::StateLocalSpace::update()
{
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(active_function_index) };
    if (solver().constants.comparators.at(fn_idx) != Comparator::EQUAL)
    {
        ++active_function_index;
        return;
    }

    if (step_coeffs.empty())
        reset_for_next_partial();

    if (column_index == solver().matrix.cols())
    {
        solver().update_matrix(gradient);

        ++active_function_index;
        reset_gradient_computation();
        gradient = Vector::Zero(solver().matrix.cols());
        return;
    }

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateLocalSpace::update(std::vector<Evaluation> const& output)
{
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(active_function_index) };
    if (fn_idx >= output.size())
        return;
    Scalar const finite_difference{
        compute_finite_difference(output.at(fn_idx).function, solver().round_constants.seed_output.at(fn_idx).function)
    };
    if (valid(finite_difference))
    {
        gradient(column_index) = finite_difference;
        step_coeffs.clear();
    }
}


SolverImpl::State SolverImpl::StateLocalSpace::transition() const
{
    if (solver().matrix.cols() == 0ULL)
        return State::FAILURE;
    if (active_function_index < solver().constants.active_function_indices.size() - 1ULL)
        return solver().state;
    return State::CONSTRAINTS;
}


void SolverImpl::StateConstraints::enter()
{
    reset_gradient_computation();
    gradients.clear();
    for (std::size_t i : solver().constants.active_function_indices)
        if (solver().constants.comparators.at(i) != Comparator::EQUAL)
            gradients.insert({ i, Vector::Zero(solver().matrix.cols()) });
    for (auto const& index_and_grad : gradients)
        partial_function_indices.insert(index_and_grad.first);
}


void SolverImpl::StateConstraints::update()
{
    if (step_coeffs.empty())
    {
        reset_for_next_partial();
        for (auto const& index_and_grad : gradients)
            partial_function_indices.insert(index_and_grad.first);
    }

    if (column_index == solver().matrix.cols())
    {
        for (auto const& index_and_grad : gradients)
            if (valid(index_and_grad.second) && index_and_grad.second.norm() >= 1e-9)
                solver().constraints.push_back({
                    index_and_grad.second.normalized(),
                    -solver().round_constants.seed_output.at(index_and_grad.first).function / index_and_grad.second.norm(),
                    solver().constants.comparators.at(index_and_grad.first)
                });
        return;
    }

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateConstraints::update(std::vector<Evaluation> const& output)
{
    for (auto it = partial_function_indices.begin(); it != partial_function_indices.end(); )
    {
        if (*it < output.size())
        {
            Scalar const finite_difference{
                compute_finite_difference(output.at(*it).function, solver().round_constants.seed_output.at(*it).function)
            };
            if (valid(finite_difference))
            {
                gradients.at(*it)(column_index) = finite_difference;
                it = partial_function_indices.erase(it);
                continue;
            }
        }
        ++it;
    }
    if (partial_function_indices.empty())
        step_coeffs.clear();
}


SolverImpl::State SolverImpl::StateConstraints::transition() const
{
    if (column_index < solver().matrix.cols())
        return solver().state;
    return State::GRADIENT;
}


void SolverImpl::StateGradient::enter()
{
    reset_gradient_computation();
    solver().gradient = Vector::Zero(solver().matrix.cols());
}


void SolverImpl::StateGradient::update()
{
    if (step_coeffs.empty())
        reset_for_next_partial();

    if (column_index == solver().matrix.cols())
        return;

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateGradient::update(std::vector<Evaluation> const& output)
{
    std::size_t const fn_idx{ solver().constants.parameter_indices.size() - 1UL };
    if (fn_idx >= output.size())
        return;
    Scalar const finite_difference{
        compute_finite_difference(output.at(fn_idx).function, solver().round_constants.seed_output.at(fn_idx).function)
    };
    if (valid(finite_difference))
    {
        solver().gradient(column_index) = finite_difference;
        step_coeffs.clear();
    }
}


SolverImpl::State SolverImpl::StateGradient::transition() const
{
    if (column_index < solver().matrix.cols())
        return solver().state;
    return State::FUZZING_GRADIENT_DESCENT;
}


void SolverImpl::StateFuzzingGradientDescent::enter()
{
    if (solver().gradient.dot(solver().gradient) < 1e-9)
        multipliers.clear();
    else
        multipliers.assign({ 1000.0, 100.0, 10.0, 0.1, 0.01, 0.001, 1.0 });
}


void SolverImpl::StateFuzzingGradientDescent::update()
{
    if (multipliers.empty())
        return;

    multiplier = multipliers.back();
    multipliers.pop_back();

    Scalar lambda;
    {
        lambda = -solver().round_constants.seed_output.back().function / solver().gradient.dot(solver().gradient);
        Scalar const epsilon = solver().epsilon_step_along_vector(solver().matrix * ((multiplier * lambda) * solver().gradient));
        switch (solver().constants.comparators.back())
        {
            case Comparator::EQUAL: break;
            case Comparator::UNEQUAL: lambda += epsilon; break;
            case Comparator::LESS: lambda -= epsilon; break;
            case Comparator::LESS_EQUAL: break;
            case Comparator::GREATER: lambda += epsilon; break;
            case Comparator::GREATER_EQUAL: break;
            default: { UNREACHABLE(); } break;
        }
    }

    Vector u{ (multiplier * lambda) * solver().gradient };
    solver().clip_by_constraints(u);

    solver().sample.vector = u;
    solver().sample.ready = true;
}


SolverImpl::State SolverImpl::StateFuzzingGradientDescent::transition() const
{
    if (!multipliers.empty())
        return solver().state;
    return State::ROUND_END;
}


void SolverImpl::StateRoundEnd::enter()
{
    if (!solver().best_io.input.empty())
    {
        solver().round_constants.seed_input = solver().best_io.input;
        solver().round_constants.seed_output = solver().best_io.output;
        improved = true;
    }
    else
        improved = false;
}


SolverImpl::State SolverImpl::StateRoundEnd::transition() const
{
    return improved ? State::ROUND_BEGIN : State::FAILURE;
}


void SolverImpl::update_matrix(Vector const& gradient)
{
    if (!valid(gradient) || gradient.norm() < 1e-9)
        return;
    Vector const g{ gradient.normalized() };
    Matrix M(matrix.rows(), 0);
    for (std::size_t i{ 0ULL }; i < matrix.cols(); ++i)
    {
        Vector w{ Vector::Unit(matrix.cols(), i) };
        w -= w.dot(g) * g;
        for (std::size_t j{ 0ULL }; j != M.cols(); ++j)
            w -= w.dot(M.col(j)) * M.col(j);
        if (valid(w) && w.norm() >= 1e-9)
        {
            M.conservativeResize(Eigen::NoChange, M.cols() + 1);
            M.col(M.cols() - 1) = (matrix * w).normalized();
        }
    }
    matrix = M;
}


Scalar SolverImpl::epsilon_step_along_vector(Vector const& u) const
{
    Scalar epsilon{ 0.0 };
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        std::visit(
            [i, &u, &epsilon](auto const x) {
                Scalar const delta{ epsilon_step(x, u(i)) };
                if (valid(delta) && delta > epsilon)
                    epsilon = delta;
            },
            round_constants.seed_input.at(constants.active_variable_indices.at(i))
        );
    return epsilon;
}


bool SolverImpl::are_constraints_satisfied(Vector const& u) const
{
    for (Constraint const& constraint : constraints)
    {
        Scalar const signed_distance{ u.dot(constraint.normal) };
        switch (constraint.comparator)
        {
            case Comparator::UNEQUAL:
                if (!(signed_distance != constraint.signed_distance))
                    return false;
                break;
            case Comparator::LESS:
                if (!(signed_distance < constraint.signed_distance))
                    return false;
                break;
            case Comparator::LESS_EQUAL:
                if (!(signed_distance <= constraint.signed_distance))
                    return false;
                break;
            case Comparator::GREATER:
                if (!(signed_distance > constraint.signed_distance))
                    return false;
                break;
            case Comparator::GREATER_EQUAL:
                if (!(signed_distance >= constraint.signed_distance))
                    return false;
                break;
            default: { UNREACHABLE(); } break;
        }
    }
    return true;
}


bool SolverImpl::clip_by_constraints(Vector& u, std::size_t const max_iterations) const
{
    for (std::size_t iteration{ 0ULL }; iteration != max_iterations; ++iteration)
    {
        bool  clipped{ false };
        for (Constraint const& constraint : constraints)
        {
            Vector direction{ constraint.normal };
            if (iteration == 0UL && gradient.dot(gradient) > 1e-9)
            {
                Vector const component_of_normal = constraint.normal - (gradient.dot(constraint.normal) / gradient.dot(gradient)) * gradient;
                if (component_of_normal.dot(component_of_normal) > 1e-9)
                    direction = component_of_normal / component_of_normal.dot(constraint.normal);
            }
            Scalar const signed_distance{ u.dot(constraint.normal) };
            if (!valid(signed_distance))
                return false;

            Scalar const epsilon{ epsilon_around(signed_distance) };
            switch (constraint.comparator)
            {
                case Comparator::UNEQUAL:
                    if (!(constraint.signed_distance != signed_distance))
                    {
                        u += ((constraint.signed_distance + epsilon) - signed_distance) * direction;
                        clipped = true;
                    }
                    break;
                case Comparator::LESS:
                    if (!(signed_distance < constraint.signed_distance))
                    {
                        u += ((constraint.signed_distance - epsilon) - signed_distance) * direction;
                        clipped = true;
                    }
                    break;
                case Comparator::LESS_EQUAL:
                    if (!(signed_distance <= constraint.signed_distance))
                    {
                        u += (constraint.signed_distance - signed_distance) * direction;
                        clipped = true;
                    }
                    break;
                case Comparator::GREATER:
                    if (!(signed_distance > constraint.signed_distance))
                    {
                        u += ((constraint.signed_distance + epsilon) - signed_distance) * direction;
                        clipped = true;
                    }
                    break;
                case Comparator::GREATER_EQUAL:
                    if (!(signed_distance >= constraint.signed_distance))
                    {
                        u += (constraint.signed_distance - signed_distance) * direction;
                        clipped = true;
                    }
                    break;
                default: { UNREACHABLE(); } break;
            }
        }
        if (!clipped)
            return true;
    }
    return false;
}


char const* SolverImpl::to_string(State const state)
{
    switch (state)
    {
        case State::ROUND_BEGIN: return "ROUND_BEGIN";
        case State::LOCAL_SPACE: return "LOCAL_SPACE";
        case State::CONSTRAINTS: return "CONSTRAINTS";
        case State::GRADIENT: return "GRADIENT";
        case State::FUZZING_GRADIENT_DESCENT: return "FUZZING_GRADIENT_DESCENT";
        case State::ROUND_END: return "ROUND_END";
        case State::SUCCESS: return "SUCCESS";
        case State::FAILURE: return "FAILURE";
        default: { UNREACHABLE(); return ""; }
    }    
}


}
