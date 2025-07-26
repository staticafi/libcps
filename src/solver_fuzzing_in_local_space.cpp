#include <cps/solver_fuzzing_in_local_space.hpp>
#include <cps/active_indices.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>
#include <limits>
#include <type_traits>

namespace cps {


SolverFuzzingInLocalSpace::SolverFuzzingInLocalSpace(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output
        )
    : Component{}
    , constants{ parameter_indices, comparators, {}, {} }
    , round_constants{ seed_input, seed_output }
    , state{ State::ROUND_BEGIN }
    , state_local_space{}
    , state_constraints{}
    , state_gradient{}
    , origin{}
    , matrix{}
{
    initialize_active_indices();
}


void SolverFuzzingInLocalSpace::initialize_active_indices()
{
    cps::compute_active_indices(constants.parameter_indices, constants.active_variable_indices, constants.active_function_indices);
}


void SolverFuzzingInLocalSpace::compute_next_input(std::vector<Variable>& input)
{
    input = round_constants.seed_input;

    sample.ready = false;
    sample.vector.resize(matrix.cols());
    sample.vector.setZero();
    while (!sample.ready)
        switch (state)
        {
            case State::ROUND_BEGIN:
                StateRoundBegin_update();
                state = StateRoundBegin_transition();
                break;
            case State::LOCAL_SPACE:
                StateLocalSpace_update();
                state = StateLocalSpace_transition();
                break;
            case State::CONSTRAINTS:
                StateConstraints_update();
                state = StateConstraints_transition();
                break;
            case State::GRADIENT:
                StateGradient_update();
                state = StateGradient_transition();
                break;
            case State::SUCCESS:
            case State::FAILURE:
                return;
            default: UNREACHABLE(); return;
        }

    Vector const u{ origin + matrix * sample.vector };
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        std::visit(
            [i, &u](auto&& x) { x = cast<std::decay_t<decltype(x)> >(u(i)); },
            input.at(constants.active_variable_indices.at(i))
        );
}


void SolverFuzzingInLocalSpace::process_output(std::vector<Evaluation> const& output_)
{
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
        std::size_t const last{ output.size() - 1ULL };
        if (output.at(last).predicate != round_constants.seed_output.at(last).predicate)
            state = State::SUCCESS;
    }

    switch (state)
    {
        case State::ROUND_BEGIN:
            break;
        case State::LOCAL_SPACE:
            StateLocalSpace_update(output);
            break;
        case State::CONSTRAINTS:
            StateConstraints_update(output);
            break;
        case State::GRADIENT:
            StateGradient_update(output);
            break;
        case State::SUCCESS:
        case State::FAILURE:
            return;
        default: UNREACHABLE(); return;
    }
}


void SolverFuzzingInLocalSpace::StateRoundBegin_update()
{
    std::size_t const n{ constants.active_variable_indices.size() };
    origin.resize(n);
    for (std::size_t i{ 0ULL }; i != n; ++i)
        std::visit(
            [i, this](auto x) { origin(i) = (double)x; },
            round_constants.seed_input.at(constants.active_variable_indices.at(i))
        );
    matrix.setIdentity(n,n);
}

SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateRoundBegin_transition()
{
    return StateLocalSpace_enter();
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateLocalSpace_enter()
{
    state_local_space.active_function_index = 0ULL;
    state_local_space.column_index = 0ULL;
    state_local_space.step = StateLocalSpace::Step::POSITIVE;
    state_local_space.epsilon = 0.0;
    state_local_space.gradient.resize(matrix.cols());
    state_local_space.gradient.setZero();
    return State::LOCAL_SPACE;
}


void SolverFuzzingInLocalSpace::StateLocalSpace_update()
{
    std::size_t const fn_idx{ constants.active_function_indices.at(state_local_space.active_function_index) };
    if (constants.comparators.at(fn_idx) != Comparator::EQUAL)
    {
        ++state_local_space.active_function_index;
        return;
    }

    if (state_local_space.step == StateLocalSpace::Step::STOP)
    {
        ++state_local_space.column_index;
        state_local_space.step = StateLocalSpace::Step::POSITIVE;
        state_local_space.epsilon = 0.0;
    }

    if (state_local_space.column_index == matrix.cols())
    {
        updateMatrix(state_local_space.gradient);

        ++state_local_space.active_function_index;
        state_local_space.column_index = 0ULL;
        state_local_space.step = StateLocalSpace::Step::POSITIVE;
        state_local_space.epsilon = 0.0;
        state_local_space.gradient.resize(matrix.cols());
        state_local_space.gradient.setZero();
        return;
    }

    Vector const u{ Vector::Unit(matrix.cols(), state_local_space.column_index) };
    state_local_space.epsilon = computeEpsilon(matrix * u) * (state_local_space.step == StateLocalSpace::Step::POSITIVE ? 1.0 : -1.0);

    sample.vector = state_local_space.epsilon * u;
    sample.ready = true;
}


void SolverFuzzingInLocalSpace::StateLocalSpace_update(std::vector<Evaluation> const& output)
{
    std::size_t const fn_idx{ constants.active_function_indices.at(state_local_space.active_function_index) };
    if (fn_idx < output.size())
    {
        Scalar const finite_difference{
            (output.at(fn_idx).function - round_constants.seed_output.at(fn_idx).function) / state_local_space.epsilon
        };
        if (!std::isnan(finite_difference) && std::isfinite(finite_difference))
        {
            state_local_space.gradient(state_local_space.column_index) = finite_difference;
            state_local_space.step = StateLocalSpace::Step::STOP;
            return;
        }
    }
    state_local_space.step = state_local_space.step == StateLocalSpace::Step::POSITIVE ?
        StateLocalSpace::Step::NEGATIVE : StateLocalSpace::Step::STOP;
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateLocalSpace_transition()
{
    if (matrix.cols() == 0ULL)
        return State::FAILURE;
    if (state_local_space.active_function_index < constants.active_function_indices.size() - 1ULL)
        return state;
    return StateConstraints_enter();
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateConstraints_enter()
{
    state_constraints.column_index = 0ULL;
    state_constraints.step = StateConstraints::Step::POSITIVE;
    state_constraints.epsilon = 0.0;
    state_constraints.gradients.clear();
    for (std::size_t i : constants.active_function_indices)
        if (constants.comparators.at(i) != Comparator::EQUAL)
            state_constraints.gradients.insert({ i, Vector::Zero(matrix.cols()) });
    for (auto const& index_and_grad : state_constraints.gradients)
        state_constraints.partial_function_indices.insert(index_and_grad.first);
    return State::CONSTRAINTS;
}


void SolverFuzzingInLocalSpace::StateConstraints_update()
{
    if (state_constraints.step == StateConstraints::Step::STOP)
    {
        ++state_constraints.column_index;
        state_constraints.step = StateConstraints::Step::POSITIVE;
        state_constraints.epsilon = 0.0;
        for (auto const& index_and_grad : state_constraints.gradients)
            state_constraints.partial_function_indices.insert(index_and_grad.first);
    }

    if (state_constraints.column_index == matrix.cols())
    {
        for (auto const& index_and_grad : state_constraints.gradients)
            if (valid(index_and_grad.second) && index_and_grad.second.norm() >= 1e-9)
                constraints.push_back({
                    index_and_grad.second.normalized(),
                    -round_constants.seed_output.at(index_and_grad.first).function / index_and_grad.second.norm(),
                    constants.comparators.at(index_and_grad.first)
                });
        return;
    }

    Vector const u{ Vector::Unit(matrix.cols(), state_constraints.column_index) };
    state_constraints.epsilon = computeEpsilon(matrix * u) * (state_constraints.step == StateConstraints::Step::POSITIVE ? 1.0 : -1.0);

    sample.vector = state_constraints.epsilon * u;
    sample.ready = true;
}


void SolverFuzzingInLocalSpace::StateConstraints_update(std::vector<Evaluation> const& output)
{
    for (auto it = state_constraints.partial_function_indices.begin(); it != state_constraints.partial_function_indices.end(); )
    {
        if (*it < output.size())
        {
            Scalar const finite_difference{
                (output.at(*it).function - round_constants.seed_output.at(*it).function) / state_constraints.epsilon
            };
            if (!std::isnan(finite_difference) && std::isfinite(finite_difference))
            {
                state_constraints.gradients.at(*it)(state_constraints.column_index) = finite_difference;
                it = state_constraints.partial_function_indices.erase(it);
                continue;
            }
        }
        ++it;
    }
    state_constraints.step = !state_constraints.partial_function_indices.empty() && state_constraints.step == StateConstraints::Step::POSITIVE ?
        StateConstraints::Step::NEGATIVE : StateConstraints::Step::STOP;
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateConstraints_transition()
{
    if (state_constraints.column_index < matrix.cols())
        return state;
    return StateGradient_enter();
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateGradient_enter()
{
    state_gradient.column_index = 0ULL;
    state_gradient.step = StateGradient::Step::POSITIVE;
    state_gradient.epsilon = 0.0;
    gradient.resize(matrix.cols());
    gradient.setZero();
    return State::GRADIENT;
}


void SolverFuzzingInLocalSpace::StateGradient_update()
{
    if (state_gradient.step == StateGradient::Step::STOP)
    {
        ++state_gradient.column_index;
        state_gradient.step = StateGradient::Step::POSITIVE;
        state_gradient.epsilon = 0.0;
    }

    if (state_gradient.column_index == matrix.cols())
        return;

    Vector const u{ Vector::Unit(matrix.cols(), state_gradient.column_index) };
    state_gradient.epsilon = computeEpsilon(matrix * u) * (state_gradient.step == StateGradient::Step::POSITIVE ? 1.0 : -1.0);

    sample.vector = state_gradient.epsilon * u;
    sample.ready = true;
}


void SolverFuzzingInLocalSpace::StateGradient_update(std::vector<Evaluation> const& output)
{
    std::size_t const fn_idx{ constants.parameter_indices.size() - 1UL };
    if (fn_idx < output.size())
    {
        Scalar const finite_difference{
            (output.at(fn_idx).function - round_constants.seed_output.at(fn_idx).function) / state_gradient.epsilon
        };
        if (!std::isnan(finite_difference) && std::isfinite(finite_difference))
        {
            gradient(state_gradient.column_index) = finite_difference;
            state_gradient.step = StateGradient::Step::STOP;
            return;
        }
    }
    state_gradient.step = state_gradient.step == StateGradient::Step::POSITIVE ?
        StateGradient::Step::NEGATIVE : StateGradient::Step::STOP;
}


SolverFuzzingInLocalSpace::State SolverFuzzingInLocalSpace::StateGradient_transition()
{
    if (state_gradient.column_index < matrix.cols())
        return state;
    // TODO:
    return State::FAILURE;
}


void SolverFuzzingInLocalSpace::updateMatrix(Vector const& gradient)
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


Scalar SolverFuzzingInLocalSpace::computeEpsilon(Vector const& u) const
{
    Scalar epsilon{ 0.0 };
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        std::visit(
            [i, &u, &epsilon](auto const x) {
                Scalar const delta{ epsilonStep(x, u(i)) };
                if (valid(delta) && delta > epsilon)
                    epsilon = delta;
            },
            round_constants.seed_input.at(constants.active_variable_indices.at(i))
        );
    return epsilon;
}


}
