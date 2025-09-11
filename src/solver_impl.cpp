#include <cps/solver_impl.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>
#include <set>
#include <limits>
#include <type_traits>

namespace cps {


static inline Variable& scalar_to_variable(Scalar const s, Variable& var)
{
    var.value.u64 = 0ULL;
    var.visit( [s]<typename T>(T&& x) { x = cast<std::decay_t<T> >(s); });
    return var;
}


static inline Scalar variable_to_scalar(Variable const& var)
{
    Scalar s;
    var.visit( [&s]<typename T>(T&& x) { s = (Scalar)x; });
    return s;
}


SolverImpl::SolverImpl(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config_
        )
    : config{ config_ }
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
    , state_fuzzing_bit_flips{ this }
    , state_fuzzing_random{ this }
    , state_round_end{ this }
    , state_success{}
    , state_failure{}
    , state_processors{
        { State::ROUND_BEGIN, &state_round_begin },
        { State::LOCAL_SPACE, &state_local_space },
        { State::CONSTRAINTS, &state_constraints },
        { State::GRADIENT, &state_gradient },
        { State::FUZZING_GRADIENT_DESCENT, &state_fuzzing_gradient_descent },
        { State::FUZZING_BIT_FLIPS, &state_fuzzing_bit_flips },
        { State::FUZZING_RANDOM, &state_fuzzing_random },
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

    statistics.insert({ "FUNCTIONS", constants.active_function_indices.size() });
    statistics.insert({ "VARIABLES", constants.active_variable_indices.size() });
    {
        auto&  count{ statistics.insert({ "EQUALITIES", 0ULL }).first->second };
        for (std::size_t i : constants.active_function_indices)
            if (constants.comparators.at(i) == Comparator::EQUAL)
                ++count;
    }
    {
        auto&  count{ statistics.insert({ "FLOATS", 0ULL }).first->second };
        for (std::size_t i : constants.active_variable_indices)
            switch (round_constants.seed_input.at(i).type)
            {
                case Variable::Type::FLOAT32:
                case Variable::Type::FLOAT64:
                    ++count;
                    break;
                default: break;
            }
    }
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
        State const next_state{ state_processors.at(state)->transition() };
        if (next_state != state)
        {
            state = next_state;
            state_processors.at(state)->enter();
        }
        else
            state_processors.at(state)->update();
    }

    Vector const u{ origin + matrix * sample.vector };
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        scalar_to_variable(u(i), input.at(constants.active_variable_indices.at(i)));
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
        else if (local::is_better_evaluation(opposite(comparator_at(last)), round_constants.seed_output.at(last).function, output.at(last).function)
                    && (best_io.input.empty() ||
                        local::is_better_evaluation(opposite(comparator_at(last)), best_io.output.at(last).function, output.at(last).function)))
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
        solver().origin(i) = variable_to_scalar(solver().round_constants.seed_input.at(solver().constants.active_variable_indices.at(i)));
    solver().matrix.setIdentity(n,n);
    solver().constraints.clear();
    solver().best_io.clear();

    ++count;
    solver().statistics["NUM_ROUNDS"] = std::min(count, solver().config.max_rounds);
}


SolverImpl::State SolverImpl::StateRoundBegin::transition() const
{
    return count <= solver().config.max_rounds ? State::LOCAL_SPACE : State::FAILURE;
}


void SolverImpl::GradientComputationBase::reset_gradient_computation(std::size_t const active_function_index_)
{
    active_function_index = active_function_index_;
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(active_function_index) };
    Scalar const fn_value{ solver().round_constants.seed_output.at(fn_idx).function };
    seed_function_value = valid(fn_value) ? fn_value : 0.0;
    column_index = -1L;
    current_step = 0.0;
    active_coordinates.clear();
    for (std::size_t const i : solver().constants.parameter_indices.at(fn_idx))
        active_coordinates.push_back(
            std::distance(
                solver().constants.active_variable_indices.begin(), 
                std::lower_bound(
                    solver().constants.active_variable_indices.begin(),
                    solver().constants.active_variable_indices.end(),
                    i
                    )
                )
            );
    step_coeffs.clear();
    left_differences.clear();
    right_differences.clear();
    gradient = Vector::Zero(solver().matrix.cols());
    reset_for_next_partial();
}


void SolverImpl::GradientComputationBase::reset_for_next_partial()
{
    current_step = 0.0;
    step_coeffs.clear();
    left_differences.clear();
    right_differences.clear();
    for (++column_index; column_index < solver().matrix.cols(); ++column_index)
    {
        bool any_affected = false;
        for (std::size_t const i : active_coordinates)
            if (std::fabs(solver().matrix.col(column_index)(i)) > 1e-9)
            {
                any_affected = true;
                break;
            }
        if (any_affected == true)
        {
            Vector const u{ solver().matrix * Vector::Unit(solver().matrix.cols(), column_index) };
            solver().epsilon_steps_along_ray(step_coeffs, solver().origin, u,  1.0, &seed_function_value);
            solver().epsilon_steps_along_ray(step_coeffs, solver().origin, u, -1.0, &seed_function_value);
            break;
        }
        gradient(column_index) = 0.0;
    }
}


Vector SolverImpl::GradientComputationBase::compute_partial_step_vector()
{
    current_step = step_coeffs.back();
    step_coeffs.pop_back();
    return current_step * Vector::Unit(solver().matrix.cols(), column_index);
}


bool SolverImpl::GradientComputationBase::update_gradient(std::vector<Evaluation> const& output)
{
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(active_function_index) };
    if (fn_idx < output.size())
    {
        Scalar const finite_difference{ (output.at(fn_idx).function - seed_function_value) / current_step };
        if (valid(finite_difference))
            (current_step < 0.0 ? left_differences : right_differences).push_back(finite_difference);
    }
    if (step_coeffs.empty())
    {
        auto const& max_abs = [](std::vector<Scalar> const& vector) {
            Scalar result = vector.front();
            for (auto const x : vector)
                if (std::fabs(x) > std::fabs(result))
                    result = x;
            return result;
        };
        if (!left_differences.empty() && !right_differences.empty())
        {
            Scalar const left_difference{ max_abs(left_differences) };
            Scalar const right_difference{ max_abs(right_differences) };
            Scalar const product = left_difference * right_difference;
            if (product >= 0.0)
                gradient(column_index) = std::fabs(left_difference) < std::fabs(right_difference) ? left_difference : right_difference;
            // else gradient(column_index) remains zero.
        }
        else if (!left_differences.empty())
            gradient(column_index) = max_abs(left_differences);
        else if (!right_differences.empty())
            gradient(column_index) = max_abs(right_differences);
        // else gradient(column_index) remains zero.

        reset_for_next_partial();
    }
    return is_gradient_ready();
}


bool SolverImpl::GradientComputationBase::is_gradient_ready() const
{
    return column_index == solver().matrix.cols();
}


void SolverImpl::StateLocalSpace::enter()
{
    reset_gradient_computation(0ULL);
}


void SolverImpl::StateLocalSpace::update()
{
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(get_active_function_index()) };
    if (solver().comparator_at(fn_idx) != Comparator::EQUAL)
    {
        reset_gradient_computation(get_active_function_index() + 1ULL);
        return;
    }

    if (is_gradient_ready())
    {
        subspace_orthogonal_to_vector(solver().matrix, get_gradient(), solver().matrix);

        if (solver().is_matrix_feasible())
            reset_gradient_computation(get_active_function_index() + 1ULL);
        else
        {
            std::size_t const n{ solver().constants.active_variable_indices.size() };
            solver().matrix.setIdentity(n,n);
            reset_gradient_computation(solver().constants.active_function_indices.size() - 1ULL);
        }

        return;
    }

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateLocalSpace::update(std::vector<Evaluation> const& output)
{
    update_gradient(output);
}


SolverImpl::State SolverImpl::StateLocalSpace::transition() const
{
    if (!solver().config.build_local_space)
        return State::CONSTRAINTS;
    if (get_active_function_index() < solver().constants.active_function_indices.size() - 1ULL)
        return solver().state;
    return State::CONSTRAINTS;
}


void SolverImpl::StateConstraints::enter()
{
    reset_gradient_computation(0ULL);
}


void SolverImpl::StateConstraints::update()
{
    std::size_t const fn_idx{ solver().constants.active_function_indices.at(get_active_function_index()) };
    if (solver().comparator_at(fn_idx) == Comparator::EQUAL)
    {
        reset_gradient_computation(get_active_function_index() + 1ULL);
        return;
    }

    if (is_gradient_ready())
    {
        if (valid(get_gradient()) && get_gradient().norm() >= 1e-9)
            solver().constraints.push_back({
                get_gradient().normalized(),
                -solver().round_constants.seed_output.at(fn_idx).function / get_gradient().norm(),
                solver().comparator_at(fn_idx)
            });

        reset_gradient_computation(get_active_function_index() + 1ULL);
        return;
    }

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateConstraints::update(std::vector<Evaluation> const& output)
{
    update_gradient(output);
}


SolverImpl::State SolverImpl::StateConstraints::transition() const
{
    if (!solver().config.build_constraints)
        return State::GRADIENT;
    if (get_active_function_index() < solver().constants.active_function_indices.size() - 1ULL)
        return solver().state;
    return State::GRADIENT;
}


void SolverImpl::StateGradient::enter()
{
    reset_gradient_computation(solver().constants.active_function_indices.size() - 1ULL);
    solver().gradient = Vector::Zero(solver().matrix.cols());
}


void SolverImpl::StateGradient::update()
{
    if (is_gradient_ready())
    {
        solver().gradient = get_gradient();
        return;
    }

    solver().sample.vector = compute_partial_step_vector();
    solver().sample.ready = true;
}


void SolverImpl::StateGradient::update(std::vector<Evaluation> const& output)
{
    if (update_gradient(output))
        solver().gradient = get_gradient();
}


SolverImpl::State SolverImpl::StateGradient::transition() const
{
    if (solver().config.use_gradient_descent || solver().config.use_random_fuzzing)
    {
        if (!is_gradient_ready())
            return solver().state;
    }
    return State::FUZZING_GRADIENT_DESCENT;
}


void SolverImpl::StateProcessorWithClipCache::push_to_clip_cache(Vector u)
{
    clip_cache.push_back(u);
    if (solver().clip_by_constraints(u, 50ULL))
        clip_cache.push_back(u);
}


void SolverImpl::StateProcessorWithClipCache::take_sample_from_clip_cache()
{
    solver().sample.vector = clip_cache.back();
    solver().sample.ready = true;

    clip_cache.pop_back();
}


void SolverImpl::StateFuzzingGradientDescent::enter()
{
    samples.clear();

    std::vector<Vector> directions{};
    directions.push_back(solver().gradient);
    if (solver().matrix.cols() > 1L)
        for (std::ptrdiff_t i{ 0L}; i != solver().matrix.cols(); ++i)
            if (std::fabs(solver().gradient(i)) >= 1e-9)
                directions.push_back(solver().gradient(i) * Vector::Unit(solver().matrix.cols(), i));
    for (Vector const& u : directions)
    {
        Scalar const lambda{ -solver().round_constants.seed_output.back().function / u.dot(u) };
        if (!valid(lambda) || !valid(lambda * u))
            continue;
        Vector const v{ solver().matrix * u };
        Vector const S{ solver().origin + lambda * v };
        std::vector<Scalar> steps;
        switch (opposite(solver().comparator_at(solver().constants.active_function_indices.back())))
        {
            case Comparator::EQUAL:
                steps.push_back(0.0);
                break;
            case Comparator::UNEQUAL:
                solver().epsilon_steps_along_ray(steps, S, v,  1.0);
                solver().epsilon_steps_along_ray(steps, S, v, -1.0);
                break;
            case Comparator::LESS:
                solver().epsilon_steps_along_ray(steps, S, v, -1.0);
                break;
            case Comparator::LESS_EQUAL:
                steps.push_back(0.0);
                solver().epsilon_steps_along_ray(steps, S, v, -1.0);
                break;
            case Comparator::GREATER:
                solver().epsilon_steps_along_ray(steps, S, v,  1.0);
                break;
            case Comparator::GREATER_EQUAL:
                steps.push_back(0.0);
                solver().epsilon_steps_along_ray(steps, S, v,  1.0);
                break;
            default: { UNREACHABLE(); } break;
        }
        for (Scalar step : steps)
            samples.push_back((lambda + step) * u);
    }                
}


void SolverImpl::StateFuzzingGradientDescent::update()
{
    if (!is_clip_cache_empty())
    {
        take_sample_from_clip_cache();
        return;
    }

    if (samples.empty())
        return;

    Vector u = samples.back();
    samples.pop_back();

    push_to_clip_cache(u);
}


SolverImpl::State SolverImpl::StateFuzzingGradientDescent::transition() const
{
    if (!solver().config.use_gradient_descent)
        return State::FUZZING_BIT_FLIPS;
    if (!is_clip_cache_empty() || !samples.empty())
        return solver().state;
    return State::FUZZING_BIT_FLIPS;
}


void SolverImpl::StateFuzzingBitFlips::enter()
{
    var = 0ULL;
    bit = 0ULL;
}


void SolverImpl::StateFuzzingBitFlips::update()
{
    if (!is_clip_cache_empty())
    {
        take_sample_from_clip_cache();
        return;
    }

    if (var == solver().constants.active_variable_indices.size())
        return;

    Variable const& var_ref{ solver().round_constants.seed_input.at(solver().constants.active_variable_indices.at(var)) };

    std::size_t num_bits;
    var_ref.visit([&num_bits]<typename T>(T) { num_bits = std::is_integral<std::decay_t<T> >::value ? 8ULL * sizeof(T) : 0ULL; });

    if (bit == num_bits)
    {
        ++var;
        bit = 0UL;
        return;
    }

    Vector u(0);
    {
        Vector n = Vector::Unit(solver().matrix.rows(), var);
        Vector W = n.transpose() * solver().matrix;
        Vector w = solver().matrix * W;

        // plane: n*(X-(solver().origin + s*n)) = 0  // s is computed later.
        // line : X=solver().origin + t*w
        // -----------------------------------------------
        // n*(solver().origin + t*w - (solver().origin + s*n)) = 0
        // n*(t*w - s*n) = 0
        // t = s/n*w

        Scalar const nw{ n.dot(w) };
        if (std::fabs(nw) > 1e-6)
        {
            Scalar s;
            var_ref.visit([this, &s]<typename T>(T const x) {
                Scalar const sign = ((std::uint64_t)x & (1ULL << bit)) != 0ULL ? -1.0 : 1.0;
                Scalar const magnitude = (Scalar)(1ULL << bit);
                s = sign * magnitude;
            });
            u = (s / nw) * W;
        }
    }

    ++bit;

    if (u.size() == 0)
        return;

    push_to_clip_cache(u);
}


SolverImpl::State SolverImpl::StateFuzzingBitFlips::transition() const
{
    if (!solver().config.use_bit_flips)
        return State::FUZZING_RANDOM;
    if (!is_clip_cache_empty() || var < solver().constants.active_variable_indices.size())
        return solver().state;
    return State::FUZZING_RANDOM;
}


void SolverImpl::StateFuzzingRandom::enter()
{
    std::size_t constexpr num_samples_per_cube{ 25ULL };

    cube_half_size = 1.0 * (std::log(std::fabs(solver().round_constants.seed_output.back().function) + 1.0) + 1.0);
    cubes.clear();
    if (solver().gradient.dot(solver().gradient) > 1e-9)
    {
        Scalar const lambda = -solver().round_constants.seed_output.back().function / solver().gradient.dot(solver().gradient);
        cubes.push_back({ .center = lambda * solver().gradient, .num_remaining = num_samples_per_cube });
    }
    cubes.push_back({ .center = Vector::Zero(solver().matrix.cols()), .num_remaining = num_samples_per_cube });
}


void SolverImpl::StateFuzzingRandom::update()
{
    if (!is_clip_cache_empty())
    {
        take_sample_from_clip_cache();
        return;
    }

    if (cubes.empty())
        return;
    if (cubes.back().num_remaining == 0ULL)
    {
        cubes.pop_back();
        return;
    }

    --cubes.back().num_remaining;

    Vector u{ cubes.back().center };
    for (std::size_t i{ 0ULL}; i != u.size(); ++i)
    {
        Scalar const sign{ get_random_natural_64_bit_in_range(0ULL, 100ULL, generator) < 50ULL ? -1.0 : 1.0 };
        Scalar const magnitude{ get_random_float_64_bit_in_range(0.0, cube_half_size, generator) };
        u(i) += sign * magnitude;
    }

    push_to_clip_cache(u);
}


SolverImpl::State SolverImpl::StateFuzzingRandom::transition() const
{
    if (!solver().config.use_random_fuzzing)
        return State::ROUND_END;
    if (!is_clip_cache_empty() || !cubes.empty())
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


Comparator SolverImpl::comparator_at(std::size_t const idx) const
{
    return round_constants.seed_output.at(idx).predicate ? constants.comparators.at(idx) : opposite(constants.comparators.at(idx));
}


void SolverImpl::epsilon_steps_along_ray(
    std::vector<Scalar>& result,
    Vector const& S,
    Vector u,
    Scalar const sign,
    Scalar const* const function_value
    ) const
{
    ASSUMPTION(S.size() == constants.active_variable_indices.size() && S.size() == u.size());

    u = sign * u;

    std::vector<Scalar> steps{ 0.0 };
    if (function_value != nullptr)
        steps.push_back(0.0);

    std::vector<std::size_t> int_coords;
    for (std::size_t i{ 0ULL }; i != constants.active_variable_indices.size(); ++i)
        if (std::fabs(u(i)) >= 1e-6)
            round_constants.seed_input.at(constants.active_variable_indices.at(i)).visit([&, i, function_value]<typename T>(T) {
                if constexpr (std::is_integral<std::decay_t<T> >::value)
                    int_coords.push_back(i);
                else
                {
                    Scalar const x{ (Scalar)cast<std::decay_t<T> >(S(i)) };
                    Scalar step{ std::fabs(epsilon_around<std::decay_t<T> >(x) / u(i)) };
                    if (step > steps.front())
                        steps.front() = step;
                    if (function_value != nullptr)
                    {
                        Scalar constexpr t{ 0.01 };
                        Scalar const y{ (1.0 - t) * x + t * (*function_value) };
                        step = std::fabs(epsilon_around<std::decay_t<T> >(y) / u(i));
                        if (step > steps.back())
                            steps.back() = step;

                    }
                }
            });

    Scalar integral_epsilon{ 0.0 };
    if (!int_coords.empty())
    {
        Vector P(int_coords.size()), v(int_coords.size());
        for (std::size_t i{ 0ULL }; i != int_coords.size(); ++i)
        {
            P(i) = S(int_coords.at(i));
            v(i) = u(int_coords.at(i));
        }
        integral_epsilon = integral_epsilon_step_along_vector(P, v);
    }

    for (Scalar& step : steps)
        step = std::max(step, integral_epsilon);

    if (function_value != nullptr && (steps.back() == 0.0 || steps.back() == steps.front()))
        steps.pop_back();

    for (Scalar step : steps)
        result.push_back(sign * step);
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
    Vector orig_u = u;
    bool any_change{ false };
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
                continue;

            Scalar const epsilon{ epsilon_around<double>(signed_distance) };
            switch (constraint.comparator)
            {
                case Comparator::UNEQUAL:
                    if (std::fabs(signed_distance - constraint.signed_distance) < 1e-9)
                    {
                        Scalar const sign = direction.dot(u - orig_u) < 0.0 ? -1.0 : 1.0;
                        u += ((constraint.signed_distance + sign * std::fabs(epsilon)) - signed_distance) * direction;
                        clipped = true;
                        any_change = true;
                    }
                    break;
                case Comparator::LESS:
                    if (!(signed_distance < constraint.signed_distance))
                    {
                        u += ((constraint.signed_distance - epsilon) - signed_distance) * direction;
                        clipped = true;
                        any_change = true;
                    }
                    break;
                case Comparator::LESS_EQUAL:
                    if (!(signed_distance <= constraint.signed_distance))
                    {
                        u += (constraint.signed_distance - signed_distance) * direction;
                        clipped = true;
                        any_change = true;
                    }
                    break;
                case Comparator::GREATER:
                    if (!(signed_distance > constraint.signed_distance))
                    {
                        u += ((constraint.signed_distance + epsilon) - signed_distance) * direction;
                        clipped = true;
                        any_change = true;
                    }
                    break;
                case Comparator::GREATER_EQUAL:
                    if (!(signed_distance >= constraint.signed_distance))
                    {
                        u += (constraint.signed_distance - signed_distance) * direction;
                        clipped = true;
                        any_change = true;
                    }
                    break;
                default: { UNREACHABLE(); } break;
            }
        }
        if (!clipped)
            return any_change;
    }
    return any_change;
}


bool SolverImpl::is_matrix_feasible()
{
    for (std::ptrdiff_t i{ 0L }; i != matrix.cols(); ++i)
    {
        Vector const& u = matrix.col(i);
        for (auto const j : constants.parameter_indices.back())
        {
            auto const k{
                    std::distance(
                        constants.active_variable_indices.begin(), 
                        std::lower_bound(constants.active_variable_indices.begin(), constants.active_variable_indices.end(), j)
                        )
                    };
            if (std::fabs(u(k)) > 1e-6f)
                return true;
        }
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
        case State::FUZZING_BIT_FLIPS: return "FUZZING_BIT_FLIPS";
        case State::FUZZING_RANDOM: return "FUZZING_RANDOM";
        case State::ROUND_END: return "ROUND_END";
        case State::SUCCESS: return "SUCCESS";
        case State::FAILURE: return "FAILURE";
        default: { UNREACHABLE(); return ""; }
    }    
}


}
