#ifndef CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED
#   define CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/config.hpp>
#   include <cps/math.hpp>
#   include <cps/statistics.hpp>
#   include <utility/random.hpp>
#   include <vector>
#   include <unordered_set>
#   include <unordered_map>
#   include <map>
#   include <memory>

namespace cps {


struct SolverImpl
{
    SolverImpl(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config_
        );

    bool is_finished() const { return success() || failure(); }
    bool success() const { return state == State::SUCCESS; }
    bool failure() const { return state == State::FAILURE; }

    std::vector<Variable> const& best_input() const { return success() ? best_io.input : round_constants.seed_input; }
    std::vector<Evaluation> const& best_output() const { return success() ? best_io.output : round_constants.seed_output; }

    void compute_next_input(std::vector<Variable>& input);
    void process_output(std::vector<Evaluation> const& output_);

    Statistics const& get_statistics() const { return statistics; }

private:

    struct Constants
    {
        std::vector<std::vector<std::size_t> > parameter_indices{};
        std::vector<Comparator> comparators{};
        std::vector<std::size_t> active_variable_indices{};
        std::vector<std::size_t> active_function_indices{};
    };

    struct RoundConstants
    {
        std::vector<Variable> seed_input{};
        std::vector<Evaluation> seed_output{};
    };

    struct Sample
    {
        bool ready{ false };
        Vector vector{};
    };

    struct BestIO
    {
        void clear() { candidate.clear(); input.clear(); output.clear(); }
        std::vector<Variable> candidate{};
        std::vector<Variable> input{};
        std::vector<Evaluation> output{};
    };

    struct Constraint
    {
        Vector normal{ Vector(0) };
        Scalar signed_distance{ 0.0 };
        Comparator comparator{ Comparator::EQUAL };
    };

    enum struct State
    {
        ROUND_BEGIN,
        LOCAL_SPACE,
        CONSTRAINTS,
        GRADIENT,
        FUZZING_GRADIENT_DESCENT,
        FUZZING_BIT_FLIPS,
        FUZZING_RANDOM,
        ROUND_END,

        SUCCESS,
        FAILURE,
    };
    static char const* to_string(State);

    struct StateProcessor
    {
        explicit StateProcessor(SolverImpl* const solver = nullptr) : solver_{ solver } {}
        virtual ~StateProcessor() {}
        SolverImpl& solver() { return *solver_; }
        SolverImpl const& solver() const { return *solver_; }
        virtual void enter() {}
        virtual void update() {}
        virtual void update(std::vector<Evaluation> const& output) {}
        virtual State transition() const = 0;
    private:
        SolverImpl* solver_;
    };

    struct StateRoundBegin : public StateProcessor
    {
        explicit StateRoundBegin(SolverImpl* const solver) : StateProcessor{ solver } {}
        void enter() override;
        State transition() const override;
    private:
        std::uint32_t count{ 0U };
    };

    struct GradientComputationBase : public StateProcessor
    {
        explicit GradientComputationBase(SolverImpl* const solver) : StateProcessor{ solver } {}
        void reset_gradient_computation(std::size_t active_function_index);
        Vector compute_partial_step_vector();
        bool update_gradient(std::vector<Evaluation> const& output);
        bool is_gradient_ready() const;
        Vector const& get_gradient() const { return gradient; }
        std::size_t get_active_function_index() const { return active_function_index; }
    private:
        void reset_for_next_partial();
        std::size_t active_function_index{ 0ULL };
        Scalar seed_function_value{ 0.0 };
        std::ptrdiff_t column_index{ -1L };
        std::vector<std::size_t> active_coordinates{};
        Scalar current_step{ 0.0 };
        std::vector<Scalar> step_coeffs{};
        std::vector<Scalar> left_differences{};
        std::vector<Scalar> right_differences{};
        Vector gradient{};
    };

    struct StateLocalSpace : public GradientComputationBase
    {
        explicit StateLocalSpace(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    };

    struct StateConstraints : public GradientComputationBase
    {
        explicit StateConstraints(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    };

    struct StateGradient : public GradientComputationBase
    {
        explicit StateGradient(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    };

    struct StateProcessorWithClipCache : public StateProcessor
    {
        explicit StateProcessorWithClipCache(SolverImpl* const solver) : StateProcessor{ solver } {}
        bool is_clip_cache_empty() const { return clip_cache.empty(); }
        void push_to_clip_cache(Vector u);
        void take_sample_from_clip_cache();
    private:
        std::vector<Vector> clip_cache{};
    };

    struct StateFuzzingGradientDescent : public StateProcessorWithClipCache
    {
        explicit StateFuzzingGradientDescent(SolverImpl* const solver) : StateProcessorWithClipCache{ solver } {}
        void enter() override;
        void update() override;
        State transition() const override;
    private:
        std::vector<Vector> samples{};
    };

    struct StateFuzzingBitFlips : public StateProcessorWithClipCache
    {
        explicit StateFuzzingBitFlips(SolverImpl* const solver) : StateProcessorWithClipCache{ solver } {}
        void enter() override;
        void update() override;
        State transition() const override;
    private:
        std::size_t var{ 0ULL };
        std::size_t bit{ 0ULL };
    };

    struct StateFuzzingRandom : public StateProcessorWithClipCache
    {
        explicit StateFuzzingRandom(SolverImpl* const solver) : StateProcessorWithClipCache{ solver } {}
        void enter() override;
        void update() override;
        State transition() const override;
    private:
        struct CubeInfo
        {
            Vector center;
            std::size_t num_remaining;
        };
        Scalar cube_half_size{ 0.0 };
        std::vector<CubeInfo> cubes{};
        random_generator_for_natural_64_bit generator{ 1ULL };
    };

    struct StateRoundEnd : public StateProcessor
    {
        explicit StateRoundEnd(SolverImpl* const solver) : StateProcessor{ solver } {}
        void enter() override;
        State transition() const override;
    private:
        bool improved{ false };
    };

    struct StateSuccess : public StateProcessor { State transition() const override { return State::SUCCESS; } };
    struct StateFailure : public StateProcessor { State transition() const override { return State::FAILURE; } };

    Comparator comparator_at(std::size_t idx) const;
    void epsilon_steps_along_ray(
        std::vector<Scalar>& result,
        Vector const& S,
        Vector u,
        Scalar sign = 1.0,
        Scalar const* function_value = nullptr
        ) const;
    bool clip_by_constraints(Vector& u, std::size_t const max_iterations) const;
    bool is_matrix_feasible();

    Config config;
    Constants constants;
    RoundConstants round_constants;
    Sample sample;
    BestIO best_io;

    State state;
    StateRoundBegin state_round_begin;
    StateLocalSpace state_local_space;
    StateConstraints state_constraints;
    StateGradient state_gradient;
    StateFuzzingGradientDescent state_fuzzing_gradient_descent;
    StateFuzzingBitFlips state_fuzzing_bit_flips;
    StateFuzzingRandom state_fuzzing_random;
    StateRoundEnd state_round_end;
    StateSuccess state_success;
    StateFailure state_failure;
    std::unordered_map<State, StateProcessor*> state_processors;

    Vector origin;
    Matrix matrix;
    std::vector<Constraint> constraints;
    Vector gradient;

    Statistics statistics;
};


}

#endif
