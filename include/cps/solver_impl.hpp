#ifndef CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED
#   define CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/config.hpp>
#   include <cps/math.hpp>
#   include <cps/statistics.hpp>
#   include <vector>
#   include <unordered_set>
#   include <unordered_map>
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

    std::vector<Variable> const& solution_input() const { return best_io.input; }
    std::vector<Evaluation> const& solution_output() const { return best_io.output; }

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
        FUZZING_BIT_MUTATIONS,
        FUZZING_RANDOM_MUTATIONS,
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
        State transition() const override { return State::LOCAL_SPACE; }
    };

    struct GradientComputationBase : public StateProcessor
    {
        explicit GradientComputationBase(SolverImpl* const solver) : StateProcessor{ solver } {}
        void reset_gradient_computation() { column_index = 0ULL; current_coeff = 0.0; step_coeffs = STEP_COEFFS; }
        void reset_for_next_partial() { ++column_index; current_coeff = 0.0; step_coeffs = STEP_COEFFS; }
        Vector compute_partial_step_vector();
        Scalar compute_finite_difference(Scalar const f_step, Scalar const f) const { return (f_step - f) / current_coeff; }
    protected:
        std::size_t column_index{ 0ULL };
        Scalar current_coeff{};
        std::vector<Scalar> step_coeffs{};
        inline static std::vector<Scalar> const STEP_COEFFS{ -1.0, 1.0 };
    };

    struct StateLocalSpace : public GradientComputationBase
    {
        explicit StateLocalSpace(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    private:
        std::size_t active_function_index{ 0ULL };
        Vector gradient{};
    };

    struct StateConstraints : public GradientComputationBase
    {
        explicit StateConstraints(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    private:
        std::unordered_map<std::size_t, Vector> gradients{};
        std::unordered_set<std::size_t> partial_function_indices{};
    };

    struct StateGradient : public GradientComputationBase
    {
        explicit StateGradient(SolverImpl* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    };

    struct StateFuzzingGradientDescent : public StateProcessor
    {
        explicit StateFuzzingGradientDescent(SolverImpl* const solver) : StateProcessor{ solver } {}
        void enter() override;
        void update() override;
        State transition() const override;
    private:
        Scalar multiplier{};
        std::vector<Scalar> multipliers{};
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

    void update_matrix(Vector const& gradient);
    Scalar epsilon_step_along_vector(Vector const& u) const;
    bool are_constraints_satisfied(Vector const& u) const;
    bool clip_by_constraints(Vector& u, std::size_t const max_iterations = 10ULL) const;

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
