#ifndef CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED
#   define CPS_SOLVER_FUZZING_IN_LOCAL_SPACE_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/component.hpp>
#   include <cps/math.hpp>
#   include <vector>
#   include <unordered_set>
#   include <unordered_map>
#   include <memory>

namespace cps {


struct SolverFuzzingInLocalSpace : public Component
{
    SolverFuzzingInLocalSpace(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output
    );
    virtual ~SolverFuzzingInLocalSpace() {}

    bool success() const override { return state == State::SUCCESS; }
    bool failure() const override { return state == State::FAILURE; }

    void compute_next_input(std::vector<Variable>& input) override;
    void process_output(std::vector<Evaluation> const& output_) override;

protected:

    virtual void initialize_active_indices();

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

    struct StateProcessor
    {
        explicit StateProcessor(SolverFuzzingInLocalSpace* const solver = nullptr) : solver_{ solver } {}
        virtual ~StateProcessor() {}
        SolverFuzzingInLocalSpace& solver() { return *solver_; }
        SolverFuzzingInLocalSpace const& solver() const { return *solver_; }
        virtual void enter() {}
        virtual void update() {}
        virtual void update(std::vector<Evaluation> const& output) {}
        virtual State transition() const = 0;
    private:
        SolverFuzzingInLocalSpace* solver_;
    };

    struct StateRoundBegin : public StateProcessor
    {
        explicit StateRoundBegin(SolverFuzzingInLocalSpace* const solver) : StateProcessor{ solver } {}
        void enter() override;
        State transition() const override { return State::LOCAL_SPACE; }
    };

    struct GradientComputationBase : public StateProcessor
    {
        explicit GradientComputationBase(SolverFuzzingInLocalSpace* const solver) : StateProcessor{ solver } {}
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
        explicit StateLocalSpace(SolverFuzzingInLocalSpace* const solver) : GradientComputationBase{ solver } {}
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
        explicit StateConstraints(SolverFuzzingInLocalSpace* const solver) : GradientComputationBase{ solver } {}
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
        explicit StateGradient(SolverFuzzingInLocalSpace* const solver) : GradientComputationBase{ solver } {}
        void enter() override;
        void update() override;
        void update(std::vector<Evaluation> const& output) override;
        State transition() const override;
    };

    struct StateSuccess : public StateProcessor { State transition() const override { return State::SUCCESS; } };
    struct StateFailure : public StateProcessor { State transition() const override { return State::FAILURE; } };

    void update_matrix(Vector const& gradient);
    Scalar compute_epsilon(Vector const& u) const;
    bool are_constraints_satisfied(Vector const& u) const;
    bool clip_by_constraints(Vector& u, std::size_t const max_iterations = 10ULL) const;

    Constants constants;
    RoundConstants round_constants;
    Sample sample;

    State state;
    StateRoundBegin state_round_begin;
    StateLocalSpace state_local_space;
    StateConstraints state_constraints;
    StateGradient state_gradient;
    StateSuccess state_success;
    StateFailure state_failure;
    std::unordered_map<State, StateProcessor*> state_processors;

    Vector origin;
    Matrix matrix;
    std::vector<Constraint> constraints;
    Vector gradient;
};


}

#endif
