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
        Scalar distance{ 0.0 };
        Comparator comparator{ Comparator::EQUAL };
    };

    void StateRoundBegin_update();
    State StateRoundBegin_transition();

    struct GradientComputationBase
    {
        void reset_gradient_computation() { column_index = 0ULL; current_coeff = 0.0; step_coeffs = STEP_COEFFS; }
        void reset_for_next_partial() { ++column_index; current_coeff = 0.0; step_coeffs = STEP_COEFFS; }
        Vector compute_partial_step_vector(SolverFuzzingInLocalSpace const* const solver);
        Scalar compute_finite_difference(Scalar const f_step, Scalar const f) const { return (f_step - f) / current_coeff; }
        std::size_t column_index{ 0ULL };
        Scalar current_coeff{};
        std::vector<Scalar> step_coeffs{};
        inline static std::vector<Scalar> const STEP_COEFFS{ -1.0, 1.0 };
    };

    struct StateLocalSpace : public GradientComputationBase
    {
        std::size_t active_function_index{ 0ULL };
        Vector gradient{};
    };
    State StateLocalSpace_enter();
    void StateLocalSpace_update();
    void StateLocalSpace_update(std::vector<Evaluation> const& output);
    State StateLocalSpace_transition();

    struct StateConstraints : public GradientComputationBase
    {
        std::unordered_map<std::size_t, Vector> gradients{};
        std::unordered_set<std::size_t> partial_function_indices{};
    };
    State StateConstraints_enter();
    void StateConstraints_update();
    void StateConstraints_update(std::vector<Evaluation> const& output);
    State StateConstraints_transition();

    struct StateGradient : public GradientComputationBase {};
    State StateGradient_enter();
    void StateGradient_update();
    void StateGradient_update(std::vector<Evaluation> const& output);
    State StateGradient_transition();

    void updateMatrix(Vector const& gradient);
    Scalar computeEpsilon(Vector const& u) const;

    Constants constants;
    RoundConstants round_constants;
    Sample sample;

    State state;
    StateLocalSpace state_local_space;
    StateConstraints state_constraints;
    StateGradient state_gradient;

    Vector origin;
    Matrix matrix;
    std::vector<Constraint> constraints;
    Vector gradient;
};


}

#endif
