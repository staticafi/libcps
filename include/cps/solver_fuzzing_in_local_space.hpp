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
        std::vector<std::vector<std::size_t> > parameter_indices;
        std::vector<Comparator> comparators;
        std::vector<std::size_t> active_variable_indices;
        std::vector<std::size_t> active_function_indices;
    };

    struct RoundConstants
    {
        std::vector<Variable> seed_input;
        std::vector<Evaluation> seed_output;
    };

    struct Sample
    {
        bool ready;
        Vector vector;
    };

    struct Constraint
    {
        Vector normal;
        Scalar distance;
        Comparator comparator;
    };

    struct GradientComputationBase
    {
        enum struct Step
        {
            POSITIVE,
            NEGATIVE,
            STOP,
        };

        std::size_t column_index{ 0ULL };
        Step step{ Step::POSITIVE };
        Scalar epsilon{ 0.0 };
    };

    void StateRoundBegin_update();
    State StateRoundBegin_transition();

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
