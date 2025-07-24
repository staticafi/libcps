#ifndef CPS_DETAIL_SOLVER_HPP_INCLUDED
#   define CPS_DETAIL_SOLVER_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/component.hpp>
#   include <cps/math.hpp>
#   include <vector>
#   include <unordered_set>
#   include <functional>
#   include <memory>
#   include <cstdint>

namespace cps {


struct SolverFuzzingInLocalSpace : public Component
{
    SolverFuzzingInLocalSpace(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output
    );

    bool is_finished() const override { return fsm.state == FiniteStateMachine::State::SUCCESS || fsm.state == FiniteStateMachine::State::FAILURE; }
    bool success() const override { return fsm.state == FiniteStateMachine::State::SUCCESS; }
    void compute_next_input(std::vector<Variable>& input) override;
    void process_output(std::vector<Evaluation> const& output) override;

private:

    struct FiniteStateMachine
    {
        enum struct State
        {
            ROUND_BEGIN,
            MATRIX,
            CONSTRAINTS,
            GRADIENT,
            FUZZING,
            ROUND_END,
            SUCCESS,
            FAILURE,
        };

        State state;
        std::unique_ptr<Component> component;
    };

    struct Constants
    {
        std::vector<std::vector<std::size_t> > parameter_indices;
        std::vector<Comparator> comparators;
        std::unordered_set<std::size_t> active_variable_indices;
        std::unordered_set<std::size_t> active_bb_function_indices;
    };

    struct RoundConstants
    {
        std::vector<Variable> seed_input;
        std::vector<Evaluation> seed_output;
    };

    FiniteStateMachine fsm;
    Constants constants;
    RoundConstants round_constants;

    Matrix from_local_to_global_space;
};


}

#endif
