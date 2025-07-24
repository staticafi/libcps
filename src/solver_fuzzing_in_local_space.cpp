#include <cps/solver_fuzzing_in_local_space.hpp>
#include <cps/active_indices.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>

namespace cps {


SolverFuzzingInLocalSpace::SolverFuzzingInLocalSpace(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output
        )
    : Component{}
    , fsm{ FiniteStateMachine::State::ROUND_BEGIN, nullptr }
    , constants{ parameter_indices, comparators, {}, {} }
    , round_constants{ seed_input, seed_output }
    , from_local_to_global_space{}
{
    compute_active_indices(constants.parameter_indices, constants.active_variable_indices, constants.active_bb_function_indices);
}


void SolverFuzzingInLocalSpace::compute_next_input(std::vector<Variable>& input)
{
    input = round_constants.seed_input;
    // while (true)
    // {
    //     if (component != nullptr && component->compute_next_valuation(valuation))
    //         return;

    //     switch (state)
    //     {
    //         case State::ROUND_BEGIN:
    //             from_local_to_global_space.setIdentity(active_variable_indices.size(), active_variable_indices.size());
    //             break;
    //         case State::MATRIX:
    //             break;
    //         case State::CONSTRAINTS:
    //             break;
    //         case State::GRADIENT:
    //             break;
    //         case State::FUZZING:
    //             break;
    //         case State::ROUND_END:
    //             break;
    //         case State::FINISHED:
    //             return;
    //         default: UNREACHABLE(); break;
    //     }
    // }
}


void SolverFuzzingInLocalSpace::process_output(std::vector<Evaluation> const& output)
{
    fsm.state = FiniteStateMachine::State::FAILURE; // TODO: This is incorrect temporary implementation!
}


}
