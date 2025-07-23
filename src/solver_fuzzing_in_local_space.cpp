#include <cps/solver_fuzzing_in_local_space.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>

namespace cps {


SolverFuzzingInLocalSpace::SolverFuzzingInLocalSpace(std::shared_ptr<CoverageProblem> const problem_ptr)
    : Component{}
    , problem{ problem_ptr }
    , state{ State::ROUND_BEGIN }
    , component{ nullptr }
    , from_local_to_global_space{}
{
}


void SolverFuzzingInLocalSpace::compute_next_input(std::vector<std::uint8_t>& input)
{
    input = problem->input;
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
    state = State::FINISHED; // TODO: This is incorrect temporary implementation!
}


}
