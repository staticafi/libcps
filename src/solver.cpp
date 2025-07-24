#include <cps/solver.hpp>
#include <cps/solver_fuzzing_in_local_space.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>

namespace cps {


Solver::Solver(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Approach const approach
        )
    : solver{ nullptr }
{
    ASSUMPTION(
        !parameter_indices.back().empty() &&
        !comparators.empty() &&
        comparators.size() == parameter_indices.size() &&
        !seed_input.empty() &&
        seed_output.size() == comparators.size()
    );
    switch (approach)
    {
        case Approach::FUZZING_IN_LOCAL_SPACE:
            solver = std::make_unique<SolverFuzzingInLocalSpace>(parameter_indices, comparators, seed_input, seed_output);
            break;
        default: UNREACHABLE(); break;
    }
}


bool solve(
    std::vector<Variable>& solution_input,
    std::vector<Evaluation>& solution_output,
    std::vector<std::vector<std::size_t> > const& parameter_indices,
    std::vector<Comparator> const& comparators,
    std::vector<Variable> const& seed_input,
    std::vector<Evaluation> const& seed_output,
    std::function<void(std::vector<Variable> const&, std::vector<bool> const&, std::vector<Evaluation>&)> const& evaluator,
    Approach const approach
    )
{
    std::vector<bool> predicates;
    for (auto const& eval : seed_output)
        predicates.push_back(eval.predicate);
    predicates.pop_back();

    Solver solver{ parameter_indices, comparators, seed_input, seed_output, approach };
    while (!solver.is_finished())
    {
        solution_input.clear();
        solver.compute_next_input(solution_input);

        solution_output.clear();
        evaluator(solution_input, predicates, solution_output);

        solver.process_output(solution_output);
    }
    return solver.success();
}


}
