#include <cps/solver.hpp>
#include <cps/solver_impl.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>

namespace cps {


Solver::Solver(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config
        )
    : Component{}
    , solver{ nullptr }
{
    ASSUMPTION(
        !parameter_indices.back().empty() &&
        !comparators.empty() &&
        comparators.size() == parameter_indices.size() &&
        !seed_input.empty() &&
        seed_output.size() == comparators.size()
    );
    solver = std::make_unique<SolverImpl>(parameter_indices, comparators, seed_input, seed_output, config);
}


bool solve(
    std::vector<Variable>& solution_input,
    std::vector<Evaluation>& solution_output,
    std::vector<std::vector<std::size_t> > const& parameter_indices,
    std::vector<Comparator> const& comparators,
    std::vector<Variable> const& seed_input,
    std::vector<Evaluation> const& seed_output,
    std::function<void(std::vector<Variable> const&, std::vector<bool> const&, std::vector<Evaluation>&)> const& evaluator,
    Config const& config,
    Statistics* statistics
    )
{
    std::vector<bool> predicates;
    for (auto const& eval : seed_output)
        predicates.push_back(eval.predicate);
    predicates.pop_back();

    Solver solver{ parameter_indices, comparators, seed_input, seed_output, config };
    while (!solver.is_finished())
    {
        solution_input.clear();
        solver.compute_next_input(solution_input);

        solution_output.clear();
        evaluator(solution_input, predicates, solution_output);

        solver.process_output(solution_output);
    }

    if (statistics != nullptr)
        *statistics = solver.get_statistics();

    return solver.success();
}


}
