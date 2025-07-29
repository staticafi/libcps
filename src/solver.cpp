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
    : solver{ nullptr }
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


bool Solver::is_finished() const { return solver->is_finished(); }
bool Solver::success() const { return solver->success(); }
bool Solver::failure() const { return solver->failure(); }

std::vector<Variable> const& Solver::solution_input() const { return solver->solution_input(); }
std::vector<Evaluation> const& Solver::solution_output() const { return solver->solution_output(); }

void Solver::compute_next_input(std::vector<Variable>& input) { return solver->compute_next_input(input); }
void Solver::process_output(std::vector<Evaluation> const& output) { return solver->process_output(output); }
Statistics const& Solver::get_statistics() const { return solver->get_statistics(); }


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
