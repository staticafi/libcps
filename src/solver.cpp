#include <cps/solver.hpp>
#include <cps/solver_impl.hpp>
#include <cps/evaluation_cache.hpp>
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
    solver = new SolverImpl(parameter_indices, comparators, seed_input, seed_output, config);
}


Solver::~Solver()
{
    delete solver;
}


bool Solver::is_finished() const { return solver->is_finished(); }
bool Solver::success() const { return solver->success(); }
bool Solver::failure() const { return solver->failure(); }

std::vector<Variable> const& Solver::best_input() const { return solver->best_input(); }
std::vector<Evaluation> const& Solver::best_output() const { return solver->best_output(); }

void Solver::compute_next_input(std::vector<Variable>& input) { return solver->compute_next_input(input); }
void Solver::process_output(std::vector<Evaluation> const& output) { return solver->process_output(output); }
Statistics const& Solver::get_statistics() const { return solver->get_statistics(); }


bool solve(
    std::vector<Variable>& best_input,
    std::vector<Evaluation>& best_output,
    std::vector<std::vector<std::size_t> > const& parameter_indices,
    std::vector<Comparator> const& comparators,
    std::vector<Variable> const& seed_input,
    std::vector<Evaluation> const& seed_output,
    std::function<void(std::vector<Variable> const&, std::function<bool(Evaluation&&)> const&)> const& evaluator,
    Config const& config,
    Statistics* statistics
    )
{
    EvaluationCache cache;
    Solver solver{ parameter_indices, comparators, seed_input, seed_output, config };
    while (!solver.is_finished())
    {
        best_input.clear();
        solver.compute_next_input(best_input);

        best_output.clear();
        EvaluationCache::Output const* const cached_output{ cache.find(best_input) };
        if (cached_output != nullptr)
            best_output = *cached_output;
        else
        {
            evaluator(best_input, [&seed_output, &best_output](Evaluation&& e) {
                best_output.emplace_back(e);
                return best_output.size() < seed_output.size() &&
                    best_output.back().predicate == seed_output.at(best_output.size() - 1ULL).predicate;
            });
            cache.insert(best_input, best_output);
        }

        solver.process_output(best_output);
    }

    best_input = solver.best_input();
    best_output = solver.best_output();

    if (statistics != nullptr)
        *statistics = solver.get_statistics();
    statistics->insert({ "CACHE_HITS", cache.hits() });

    return solver.success();
}


}
