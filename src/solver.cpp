#include <cps/solver.hpp>
#include <utility/assumptions.hpp>

namespace cps {


Solver::Solver(Setup const& setup_, Valuation const& seed_valuation_, Result const& seed_result_)
    : setup{ setup_ }
    , seed_valuation{ seed_valuation_ }
    , seed_result{ seed_result_ }
    , active_variable_indices{}
    , active_bb_function_indices{}
    , finished{ false }
{
    ASSUMPTION(
        !setup.variables.empty() &&
        !setup.comparators.empty() &&
        setup.comparators.size() == setup.parameter_indices.size() &&
        !setup.parameter_indices.back().empty() &&
        seed_result.bb_functions_values.size() == setup.comparators.size() &&
        seed_result.predicate_values.size() == setup.comparators.size()
    );

    compute_active_indices();
}


void Solver::compute_active_indices()
{
    active_variable_indices.insert(setup.parameter_indices.back().begin(), setup.parameter_indices.back().end());
    active_bb_function_indices.insert(setup.parameter_indices.size() - 1ULL);
    std::unordered_set<std::size_t>  work_set{};
    for (std::size_t  i = 1UL; i < setup.parameter_indices.size(); ++i)
        work_set.insert(i - 1UL);
    while (true)
    {
        bool changed{ false };
        for (auto it = work_set.begin(); it != work_set.end(); )
        {
            auto const& params{ setup.parameter_indices.at(*it) };
            bool intersection{ false };
            for (auto i : params)
                if (active_variable_indices.contains(i))
                {
                    intersection = true;
                    break;
                }
            if (intersection)
            {
                active_variable_indices.insert(params.begin(), params.end());
                active_bb_function_indices.insert(*it);

                changed = true;

                it = work_set.erase(it);
            }
            else
                ++it;
        }
        if (changed == false)
            break;
    }
}


bool Solver::is_solution(Result const& result) const
{
    if (!result.valid() || result.predicate_values.size() < seed_result.predicate_values.size())
        return false;
    std::size_t const last{ seed_result.predicate_values.size() - 1ULL };
    for (std::size_t i{ 0ULL }; i != last; ++i)
        if (result.predicate_values.at(i) != seed_result.predicate_values.at(i))
            return false;
    return result.predicate_values.at(last) != seed_result.predicate_values.at(last);
}


void Solver::compute_next_valuation(Valuation& valuation)
{
    valuation = seed_valuation;
}


void Solver::process_result(Result const& result)
{
    finished = true; // TODO: This is incorrect temporary implementation!
}


Result solve(
    Setup const& setup,
    Valuation const& seed_valuation,
    Result const& seed_result,
    std::function<Result(Valuation const&, Setup const&)> const& evaluator
    )
{
    Solver solver{ setup, seed_valuation, seed_result };
    Valuation valuation;
    while (!solver.is_finished())
    {
        solver.compute_next_valuation(valuation);
        Result result{ evaluator(valuation, solver.get_setup()) };
        if (solver.is_solution(result))
            return result;
        solver.process_result(result);
    }
    return {};
}


}
