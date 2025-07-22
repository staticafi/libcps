#ifndef CPS_SOLVER_HPP_INCLUDED
#   define CPS_SOLVER_HPP_INCLUDED

#   include <cps/setup.hpp>
#   include <cps/valuation.hpp>
#   include <cps/result.hpp>
#   include <vector>
#   include <unordered_set>
#   include <functional>
#   include <cstdint>

namespace cps {


struct Solver
{
    Solver(Setup const& setup_, Valuation const& seed_valuation_, Result const& seed_result_);

    Setup const& get_setup() const { return setup; }
    Valuation const& get_seed_valuation() const { return seed_valuation; }
    Result const& get_seed_result() const { return seed_result; }
    std::unordered_set<std::size_t> const& get_active_variable_indices() const { return active_variable_indices; }
    std::unordered_set<std::size_t> const& get_active_bb_function_indices() const { return active_bb_function_indices; }

    bool is_finished() const { return finished; }
    bool is_solution(Result const& result) const;

    void compute_next_valuation(Valuation& valuation);
    void process_result(Result const& result);

private:
    void compute_active_indices();

    Setup setup;
    Valuation seed_valuation;
    Result seed_result;
    std::unordered_set<std::size_t> active_variable_indices;
    std::unordered_set<std::size_t> active_bb_function_indices;
    bool finished;
};


Result solve(
    Setup const& setup,
    Valuation const& seed_valuation,
    Result const& seed_result,
    std::function<Result(Valuation const&, Setup const&)> const& evaluator
    );


}

#endif
