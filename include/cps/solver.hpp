#ifndef CPS_SOLVER_HPP_INCLUDED
#   define CPS_SOLVER_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/config.hpp>
#   include <cps/statistics.hpp>
#   include <vector>
#   include <functional>
#   include <memory>
#   include <cstdint>

namespace cps {


struct SolverImpl;


struct Solver
{
    Solver(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config = Config{}
    );

    bool is_finished() const;
    bool success() const;
    bool failure() const;

    std::vector<Variable> const& solution_input() const;
    std::vector<Evaluation> const& solution_output() const;

    void compute_next_input(std::vector<Variable>& input);
    void process_output(std::vector<Evaluation> const& output);

    Statistics const& get_statistics() const;

private:
    std::unique_ptr<SolverImpl> solver;
};


bool solve(
    std::vector<Variable>& solution_input,
    std::vector<Evaluation>& solution_output,
    std::vector<std::vector<std::size_t> > const& parameter_indices,
    std::vector<Comparator> const& comparators,
    std::vector<Variable> const& seed_input,
    std::vector<Evaluation> const& seed_output,
    std::function<void(std::vector<Variable> const&, std::vector<bool> const&, std::vector<Evaluation>&)> const& evaluator,
    Config const& config = Config{},
    Statistics* statistics = nullptr
    );


}

#endif
