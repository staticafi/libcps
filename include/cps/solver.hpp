#ifndef CPS_SOLVER_HPP_INCLUDED
#   define CPS_SOLVER_HPP_INCLUDED

#   include <cps/coverage_problem.hpp>
#   include <cps/component.hpp>
#   include <vector>
#   include <functional>
#   include <memory>
#   include <cstdint>

namespace cps {


enum struct Approach : std::uint8_t
{
    FUZZING_IN_LOCAL_SPACE = 0,
};


struct Solver
{
    Solver(CoverageProblem const& problem_, Approach approach = Approach::FUZZING_IN_LOCAL_SPACE);

    CoverageProblem const& coverage_problem() const { return *problem; }

    bool is_finished() const { return solver->is_finished(); }
    void compute_next_input(std::vector<std::uint8_t>& input) { return solver->compute_next_input(input); }
    void process_output(std::vector<Evaluation> const& output) { return solver->process_output(output); }

private:
    std::shared_ptr<CoverageProblem> problem;
    std::unique_ptr<Component> solver;
};


bool solve(
    std::vector<std::uint8_t>& solution_input,
    std::vector<Evaluation>& solution_output,
    CoverageProblem const& problem,
    std::function<void(std::vector<std::uint8_t> const&, CoverageProblem const&, std::vector<Evaluation>&)> const& evaluator,
    Approach approach = Approach::FUZZING_IN_LOCAL_SPACE
    );


}

#endif
