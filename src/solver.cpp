#include <cps/solver.hpp>
#include <cps/solver_fuzzing_in_local_space.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>

namespace cps {


Solver::Solver(CoverageProblem const& problem_, Approach const approach)
    : problem{ std::make_shared<CoverageProblem>(problem_) }
    , solver{ nullptr }
{
    ASSUMPTION(problem->valid());
    switch (approach)
    {
        case Approach::FUZZING_IN_LOCAL_SPACE:
            solver = std::make_unique<SolverFuzzingInLocalSpace>(problem);
            break;
        default: UNREACHABLE(); break;
    }
}


bool solve(
    std::vector<std::uint8_t>& solution_input,
    std::vector<Evaluation>& solution_output,
    CoverageProblem const& problem,
    std::function<void(std::vector<std::uint8_t> const&, CoverageProblem const&, std::vector<Evaluation>&)> const& evaluator,
    Approach const approach
    )
{
    Solver solver{ problem, approach };
    while (!solver.is_finished())
    {
        solution_input.clear();
        solver.compute_next_input(solution_input);

        solution_output.clear();
        evaluator(solution_input, solver.coverage_problem(), solution_output);

        if (solver.coverage_problem().is_solution(solution_output))
            return true;

        solver.process_output(solution_output);
    }
    return false;
}


}
