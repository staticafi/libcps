#ifndef CPS_DETAIL_SOLVER_HPP_INCLUDED
#   define CPS_DETAIL_SOLVER_HPP_INCLUDED

#   include <cps/coverage_problem.hpp>
#   include <cps/component.hpp>
#   include <cps/math.hpp>
#   include <cps/component.hpp>
#   include <vector>
#   include <unordered_set>
#   include <functional>
#   include <memory>
#   include <cstdint>

namespace cps {


struct SolverFuzzingInLocalSpace : public Component
{
    SolverFuzzingInLocalSpace(std::shared_ptr<CoverageProblem> const problem_ptr);

    bool is_finished() const override { return state == State::FINISHED; }
    void compute_next_input(std::vector<std::uint8_t>& input) override;
    void process_output(std::vector<Evaluation> const& output) override;

private:

    enum struct State
    {
        ROUND_BEGIN,
        MATRIX,
        CONSTRAINTS,
        GRADIENT,
        FUZZING,
        ROUND_END,
        FINISHED
    };

    void compute_active_indices();

    std::shared_ptr<CoverageProblem> problem;
    State state;
    std::unique_ptr<Component> component;
    Matrix from_local_to_global_space;
};


}

#endif
