#ifndef CPS_SOLVER_HPP_INCLUDED
#   define CPS_SOLVER_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <cps/comparator.hpp>
#   include <cps/component.hpp>
#   include <cps/config.hpp>
#   include <vector>
#   include <functional>
#   include <memory>
#   include <cstdint>

namespace cps {


struct Solver : public Component
{
    Solver(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output,
        Config const& config = Config{}
    );

    bool success() const override { return solver->success(); }
    bool failure() const override { return solver->failure(); }
    void compute_next_input(std::vector<Variable>& input) override { return solver->compute_next_input(input); }
    void process_output(std::vector<Evaluation> const& output) override { return solver->process_output(output); }

private:
    std::unique_ptr<Component> solver;
};


bool solve(
    std::vector<Variable>& solution_input,
    std::vector<Evaluation>& solution_output,
    std::vector<std::vector<std::size_t> > const& parameter_indices,
    std::vector<Comparator> const& comparators,
    std::vector<Variable> const& seed_input,
    std::vector<Evaluation> const& seed_output,
    std::function<void(std::vector<Variable> const&, std::vector<bool> const&, std::vector<Evaluation>&)> const& evaluator,
    Config const& config = Config{}
    );


}

#endif
