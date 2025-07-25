#ifndef CPS_SOLVER_FUZZING_IN_GLOBAL_SPACE_HPP_INCLUDED
#   define CPS_SOLVER_FUZZING_IN_GLOBAL_SPACE_HPP_INCLUDED

#   include <cps/solver_fuzzing_in_local_space.hpp>
#   include <vector>
#   include <algorithm>

namespace cps {


struct SolverFuzzingInGlobalSpace : public SolverFuzzingInLocalSpace
{
    SolverFuzzingInGlobalSpace(
        std::vector<std::vector<std::size_t> > const& parameter_indices,
        std::vector<Comparator> const& comparators,
        std::vector<Variable> const& seed_input,
        std::vector<Evaluation> const& seed_output
        )
        : SolverFuzzingInLocalSpace(parameter_indices, comparators, seed_input, seed_output)
    {}
    virtual ~SolverFuzzingInGlobalSpace() {}

protected:

    void initialize_active_indices() override
    {
        constants.active_variable_indices = constants.parameter_indices.back();
        std::sort(constants.active_variable_indices.begin(), constants.active_variable_indices.end());
        constants.active_function_indices.push_back(constants.parameter_indices.size() - 1ULL);
    }
};


}

#endif
