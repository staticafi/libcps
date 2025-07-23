#ifndef CPS_COVERAGE_PROBLEM_HPP_INCLUDED
#   define CPS_COVERAGE_PROBLEM_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/comparator.hpp>
#   include <cps/evaluation.hpp>
#   include <vector>
#   include <unordered_set>
#   include <cstdint>

namespace cps {


struct CoverageProblem
{
    bool valid() const;
    bool is_solution(std::vector<Evaluation> const& output) const;

    template<typename T>
    T get_input(std::size_t const var_index) const { return *(T*)(input.data() + variables.at(var_index).start_byte_index); }

    template<typename T>
    void set_input(std::size_t const var_index, T const value) { *(T*)(input.data() + variables.at(var_index).start_byte_index) = value; }

    void compute_active_indices();

    std::vector<Variable> variables;
    std::vector<std::vector<std::size_t> > parameter_indices;
    std::vector<Comparator> comparators;
    std::vector<std::uint8_t> input;
    std::vector<Evaluation> output;

    // The following sets should be computed by 'compute_active_indices()' before the problem is passed to a solver.
    std::unordered_set<std::size_t> active_variable_indices;
    std::unordered_set<std::size_t> active_bb_function_indices;
};


}

#endif
