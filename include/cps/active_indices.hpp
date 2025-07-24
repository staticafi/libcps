#ifndef CPS_ACTIVE_INDICES_HPP_INCLUDED
#   define CPS_ACTIVE_INDICES_HPP_INCLUDED

#   include <cps/solver_fuzzing_in_local_space.hpp>
#   include <vector>
#   include <unordered_set>

namespace cps {


void compute_active_indices(
    std::vector<std::vector<std::size_t> > const & parameter_indices,
    std::unordered_set<std::size_t>& active_variable_indices,
    std::unordered_set<std::size_t>& active_bb_function_indices
    );


}

#endif
