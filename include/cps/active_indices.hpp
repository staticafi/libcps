#ifndef CPS_ACTIVE_INDICES_HPP_INCLUDED
#   define CPS_ACTIVE_INDICES_HPP_INCLUDED

#   include <cps/solver_fuzzing_in_local_space.hpp>
#   include <vector>

namespace cps {


void compute_active_indices(
    std::vector<std::vector<std::size_t> > const & parameter_indices,
    std::vector<std::size_t>& active_variable_indices,
    std::vector<std::size_t>& active_function_indices
    );


}

#endif
