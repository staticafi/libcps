#ifndef CPS_SETUP_HPP_INCLUDED
#   define CPS_SETUP_HPP_INCLUDED

#   include <cps/comparator.hpp>
#   include <cps/variable.hpp>
#   include <vector>
#   include <cstdint>

namespace cps {


struct Setup
{
    std::vector<Variable> variables;
    std::vector<std::vector<std::size_t> > parameter_indices;
    std::vector<Comparator> comparators;
};


}

#endif
