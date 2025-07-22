#ifndef CPS_RESULT_HPP_INCLUDED
#   define CPS_RESULT_HPP_INCLUDED

#   include <vector>

namespace cps {


struct Result
{
    bool valid() const { return !predicate_values.empty() && predicate_values.size() == bb_functions_values.size(); }
    std::vector<double> bb_functions_values{};
    std::vector<bool> predicate_values{};
};


}

#endif
