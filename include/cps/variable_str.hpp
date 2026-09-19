#ifndef CPS_VARIABLE_STR_HPP_INCLUDED
#   define CPS_VARIABLE_STR_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <string>

namespace cps {


std::string type_as_str(Variable const& var);
std::string value_as_str(Variable const& var);


}

#endif
