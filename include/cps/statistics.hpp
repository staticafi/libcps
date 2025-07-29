#ifndef CPS_STATISTICS_HPP_INCLUDED
#   define CPS_STATISTICS_HPP_INCLUDED

#   include <unordered_map>
#   include <string>
#   include <cstdint>

namespace cps {


using Statistics = std::unordered_map<std::string, std::uint64_t>;


}

#endif
