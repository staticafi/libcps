#ifndef CPS_VARIABLE_HPP_INCLUDED
#   define CPS_VARIABLE_HPP_INCLUDED

#   include <variant>
#   include <cstdint>

namespace cps {


using Variable = std::variant<
    bool,
    std::uint8_t,
    std::int8_t,
    std::uint16_t,
    std::int16_t,
    std::uint32_t,
    std::int32_t,
    std::uint64_t,
    std::int64_t,
    float,
    double
>;


}

#endif
