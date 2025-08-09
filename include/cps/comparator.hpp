#ifndef CPS_COMPARATOR_HPP_INCLUDED
#   define CPS_COMPARATOR_HPP_INCLUDED

#   include <cstdint>

namespace cps {


enum struct Comparator : std::uint8_t
{
    EQUAL           = 0,
    UNEQUAL         = 1,
    LESS            = 2,
    LESS_EQUAL      = 3,
    GREATER         = 4,
    GREATER_EQUAL   = 5
};


Comparator opposite(Comparator comparator);


}

#endif
