#include <cps/comparator_str.hpp>
#include <utility/invariants.hpp>

namespace cps {


std::string as_str(Comparator const comparator)
{
    switch (comparator)
    {
        case Comparator::EQUAL:         return "EQUAL";
        case Comparator::UNEQUAL:       return "UNEQUAL";
        case Comparator::LESS:          return "LESS";
        case Comparator::LESS_EQUAL:    return "LESS_EQUAL";
        case Comparator::GREATER:       return "GREATER";
        case Comparator::GREATER_EQUAL: return "GREATER_EQUAL";
        default: UNREACHABLE();         return "UNKNOWN";
    }
}


}
