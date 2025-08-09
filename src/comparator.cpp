#include <cps/comparator.hpp>
#include <utility/invariants.hpp>

namespace cps {


Comparator opposite(Comparator const comparator)
{
    switch (comparator)
    {
        case Comparator::EQUAL: return Comparator::UNEQUAL;
        case Comparator::UNEQUAL: return Comparator::EQUAL;
        case Comparator::LESS: return Comparator::GREATER_EQUAL;
        case Comparator::LESS_EQUAL: return Comparator::GREATER;
        case Comparator::GREATER: return Comparator::LESS_EQUAL;
        case Comparator::GREATER_EQUAL: return Comparator::LESS;
        default: UNREACHABLE();
    }
}

}
