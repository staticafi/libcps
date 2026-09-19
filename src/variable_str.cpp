#include <cps/variable_str.hpp>
#include <utility/invariants.hpp>

namespace cps {


std::string type_as_str(Variable const& var)
{
    switch (var.type)
    {
        case Variable::Type::BOOLEAN:   return "BOOLEAN";
        case Variable::Type::UINT8:     return "UINT8";
        case Variable::Type::SINT8:     return "SINT8";
        case Variable::Type::UINT16:    return "UINT16";
        case Variable::Type::SINT16:    return "SINT16";
        case Variable::Type::UINT32:    return "UINT32";
        case Variable::Type::SINT32:    return "SINT32";
        case Variable::Type::UINT64:    return "UINT64";
        case Variable::Type::SINT64:    return "SINT64";
        case Variable::Type::FLOAT32:   return "FLOAT32";
        case Variable::Type::FLOAT64:   return "FLOAT64";
        default: UNREACHABLE();         return "UNKNOWN";
    }
}


std::string value_as_str(Variable const& var)
{
    std::string s;
    var.visit([&s](auto const v) { s = std::to_string(v); });
    return s;
}


}
