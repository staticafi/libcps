#include <cps/variable.hpp>
#include <utility/hash_combine.hpp>

namespace cps {


bool operator==(Variable const& lhs, Variable const& rhs)
{
    if (lhs.type == rhs.type)
        switch (lhs.type)
        {
            case Variable::Type::BOOLEAN:   return lhs.value.b   == rhs.value.b;
            case Variable::Type::UINT8:     return lhs.value.u8  == rhs.value.u8;
            case Variable::Type::SINT8:     return lhs.value.s8  == rhs.value.s8;
            case Variable::Type::UINT16:    return lhs.value.u16 == rhs.value.u16;
            case Variable::Type::SINT16:    return lhs.value.s16 == rhs.value.s16;
            case Variable::Type::UINT32:    return lhs.value.u32 == rhs.value.u32;
            case Variable::Type::SINT32:    return lhs.value.s32 == rhs.value.s32;
            case Variable::Type::UINT64:    return lhs.value.u64 == rhs.value.u64;
            case Variable::Type::SINT64:    return lhs.value.s64 == rhs.value.s64;
            case Variable::Type::FLOAT32:   return lhs.value.f32 == rhs.value.f32;
            case Variable::Type::FLOAT64:   return lhs.value.f64 == rhs.value.f64;
        }
    return false;
}


std::size_t VariableHasher::operator()(Variable const& var) const
{
    std::size_t seed{ std::hash<std::uint8_t>{}(static_cast<std::uint8_t>(var.type)) };
    var.visit([&seed]<typename T>(T const x) { hash_combine(seed, std::hash<std::decay_t<T> >{}(x)); });
    return seed;
}


}
