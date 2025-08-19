#ifndef CPS_VARIABLE_HPP_INCLUDED
#   define CPS_VARIABLE_HPP_INCLUDED

#   include <cstdint>

namespace cps {


struct Variable
{
    enum struct Type : std::uint8_t
    {
        BOOLEAN     = 0U,
        UINT8       = 1U,
        SINT8       = 2U,
        UINT16      = 3U,
        SINT16      = 4U,
        UINT32      = 5U,
        SINT32      = 6U,
        UINT64      = 7U,
        SINT64      = 8U,
        FLOAT32     = 9U,
        FLOAT64     = 10U,
    };

    union  Value
    {
        bool            b;
        std::uint8_t    u8;
        std::int8_t     s8;
        std::uint16_t   u16;
        std::int16_t    s16;
        std::uint32_t   u32;
        std::int32_t    s32;
        std::uint64_t   u64;
        std::int64_t    s64;
        float           f32;
        double          f64;
    };

    template<typename T> inline Variable(T value_);

    template<typename Visitor> void visit(Visitor visitor);
    template<typename Visitor> void visit(Visitor visitor) const;

    template<typename T> inline void write_value_to(T&& x) const { visit([&x](auto const v) { x = (T)v; }); }

    Type type;
    Value value;
};

template<> inline Variable::Variable(bool const value_)            : type{ Type::BOOLEAN },  value{ .u64 = 0ULL } { value.b   = value_; }
template<> inline Variable::Variable(std::uint8_t const value_)    : type{ Type::UINT8 },    value{ .u64 = 0ULL } { value.u8  = value_; }
template<> inline Variable::Variable(std::int8_t const value_)     : type{ Type::SINT8 },    value{ .u64 = 0ULL } { value.s8  = value_; }
template<> inline Variable::Variable(std::uint16_t const value_)   : type{ Type::UINT16 },   value{ .u64 = 0ULL } { value.u16 = value_; }
template<> inline Variable::Variable(std::int16_t const value_)    : type{ Type::SINT16 },   value{ .u64 = 0ULL } { value.s16 = value_; }
template<> inline Variable::Variable(std::uint32_t const value_)   : type{ Type::UINT32 },   value{ .u64 = 0ULL } { value.u32 = value_; }
template<> inline Variable::Variable(std::int32_t const value_)    : type{ Type::SINT32 },   value{ .u64 = 0ULL } { value.s32 = value_; }
template<> inline Variable::Variable(std::uint64_t const value_)   : type{ Type::UINT64 },   value{ .u64 = 0ULL } { value.u64 = value_; }
template<> inline Variable::Variable(std::int64_t const value_)    : type{ Type::SINT64 },   value{ .u64 = 0ULL } { value.s64 = value_; }
template<> inline Variable::Variable(float const value_)           : type{ Type::FLOAT32 },  value{ .u64 = 0ULL } { value.f32 = value_; }
template<> inline Variable::Variable(double const value_)          : type{ Type::FLOAT64 },  value{ .u64 = 0ULL } { value.f64 = value_; }


template<typename Visitor>
void Variable::visit(Visitor visitor)
{
    switch (type)
    {
        case Type::BOOLEAN:   visitor(value.b); break;
        case Type::UINT8:     visitor(value.u8); break;
        case Type::SINT8:     visitor(value.s8); break;
        case Type::UINT16:    visitor(value.u16); break;
        case Type::SINT16:    visitor(value.s16); break;
        case Type::UINT32:    visitor(value.u32); break;
        case Type::SINT32:    visitor(value.s32); break;
        case Type::UINT64:    visitor(value.u64); break;
        case Type::SINT64:    visitor(value.s64); break;
        case Type::FLOAT32:   visitor(value.f32); break;
        case Type::FLOAT64:   visitor(value.f64); break;
    }
}

template<typename Visitor>
void Variable::visit(Visitor visitor) const
{
    switch (type)
    {
        case Type::BOOLEAN:   visitor(value.b); break;
        case Type::UINT8:     visitor(value.u8); break;
        case Type::SINT8:     visitor(value.s8); break;
        case Type::UINT16:    visitor(value.u16); break;
        case Type::SINT16:    visitor(value.s16); break;
        case Type::UINT32:    visitor(value.u32); break;
        case Type::SINT32:    visitor(value.s32); break;
        case Type::UINT64:    visitor(value.u64); break;
        case Type::SINT64:    visitor(value.s64); break;
        case Type::FLOAT32:   visitor(value.f32); break;
        case Type::FLOAT64:   visitor(value.f64); break;
    }
}


bool operator==(Variable const& lhs, Variable const& rhs);
struct VariableHasher { std::size_t operator()(Variable const& var) const; };


}

#endif
