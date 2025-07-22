#ifndef CPS_VARIABLE_HPP_INCLUDED
#   define CPS_VARIABLE_HPP_INCLUDED

#   include <cstdint>
#   include <utility>

namespace cps {


struct Variable
{
    enum struct Type : std::uint8_t
    {
        BOOLEAN = 0U,

        UINT8 = 1U,
        SINT8 = 2U,

        UINT16 = 3U,
        SINT16 = 4U,

        UINT32 = 5U,
        SINT32 = 6U,

        UINT64 = 7U,
        SINT64 = 8U,

        FLOAT32 = 9U,
        FLOAT64 = 10U,
    };

    std::size_t start_byte_index;
    Type type;
};


}

#endif
