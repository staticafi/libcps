#ifndef CPS_CONFIG_HPP_INCLUDED
#   define CPS_CONFIG_HPP_INCLUDED

namespace cps {


struct Config
{
    bool build_local_space{ true };
    bool build_constraints{ true };
    bool use_gradient_descent{ true };
    bool use_bit_flips{ true };
    bool use_random_fuzzing{ true };
};


}

#endif
