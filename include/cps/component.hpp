#ifndef CPS_COMPONENT_HPP_INCLUDED
#   define CPS_COMPONENT_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <vector>

namespace cps {


struct Component
{
    virtual ~Component() {}
    virtual bool is_finished() const = 0; 
    virtual bool success() const = 0; 
    virtual void compute_next_input(std::vector<Variable>& input) = 0;
    virtual void process_output(std::vector<Evaluation> const& output) = 0;
};


}

#endif
