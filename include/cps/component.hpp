#ifndef CPS_COMPONENT_HPP_INCLUDED
#   define CPS_COMPONENT_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <vector>

namespace cps {


struct Component
{
    virtual ~Component() {}
    bool is_finished() const { return success() || failure(); }
    virtual bool success() const = 0; 
    virtual bool failure() const = 0; 
    virtual std::vector<Variable> const& solution_input() const = 0;
    virtual std::vector<Evaluation> const& solution_output() const = 0;
    virtual void compute_next_input(std::vector<Variable>& input) = 0;
    virtual void process_output(std::vector<Evaluation> const& output) = 0;
};


}

#endif
