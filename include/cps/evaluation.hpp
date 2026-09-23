#ifndef CPS_EVALUATION_HPP_INCLUDED
#   define CPS_EVALUATION_HPP_INCLUDED

namespace cps {


struct Evaluation
{
    double function;
    bool predicate;
};


inline bool operator==(Evaluation const& lhs, Evaluation const& rhs)
{ return lhs.predicate == rhs.predicate && lhs.function == rhs.function; }

inline bool operator!=(Evaluation const& lhs, Evaluation const& rhs) { return !(lhs == rhs); }

inline bool operator<(Evaluation const& lhs, Evaluation const& rhs)
{ return lhs.predicate == rhs.predicate ? lhs.function < rhs.function : lhs.predicate < rhs.predicate; }


}

#endif
