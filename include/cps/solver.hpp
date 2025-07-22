#ifndef CPS_SOLVER_HPP_INCLUDED
#   define CPS_SOLVER_HPP_INCLUDED

#   include <vector>
#   include <functional>
#   include <cstdint>

namespace cps {


using Valuation = std::vector<std::uint8_t>;


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


enum struct Comparator : std::uint8_t
{
    EQUAL           = 0,
    UNEQUAL         = 1,
    LESS            = 2,
    LESS_EQUAL      = 3,
    GREATER         = 4,
    GREATER_EQUAL   = 5
};


struct Setup
{
    std::vector<Variable> variables;
    std::vector<std::vector<std::size_t> > parameter_indices;
    std::vector<Comparator> comparators;
};


using Real = double;


struct Result
{
    bool valid() const { return !predicate_values.empty() && predicate_values.size() == bb_functions_values.size(); }
    std::vector<Real> bb_functions_values{};
    std::vector<bool> predicate_values{};
};


struct Solver
{
    Solver(Setup const& setup_, Valuation const& seed_valuation_, Result const& seed_result_);

    Setup const& get_setup() const { return setup; }
    Valuation const& get_seed_valuation() const { return seed_valuation; }
    Result const& get_seed_result() const { return seed_result; }

    bool is_finished() const { return finished; }
    bool is_solution(Result const& result) const;

    void compute_next_valuation(Valuation& valuation);
    void process_result(Result const& result);

private:
    Setup setup;
    Valuation seed_valuation;
    Result seed_result;
    bool finished;
};


Result solve(
    Setup const& setup,
    Valuation const& seed_valuation,
    Result const& seed_result,
    std::function<Result(Valuation const&, Setup const&)> const& evaluator
    );


}

#endif
