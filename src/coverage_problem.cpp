#include <cps/coverage_problem.hpp>

namespace cps {


bool CoverageProblem::valid() const
{
    // TODO: check also parameter indices and also start byte indices of variables.
    return
        !variables.empty() &&
        !comparators.empty() &&
        comparators.size() == parameter_indices.size() &&
        !parameter_indices.back().empty() &&
        output.size() == comparators.size();
}


bool CoverageProblem::is_solution(std::vector<Evaluation> const& candidate) const
{
    if (candidate.size() < output.size())
        return false;
    std::size_t const last{ output.size() - 1ULL };
    for (std::size_t i{ 0ULL }; i != last; ++i)
        if (candidate.at(i).predicate != output.at(i).predicate)
            return false;
    return candidate.at(last).predicate != output.at(last).predicate;
}


void CoverageProblem::compute_active_indices()
{
    active_variable_indices.insert(parameter_indices.back().begin(), parameter_indices.back().end());
    active_bb_function_indices.insert(parameter_indices.size() - 1ULL);
    std::unordered_set<std::size_t>  work_set{};
    for (std::size_t  i = 1UL; i < parameter_indices.size(); ++i)
        work_set.insert(i - 1UL);
    while (true)
    {
        bool changed{ false };
        for (auto it = work_set.begin(); it != work_set.end(); )
        {
            auto const& params{ parameter_indices.at(*it) };
            bool intersection{ false };
            for (auto i : params)
                if (active_variable_indices.contains(i))
                {
                    intersection = true;
                    break;
                }
            if (intersection)
            {
                active_variable_indices.insert(params.begin(), params.end());
                active_bb_function_indices.insert(*it);

                changed = true;

                it = work_set.erase(it);
            }
            else
                ++it;
        }
        if (changed == false)
            break;
    }
}


}
