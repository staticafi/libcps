#include <cps/active_indices.hpp>
#include <set>
#include <unordered_set>

namespace cps {


void compute_active_indices(
    std::vector<std::vector<std::size_t> > const & parameter_indices,
    std::vector<std::size_t>& active_variable_indices,
    std::vector<std::size_t>& active_function_indices
    )
{
    std::set<std::size_t> variable_indices;
    std::set<std::size_t> function_indices;
    variable_indices.insert(parameter_indices.back().begin(), parameter_indices.back().end());
    function_indices.insert(parameter_indices.size() - 1ULL);
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
                if (variable_indices.contains(i))
                {
                    intersection = true;
                    break;
                }
            if (intersection)
            {
                variable_indices.insert(params.begin(), params.end());
                function_indices.insert(*it);

                changed = true;

                it = work_set.erase(it);
            }
            else
                ++it;
        }
        if (changed == false)
            break;
    }
    active_variable_indices.assign(variable_indices.begin(), variable_indices.end());
    active_function_indices.assign(function_indices.begin(), function_indices.end());
}


}
