#ifndef CPS_EVALUATION_CACHE_HPP_INCLUDED
#   define CPS_EVALUATION_CACHE_HPP_INCLUDED

#   include <cps/variable.hpp>
#   include <cps/evaluation.hpp>
#   include <unordered_map>
#   include <vector>
#   include <memory>
#   include <deque>

namespace cps {


struct EvaluationCache
{
    using Input = std::vector<Variable>;
    using Output = std::vector<Evaluation>;

    explicit EvaluationCache(std::size_t max_size_ = 10000UL);

    bool insert(Input const& input, Output const& output);
    Output const* find(Input const& input) const;

    std::size_t size() const { return cache.size(); }
    std::size_t hits() const { return num_hits; }

private:
    struct Key
    {
        struct Hasher { std::size_t operator()(Key key) const; };
        bool operator==(Key other) const;
        Input const* input{ nullptr };
    };

    std::unordered_map<Key, Output, Key::Hasher> cache;
    std::deque<std::unique_ptr<Input> > history;
    std::size_t max_size;
    mutable std::size_t num_hits;
};


}

#endif
