#include <cps/evaluation_cache.hpp>
#include <utility/invariants.hpp>
#include <utility/hash_combine.hpp>

namespace cps {


EvaluationCache::EvaluationCache(std::size_t const max_size_)
    : cache{}
    , history{}
    , max_size{ max_size_ }
    , num_hits{ 0UL }
{}


bool EvaluationCache::insert(Input const& input, Output const& output)
{
    while (max_size > 0UL && cache.size() + 1UL > max_size)
    {
        cache.erase({ history.front().get() });
        history.pop_front();
    }
    history.emplace_back(std::make_unique<Input>(input));
    return cache.insert({ { history.back().get() }, output }).second;
}


EvaluationCache::Output const* EvaluationCache::find(Input const& input) const
{
    auto const it{ cache.find({ &input }) };
    if (it == cache.end())
        return nullptr;
    ++num_hits;
    return &it->second;
}


bool EvaluationCache::Key::operator==(Key other) const
{
    return input == other.input || *input == *other.input;
}


std::size_t EvaluationCache::Key::Hasher::operator()(Key key) const
{
    std::size_t seed{ 0UL };
    for (auto const& var : *key.input)
        hash_combine(seed, VariableHasher{}(var));
    return seed;
}


}
