#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

#include <array>
#include <cassert>
#include <iostream>

template <std::size_t, typename T> using alwaysT_t = T;

template<typename T, std::size_t ... Dims>
class MultiArray
{
public:
    const T& operator() (alwaysT_t<Dims, std::size_t>... indexes) const
    {
        return values[computeIndex(indexes...)];
    }
    T& operator() (alwaysT_t<Dims, std::size_t>... indexes)
    {
        return values[computeIndex(indexes...)];
    }

private:
    size_t computeIndex(alwaysT_t<Dims, std::size_t>... indexes_args) const
    {
        constexpr std::size_t dimensions[] = {Dims...};
        std::size_t indexes[] = {indexes_args...};

        size_t index = 0;
        size_t mul = 1;

        for (size_t i = 0; i != sizeof...(Dims); ++i) {
            assert(indexes[i] < dimensions[i]);
            index += indexes[i] * mul;
            mul *= dimensions[i];
        }
        assert(index < (Dims * ...));
        return index;
    }

private:
    std::array<T, (Dims * ...)> values;
};

#endif