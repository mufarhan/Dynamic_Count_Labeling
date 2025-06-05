#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <barrier>
#include <cassert>


#include "road_network.h"

namespace util {

// start new time measurement
void start_timer();
// returns time in seconds since last unconsumed start_timer call and consumes it
double stop_timer();

// sort vector and remove duplicate elements
template<typename T>
void make_set(std::vector<T> &v)
{
    size_t v_size = v.size();
    if (v_size == 0)
        return;
    std::sort(v.begin(), v.end());
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next] != v[last_distinct])
        {
            last_distinct++;
            std::swap(v[next], v[last_distinct]);
        }
    v.resize(last_distinct + 1);
}

template<typename T, class Compare>
void make_set(std::vector<T> &v, Compare comp)
{
    size_t v_size = v.size();
    if (v_size < 2)
        return;
    std::sort(v.begin(), v.end(), comp);
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next].node != v[last_distinct].node)
        {
            last_distinct++;
	    road_network::Neighbor temp = v[next];
	    v[next] = v[last_distinct];
	    v[last_distinct] = temp;
        }
    v.erase(v.begin() + (last_distinct + 1), v.end());
}

// remove elements in set from v
// set must be sorted
template<typename T>
void remove_set(std::vector<T> &v, const std::vector<T> set)
{
    assert(is_sorted(set.cbegin(), set.cend()));
    if (v.empty() || set.empty())
        return;
    std::erase_if(v, [&set](T value) { return std::binary_search(set.cbegin(), set.cend(), value); });
}

struct Summary
{
    double min;
    double max;
    double avg;
    Summary operator*(double x) const;
};

template<typename T, class Map>
Summary summarize(const std::vector<T> &v, Map f)
{
    Summary summary = {};
    if (!v.empty())
        summary.min = f(v[0]);
    for (const T& e : v)
    {
        double x = f(e);
        summary.avg += x;
        if (x < summary.min)
            summary.min = x;
        if (x > summary.max)
            summary.max = x;
    }
    if (!v.empty())
        summary.avg /= v.size();
    return summary;
}

// compute total number of elements in vector of collections
template<typename T>
size_t size_sum(const std::vector<T> &v)
{
    size_t sum = 0;
    for (const T &x : v)
        sum += x.size();
    return sum;
}

// extract size values of vector of collection
template<typename T>
std::vector<size_t> sizes(const std::vector<T> &v)
{
    std::vector<size_t> s;
    for (const T &x : v)
        s.push_back(x.size());
    return s;
}

template<typename T>
T random(const std::vector<T> &v)
{
    assert(v.size() > 0);
    return v[rand() % v.size()];
}

template<typename T>
class min_bucket_queue
{
    std::vector<std::vector<T> > buckets;
    size_t min_bucket; // minimum non-empty bucket
    // mutex for thread synchronization
    //std::mutex m_mutex;
public:
    min_bucket_queue() : min_bucket(0) {}
    void push(T value, size_t bucket)
    {
        //std::lock_guard<std::mutex> lock(m_mutex);
        if (empty() || min_bucket > bucket)
            min_bucket = bucket;
        if (buckets.size() <= bucket)
            buckets.resize(bucket + 1);
        buckets[bucket].push_back(value);
    }
    bool empty() const
    {
        return min_bucket >= buckets.size();
    }
    T pop() {
        assert(min_bucket < buckets.size() && !buckets[min_bucket].empty());
        T top = buckets[min_bucket].back();
        buckets[min_bucket].pop_back();
        // skip empty buckets
        while (min_bucket < buckets.size() && buckets[min_bucket].empty())
            min_bucket++;
        return top;
    }
    size_t size() const
    {
        return buckets.size();
    }
};

template<typename T, size_t threads>
class par_max_bucket_list
{
    std::vector<std::vector<T> > buckets;
    size_t current_bucket = 0, next_in_bucket = 0;
    std::barrier<> sync_point;
    std::mutex m_mutex;
    bool is_empty = true;
    void on_complete()
    {
        assert(!is_empty);
        if (current_bucket == buckets.size() - 1)
        {
            is_empty = true;
            return;
        }
        next_in_bucket = 0;
        // advance to next non-empty bucket
        while (buckets[++current_bucket].empty())
            if (current_bucket == buckets.size() - 1)
            {
                is_empty = true;
                return;
            }
    }
public:
    par_max_bucket_list(size_t max_bucket) : sync_point(threads)
    {
        // prevent re-allocation which could cause syncronization issues
        buckets.reserve(max_bucket + 1);
    }
    void push(T value, size_t bucket)
    {
        if (buckets.size() <= bucket)
        {
            buckets.resize(bucket + 1);
        }
        buckets[bucket].push_back(value);
        is_empty = false;
    }

    bool next(T& value, size_t thread)
    {
        while (true)
        {
            if (is_empty)
                return false;
            {
                // if there's a value available in current bucket, simply return it
                std::lock_guard<std::mutex> lock(m_mutex);
                if (next_in_bucket < buckets[current_bucket].size()) {
                    value = buckets[current_bucket][next_in_bucket++];
                    return true;
                }
            }
            sync_point.arrive_and_wait();
            // ugly workaround for barrier template fuckery
            if (thread == 0)
                on_complete();
            sync_point.arrive_and_wait();
        }
    }
    void reset()
    {
        current_bucket = buckets.empty() ? 0 : buckets.size() - 1;
        next_in_bucket = 0;
        is_empty = false;
    }
};

// thread-safe queue of buckets
template <typename T>
class TSBucketQueue {
private:
    // underlying data structure
    std::vector<std::vector<T> > buckets;
    // minimum non-empty bucket
    size_t min_bucket;
    // mutex for thread synchronization
    std::mutex m_mutex;

    bool empty() const
    {
        return min_bucket >= buckets.size();
    }
public:
    // pushes an element to a given bucket in the queue (NOT thread-safe)
    void push(T item, size_t bucket)
    {
        //std::lock_guard<std::mutex> lock(m_mutex);
        if (empty() || min_bucket > bucket)
            min_bucket = bucket;
        if (buckets.size() <= bucket)
            buckets.resize(bucket + 1);
        buckets[bucket].push_back(item);
    }
    // pops next non-empty bucket off the queue, along with its bucket number; returns whether successful (thread-safe)
    bool next_bucket(std::vector<T>& items, size_t &bucket)
    {
        // restrict mutex lock to critical block (delay copy)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (empty())
                return false;
            bucket = min_bucket++;
            // skip empty buckets
            while (min_bucket < buckets.size() && buckets[min_bucket].empty())
                min_bucket++;
        }
        // copy bucket content
        items = buckets[bucket];
        buckets[bucket].clear();
        return true;
    }
};

} // util

namespace std {

enum class ListFormat { plain, indexed };

void set_list_format(ListFormat format);
ListFormat get_list_format();

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
    if (v.empty())
        return os << "[]";
    if (get_list_format() == ListFormat::indexed)
    {
        os << "[0:" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << i << ":" << v[i];
    }
    else
    {
        os << "[" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << v[i];
    }
    return os << ']';
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A,B> &p)
{
    return os << "(" << p.first << "," << p.second << ")";
}

std::ostream& operator<<(std::ostream& os, util::Summary s);

} // std
