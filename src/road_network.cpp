#include "road_network.h"
#include "util.h"

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <bitset>
#include <unordered_set>
#include <thread>
#include <atomic>
#include <cstring>
#include <random>

using namespace std;

#define DEBUG(X) //cerr << X << endl

// algorithm config
//#define CUT_REPEAT 3 // repeat whole cut computation multiple times (with random starting points for rough partition) and pick best result
#define MULTI_CUT // extract two different min-cuts from max-flow and pick more balanced result
static const bool weighted_furthest = false; // use edge weights for finding distant nodes during rough partitioning
static const bool weighted_diff = false; // use edge weights for computing rough partition

namespace road_network {

static const NodeID NO_NODE = 0; // null value equivalent for integers identifying nodes
static const SubgraphID NO_SUBGRAPH = 0; // used to indicate that node does not belong to any active subgraph
static const uint16_t MAX_CUT_LEVEL = 58; // maximum height of decomposition tree; 58 bits to store binary path, plus 6 bits to store path length = 64 bit integer

// profiling
#ifndef NPROFILE
    static atomic<double> t_partition, t_label, t_shortcut;
    #define START_TIMER util::start_timer()
    #define STOP_TIMER(var) var += util::stop_timer()
#else
    #define START_TIMER
    #define STOP_TIMER(var)
#endif

// progress of 0 resets counter
static bool log_progress_on = false;
void log_progress(size_t p, ostream &os = cout)
{
    static const size_t P_DIFF = 1000000L;
    static size_t progress = 0;
    if (log_progress_on)
    {
        size_t old_log = progress / P_DIFF;
        if (p == 0)
        {
            // terminate progress line & reset
            if (old_log > 0)
                os << endl;
            progress = 0;
            return;
        }
        progress += p;
        size_t new_log = progress / P_DIFF;
        if (old_log < new_log)
        {
            for (size_t i = old_log; i < new_log; i++)
                os << '.';
            os << flush;
        }
    }
    else
        progress += p;
}

// half-matrix index for storing half-matrix in flat vector
static size_t hmi(size_t a, size_t b)
{
    assert(a != b);
    return a < b ? (b * (b - 1) >> 1) + a : (a * (a - 1) >> 1) + b;
}

// offset by cut level
uint16_t get_offset(const uint16_t *dist_index, size_t cut_level)
{
    return cut_level ? dist_index[cut_level - 1] : 0;
}

//--------------------------- CutIndex ------------------------------

CutIndex::CutIndex() : partition(0), cut_level(0)
{
#ifdef PRUNING
    pruning_2hop = pruning_3hop = pruning_tail = 0;
#endif
}

#ifdef PRUNING
void CutIndex::prune_tail()
{
    assert(is_consistent(true));
    assert(dist_index.back() == distances.size());
    assert(distances.size() > 0);
    // cut_level may not be set yet
    size_t cl = dist_index.size() - 1;
    // only prune latest cut
    size_t last_unpruned = get_offset(&dist_index[0], cl);
    // nothing to prune for empty cuts
    if (last_unpruned == distances.size())
        return;
    // first node must never be pruned
    assert(distances[last_unpruned] & 1);
    // fix distances and recall last unprunded label
    for (size_t i = last_unpruned; i < distances.size(); i++)
    {
        if (distances[i] & 1)
            last_unpruned = i;
        distances[i] >>= 1;
    }
    size_t new_size = last_unpruned + 1;
    assert(new_size <= distances.size());
    if (new_size < distances.size())
    {
        pruning_tail += distances.size() - new_size;
        distances.resize(new_size);
        dist_index.back() = new_size;
        DEBUG("pruned tail: " << *this);
    }
}
#endif

bool CutIndex::is_consistent(bool partial) const
{
    const uint64_t one = 1;
    if (cut_level > MAX_CUT_LEVEL)
    {
        cerr << "cut_level=" << (int)cut_level << endl;
        return false;
    }
    if (!partial && partition >= (one << cut_level))
    {
        cerr << "partition=" << partition << " for cut_level=" << (int)cut_level << endl;
        return false;
    }
    if (!partial && dist_index.size() != cut_level + one)
    {
        cerr << "dist_index.size()=" << dist_index.size() << " for cut_level=" << (int)cut_level << endl;
        return false;
    }
    if (!is_sorted(dist_index.cbegin(), dist_index.cend()))
    {
        cerr << "unsorted dist_index: " << dist_index << endl;
        return false;
    }
    return true;
}

bool CutIndex::empty() const
{
    return dist_index.empty();
}

// need to implement distance calculation (for given cut level) using CutIndex as it's used to identify redundant shortcuts
static distance_t get_cut_level_distance(const CutIndex &a, const CutIndex &b, size_t cut_level)
{
    distance_t min_dist = infinity;
    uint16_t a_offset = get_offset(&a.dist_index[0], cut_level);
    uint16_t b_offset = get_offset(&b.dist_index[0], cut_level);
    const distance_t* a_ptr = &a.distances[0] + a_offset;
    const distance_t* b_ptr = &b.distances[0] + b_offset;
    const distance_t* a_end = a_ptr + min(a.dist_index[cut_level] - a_offset, b.dist_index[cut_level] - b_offset);
    // find min 2-hop distance within partition
    while (a_ptr != a_end)
    {
        distance_t dist = *a_ptr + *b_ptr;
        if (dist < min_dist)
            min_dist = dist;
        a_ptr++;
        b_ptr++;
    }
    return min_dist;
}

//--------------------------- PBV -----------------------------------

namespace PBV
{

uint64_t from(uint64_t bits, uint16_t length)
{
    if (length == 0)
        return 0;
    return (bits << (64 - length) >> (58 - length)) | length;
}

uint64_t partition(uint64_t bv)
{
    // cutlevel is stored in lowest 6 bits
    return bv >> 6;
}

uint16_t cut_level(uint64_t bv)
{
    // cutlevel is stored in lowest 6 bits
    return bv & 63ul;
}

uint16_t lca_level(uint64_t bv1, uint64_t bv2)
{
    // find lowest level at which partitions differ
    uint16_t lca_level = min(cut_level(bv1), cut_level(bv2));
    uint64_t p1 = partition(bv1), p2 = partition(bv2);
    if (p1 != p2)
    {
        uint16_t diff_level = __builtin_ctzll(p1 ^ p2); // count trailing zeros
        if (diff_level < lca_level)
            lca_level = diff_level;
    }
    return lca_level;
}

uint64_t lca(uint64_t bv1, uint64_t bv2)
{
    uint64_t cut_level = lca_level(bv1, bv2);
    // shifting by 64 does not work
    if (cut_level == 0)
        return 0;
    return (bv1 >> 6) << (64 - cut_level) >> (58 - cut_level) | cut_level;
}

bool is_ancestor(uint64_t bv_ancestor, uint64_t bv_descendant)
{
    uint16_t cla = cut_level(bv_ancestor), cld = cut_level(bv_descendant);
    // shifting by 64 does not work, so need to check for cla == 0
    return cla == 0 || (cla <= cld && (bv_ancestor ^ bv_descendant) >> 6 << (64 - cla) == 0);
}

}

//--------------------------- FlatCutIndex --------------------------

// helper function for memory alignment
template<typename T>
size_t aligned(size_t size);

template<>
size_t aligned<uint32_t>(size_t size)
{
    size_t mod = size & 3ul;
    return mod ? size + (4 - mod) : size;
}

FlatCutIndex::FlatCutIndex() : data(nullptr)
{
}

FlatCutIndex::FlatCutIndex(const CutIndex &ci)
{
    assert(ci.is_consistent());
    // allocate memory for partition bitvector, dist_index, paths count and distances
    size_t data_size = sizeof(uint64_t) + aligned<distance_t>(ci.dist_index.size() * sizeof(uint16_t)) + ci.distances.size() * sizeof(distance_t) + ci.paths.size() * sizeof(uint16_t);
    data = (char*)calloc(data_size, 1);
    // copy partition bitvector, dist_index, paths count and distances into data
    *partition_bitvector() = PBV::from(ci.partition, ci.cut_level); 
    memcpy(dist_index(), &ci.dist_index[0], ci.dist_index.size() * sizeof(uint16_t));  
    memcpy(distances(), &ci.distances[0], ci.distances.size() * sizeof(distance_t)); 
    memcpy(paths(), &ci.paths[0], ci.paths.size() * sizeof(uint16_t));
}

bool FlatCutIndex::operator==(FlatCutIndex other) const
{
    return data == other.data;
}

uint64_t* FlatCutIndex::partition_bitvector()
{
    assert(!empty());
    return (uint64_t*)data;
}

const uint64_t* FlatCutIndex::partition_bitvector() const
{
    assert(!empty());
    return (uint64_t*)data;
}

uint16_t* FlatCutIndex::dist_index()
{
    assert(!empty());
    return (uint16_t*)(data + sizeof(uint64_t));
}

const uint16_t* FlatCutIndex::dist_index() const
{
    assert(!empty());
    return (uint16_t*)(data + sizeof(uint64_t));
}

uint16_t* FlatCutIndex::paths()
{
    assert(!empty());
    return (uint16_t*) (data + sizeof(uint64_t) + aligned<distance_t>((cut_level() + 1) * sizeof(uint16_t)) + label_count() * sizeof(distance_t));
}

const uint16_t* FlatCutIndex::paths() const
{
    assert(!empty());
    return (uint16_t*)(data + sizeof(uint64_t) + aligned<distance_t>((cut_level() + 1) * sizeof(uint16_t)) + label_count() * sizeof(distance_t));
}

distance_t* FlatCutIndex::distances()
{
    assert(!empty());
    return (distance_t*)(data + sizeof(uint64_t) + aligned<distance_t>((cut_level() + 1) * sizeof(uint16_t)));
}

const distance_t* FlatCutIndex::distances() const
{
    assert(!empty());
    return (distance_t*)(data + sizeof(uint64_t) + aligned<distance_t>((cut_level() + 1) * sizeof(uint16_t)));
}

uint64_t FlatCutIndex::partition() const
{
    return PBV::partition(*partition_bitvector());
}

uint16_t FlatCutIndex::cut_level() const
{
    return PBV::cut_level(*partition_bitvector());
}

size_t FlatCutIndex::size() const
{
    size_t total = sizeof(uint64_t);
    total += aligned<distance_t>((cut_level() + 1) * sizeof(uint16_t));
    total += dist_index()[cut_level()] * sizeof(distance_t);
    total += dist_index()[cut_level()] * sizeof(uint16_t);
    return total;
}

size_t FlatCutIndex::label_count() const
{
    return dist_index()[cut_level()];
}

size_t FlatCutIndex::cut_size(size_t cl) const
{
    return cl == 0 ? *dist_index() : dist_index()[cl] - dist_index()[cl - 1];
}

size_t FlatCutIndex::bottom_cut_size() const
{
    return cut_size(cut_level());
}

bool FlatCutIndex::empty() const
{
    return data == nullptr;
}

const distance_t* FlatCutIndex::cl_begin(size_t cl) const
{
    return distances() + get_offset(dist_index(), cl);
}

const distance_t* FlatCutIndex::cl_end(size_t cl) const
{
    return distances() + dist_index()[cl];
}

const uint16_t* FlatCutIndex::pl_begin(size_t cl) const
{
    return paths() + get_offset(dist_index(), cl);
}

const uint16_t* FlatCutIndex::pl_end(size_t cl) const
{
    return paths() + dist_index()[cl];
}

vector<vector<distance_t> > FlatCutIndex::unflatten() const
{
    vector<vector<distance_t>> labels;
    for (size_t cl = 0; cl <= cut_level(); cl++)
    {
        vector<distance_t> cut_labels;
        for (const distance_t *l = cl_begin(cl); l != cl_end(cl); l++)
            cut_labels.push_back(*l);
        labels.push_back(cut_labels);
    }
    return labels;
}

vector<vector<pair<distance_t, uint16_t> > > FlatCutIndex::unflatten_spc() const
{
    vector<vector<pair<distance_t, uint16_t> > > labels;
    for (size_t cl = 0; cl <= cut_level(); cl++)
    {
        vector<uint16_t> path_labels;
        for (const uint16_t *l = pl_begin(cl); l != pl_end(cl); l++)
            path_labels.push_back(*l);
	vector<distance_t> distance_labels;
        for (const distance_t *l = cl_begin(cl); l != cl_end(cl); l++)
            distance_labels.push_back(*l);

	vector<pair<distance_t,uint16_t> > label_pairs;
	for(size_t i = 0; i < distance_labels.size(); i++)
            label_pairs.push_back(make_pair(distance_labels[i], path_labels[i]));
	labels.push_back(label_pairs);
    }
    return labels;
}

//--------------------------- ContractionLabel ----------------------

ContractionLabel::ContractionLabel() : distance_offset(0), parent(NO_NODE)
{
}

size_t ContractionLabel::size() const
{
    size_t total = sizeof(ContractionLabel);
    // only count index data if owned
    if (distance_offset == 0)
        total += cut_index.size();
    return total;
}

//--------------------------- ContractionIndex ----------------------

template<typename T>
static void clear_and_shrink(vector<T> &v)
{
    v.clear();
    v.shrink_to_fit();
}

ContractionIndex::ContractionIndex(vector<CutIndex> &ci, vector<Neighbor> &closest)
{
    assert(ci.size() == closest.size());
    labels.resize(ci.size());
    // handle core nodes
    for (NodeID node = 1; node < closest.size(); node++)
    {
        if (closest[node].node == node)
        {
            assert(closest[node].distance == 0);
            labels[node].cut_index = FlatCutIndex(ci[node]);
        }
        // conserve memory
        clear_and_shrink(ci[node].dist_index);
	clear_and_shrink(ci[node].paths);
        clear_and_shrink(ci[node].distances);
    }
    // handle periferal nodes
    for (NodeID node = 1; node < closest.size(); node++)
    {
        Neighbor n = closest[node];
        // isolated nodes got removed (n.node == NO_NODE)
        if (n.node != node && n.node != NO_NODE)
        {
            assert(n.distance > 0);
            // find root & distance
            NodeID root = n.node;
            distance_t root_dist = n.distance;
            while (closest[root].node != root)
            {
                root_dist += closest[root].distance;
                root = closest[root].node;
            }
            // copy index
            assert(!labels[root].cut_index.empty());
            labels[node].cut_index = labels[root].cut_index;
            labels[node].distance_offset = root_dist;
            labels[node].parent = n.node;
        }
    }
    clear_and_shrink(ci);
    clear_and_shrink(closest);
}

ContractionIndex::ContractionIndex(std::vector<CutIndex> &ci)
{
    labels.resize(ci.size());
    for (NodeID node = 1; node < ci.size(); node++)
        if (!ci[node].empty())
        {
            labels[node].cut_index = FlatCutIndex(ci[node]);
            // conserve memory
            clear_and_shrink(ci[node].dist_index);
	    clear_and_shrink(ci[node].paths);
            clear_and_shrink(ci[node].distances);
        }
    clear_and_shrink(ci);
}

ContractionIndex::~ContractionIndex()
{
    for (NodeID node = 1; node < labels.size(); node++)
        // not all labels own their cut index data
        if (!labels[node].cut_index.empty() && labels[node].distance_offset == 0)
            free(labels[node].cut_index.data);
}

distance_t ContractionIndex::get_distance(NodeID v, NodeID w) const
{
    ContractionLabel cv = labels[v], cw = labels[w];
    assert(!cv.cut_index.empty() && !cw.cut_index.empty());
    if (cv.cut_index == cw.cut_index)
    {
        if (v == w)
            return 0;
        if (cv.distance_offset == 0)
            return cw.distance_offset;
        if (cw.distance_offset == 0)
            return cv.distance_offset;
        if (cv.parent == w)
            return cv.distance_offset - cw.distance_offset;
        if (cw.parent == v)
            return cw.distance_offset - cv.distance_offset;
        // find lowest common ancestor
        NodeID v_anc = v, w_anc = w;
        ContractionLabel cv_anc = cv, cw_anc = cw;
        while (v_anc != w_anc)
        {
            if (cv_anc.distance_offset < cw_anc.distance_offset)
            {
                w_anc = cw_anc.parent;
                cw_anc = labels[w_anc];
            }
            else if (cv_anc.distance_offset > cw_anc.distance_offset)
            {
                v_anc = cv_anc.parent;
                cv_anc = labels[v_anc];
            }
            else
            {
                v_anc = cv_anc.parent;
                w_anc = cw_anc.parent;
                cv_anc = labels[v_anc];
                cw_anc = labels[w_anc];
            }
        }
        return cv.distance_offset + cw.distance_offset - 2 * cv_anc.distance_offset;
    }
    return cv.distance_offset + cw.distance_offset + get_distance(cv.cut_index, cw.cut_index);
}

uint16_t ContractionIndex::get_spc(NodeID v, NodeID w) const
{
    ContractionLabel cv = labels[v], cw = labels[w];
    assert(!cv.cut_index.empty() && !cw.cut_index.empty());
    if (cv.cut_index == cw.cut_index)
        return 1;
    return get_paths(cv.cut_index, cw.cut_index);
}

size_t ContractionIndex::get_hoplinks(NodeID v, NodeID w) const
{
    FlatCutIndex cv = labels[v].cut_index, cw = labels[w].cut_index;
    if (cv == cw)
        return 0;
    return get_hoplinks(cv, cw);
}

double ContractionIndex::avg_hoplinks(const std::vector<std::pair<NodeID,NodeID>> &queries) const
{
    size_t hop_count = 0;
    for (pair<NodeID,NodeID> q : queries)
        hop_count += get_hoplinks(q.first, q.second);
    return hop_count / (double)queries.size();
}

distance_t ContractionIndex::get_cut_level_distance(FlatCutIndex a, FlatCutIndex b, size_t cut_level)
{
    distance_t min_dist = infinity;
    uint16_t a_offset = get_offset(a.dist_index(), cut_level);
    uint16_t b_offset = get_offset(b.dist_index(), cut_level);
    const distance_t* a_ptr = a.distances() + a_offset;
    const distance_t* b_ptr = b.distances() + b_offset;
    const distance_t* a_end = a_ptr + min(a.dist_index()[cut_level] - a_offset, b.dist_index()[cut_level] - b_offset);
    // find min 2-hop distance within partition
    while (a_ptr != a_end)
    {
        distance_t dist = *a_ptr + *b_ptr;
        if (dist < min_dist)
            min_dist = dist;
        a_ptr++;
        b_ptr++;
    }
    return min_dist;
}

size_t ContractionIndex::get_cut_level_hoplinks(FlatCutIndex a, FlatCutIndex b, size_t cut_level)
{
    return min(a.cut_size(cut_level), b.cut_size(cut_level));
}

distance_t ContractionIndex::get_distance(FlatCutIndex a, FlatCutIndex b)
{
    // find lowest level at which partitions differ
    size_t cut_level = PBV::lca_level(*a.partition_bitvector(), *b.partition_bitvector());
#ifdef NO_SHORTCUTS
    distance_t min_dist = infinity;
#ifdef PRUNING
    for (size_t cl = 0; cl <= cut_level; cl++)
        min_dist = min(min_dist, get_cut_level_distance(a, b, cl));
#else
    // no pruning means we have a continuous block to check
    const distance_t* a_ptr = a.distances();
    const distance_t* b_ptr = b.distances();
    const distance_t* a_end = a_ptr + min(a.dist_index()[cut_level], b.dist_index()[cut_level]);
    while (a_ptr != a_end)
    {
        distance_t dist = *a_ptr + *b_ptr;
        if (dist < min_dist)
            min_dist = dist;
        a_ptr++;
        b_ptr++;
    }
#endif
    return min_dist;
#else
    return get_cut_level_distance(a, b, cut_level);
#endif
}

uint16_t ContractionIndex::get_paths(FlatCutIndex a, FlatCutIndex b)
{
    // find lowest level at which partitions differ
    size_t cut_level = PBV::lca_level(*a.partition_bitvector(), *b.partition_bitvector());

    distance_t min_dist = infinity; uint16_t spc = 0;
    const distance_t* a_ptr = a.distances();
    const distance_t* b_ptr = b.distances();
    const uint16_t* x_ptr = a.paths();
    const uint16_t* y_ptr = b.paths();
    const distance_t* a_end = a_ptr + min(a.dist_index()[cut_level], b.dist_index()[cut_level]);
    while (a_ptr != a_end)
    {
        distance_t d = *a_ptr + *b_ptr;
	uint16_t c = *x_ptr * *y_ptr;
        if (d < min_dist) {
            min_dist = d;
	    spc = c;
	} else if(d == min_dist) {
	    spc = spc + c;
	}
        a_ptr++; b_ptr++;
	x_ptr++; y_ptr++;
    }

    return spc;
}


bool ContractionIndex::is_contracted(NodeID node) const
{
    return labels[node].parent != NO_NODE;
}

size_t ContractionIndex::uncontracted_count() const
{
    size_t total = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!is_contracted(node))
            total++;
    return total;
}

bool ContractionIndex::in_partition_subgraph(NodeID node, uint64_t partition_bitvector) const
{
    return !is_contracted(node) && PBV::is_ancestor(partition_bitvector, *labels[node].cut_index.partition_bitvector());
}

uint16_t ContractionIndex::dist_index(NodeID node) const
{
    FlatCutIndex const& ci = labels[node].cut_index;
    uint16_t index = get_offset(ci.dist_index(), ci.cut_level());
    while (ci.distances()[index] != 0)
        index++;
    return index;
}

ContractionLabel ContractionIndex::get_contraction_label(NodeID v) const
{
    return labels[v];
}

void ContractionIndex::update_distance_offset(NodeID n, distance_t d)
{
    labels[n].distance_offset = d;
}

size_t ContractionIndex::get_hoplinks(FlatCutIndex a, FlatCutIndex b)
{
    // find lowest level at which partitions differ
    size_t cut_level = min(a.cut_level(), b.cut_level());
    uint64_t pa = a.partition(), pb = b.partition();
    if (pa != pb)
    {
        size_t diff_level = __builtin_ctzll(pa ^ pb); // count trailing zeros
        if (diff_level < cut_level)
            cut_level = diff_level;
    }
#ifdef NO_SHORTCUTS
    size_t hoplinks = 0;
    for (size_t cl = 0; cl <= cut_level; cl++)
        hoplinks += get_cut_level_hoplinks(a, b, cl);
    return hoplinks;
#else
    return get_cut_level_hoplinks(a, b, cut_level);
#endif
}

size_t ContractionIndex::size() const
{
    size_t total = 0;
    for (NodeID node = 1; node < labels.size(); node++)
    {
        // skip isolated nodes (subgraph)
        if (!labels[node].cut_index.empty())
            total += labels[node].size();
    }
    return total;
}

double ContractionIndex::avg_cut_size() const
{
    double cut_sum = 0, label_count = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!labels[node].cut_index.empty())
        {
            cut_sum += labels[node].cut_index.cut_level() + 1;
            label_count += labels[node].cut_index.label_count();
            // adjust for label pruning
            label_count += labels[node].cut_index.bottom_cut_size() + 1;
        }
    return label_count / max(1.0, cut_sum);
}

size_t ContractionIndex::max_cut_size() const
{
    size_t max_cut = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!labels[node].cut_index.empty())
            max_cut = max(max_cut, 1 + labels[node].cut_index.bottom_cut_size());
    return max_cut;
}

size_t ContractionIndex::height() const
{
    uint16_t max_cut_level = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!labels[node].cut_index.empty())
            max_cut_level = max(max_cut_level, labels[node].cut_index.cut_level());
    return max_cut_level;
}

size_t ContractionIndex::max_label_count() const
{
    size_t max_label_count = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!labels[node].cut_index.empty())
            max_label_count = max(max_label_count, labels[node].cut_index.label_count());
    return max_label_count;
}

size_t ContractionIndex::label_count() const
{
    size_t total = 0;
    for (NodeID node = 1; node < labels.size(); node++)
        if (!labels[node].cut_index.empty() && labels[node].distance_offset == 0)
            total += labels[node].cut_index.label_count();
    return total;
}

size_t ContractionIndex::non_empty_cuts() const
{
    size_t total = 0;
    for (NodeID node = 1; node < labels.size(); node++)
    {
        if (is_contracted(node))
            continue;
        // count nodes that come first within their cut
        FlatCutIndex const& ci = labels[node].cut_index;
        if (ci.distances()[get_offset(ci.dist_index(), ci.cut_level())] == 0)
            total++;
    }
    return total;
}

bool ContractionIndex::check_query(std::pair<NodeID,NodeID> query, Graph &g) const
{
    distance_t d_index = get_distance(query.first, query.second);
    uint16_t p_index = get_spc(query.first, query.second);
    distance_t d_dijkstra = g.get_distance(query.first, query.second, true);
    uint16_t p_dijkstra = g.get_path_count(query.first, query.second, true);
    if (d_index != d_dijkstra)
    {
        cerr << "BUG: d_index=" << d_index << ", d_dijkstra=" << d_dijkstra << endl;
        cerr << "index[" << query.first << "]=" << labels[query.first] << endl;
        cerr << "index[" << query.second << "]=" << labels[query.second] << endl;
    }

    if (p_index != p_dijkstra)
    {
        cerr << "BUG: p_index=" << p_index << ", p_dijkstra=" << p_dijkstra << endl;
        cerr << "index[" << query.first << "]=" << labels[query.first] << endl;
        cerr << "index[" << query.second << "]=" << labels[query.second] << endl;
    }
    return d_index == d_dijkstra;
}

pair<NodeID,NodeID> ContractionIndex::random_query() const
{
    assert(labels.size() > 1);
    NodeID node_count = labels.size() - 1;
    NodeID a = 1 + rand() % node_count;
    NodeID b = 1 + rand() % node_count;
    return make_pair(a, b);
}

void ContractionIndex::write(ostream& os) const
{
    size_t node_count = labels.size() - 1;
    os.write((char*)&node_count, sizeof(size_t));
    for (NodeID node = 1; node < labels.size(); node++)
    {
        ContractionLabel cl = labels[node];
        os.write((char*)&cl.distance_offset, sizeof(distance_t));
        if (cl.distance_offset == 0)
        {
            size_t data_size = cl.cut_index.size();
            os.write((char*)&data_size, sizeof(size_t));
            os.write(cl.cut_index.data, data_size);
        }
        else
            os.write((char*)&cl.parent, sizeof(NodeID));
    }
}

void ContractionIndex::write_json(std::ostream& os) const
{
    ListFormat lf = get_list_format();
    set_list_format(ListFormat::plain);
    // print json
    /*os << '{' << endl;
    for (NodeID node = 1; node < labels.size(); node++)
    {
        os << node << ":";
        ContractionLabel cl = labels[node];
        if (cl.distance_offset == 0) 
            os << cl.cut_index.unflatten();    
	else
            os << "{\"p\":" << cl.parent << ",\"d\":" << cl.distance_offset << "}";
        os << (node == labels.size() - 1 ? "" : ",") << endl;
    }
    os << '}' << endl;
    // reset formatting
    set_list_format(lf);*/

    os << '{' << endl;
    for (NodeID node = 1; node < labels.size(); node++)
    {
        os << node << ":";
        ContractionLabel cl = labels[node];
        if (cl.distance_offset == 0)
            os << cl.cut_index.unflatten_spc();
        else
            os << "{\"p\":" << cl.parent << ",\"d\":" << cl.distance_offset << "}";
        os << (node == labels.size() - 1 ? "" : ",") << endl;
    }
    os << '}' << endl;
    // reset formatting
    set_list_format(lf);
}

ContractionIndex::ContractionIndex(istream& is)
{
    // read index data
    size_t node_count = 0;
    is.read((char*)&node_count, sizeof(size_t));
    labels.resize(node_count + 1);
    for (NodeID node = 1; node < labels.size(); node++)
    {
        ContractionLabel &cl = labels[node];
        is.read((char*)&cl.distance_offset, sizeof(distance_t));
        if (cl.distance_offset == 0)
        {
            size_t data_size = 0;
            is.read((char*)&data_size, sizeof(size_t));
            cl.cut_index.data = (char*)malloc(data_size);
            is.read(cl.cut_index.data, data_size);
        }
        else
            is.read((char*)&cl.parent, sizeof(NodeID));
    }
    // fix data references
    for (NodeID node = 1; node < labels.size(); node++)
    {
        ContractionLabel &cl = labels[node];
        if (cl.distance_offset != 0)
        {
            NodeID root = cl.parent;
            while (labels[root].distance_offset != 0)
                root = labels[root].parent;
            cl.cut_index = labels[root].cut_index;
        }
    }
}

//--------------------------- Graph ---------------------------------

SubgraphID next_subgraph_id(bool reset)
{
    static atomic<SubgraphID> next_id = 1;
    if (reset)
        next_id = 1;
    return next_id++;
}

Neighbor::Neighbor(NodeID node, distance_t distance) : node(node), distance(distance)
{
}

Neighbor::Neighbor(NodeID node, distance_t distance, size_t path_count) : node(node), distance(distance), path_count(path_count)
{
}

bool Neighbor::operator<(const Neighbor &other) const
{
    return node < other.node;
}

Node::Node(SubgraphID subgraph_id) : subgraph_id(subgraph_id)
{
    distance = outcopy_distance = 0;
    inflow = outflow = NO_NODE;
    landmark_level = 0;
}

Node& MultiThreadNodeData::operator[](size_type pos)
{
    if (pos == Graph::s)
        return s_data;
    if (pos == Graph::t)
        return t_data;
    return vector::operator[](pos);
}

const Node& MultiThreadNodeData::operator[](size_type pos) const
{
    if (pos == Graph::s)
        return s_data;
    if (pos == Graph::t)
        return t_data;
    return vector::operator[](pos);
}

void MultiThreadNodeData::normalize()
{
    vector::operator[](Graph::s) = s_data;
    vector::operator[](Graph::t) = t_data;
}

double Partition::rating() const
{
    size_t l = left.size(), r = right.size(), c = cut.size();
    return min(l, r) / (c * c + 1.0);
}

Edge::Edge(NodeID a, NodeID b, distance_t d) : a(a), b(b), d(d)
{
}

bool Edge::operator<(Edge other) const
{
    return a < other.a
        || (a == other.a && b < other.b)
        || (a == other.a && b == other.b && d < other.d);
}

int32_t DiffData::diff() const
{
    return static_cast<int32_t>(dist_a) - static_cast<int32_t>(dist_b);
}

distance_t DiffData::min() const
{
    return std::min(dist_a, dist_b);
}

DiffData::DiffData(NodeID node, distance_t dist_a, distance_t dist_b) : node(node), dist_a(dist_a), dist_b(dist_b)
{
}

bool DiffData::cmp_diff(DiffData x, DiffData y)
{
    return x.diff() < y.diff();
}

// definition of static members
thread_local Node MultiThreadNodeData::s_data(NO_SUBGRAPH), MultiThreadNodeData::t_data(NO_SUBGRAPH);
#ifdef MULTI_THREAD
MultiThreadNodeData Graph::node_data;
size_t Graph::thread_threshold;
#else
vector<Node> Graph::node_data;
#endif
NodeID Graph::s, Graph::t;

void Graph::show_progress(bool state)
{
    log_progress_on = state;
}

bool Graph::contains(NodeID node) const
{
    return node_data[node].subgraph_id == subgraph_id;
}

Graph::Graph(size_t node_count)
{
    subgraph_id = next_subgraph_id(true);
    node_data.clear();
    resize(node_count);
    CHECK_CONSISTENT;
}

Graph::Graph(size_t node_count, const vector<Edge> &edges) : Graph(node_count)
{
    for (Edge e : edges)
        add_edge(e.a, e.b, e.d, true);
}

void Graph::resize(size_t node_count)
{
    assert(nodes.empty());
    // node numbering starts from 1, and we reserve two additional nodes for s & t
    node_data.clear();
    node_data.resize(node_count + 3, Node(subgraph_id));
    s = node_count + 1;
    t = node_count + 2;
    node_data[0].subgraph_id = node_data[s].subgraph_id = node_data[t].subgraph_id = NO_SUBGRAPH;
    nodes.reserve(node_count);
    for (NodeID node = 1; node <= node_count; node++)
        nodes.push_back(node);
#ifdef MULTI_THREAD
    thread_threshold = max(node_count / MULTI_THREAD, static_cast<size_t>(1000));
#endif
}

void Graph::add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse)
{
    assert(v < node_data.size());
    assert(w < node_data.size());
    assert(distance > 0);
    // check for existing edge
    bool exists = false;
    for (Neighbor &n : node_data[v].neighbors)
        if (n.node == w)
        {
            exists = true;
            n.distance = min(n.distance, distance);
            break;
        }
    if (!exists)
        node_data[v].neighbors.push_back(Neighbor(w, distance));
    if (add_reverse)
        add_edge(w, v, distance, false);
}

void Graph::remove_edge(NodeID v, NodeID w)
{
    std::erase_if(node_data[v].neighbors, [w](const Neighbor &n) { return n.node == w; });
    std::erase_if(node_data[w].neighbors, [v](const Neighbor &n) { return n.node == v; });
}

pair<distance_t, pair<NodeID, NodeID> > Graph::random_update()
{
        NodeID a = random_node();
        NodeID b = rand() % node_data[a].neighbors.size();

        return make_pair(node_data[a].neighbors[b].distance, make_pair(a, node_data[a].neighbors[b].node));
}

void Graph::update_edge(NodeID v, NodeID w, distance_t d)
{
        for (Neighbor &n : node_data[v].neighbors)
                if (n.node == w) {
                        n.distance = d;
                        break;
                }
}

void Graph::remove_isolated()
{
    unordered_set<NodeID> isolated;
    for (NodeID node : nodes)
        if (degree(node) == 0)
        {
            isolated.insert(node);
            node_data[node].subgraph_id = NO_SUBGRAPH;
        }
    std::erase_if(nodes, [&isolated](NodeID node) { return isolated.contains(node); });
}

void Graph::reset()
{
    nodes.clear();
    for (NodeID node = 1; node < node_data.size() - 2; node++)
    {
        if (!node_data[node].neighbors.empty())
        {
            nodes.push_back(node);
            node_data[node].subgraph_id = subgraph_id;
        }
    }
    node_data[s].subgraph_id = NO_SUBGRAPH;
    node_data[t].subgraph_id = NO_SUBGRAPH;
}

void Graph::add_node(NodeID v)
{
    assert(v < node_data.size());
    nodes.push_back(v);
    node_data[v].subgraph_id = subgraph_id;
}

void Graph::remove_nodes(const vector<NodeID> &node_set)
{
    util::remove_set(nodes, node_set);
    for (NodeID node : node_set)
        node_data[node].subgraph_id = NO_NODE;
}

size_t Graph::node_count() const
{
    return nodes.size();
}

size_t Graph::edge_count() const
{
    size_t ecount = 0;
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (contains(n.node))
                ecount++;
    return ecount / 2;
}

size_t Graph::degree(NodeID v) const
{
    assert(contains(v));
    size_t deg = 0;
    for (Neighbor n : node_data[v].neighbors)
        if (contains(n.node))
            deg++;
    return deg;
}

Neighbor Graph::single_neighbor(NodeID v) const
{
    assert(contains(v));
    Neighbor neighbor(NO_NODE, 0);
    for (Neighbor n : node_data[v].neighbors)
        if (contains(n.node))
        {
            if (neighbor.node == NO_NODE)
                neighbor = n;
            else
                return Neighbor(NO_NODE, 0);
        }
    return neighbor;
}

size_t Graph::super_node_count()
{
    return node_data.size() - 3;
}

const vector<NodeID>& Graph::get_nodes() const
{
    return nodes;
}

void Graph::get_edges(vector<Edge> &edges) const
{
    edges.clear();
    for (NodeID a : nodes)
        for (const Neighbor &n : node_data[a].neighbors)
            if (n.node > a && contains(n.node))
                edges.push_back(Edge(a, n.node, n.distance));
}

void Graph::assign_nodes()
{
    for (NodeID node : nodes)
        node_data[node].subgraph_id = subgraph_id;
}

//--------------------------- Synchronous queue ----------------------

// Pushes an element to the queue
template <typename T>
void TSQueue<T>::push(T item)
{
        std::lock_guard<std::mutex> lock(m_mutex);
        m_queue.push(item);
}

template <typename T>
bool TSQueue<T>::next(T& item)
{
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_queue.empty()) {
                return false;
        }

        item = m_queue.front();
        m_queue.pop();
        return true;
}

//--------------------------- Graph algorithms ----------------------

// helper struct to enque nodes by distance
struct SearchNode
{
    distance_t distance;
    NodeID node;
    // reversed for min-heap ordering
    bool operator<(const SearchNode &other) const { return distance > other.distance; }
    SearchNode(distance_t distance, NodeID node) : distance(distance), node(node) {}
};

void Graph::run_dijkstra(NodeID v)
{
    CHECK_CONSISTENT;
    assert(contains(v));
    // init distances
    for (NodeID node : nodes) {
        node_data[node].distance = infinity;
	node_data[node].path_count = 0;
    }
    node_data[v].distance = 0;
    node_data[v].path_count = 1;
    // init queue
    priority_queue<SearchNode> q;
    q.push(SearchNode(0, v));
    // dijkstra
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();

        for (Neighbor n : node_data[next.node].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph
            if (!contains(n.node))
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (new_dist < node_data[n.node].distance)
            {
                node_data[n.node].distance = new_dist;
		node_data[n.node].path_count = node_data[next.node].path_count;
                q.push(SearchNode(new_dist, n.node));
            } else if(new_dist == node_data[n.node].distance) {
	        node_data[n.node].path_count += node_data[next.node].path_count;
	    }
        }
    }
}

void Graph::run_dijkstra_llsub(NodeID v)
{
    CHECK_CONSISTENT;
    assert(contains(v));
    const uint16_t pruning_level = node_data[v].landmark_level;
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    priority_queue<SearchNode> q;
    q.push(SearchNode(0, v));
    // dijkstra
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();

        for (Neighbor n : node_data[next.node].neighbors)
        {
            Node &n_data = node_data[n.node];
            // filter neighbors nodes not belonging to subgraph or having higher landmark level
            if (!contains(n.node) || n_data.landmark_level >= pruning_level)
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (new_dist < n_data.distance)
            {
                n_data.distance = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
}

#ifdef PRUNING
void Graph::run_dijkstra_ll(NodeID v)
{
    CHECK_CONSISTENT;
    assert(contains(v));
    const uint16_t pruning_level = node_data[v].landmark_level;
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 1;
    // init queue
    priority_queue<SearchNode> q;
    for (Neighbor n : node_data[v].neighbors)
    {
        distance_t n_dist = (n.distance << 1) | 1;
        node_data[n.node].distance = n_dist;
        q.push(SearchNode(n_dist, n.node));
    }
    // dijkstra
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();

        const Node &next_data = node_data[next.node];
        distance_t current_dist = next_data.landmark_level >= pruning_level ? next.distance & ~static_cast<distance_t>(1) : next.distance;
        for (Neighbor n : next_data.neighbors)
        {
            // filter neighbors nodes not belonging to subgraph
            if (!contains(n.node))
                continue;
            // update distance and enque
            distance_t new_dist = current_dist + (n.distance << 1);
            if (new_dist < node_data[n.node].distance)
            {
                node_data[n.node].distance = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
}
#endif

#ifdef MULTI_THREAD_DISTANCES
void Graph::run_dijkstra_par(const vector<NodeID> &vertices)
{
    CHECK_CONSISTENT;
    vector<thread> threads;
    auto dijkstra = [this](NodeID v, size_t distance_id) {
        assert(contains(v));
        assert(distance_id < MULTI_THREAD_DISTANCES);
        // init distances
        for (NodeID node : nodes)
            node_data[node].distances[distance_id] = infinity;
        node_data[v].distances[distance_id] = 0;
        // init queue
        priority_queue<SearchNode> q;
        q.push(SearchNode(0, v));
        // dijkstra
        while (!q.empty())
        {
            SearchNode next = q.top();
            q.pop();

            for (Neighbor n : node_data[next.node].neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // update distance and enque
                distance_t new_dist = next.distance + n.distance;
                if (new_dist < node_data[n.node].distances[distance_id])
                {
                    node_data[n.node].distances[distance_id] = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                }
            }
        }
    };
    for (size_t i = 0; i < vertices.size(); i++)
        threads.push_back(thread(dijkstra, vertices[i], i));
    for (size_t i = 0; i < vertices.size(); i++)
        threads[i].join();
}

void Graph::run_dijkstra_llsub_par(const std::vector<NodeID> &vertices)
{
    CHECK_CONSISTENT;
    vector<thread> threads;
    auto dijkstra = [this](NodeID v, size_t distance_id) {
        assert(contains(v));
        assert(distance_id < MULTI_THREAD_DISTANCES);
        const uint16_t pruning_level = node_data[v].landmark_level;
        // init distances
        for (NodeID node : nodes)
            node_data[node].distances[distance_id] = infinity;
        node_data[v].distances[distance_id] = 0;
        // init queue
        priority_queue<SearchNode> q;
        q.push(SearchNode(0, v));
        // dijkstra
        while (!q.empty())
        {
            SearchNode next = q.top();
            q.pop();

            for (Neighbor n : node_data[next.node].neighbors)
            {
                Node &n_data = node_data[n.node];
                // filter neighbors nodes not belonging to subgraph or having higher landmark level
                if (!contains(n.node) || n_data.landmark_level >= pruning_level)
                    continue;
                // update distance and enque
                distance_t new_dist = next.distance + n.distance;
                if (new_dist < n_data.distances[distance_id])
                {
                    n_data.distances[distance_id] = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                }
            }
        }
    };
    for (size_t i = 0; i < vertices.size(); i++)
        threads.push_back(thread(dijkstra, vertices[i], i));
    for (size_t i = 0; i < vertices.size(); i++)
        threads[i].join();
}

#ifdef PRUNING
void Graph::run_dijkstra_ll_par(const vector<NodeID> &vertices)
{
    CHECK_CONSISTENT;
    vector<thread> threads;
    auto dijkstra = [this](NodeID v, size_t distance_id) {
        assert(contains(v));
        assert(distance_id < MULTI_THREAD_DISTANCES);
        const uint16_t pruning_level = node_data[v].landmark_level;
        // init distances
        for (NodeID node : nodes)
            node_data[node].distances[distance_id] = infinity;
        node_data[v].distances[distance_id] = 1;
        // init queue
        priority_queue<SearchNode> q;
        for (Neighbor n : node_data[v].neighbors)
        {
            distance_t n_dist = (n.distance << 1) | 1;
            node_data[n.node].distances[distance_id] = n_dist;
            q.push(SearchNode(n_dist, n.node));
        }
        // dijkstra
        while (!q.empty())
        {
            SearchNode next = q.top();
            q.pop();

            const Node &next_data = node_data[next.node];
            distance_t current_dist = next_data.landmark_level >= pruning_level ? next.distance & ~static_cast<distance_t>(1) : next.distance;
            for (Neighbor n : next_data.neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // update distance and enque
                distance_t new_dist = current_dist + (n.distance << 1);
                if (new_dist < node_data[n.node].distances[distance_id])
                {
                    node_data[n.node].distances[distance_id] = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                }
            }
        }
    };
    for (size_t i = 0; i < vertices.size(); i++)
        threads.push_back(thread(dijkstra, vertices[i], i));
    for (size_t i = 0; i < vertices.size(); i++)
        threads[i].join();
}
#endif
#endif

void Graph::run_bfs(NodeID v)
{
    CHECK_CONSISTENT;
    assert(contains(v));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    queue<NodeID> q;
    q.push(v);
    // BFS
    while (!q.empty())
    {
        NodeID next = q.front();
        q.pop();

        distance_t new_dist = node_data[next].distance + 1;
        for (Neighbor n : node_data[next].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph or already visited
            if (contains(n.node) && node_data[n.node].distance == infinity)
            {
                // update distance and enque
                node_data[n.node].distance = new_dist;
                q.push(n.node);
            }
        }
    }
}

// node in flow graph which splits nodes into incoming and outgoing copies
struct FlowNode
{
    NodeID node;
    bool outcopy; // outgoing copy of node?
    FlowNode(NodeID node, bool outcopy) : node(node), outcopy(outcopy) {}
};
ostream& operator<<(ostream &os, FlowNode fn)
{
    return os << "(" << fn.node << "," << (fn.outcopy ? "T" : "F") << ")";
}

// helper function
bool update_distance(distance_t &d, distance_t d_new)
{
    if (d > d_new)
    {
        d = d_new;
        return true;
    }
    return false;
}

void Graph::run_flow_bfs_from_s()
{
    CHECK_CONSISTENT;
    assert(contains(s) && contains(t));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = node_data[node].outcopy_distance = infinity;
    node_data[t].distance = node_data[t].outcopy_distance = 0;
    // init queue - start with neighbors of s as s requires special flow handling
    queue<FlowNode> q;
    for (Neighbor n : node_data[s].neighbors)
        if (contains(n.node) && node_data[n.node].inflow != s)
        {
            assert(node_data[n.node].inflow == NO_NODE);
            node_data[n.node].distance = 1;
            node_data[n.node].outcopy_distance = 1; // treat inner-node edges as length 0
            q.push(FlowNode(n.node, false));
        }
    // BFS
    while (!q.empty())
    {
        FlowNode fn = q.front();
        q.pop();

        distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
        NodeID inflow = node_data[fn.node].inflow;
        // special treatment is needed for node with flow through it
        if (inflow != NO_NODE && !fn.outcopy)
        {
            // inflow is only valid neighbor
            if (update_distance(node_data[inflow].outcopy_distance, fn_dist + 1))
            {
                // need to set distance for 0-distance nodes immediately
                // otherwise a longer path may set wrong distance value first
                update_distance(node_data[inflow].distance, fn_dist + 1);
                q.push(FlowNode(inflow, true));
            }
        }
        else
        {
            // when arriving at the outgoing copy of flow node, all neighbors except outflow are valid
            // outflow must have been already visited in this case, so checking all neighbors is fine
            for (Neighbor n : node_data[fn.node].neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // following inflow by inverting flow requires special handling
                if (n.node == inflow)
                {
                    if (update_distance(node_data[n.node].outcopy_distance, fn_dist + 1))
                    {
                        // neighbor must be a flow node
                        update_distance(node_data[n.node].distance, fn_dist + 1);
                        q.push(FlowNode(n.node, true));
                    }
                }
                else
                {
                    if (update_distance(node_data[n.node].distance, fn_dist + 1))
                    {
                        // neighbor may be a flow node
                        if (node_data[n.node].inflow == NO_NODE)
                            update_distance(node_data[n.node].outcopy_distance, fn_dist + 1);
                        q.push(FlowNode(n.node, false));
                    }
                }
            }
        }
    }
}

void Graph::run_flow_bfs_from_t()
{
    CHECK_CONSISTENT;
    assert(contains(s) && contains(t));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = node_data[node].outcopy_distance = infinity;
    node_data[t].distance = node_data[t].outcopy_distance = 0;
    // init queue - start with neighbors of t as t requires special flow handling
    queue<FlowNode> q;
    for (Neighbor n : node_data[t].neighbors)
        if (contains(n.node) && node_data[n.node].outflow != t)
        {
            assert(node_data[n.node].outflow == NO_NODE);
            node_data[n.node].outcopy_distance = 1;
            node_data[n.node].distance = 1; // treat inner-node edges as length 0
            q.push(FlowNode(n.node, true));
        }
    // BFS
    while (!q.empty())
    {
        FlowNode fn = q.front();
        q.pop();

        distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
        NodeID outflow = node_data[fn.node].outflow;
        // special treatment is needed for node with flow through it
        if (outflow != NO_NODE && fn.outcopy)
        {
            // outflow is only valid neighbor
            if (update_distance(node_data[outflow].distance, fn_dist + 1))
            {
                // need to set distance for 0-distance nodes immediately
                // otherwise a longer path may set wrong distance value first
                update_distance(node_data[outflow].outcopy_distance, fn_dist + 1);
                q.push(FlowNode(outflow, false));
            }
        }
        else
        {
            // when arriving at the incoming copy of flow node, all neighbors except inflow are valid
            // inflow must have been already visited in this case, so checking all neighbors is fine
            for (Neighbor n : node_data[fn.node].neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // following outflow by inverting flow requires special handling
                if (n.node == outflow)
                {
                    if (update_distance(node_data[n.node].distance, fn_dist + 1))
                    {
                        // neighbor must be a flow node
                        update_distance(node_data[n.node].outcopy_distance, fn_dist + 1);
                        q.push(FlowNode(n.node, false));
                    }
                }
                else
                {
                    if (update_distance(node_data[n.node].outcopy_distance, fn_dist + 1))
                    {
                        // neighbor may be a flow node
                        if (node_data[n.node].outflow == NO_NODE)
                            update_distance(node_data[n.node].distance, fn_dist + 1);
                        q.push(FlowNode(n.node, true));
                    }
                }
            }
        }
    }
}

distance_t Graph::get_distance(NodeID v, NodeID w, bool weighted)
{
    assert(contains(v) && contains(w));
    weighted ? run_dijkstra(v) : run_bfs(v);
    return node_data[w].distance;
}

distance_t Graph::get_path_count(NodeID v, NodeID w, bool weighted)
{
    assert(contains(v) && contains(w));
    weighted ? run_dijkstra(v) : run_bfs(v);
    return node_data[w].path_count;
}

pair<NodeID,distance_t> Graph::get_furthest(NodeID v, bool weighted)
{
    NodeID furthest = v;

    weighted ? run_dijkstra(v) : run_bfs(v);
    for (NodeID node : nodes)
        if (node_data[node].distance > node_data[furthest].distance)
            furthest = node;
    return make_pair(furthest, node_data[furthest].distance);
}

Edge Graph::get_furthest_pair(bool weighted)
{
    assert(nodes.size() > 1);
    distance_t max_dist = 0;
    NodeID start = nodes[0];
    pair<NodeID,distance_t> furthest = get_furthest(start, weighted);
    while (furthest.second > max_dist)
    {
        max_dist = furthest.second;
        start = furthest.first;
        furthest = get_furthest(start, weighted);
    }
    return Edge(start, furthest.first, max_dist);
}

distance_t Graph::diameter(bool weighted)
{
    if (nodes.size() < 2)
        return 0;
    return get_furthest_pair(weighted).d;
}

void Graph::get_diff_data(std::vector<DiffData> &diff, NodeID a, NodeID b, bool weighted, bool pre_computed)
{
    CHECK_CONSISTENT;
    assert(diff.empty());
    assert(!pre_computed || node_data[a].distance == 0);
    diff.reserve(nodes.size());
    // init with distances to a
    if (!pre_computed)
        weighted ? run_dijkstra(a) : run_bfs(a);
    for (NodeID node : nodes)
        diff.push_back(DiffData(node, node_data[node].distance, 0));
    // add distances to b
    weighted ? run_dijkstra(b) : run_bfs(b);
    for (DiffData &dd : diff)
        dd.dist_b = node_data[dd.node].distance;
}

// helper function for sorting connected components by size
static bool cmp_size_desc(const vector<NodeID> &a, const vector<NodeID> &b)
{
    return a.size() > b.size();
};

// helper function for adding nodes to smaller of two sets
static void add_to_smaller(vector<NodeID> &pa, vector<NodeID> &pb, const vector<NodeID> &cc)
{
    vector<NodeID> &smaller = pa.size() <= pb.size() ? pa : pb;
    smaller.insert(smaller.begin(), cc.cbegin(), cc.cend());
}

bool Graph::get_rough_partition(Partition &p, double balance, bool disconnected)
{
    DEBUG("get_rough_partition, p=" << p << ", disconnected=" << disconnected << " on " << *this);
    CHECK_CONSISTENT;
    assert(p.left.empty() && p.cut.empty() && p.right.empty());
    if (disconnected)
    {
        vector<vector<NodeID>> cc;
        get_connected_components(cc);
        if (cc.size() > 1)
        {
            DEBUG("found multiple connected components: " << cc);
            sort(cc.begin(), cc.end(), cmp_size_desc);
            // for size zero cuts we loosen the balance requirement
            if (cc[0].size() < nodes.size() * (1 - balance/2))
            {
                for (vector<NodeID> &c : cc)
                    add_to_smaller(p.left, p.right, c);
                return true;
            }
            // get rough partion over main component
            Graph main_cc(cc[0].begin(), cc[0].end());
            bool is_fine = main_cc.get_rough_partition(p, balance, false);
            // reset subgraph ids
            for (NodeID node : main_cc.nodes)
                node_data[node].subgraph_id = subgraph_id;
            if (is_fine)
            {
                // distribute remaining components
                for (size_t i = 1; i < cc.size(); i++)
                    add_to_smaller(p.left, p.right, cc[i]);
            }
            return is_fine;
        }
    }
    // graph is connected - find two extreme points
#ifdef NDEBUG
    NodeID a = get_furthest(random_node(), weighted_furthest).first;
#else
    NodeID a = get_furthest(nodes[0], weighted_furthest).first;
#endif
    NodeID b = get_furthest(a, weighted_furthest).first;
    DEBUG("furthest nodes: a=" << a << ", b=" << b);
    // get distances from a and b and sort by difference
    vector<DiffData> diff;
    get_diff_data(diff, a, b, weighted_diff, weighted_furthest);
    sort(diff.begin(), diff.end(), DiffData::cmp_diff);
    DEBUG("diff=" << diff);
    // get parition bounds based on balance; round up if possible
    size_t max_left = min(nodes.size() / 2, static_cast<size_t>(ceil(nodes.size() * balance)));
    size_t min_right = nodes.size() - max_left;
    DEBUG("max_left=" << max_left << ", min_right=" << min_right);
    assert(max_left <= min_right);
    // check for corner case where most nodes have same distance difference
    if (diff[max_left - 1].diff() == diff[min_right].diff())
    {
        // find bottleneck(s)
        const int32_t center_diff_value = diff[min_right].diff();
        distance_t min_dist = infinity;
        vector<NodeID> bottlenecks;
        for (DiffData dd : diff)
            if (dd.diff() == center_diff_value)
            {
                if (dd.min() < min_dist)
                {
                    min_dist = dd.min();
                    bottlenecks.clear();
                }
                if (dd.min() == min_dist)
                    bottlenecks.push_back(dd.node);
            }
        sort(bottlenecks.begin(), bottlenecks.end());
        DEBUG("bottlenecks=" << bottlenecks);
        // try again with bottlenecks removed
        remove_nodes(bottlenecks);
        bool is_fine = get_rough_partition(p, balance, true);
        // add bottlenecks back to graph and to center partition
        for (NodeID bn : bottlenecks)
        {
            add_node(bn);
            p.cut.push_back(bn);
        }
        // if bottlenecks are the only cut vertices, they must form a minimal cut
        return is_fine && p.cut.size() == bottlenecks.size();
    }
    // ensure left and right pre-partitions are connected
    while (diff[max_left - 1].diff() == diff[max_left].diff())
        max_left++;
    while (diff[min_right - 1].diff() == diff[min_right].diff())
        min_right--;
    // assign nodes to left/cut/right
    for (size_t i = 0; i < diff.size(); i++)
    {
        if (i < max_left)
            p.left.push_back(diff[i].node);
        else if (i < min_right)
            p.cut.push_back(diff[i].node);
        else
            p.right.push_back(diff[i].node);
    }
    return false;
}

void Graph::min_vertex_cuts(vector<vector<NodeID>> &cuts)
{
    DEBUG("min_vertex_cut over " << *this);
    CHECK_CONSISTENT;
    assert(contains(s) && contains(t));
    // set flow to empty
    for (NodeID node : nodes)
        node_data[node].inflow = node_data[node].outflow = NO_NODE;
#ifndef NDEBUG
    size_t last_s_distance = 1; // min s_distance is 2
#endif
    // find max s-t flow using Dinitz' algorithm
    while (true)
    {
        // construct BFS tree from t
        run_flow_bfs_from_t();
        DEBUG("BFS-tree: " << distances());
        const distance_t s_distance = node_data[s].outcopy_distance;
        if (s_distance == infinity)
            break;
        assert(s_distance > last_s_distance && (last_s_distance = s_distance));
        // run DFS from s along inverse BFS tree edges
        vector<NodeID> path;
        vector<FlowNode> stack;
        // iterating over neighbors of s directly simplifies stack cleanup after new s-t path is found
        for (Neighbor sn : node_data[s].neighbors)
        {
            if (!contains(sn.node) || node_data[sn.node].distance != s_distance - 1)
                continue;
            // ensure edge from s to neighbor exists in residual graph
            if (node_data[sn.node].inflow != NO_NODE)
            {
                assert(node_data[sn.node].inflow == s);
                continue;
            }
            stack.push_back(FlowNode(sn.node, false));
            while (!stack.empty())
            {
                FlowNode fn = stack.back();
                stack.pop_back();
                DEBUG("fn=" << fn);
                // clean up path (back tracking)
                distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
                // safeguard against re-visiting node during DFS (may have been enqueued before first visit)
                if (fn_dist == infinity)
                    continue;
                assert(fn_dist < s_distance && s_distance - fn_dist - 1 <= path.size());
                path.resize(s_distance - fn_dist - 1);
                // increase flow when s-t path is found
                if (fn.node == t)
                {
                    DEBUG("flow path=" << path);
                    assert(node_data[path.front()].inflow == NO_NODE);
                    node_data[path.front()].inflow = s;
                    for (size_t path_pos = 1; path_pos < path.size(); path_pos++)
                    {
                        NodeID from = path[path_pos - 1];
                        NodeID to = path[path_pos];
                        // we might be reverting existing flow
                        // from.inflow may have been changed already => check outflow
                        if (node_data[to].outflow == from)
                        {
                            node_data[to].outflow = NO_NODE;
                            if (node_data[from].inflow == to)
                                node_data[from].inflow = NO_NODE;
                        }
                        else
                        {
                            node_data[from].outflow = to;
                            node_data[to].inflow = from;
                        }
                    }
                    assert(node_data[path.back()].outflow == NO_NODE);
                    node_data[path.back()].outflow = t;
                    // skip to next neighbor of s
                    stack.clear();
                    path.clear();
                    DEBUG("new flow=" << flow());
                    break;
                }
                // ensure vertex is not re-visited during current DFS iteration
                if (fn.outcopy)
                    node_data[fn.node].outcopy_distance = infinity;
                else
                    node_data[fn.node].distance = infinity;
                // continue DFS from node
                path.push_back(fn.node);
                distance_t next_distance = fn_dist - 1;
                // when arriving at outgoing copy of a node with flow through it,
                // we are inverting outflow, so all neighbors are valid (except outflow)
                // otherwise inverting the inflow is the only possible option
                NodeID inflow = node_data[fn.node].inflow;
                if (inflow != NO_NODE && !fn.outcopy)
                {
                    if (node_data[inflow].outcopy_distance == next_distance)
                        stack.push_back(FlowNode(inflow, true));
                }
                else
                {
                    for (Neighbor n : node_data[fn.node].neighbors)
                    {
                        if (!contains(n.node))
                            continue;
                        // inflow inversion requires special handling
                        if (n.node == inflow)
                        {
                            if (node_data[inflow].outcopy_distance == next_distance)
                                stack.push_back(FlowNode(inflow, true));
                        }
                        else
                        {
                            if (node_data[n.node].distance == next_distance)
                                stack.push_back(FlowNode(n.node, false));
                        }
                    }
                }
            }
        }
    }
    // find min cut
    assert(cuts.empty());
    cuts.resize(1);
    // node-internal edge appears in cut iff outgoing copy is reachable from t in inverse residual graph and incoming copy is not
    // for node-external edges reachability of endpoint but unreachability of starting point is only possible if endpoint is t
    // in that case, starting point must become the cut vertex
    for (NodeID node : nodes)
    {
        NodeID outflow = node_data[node].outflow;
        // distance already stores distance from t in inverse residual graph
        if (outflow != NO_NODE)
        {
            assert(node_data[node].inflow != NO_NODE);
            if (node_data[node].outcopy_distance < infinity)
            {
                // check inner edge
                if (node_data[node].distance == infinity)
                    cuts[0].push_back(node);
            }
            else
            {
                // check outer edge
                if (outflow == t)
                    cuts[0].push_back(node);
            }
        }
    }
#ifdef MULTI_CUT
    // same thing but w.r.t. reachability from s in residual graph
    run_flow_bfs_from_s();
    cuts.resize(2);
    // distance now stores distance from s in residual graph
    for (NodeID node : nodes)
    {
        NodeID inflow = node_data[node].inflow;
        if (inflow != NO_NODE)
        {
            assert(node_data[node].outflow != NO_NODE);
            if (node_data[node].distance < infinity)
            {
                // check inner edge
                if (node_data[node].outcopy_distance == infinity)
                    cuts[1].push_back(node);
            }
            else
            {
                // check outer edge
                if (inflow == s)
                    cuts[1].push_back(node);
            }
        }
    }
    // eliminate potential duplicate
    if (cuts[0] == cuts[1])
        cuts.resize(1);
#endif
    DEBUG("cuts=" << cuts);
}

void Graph::get_connected_components(vector<vector<NodeID>> &components)
{
    CHECK_CONSISTENT;
    components.clear();
    for (NodeID start_node : nodes)
    {
        // visited nodes are temporarily removed
        if (!contains(start_node))
            continue;
        node_data[start_node].subgraph_id = NO_SUBGRAPH;
        // create new connected component
        components.push_back(vector<NodeID>());
        vector<NodeID> &cc = components.back();
        vector<NodeID> stack;
        stack.push_back(start_node);
        while (!stack.empty())
        {
            NodeID node = stack.back();
            stack.pop_back();
            cc.push_back(node);
            for (Neighbor n : node_data[node].neighbors)
                if (contains(n.node))
                {
                    node_data[n.node].subgraph_id = NO_SUBGRAPH;
                    stack.push_back(n.node);
                }
        }
    }
    // reset subgraph IDs
    assign_nodes();
    DEBUG("components=" << components);
    assert(util::size_sum(components) == nodes.size());
}

void Graph::rough_partition_to_cuts(vector<vector<NodeID>> &cuts, const Partition &p)
{
    // build subgraphs for rough partitions
    Graph left(p.left.cbegin(), p.left.cend());
    Graph center(p.cut.cbegin(), p.cut.cend());
    Graph right(p.right.cbegin(), p.right.cend());
    // construct s-t flow graph
    center.add_node(s);
    center.add_node(t);
    // handle corner case of edges between left and right partition
    // do this first as it can eliminate other s/t neighbors
    vector<NodeID> s_neighbors, t_neighbors;
    for (NodeID node : left.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (right.contains(n.node))
            {
                s_neighbors.push_back(node);
                t_neighbors.push_back(n.node);
            }
    util::make_set(s_neighbors);
    util::make_set(t_neighbors);
    // update pre-partition
    DEBUG("moving " << s_neighbors << " and " << t_neighbors << " to center");
    left.remove_nodes(s_neighbors);
    for (NodeID node : s_neighbors)
        center.add_node(node);
    right.remove_nodes(t_neighbors);
    for (NodeID node : t_neighbors)
        center.add_node(node);
    DEBUG("pre-partition=" << left.nodes << "|" << center.nodes << "|" << right.nodes);
    // identify additional neighbors of s and t
    for (NodeID node : left.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (center.contains(n.node))
                s_neighbors.push_back(n.node);
    for (NodeID node : right.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (center.contains(n.node))
                t_neighbors.push_back(n.node);
    util::make_set(s_neighbors);
    util::make_set(t_neighbors);
    // add edges incident to s and t
    for (NodeID node : s_neighbors)
        center.add_edge(s, node, 1, true);
    for (NodeID node : t_neighbors)
        center.add_edge(t, node, 1, true);
    // find minimum cut
    center.min_vertex_cuts(cuts);
    // revert s-t addition
    for (NodeID node : t_neighbors)
    {
        assert(node_data[node].neighbors.back().node == t);
        node_data[node].neighbors.pop_back();
    }
    node_data[t].neighbors.clear();
    for (NodeID node : s_neighbors)
    {
        assert(node_data[node].neighbors.back().node == s);
        node_data[node].neighbors.pop_back();
    }
    node_data[s].neighbors.clear();
    // repair subgraph IDs
    assign_nodes();
}

void Graph::complete_partition(Partition &p)
{
    CHECK_CONSISTENT;
    util::make_set(p.cut);
    remove_nodes(p.cut);
    // create left/right partitions
    p.left.clear(); p.right.clear();
    vector<vector<NodeID>> components;
    get_connected_components(components);
    sort(components.begin(), components.end(), cmp_size_desc);
    for (const vector<NodeID> &cc : components)
        add_to_smaller(p.left, p.right, cc);
    // add cut vertices back to graph
    for (NodeID node : p.cut)
        add_node(node);
    assert(p.left.size() + p.right.size() + p.cut.size() == nodes.size());
}

void Graph::create_partition(Partition &p, double balance)
{
    CHECK_CONSISTENT;
    assert(nodes.size() > 1);
    DEBUG("create_partition, p=" << p << " on " << *this);
    // find initial rough partition
#ifdef NO_SHORTCUTS
    bool is_fine = get_rough_partition(p, balance, true);
#else
    bool is_fine = get_rough_partition(p, balance, false);
#endif
    if (is_fine)
    {
        DEBUG("get_rough_partition found partition=" << p);
        return;
    }
    // find minimum cut
    vector<vector<NodeID>> cuts;
    rough_partition_to_cuts(cuts, p);
    assert(cuts.size() > 0);
    // create partition
    p.cut = cuts[0];
    complete_partition(p);
    for (size_t i = 1; i < cuts.size(); i++)
    {
        Partition p_alt;
        p_alt.cut = cuts[i];
        complete_partition(p_alt);
        if (p.rating() < p_alt.rating())
            p = p_alt;
    }
    DEBUG("partition=" << p);
}

void Graph::add_shortcuts(const vector<NodeID> &cut, const vector<CutIndex> &ci)
{
    CHECK_CONSISTENT;
    DEBUG("adding shortscuts on g=" << *this << ", cut=" << cut);
    // compute border nodes
    vector<NodeID> border;
    for (NodeID cut_node : cut)
        for (Neighbor n : node_data[cut_node].neighbors)
            if (contains(n.node))
                border.push_back(n.node);
    util::make_set(border);
    assert(!border.empty());
    // for distance in parent graph we use distances to cut nodes, which must already be in index
    size_t cut_level = ci[cut[0]].cut_level;
    // compute distances between border nodes within subgraph and parent graph
    vector<distance_t> d_partition, d_graph;
#ifdef MULTI_THREAD_DISTANCES
    if (nodes.size() > thread_threshold)
    {
        size_t next_offset;
        for (size_t offset = 0; offset < border.size(); offset = next_offset)
        {
            next_offset = min(offset + MULTI_THREAD_DISTANCES, border.size());
            const vector<NodeID> partial_cut(border.begin() + offset, border.begin() + next_offset);
            run_dijkstra_par(partial_cut);
            for (size_t distance_id = 0; distance_id < partial_cut.size(); distance_id++)
            {
                NodeID n_i = border[distance_id + offset];
                for (size_t j = 0; j < distance_id + offset; j++)
                {
                    NodeID n_j = border[j];
                    distance_t d_ij = node_data[n_j].distances[distance_id];
                    d_partition.push_back(d_ij);
                    distance_t d_cut = get_cut_level_distance(ci[n_i], ci[n_j], cut_level);
                    d_graph.push_back(min(d_ij, d_cut));
                }
            }
        }
    }
    else
#endif
    for (size_t i = 1; i < border.size(); i++)
    {
        NodeID n_i = border[i];
        run_dijkstra(n_i);
        for (size_t j = 0; j < i; j++)
        {
            assert(d_partition.size() == hmi(i, j));
            NodeID n_j = border[j];
            distance_t d_ij = node_data[n_j].distance;
            d_partition.push_back(d_ij);
            distance_t d_cut = get_cut_level_distance(ci[n_i], ci[n_j], cut_level);
            d_graph.push_back(min(d_ij, d_cut));
        }
    }
    // find & add non-redundant shortcuts
    // separate loop as d_graph must be fully computed for redundancy check
    size_t idx_ij = 0;
    for (size_t i = 1; i < border.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            assert(idx_ij == hmi(i, j));
            distance_t dg_ij = d_graph[idx_ij];
            if (d_partition[idx_ij] > dg_ij)
            {
                bool redundant = false;
                // check for redundancy due to shortest path through third border node k
                for (size_t k = 0; k < border.size(); k++)
                {
                    if (k == i || k == j)
                        continue;
                    if (d_graph[hmi(i, k)] + d_graph[hmi(k, j)] == dg_ij)
                    {
                        redundant = true;
                        break;
                    }
                }
                if (!redundant)
                {
                    DEBUG("shortcut: " << border[i] << "-[" << dg_ij << "]-" << border[j]);
                    add_edge(border[i], border[j], dg_ij, true);
                }
            }
            idx_ij++;
        }
    }
}

void Graph::sort_cut_for_pruning(vector<NodeID> &cut, [[maybe_unused]] vector<CutIndex> &ci)
{
    // compute pruning potential for each cut node
    vector<pair<size_t,NodeID>> pruning_potential;
    for (NodeID node : cut)
        pruning_potential.push_back(make_pair(0, node));
#ifdef PRUNING
    // mimics code in extend_on_partition
    for (size_t c = 0; c < cut.size(); c++)
        node_data[cut[c]].landmark_level = 1;
    #ifdef MULTI_THREAD_DISTANCES
    if (nodes.size() > thread_threshold)
    {
        size_t next_offset;
        for (size_t offset = 0; offset < cut.size(); offset = next_offset)
        {
            next_offset = min(offset + MULTI_THREAD_DISTANCES, cut.size());
            const vector<NodeID> partial_cut(cut.begin() + offset, cut.begin() + next_offset);
            run_dijkstra_ll_par(partial_cut);
            for (size_t distance_id = 0; distance_id < partial_cut.size(); distance_id++)
                for (NodeID node : nodes)
                {
                    distance_t dist_and_flag = node_data[node].distances[distance_id];
                    if ((dist_and_flag & 1) == 0)
                    {
                        pruning_potential[offset + distance_id].first++;
                        ci[node].pruning_3hop++;
                    }
                }
        }
    }
    else
    #endif
    for (size_t c = 0; c < cut.size(); c++)
    {
        run_dijkstra_ll(cut[c]);
        for (NodeID node : nodes)
        {
            distance_t dist_and_flag = node_data[node].distance;
            if ((dist_and_flag & 1) == 0)
            {
                pruning_potential[c].first++;
                ci[node].pruning_3hop++;
            }
        }
    }
#endif
    // sort cut
    sort(pruning_potential.begin(), pruning_potential.end());
    for (size_t c = 0; c < cut.size(); c++)
        cut[c] = pruning_potential[c].second;
}

void Graph::extend_on_partition(vector<CutIndex> &ci, double balance, uint8_t cut_level, const vector<NodeID> &p, [[maybe_unused]] const vector<NodeID> &cut)
{
    DEBUG("extend_on_partition, p=" << p << ", cut=" << cut);
    if (p.size() > 1)
    {
        Graph g(p.begin(), p.end());
#ifndef NO_SHORTCUTS
        START_TIMER;
        g.add_shortcuts(cut, ci);
        STOP_TIMER(t_shortcut);
#endif
        g.extend_cut_index(ci, balance, cut_level + 1);
    }
    else if (p.size() == 1)
    {
        ci[p[0]].cut_level = cut_level + 1;
        ci[p[0]].dist_index.push_back(ci[p[0]].dist_index[cut_level] + 1);
        assert(ci[p[0]].is_consistent());
    }
}

void Graph::extend_cut_index(vector<CutIndex> &ci, double balance, uint8_t cut_level)
{
    //cout << (int)cut_level << flush;
    DEBUG("extend_cut_index at level " << (int)cut_level << " on " << *this);
    DEBUG("cut index=" << ci);
    CHECK_CONSISTENT;
    assert(cut_level <= MAX_CUT_LEVEL);
    if (node_count() < 2)
    {
        assert(cut_level == 0);
        for (NodeID node : nodes)
        {
            ci[node].cut_level = 0;
            ci[node].dist_index.push_back(0);
        }
        return;
    }
    // find balanced cut
    Partition p;
    if (cut_level < MAX_CUT_LEVEL)
    {
        START_TIMER;
        create_partition(p, balance);
#ifdef CUT_REPEAT
        for (size_t i = 1; i < CUT_REPEAT; i++)
        {
            Partition p_new;
            create_partition(p_new, balance);
            if (p_new.rating() > p.rating())
                p = p_new;
        }
#endif
        STOP_TIMER(t_partition);
    }
    else
        p.cut = nodes;

    for (size_t c = 0; c < p.cut.size(); c++)
        node_data[p.cut[c]].landmark_level = p.cut.size() - c;

    // update dist_index
    for (NodeID node : nodes)
    {
        assert(ci[node].dist_index.size() == cut_level);
        if(node_data[node].landmark_level == 0)
            ci[node].dist_index.push_back(ci[node].dist_index[cut_level - 1] + p.cut.size());
	else 
	    ci[node].dist_index.push_back(ci[node].dist_index[cut_level - 1] + (p.cut.size() - node_data[node].landmark_level + 1));
    }

    // set cut_level
    for (NodeID c : p.cut)
    {
        ci[c].cut_level = cut_level;
        assert(ci[c].is_consistent());
    }
    // update partition bitstring
    for (NodeID node : p.right)
        ci[node].partition |= (static_cast<uint64_t>(1) << cut_level);
    DEBUG("cut index extended to " << ci);

    // reset landmark flags
    for (NodeID c : p.cut)
        node_data[c].landmark_level = 0;
    STOP_TIMER(t_label);

    // add shortcuts and recurse
#ifdef MULTI_THREAD
    if (nodes.size() > thread_threshold)
    {
        std::thread t_left(extend_on_partition, std::ref(ci), balance, cut_level, std::cref(p.left), std::cref(p.cut));
        extend_on_partition(ci, balance, cut_level, p.right, p.cut);
        t_left.join();
    }
    else
#endif
    {
        extend_on_partition(ci, balance, cut_level, p.left, p.cut);
        extend_on_partition(ci, balance, cut_level, p.right, p.cut);
    }
}

size_t Graph::create_cut_index(std::vector<CutIndex> &ci, double balance)
{
#ifndef NPROFILE
    t_partition = t_label = t_shortcut = 0;
#endif
    assert(is_undirected());
#ifndef NDEBUG
    // sort neighbors to make algorithms deterministic
    for (NodeID node : nodes)
        sort(node_data[node].neighbors.begin(), node_data[node].neighbors.end());
#endif
    // store original neighbor counts
    vector<NodeID> original_nodes = nodes;
    vector<size_t> original_neighbors(node_data.size());
    for (NodeID node : nodes)
        original_neighbors[node] = node_data[node].neighbors.size();
    // create index
    ci.clear();
    ci.resize(node_data.size() - 2);
    // reduce memory fragmentation by pre-allocating sensible values
    for (NodeID node : nodes)
    {
        ci[node].dist_index.reserve(32);
    }
    extend_cut_index(ci, balance, 0);
    log_progress(0);
    // reset nodes (top-level cut vertices got removed)
    nodes = original_nodes;
#ifndef NDEBUG
    for (NodeID node : nodes)
        if (!ci[node].is_consistent())
            cerr << "inconsistent cut index for node " << node << ": "<< ci[node] << endl;
#endif
#ifndef NPROFILE
    cerr << "partitioning took " << t_partition << "s" << endl;
    cerr << "labeling took " << t_label << "s" << endl;
    cerr << "shortcuts took " << t_shortcut << "s" << endl;
#endif
    //return shortcuts / 2;
    return 0;
}

void Graph::get_redundant_edges(std::vector<Edge> &edges)
{
    CHECK_CONSISTENT;
    assert(edges.empty());
    // reset distances for all nodes
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    // run localized Dijkstra from each node
    vector<NodeID> visited;
    priority_queue<SearchNode> q;
    for (NodeID v : nodes)
    {
        node_data[v].distance = 0;
        visited.push_back(v);
        distance_t max_dist = 0;
        // init queue - starting from neighbors ensures that only paths of length 2+ are considered
        for (Neighbor n : node_data[v].neighbors)
            if (contains(n.node))
            {
                q.push(SearchNode(n.distance, n.node));
                if (v < n.node)
                    max_dist = max(max_dist, n.distance);
            }
        // dijkstra
        while (!q.empty())
        {
            SearchNode next = q.top();
            q.pop();

            for (Neighbor n : node_data[next.node].neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // update distance and enque
                distance_t new_dist = next.distance + n.distance;
                if (new_dist <= max_dist && new_dist < node_data[n.node].distance)
                {
                    node_data[n.node].distance = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                    visited.push_back(n.node);
                }
            }
        }
        // identify redundant edges
        for (Neighbor n : node_data[v].neighbors)
            // only add redundant edges once
            if (v < n.node && contains(n.node) && node_data[n.node].distance <= n.distance)
                edges.push_back(Edge(v, n.node, n.distance));
        // cleanup
        for (NodeID w : visited)
            node_data[w].distance = infinity;
        visited.clear();
    }
}

void Graph::contract(vector<Neighbor> &closest)
{
    closest.resize(node_data.size() - 2, Neighbor(NO_NODE, 0));
    for (NodeID node : nodes)
        closest[node] = Neighbor(node, 0);
    // helper function to identify degree one nodes and associated neighbors
    auto find_degree_one = [this, &closest](const vector<NodeID> &nodes, vector<NodeID> &degree_one, vector<NodeID> &neighbors) {
        degree_one.clear();
        neighbors.clear();
        for (NodeID node : nodes)
        {
            Neighbor neighbor = single_neighbor(node);
            if (neighbor.node != NO_NODE)
            {
                // avoid complete contraction (screws with testing)
                if (single_neighbor(neighbor.node).node == NO_NODE)
                {
                    closest[node] = neighbor;
                    degree_one.push_back(node);
                    neighbors.push_back(neighbor.node);
                }
            }
        }
    };
    // remove nodes
    vector<NodeID> degree_one, neighbors;
    find_degree_one(nodes, degree_one, neighbors);
    while (!degree_one.empty())
    {
        sort(degree_one.begin(), degree_one.end());
        remove_nodes(degree_one);
        vector<NodeID> old_neighbors = neighbors;
        find_degree_one(old_neighbors, degree_one, neighbors);
    }
}

//--------------------------- ContractionHierarchy ------------------

ContractionHierarchy::ContractionHierarchy() {

}

ContractionHierarchy::ContractionHierarchy(istream &is) {

    size_t count, node_count, total_size = 0;
    is.read((char*)&node_count, sizeof(size_t));
    nodes.resize(node_count);
    for(NodeID i = 1; i < node_count; i++) {
        is.read((char*)&nodes[i].dist_index, sizeof(uint16_t));
        if(nodes[i].dist_index == 65535)
            continue;
        is.read((char*)&count, sizeof(size_t));
        nodes[i].up_neighbors.reserve(count); total_size += count;
        for(size_t j = 0; j < count; j++) {
            Neighbor n(NO_NODE, 0, 0);
            is.read((char*)&n.node, sizeof(NodeID));
            is.read((char*)&n.distance, sizeof(distance_t));
            is.read((char*)&n.path_count, sizeof(uint16_t));
            nodes[i].up_neighbors.push_back(n);
        }
        is.read((char*)&count, sizeof(size_t));
        nodes[i].down_neighbors.reserve(count);
        for(size_t j = 0; j < count; j++) {
            NodeID n;
            is.read((char*)&n, sizeof(NodeID));
            nodes[i].down_neighbors.push_back(n);
        }
    }
    cout << total_size << endl;
}

void ContractionHierarchy::write(ostream &os) {

    size_t count = nodes.size();
    os.write((char*)&count, sizeof(size_t));
    for(NodeID i = 1; i < nodes.size() ; i++) {
        os.write((char*)&nodes[i].dist_index, sizeof(uint16_t));
        if(nodes[i].dist_index == 65535)
            continue;
        count = nodes[i].up_neighbors.size();
        os.write((char*)&count, sizeof(size_t));
        for(Neighbor n: nodes[i].up_neighbors) {
            os.write((char*)&n.node, sizeof(NodeID));
            os.write((char*)&n.distance, sizeof(distance_t));
            os.write((char*)&n.path_count, sizeof(uint16_t));
        }
        count = nodes[i].down_neighbors.size();
        os.write((char*)&count, sizeof(size_t));
        for(NodeID n: nodes[i].down_neighbors) {
            os.write((char*)&n, sizeof(NodeID));
        }
    }
}

size_t ContractionHierarchy::size() const
{
    size_t total = 0;
    for(NodeID i = 1; i < nodes.size() ; i++) {
	if(nodes[i].dist_index == 65535)
            continue;
        total += sizeof(uint64_t);
	total += nodes[i].up_neighbors.size() * sizeof(NodeID);
	total += nodes[i].up_neighbors.size() * sizeof(distance_t);
	total += nodes[i].up_neighbors.size() * sizeof(uint16_t);
	total += nodes[i].down_neighbors.size() * sizeof(NodeID);
    }       
    return total;
}

size_t ContractionHierarchy::edge_count() const
{
    size_t total = 0;
    for (CHNode const& node : nodes)
        total += node.up_neighbors.size();
    return total;
}

void Graph::create_sc_graph(ContractionHierarchy &ch, vector<CutIndex> &ci)
{
    vector<thread> threads;
    auto dclp = [this](ContractionHierarchy &ch, vector<CutIndex> &ci, util::par_max_bucket_list<NodeID, MULTI_THREAD_DISTANCES> &que, size_t id) {

        NodeID x;
        while (que.next(x, id))
        {
            for(Neighbor &n: ch.nodes[x].up_neighbors) {
                for(size_t anc = 0; anc < ch.nodes[n.node].dist_index; anc++) {
                    distance_t dist = n.distance + ci[n.node].distances[anc];
                    uint16_t path_count = n.path_count * ci[n.node].paths[anc];
                    if(dist < ci[x].distances[anc]) {
                        ci[x].distances[anc] = dist;
                        ci[x].paths[anc] = path_count;
                    } else if (dist == ci[x].distances[anc]) {
                        ci[x].paths[anc] += path_count;
                    }
                }
            }
            ci[x].distances.push_back(0);
            ci[x].paths.push_back(1);
        }
    };

    vector<NodeID> bottom_up_nodes;
    bottom_up_nodes.reserve(nodes.size() + 1);
    // initialize distance index to determine edge direction
    ch.nodes.resize(nodes.size() + 1);
    for (NodeID node : nodes) {
        ch.nodes[node].dist_index = ci[node].dist_index[ci[node].cut_level] - 1;
	ci[node].distances.resize(ch.nodes[node].dist_index, infinity);
	ci[node].paths.resize(ch.nodes[node].dist_index, 0);
    }

    // initialize with upwards graph edges
    for (NodeID node : nodes)
    {
        bottom_up_nodes.push_back(node);
        for (Neighbor &n : node_data[node].neighbors)
	    if (ch.nodes[n.node].dist_index < ch.nodes[node].dist_index) {
		ch.nodes[node].up_neighbors.push_back(Neighbor(n.node, n.distance, 1));
                ci[node].distances[ch.nodes[n.node].dist_index] = n.distance;
                ci[node].paths[ch.nodes[n.node].dist_index] = 1;
	    }
    };

    // add shortcuts bottom-up
    auto di_order = [&ch](NodeID a, NodeID b) -> bool
    {
        return ch.nodes[a].dist_index > ch.nodes[b].dist_index;
    };
    auto di_order1 = [&ch](Neighbor a, Neighbor b) -> bool
    {
        if(ch.nodes[a.node].dist_index > ch.nodes[b.node].dist_index) return true;
	if(ch.nodes[a.node].dist_index == ch.nodes[b.node].dist_index && a.distance < b.distance) return true;
	if(ch.nodes[a.node].dist_index == ch.nodes[b.node].dist_index && a.distance == b.distance && a.path_count > b.path_count) return true;
	return false;
    };

    std::sort(bottom_up_nodes.begin(), bottom_up_nodes.end(), di_order);

    for (NodeID node : bottom_up_nodes)
    {
        vector<Neighbor> &up = ch.nodes[node].up_neighbors;
        util::make_set(up, di_order1);

	for (size_t i = 0; i + 1 < up.size(); i++) {
            for (size_t j = i + 1; j < up.size(); j++) {

                distance_t weight = up[i].distance + up[j].distance;
                uint16_t path_count = up[i].path_count * up[j].path_count;
                if(weight < ci[up[i].node].distances[ch.nodes[up[j].node].dist_index]) {
                    ch.nodes[up[i].node].up_neighbors.push_back(Neighbor(up[j].node, weight, path_count));
                    ci[up[i].node].distances[ch.nodes[up[j].node].dist_index] = weight;
                    ci[up[i].node].paths[ch.nodes[up[j].node].dist_index] = path_count;
                } else if(weight == ci[up[i].node].distances[ch.nodes[up[j].node].dist_index]) {
		    ci[up[i].node].paths[ch.nodes[up[j].node].dist_index] += path_count;
                    ch.nodes[up[i].node].up_neighbors.push_back(Neighbor(up[j].node, weight, ci[up[i].node].paths[ch.nodes[up[j].node].dist_index]));
                }
            }
        }

        // create downward neighbors from upward ones
        for (Neighbor upn: up)
            ch.nodes[upn.node].down_neighbors.push_back(node);
    }

#ifdef MULTI_THREAD_DISTANCES
    util::par_max_bucket_list<NodeID, MULTI_THREAD_DISTANCES> que(ch.nodes[bottom_up_nodes[0]].dist_index);
    for (NodeID node : bottom_up_nodes)
        que.push(node, ch.nodes[node].dist_index);

    // add nodes to queue (pre-compute)
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads.push_back(thread(dclp, std::ref(ch), std::ref(ci), std::ref(que), i));
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads[i].join();
#else
    // compute DHCL distances
    for(auto it = bottom_up_nodes.rbegin(); it != bottom_up_nodes.rend(); it++) {
	sort(ch.nodes[*it].down_neighbors.begin(), ch.nodes[*it].down_neighbors.end());
    	for(Neighbor n: ch.nodes[*it].up_neighbors) {
            for(size_t anc = 0; anc < ch.nodes[n.node].dist_index; anc++) {
                distance_t dist = n.distance + ci[n.node].distances[anc];
		uint16_t path_count = n.path_count * ci[n.node].paths[anc];
                if(dist < ci[*it].distances[anc]) {
                    ci[*it].distances[anc] = dist;
                    ci[*it].paths[anc] = path_count;
                } else if (dist == ci[*it].distances[anc]) {
                    ci[*it].paths[anc] += path_count;
                }
            }
        }
        ci[*it].distances.push_back(0);
        ci[*it].paths.push_back(1);
    }
#endif
}

void Graph::create_sc_graph(ContractionHierarchy &ch, vector<CutIndex> &ci, vector<Neighbor> &closest)
{
    vector<thread> threads;
    auto dclp = [this](ContractionHierarchy &ch, vector<CutIndex> &ci, util::par_max_bucket_list<NodeID, MULTI_THREAD_DISTANCES> &que, size_t id) {

        NodeID x;
        while (que.next(x, id))
        {
            for(Neighbor &n: ch.nodes[x].up_neighbors) {
                for(size_t anc = 0; anc < ch.nodes[n.node].dist_index; anc++) {
                    distance_t dist = n.distance + ci[n.node].distances[anc];
                    uint16_t path_count = n.path_count * ci[n.node].paths[anc];
                    if(dist < ci[x].distances[anc]) {
                        ci[x].distances[anc] = dist;
                        ci[x].paths[anc] = path_count;
                    } else if (dist == ci[x].distances[anc]) {
                        ci[x].paths[anc] += path_count;
                    }
                }
            }
            ci[x].distances.push_back(0);
            ci[x].paths.push_back(1);
        }
    };

    vector<NodeID> bottom_up_nodes;
    bottom_up_nodes.reserve(nodes.size() + 1); 
    // initialize distance index to determine edge direction
    ch.nodes.resize(nodes.size() + 1);
    for (NodeID node : nodes) {
        if(closest[node].node == node) {

            bottom_up_nodes.push_back(node);
            ch.nodes[node].dist_index = ci[node].dist_index[ci[node].cut_level] - 1;
            ci[node].distances.resize(ch.nodes[node].dist_index, infinity);
	    ci[node].paths.resize(ch.nodes[node].dist_index, 0);
        } else
            ch.nodes[node].dist_index = 65535;
    } 

    // initialize with upwards graph edges
    for (NodeID node : bottom_up_nodes)
    {
        for (Neighbor &n : node_data[node].neighbors)
            if (closest[n.node].node == n.node && ch.nodes[n.node].dist_index < ch.nodes[node].dist_index) {
                ch.nodes[node].up_neighbors.push_back(Neighbor(n.node, n.distance, 1));
                ci[node].distances[ch.nodes[n.node].dist_index] = n.distance;
		ci[node].paths[ch.nodes[n.node].dist_index] = 1;
            }
    }

    // add shortcuts bottom-up
    auto di_order = [&ch](NodeID a, NodeID b) -> bool
    {
        return ch.nodes[a].dist_index > ch.nodes[b].dist_index;
    };
    auto di_order1 = [&ch](Neighbor a, Neighbor b) -> bool
    {
        if(ch.nodes[a.node].dist_index > ch.nodes[b.node].dist_index) return true;
        if(ch.nodes[a.node].dist_index == ch.nodes[b.node].dist_index && a.distance < b.distance) return true;
	if(ch.nodes[a.node].dist_index == ch.nodes[b.node].dist_index && a.distance == b.distance && a.path_count > b.path_count) return true;
        return false;
    };

    std::sort(bottom_up_nodes.begin(), bottom_up_nodes.end(), di_order);
    for (NodeID node : bottom_up_nodes)
    {
        vector<Neighbor> &up = ch.nodes[node].up_neighbors;
        util::make_set(up, di_order1);

        for (size_t i = 0; i + 1 < up.size(); i++) {
            for (size_t j = i + 1; j < up.size(); j++) {
                distance_t weight = up[i].distance + up[j].distance;
		size_t path_count = up[i].path_count * up[j].path_count;
                if(weight < ci[up[i].node].distances[ch.nodes[up[j].node].dist_index]) {
                    ch.nodes[up[i].node].up_neighbors.push_back(Neighbor(up[j].node, weight, path_count));
                    ci[up[i].node].distances[ch.nodes[up[j].node].dist_index] = weight;
		    ci[up[i].node].paths[ch.nodes[up[j].node].dist_index] = path_count;
                } else if(weight == ci[up[i].node].distances[ch.nodes[up[j].node].dist_index]) {
		    ci[up[i].node].paths[ch.nodes[up[j].node].dist_index] += path_count;
		    ch.nodes[up[i].node].up_neighbors.push_back(Neighbor(up[j].node, weight, ci[up[i].node].paths[ch.nodes[up[j].node].dist_index]));
		}
            }
        }

        // create downward neighbors from upward ones
        for (Neighbor upn: up)
            ch.nodes[upn.node].down_neighbors.push_back(node);
    }

#ifdef MULTI_THREAD_DISTANCES
    util::par_max_bucket_list<NodeID, MULTI_THREAD_DISTANCES> que(ch.nodes[bottom_up_nodes[0]].dist_index);
    for (NodeID node : bottom_up_nodes)
        que.push(node, ch.nodes[node].dist_index);

    // add nodes to queue (pre-compute)
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads.push_back(thread(dclp, std::ref(ch), std::ref(ci), std::ref(que), i));
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads[i].join();
#else
    // compute DHCL distances
    for(auto it = bottom_up_nodes.rbegin(); it != bottom_up_nodes.rend(); it++) {
	sort(ch.nodes[*it].down_neighbors.begin(), ch.nodes[*it].down_neighbors.end());
        for(Neighbor n: ch.nodes[*it].up_neighbors) {
            for(size_t anc = 0; anc < ch.nodes[n.node].dist_index; anc++) {
	        distance_t dist = n.distance + ci[n.node].distances[anc];
		uint16_t path_count = n.path_count * ci[n.node].paths[anc];
	        if(dist < ci[*it].distances[anc]) {
		    ci[*it].distances[anc] = dist;
		    ci[*it].paths[anc] = path_count;
		} else if (dist == ci[*it].distances[anc]) {
		    ci[*it].paths[anc] += path_count;
		}
	    }
        }
        ci[*it].distances.push_back(0);
	ci[*it].paths.push_back(1);
    }
#endif
}

Neighbor& Graph::UpNeighbor(ContractionHierarchy &ch, NodeID v, NodeID w) {

    for(Neighbor &n: ch.nodes[v].up_neighbors) {
        if(n.node == w) {
            return n;
        }
    }
}

struct DCHSearchNode
{
    uint16_t dist_index;
    NodeID v;
    NodeID w;
    distance_t distance;
    uint16_t path_count;
    bool operator<(const DCHSearchNode &other) const { return dist_index < other.dist_index; }
    DCHSearchNode(uint16_t dist_index, NodeID v, NodeID w, distance_t distance, uint16_t path_count) : dist_index(dist_index), v(v), w(w), distance(distance), path_count(path_count) {}
};

struct ICHSearchNode
{
    NodeID v;
    uint16_t i;
    distance_t distance;
    uint16_t path_count;
    ICHSearchNode(NodeID v, uint16_t i, distance_t distance, uint16_t path_count) : v(v), i(i), distance(distance), path_count(path_count) {}
};

struct ICHSearchNode_P
{
    NodeID v;
    distance_t distance;
    uint16_t path_count;
    ICHSearchNode_P(NodeID v, distance_t distance, uint16_t path_count) : v(v), distance(distance), path_count(path_count) {}
};

////////////////////// Shortcut Count Graph Maintenance

void Graph::merge_edges(vector<pair<edge_t,edata_t> > &v)
{
    size_t v_size = v.size();
    if (v_size < 2)
        return;
    sort(v.begin(), v.end());

    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
    {
        if (v[next].first == v[last_distinct].first)
        {
            // combine with last distinct edge: replace or increase count
            if (v[next].second.first < v[last_distinct].second.first) {
                v[last_distinct].second = v[next].second;
            } else if (v[next].second.first == v[last_distinct].second.first)
                v[last_distinct].second.second += v[next].second.second;
        }
        else
            v[++last_distinct] = v[next];
    }
    v.erase(v.begin() + (last_distinct + 1), v.end());
}

void Graph::GS_Dec(ContractionHierarchy &ch, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates, vector<pair<edge_t, edata_t> > &C) 
{
    priority_queue<DCHSearchNode> q; NodeID a, b;
    for(pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > iter: updates) {

        a = iter.second.first, b = iter.second.second;
        if(ch.nodes[a].dist_index < ch.nodes[b].dist_index) swap(a, b);
        if(UpNeighbor(ch, a, b).distance >= iter.first.second)
            q.push(DCHSearchNode(ch.nodes[a].dist_index, a, b, iter.first.second, 1));
    }

    vector<pair<edge_t, edata_t> > temp;
    while(!q.empty()) {
        DCHSearchNode next = q.top(); q.pop();

        Neighbor &x = UpNeighbor(ch, next.v, next.w);
        if(next.distance < x.distance) {
            x.distance = next.distance;
            x.path_count = next.path_count;
        } else if(next.distance == x.distance) {
            x.path_count = x.path_count + next.path_count;
        } else
            continue;

        for(Neighbor n: ch.nodes[next.v].up_neighbors) {
            if(n.node != next.w) {
                distance_t dist = next.distance + n.distance;
                uint16_t path_count = next.path_count * n.path_count;

                a = next.w, b = n.node;
                if(ch.nodes[a].dist_index < ch.nodes[b].dist_index) swap(a, b);
                if(UpNeighbor(ch, a, b).distance >= dist)
                    q.push(DCHSearchNode(ch.nodes[a].dist_index, a, b, dist, path_count));
            }
        }
        C.push_back(make_pair(make_pair(next.v, next.w), make_pair(next.distance, next.path_count)));
    }
    merge_edges(C);
}

void Graph::GS_Inc(ContractionHierarchy &ch, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates, vector<pair<edge_t, edata_t> > &C) 
{
    priority_queue<DCHSearchNode> q; NodeID a, b;
    for(pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > iter: updates) {

        a = iter.second.first, b = iter.second.second;
        if(ch.nodes[a].dist_index < ch.nodes[b].dist_index) swap(a, b);

        if(UpNeighbor(ch, a, b).distance == iter.first.first)
            q.push(DCHSearchNode(ch.nodes[a].dist_index, a, b, iter.first.first, 1));
    }

    while(!q.empty()) {
        DCHSearchNode next = q.top(); q.pop();

        for(Neighbor &n: ch.nodes[next.v].up_neighbors) {
            if(n.node != next.w) {
                distance_t dist = next.distance + n.distance;
                uint16_t path_count = next.path_count * n.path_count;

                a = next.w, b = n.node;
                if(ch.nodes[a].dist_index < ch.nodes[b].dist_index) swap(a, b);
                if(UpNeighbor(ch, a, b).distance == dist)
                    q.push(DCHSearchNode(ch.nodes[a].dist_index, a, b, dist, path_count));
            }
        }

        Neighbor &x = UpNeighbor(ch, next.v, next.w);
        if(x.path_count > next.path_count) {
            x.path_count = x.path_count - next.path_count;
        } else {
            x.distance = infinity; x.path_count = 1;
            for (Neighbor &n : node_data[next.v].neighbors) {
                if (n.node == next.w) {
                    x.distance = n.distance;
                    break;
                }
            }

            size_t i = 0, j = 0;
            while (i < ch.nodes[next.v].down_neighbors.size() && j < ch.nodes[next.w].down_neighbors.size()) {
                a = ch.nodes[next.v].down_neighbors[i]; b = ch.nodes[next.w].down_neighbors[j];
                if (a < b) i++;
                else if (b < a) j++;
                else {
                    Neighbor &av = UpNeighbor(ch, a, next.v);
                    Neighbor &aw = UpNeighbor(ch, a, next.w);
                    distance_t dist = av.distance + aw.distance;
                    uint16_t path_count = av.path_count * aw.path_count;
                    if(dist < x.distance) {
                        x.distance = dist;
                        x.path_count = path_count;
                    } else if(dist == x.distance) {
                        x.path_count = x.path_count + path_count;
                    }
                    i++; j++;
                }
            }
        }
	C.push_back(make_pair(make_pair(next.v, next.w), make_pair(next.distance, next.path_count)));
    }
    merge_edges(C);
}

////////////////////// 2-Hop Count Labeling Maintenance

void Graph::DCL_Dec(ContractionHierarchy &ch, ContractionIndex &ci, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates) 
{
    vector<pair<edge_t, edata_t> > C;
    GS_Dec(ch, updates, C);

    //update distances involving ancestors
    util::min_bucket_queue<ICHSearchNode> q;
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first <= a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];

                if(a.distances()[i] >= dist) {
                    uint16_t path_count = iter.second.second * b.paths()[i];
                    q.push(ICHSearchNode(iter.first.first, i, dist, path_count), ch.nodes[iter.first.first].dist_index);
                }
            }
        }
    }

    // update distances involving descendants
    while(!q.empty()) {
        ICHSearchNode next = q.pop();

        FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
        if(cv.distances()[next.i] > next.distance) {
            cv.distances()[next.i] = next.distance;
            cv.paths()[next.i] = next.path_count;
        } else if(cv.distances()[next.i] == next.distance) {
            cv.paths()[next.i] = cv.paths()[next.i] + next.path_count;
        } else
            continue;

        // queue updates for descendants
        for(NodeID u: ch.nodes[next.v].down_neighbors) {
            Neighbor &x = UpNeighbor(ch, u, next.v);
            distance_t dist = x.distance + next.distance;

            FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
            if(cu.distances()[next.i] >= dist) {
                uint16_t path_count = x.path_count * next.path_count;
                q.push(ICHSearchNode(u, next.i, dist, path_count), ch.nodes[u].dist_index);
            }
        }
    }
}

void Graph::DCL_Inc(ContractionHierarchy &ch, ContractionIndex &ci, std::vector<std::pair<std::pair<distance_t, distance_t>, std::pair<NodeID, NodeID> > >& updates) 
{
    vector<pair<edge_t, edata_t> > C;
    GS_Inc(ch, updates, C);

    //update distances involving ancestors
    util::min_bucket_queue<ICHSearchNode> q;
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first == a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];
                uint16_t path_count = iter.second.second * b.paths()[i];

                if(dist == a.distances()[i])
                    q.push(ICHSearchNode(iter.first.first, i, dist, path_count), ch.nodes[iter.first.first].dist_index);
            }
        }
    }

    // update distances involving descendants
    while(!q.empty()) {
        ICHSearchNode next = q.pop();

        // update descendants
        FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
        for(NodeID u: ch.nodes[next.v].down_neighbors) {
            Neighbor &x = UpNeighbor(ch, u, next.v);
            FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
            distance_t dist = x.distance + cv.distances()[next.i];
            uint16_t path_count = x.path_count * next.path_count;

            if(dist == cu.distances()[next.i])
                q.push(ICHSearchNode(u, next.i, dist, path_count), ch.nodes[u].dist_index);
        }

        if(cv.paths()[next.i] > next.path_count) { // update path count, distance does not change
            cv.paths()[next.i] = cv.paths()[next.i] - next.path_count;
        } else { // recompute distance and path count
            cv.distances()[next.i] = infinity;
            for(Neighbor &u: ch.nodes[next.v].up_neighbors) {
                if(ch.nodes[u.node].dist_index >= next.i) {
                    Neighbor &x = UpNeighbor(ch, next.v, u.node);
                    FlatCutIndex cu = ci.get_contraction_label(u.node).cut_index;
                    distance_t dist = x.distance + cu.distances()[next.i];
                    uint16_t path_count = x.path_count * cu.paths()[next.i];

                    if(dist < cv.distances()[next.i]) {
                        cv.distances()[next.i] = dist;
                        cv.paths()[next.i] = path_count;
                    } else if(dist == cv.distances()[next.i]) {
                        cv.paths()[next.i] = cv.paths()[next.i] + path_count;
                    }
                }
            }
        }
    }
}

////////////////////// Parallel Maintenance

void Graph::DCL_Dec_Par(ContractionHierarchy &ch, ContractionIndex &ci, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates) {

    vector<thread> threads;
    auto dhcldec = [this](ContractionHierarchy &ch, ContractionIndex& ci, util::TSBucketQueue<ICHSearchNode_P>& que) {

        util::min_bucket_queue<ICHSearchNode_P> bq;
        vector<ICHSearchNode_P> bucket; size_t label_index;
        while (que.next_bucket(bucket, label_index))
        {
            for (ICHSearchNode_P obj : bucket)
                bq.push(obj, label_index);

            // update distances involving descendants
            while(!bq.empty()) {
                ICHSearchNode_P next = bq.pop();

                FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
                if(cv.distances()[label_index] > next.distance) {
                    cv.distances()[label_index] = next.distance;
                    cv.paths()[label_index] = next.path_count;
                } else if(cv.distances()[label_index] == next.distance) {
                    cv.paths()[label_index] = cv.paths()[label_index] + next.path_count;
                } else
                    continue;

                // queue updates for descendants
                for(NodeID u: ch.nodes[next.v].down_neighbors) {
                    Neighbor &x = UpNeighbor(ch, u, next.v);
                    distance_t dist = x.distance + next.distance;

                    FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
                    if(cu.distances()[label_index] >= dist) {
                        uint16_t path_count = x.path_count * next.path_count;
                        bq.push(ICHSearchNode_P(u, dist, path_count), label_index);
                    }
                }
            }
        }
    };

    vector<pair<edge_t, edata_t> > C;
    GS_Dec(ch, updates, C);

    //update distances involving ancestors
    util::TSBucketQueue<ICHSearchNode_P> grouping;
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first <= a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];

                if(a.distances()[i] >= dist) {
                    uint16_t path_count = iter.second.second * b.paths()[i];
                    grouping.push(ICHSearchNode_P(iter.first.first, dist, path_count), i);
                }
            }
        }
    }

    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads.push_back(thread(dhcldec, std::ref(ch), std::ref(ci), std::ref(grouping)));
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads[i].join();
}

void Graph::DCL_Inc_Par(ContractionHierarchy &ch, ContractionIndex &ci, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates) {

    vector<thread> threads;
    auto dhclinc = [this](ContractionHierarchy &ch, ContractionIndex& ci, util::TSBucketQueue<ICHSearchNode_P>& que) {

        util::min_bucket_queue<ICHSearchNode_P> bq;
        vector<ICHSearchNode_P> bucket; size_t label_index;
        while (que.next_bucket(bucket, label_index))
        {
            // move items to bucket queue
            for (ICHSearchNode_P obj: bucket)
                bq.push(obj, label_index);

            while(!bq.empty()) {
                ICHSearchNode_P next = bq.pop();

                // update descendants
                FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
                for(NodeID u: ch.nodes[next.v].down_neighbors) {
                    Neighbor &x = UpNeighbor(ch, u, next.v);
                    FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
                    distance_t dist = x.distance + cv.distances()[label_index];

                    if(dist == cu.distances()[label_index]) {
                        uint16_t path_count = x.path_count * next.path_count;
                        bq.push(ICHSearchNode_P(u, dist, path_count), label_index);
                    }
                }

                if(cv.paths()[label_index] > next.path_count) { // update path count, distance does not change
                    cv.paths()[label_index] = cv.paths()[label_index] - next.path_count;
                } else { // recompute distance and path count
                    cv.distances()[label_index] = infinity;
                    for(Neighbor &u: ch.nodes[next.v].up_neighbors) {
                        if(ch.nodes[u.node].dist_index >= label_index) {
                            Neighbor &x = UpNeighbor(ch, next.v, u.node);
                            FlatCutIndex cu = ci.get_contraction_label(u.node).cut_index;
                            distance_t dist = x.distance + cu.distances()[label_index];
                            uint16_t path_count = x.path_count * cu.paths()[label_index];

                            if(dist < cv.distances()[label_index]) {
                                cv.distances()[label_index] = dist;
                                cv.paths()[label_index] = path_count;
                            } else if(dist == cv.distances()[label_index]) {
                                cv.paths()[label_index] = cv.paths()[label_index] + path_count;
                            }
                        }
                    }
                }
            }
        }
    };

    vector<pair<edge_t, edata_t> > C;
    GS_Inc(ch, updates, C);

    //update distances involving ancestors
    util::TSBucketQueue<ICHSearchNode_P> grouping;
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first == a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];

		if(dist == a.distances()[i]) {
                    uint16_t path_count = iter.second.second * b.paths()[i];
                    grouping.push(ICHSearchNode_P(iter.first.first, dist, path_count), i);
                }
            }
        }
    }

    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads.push_back(thread(dhclinc, std::ref(ch), std::ref(ci), std::ref(grouping)));
    for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
        threads[i].join();
}


////////////////////// Optimized Maintenance (Avoiding Propagation Overhead in Count Updates)

util::min_bucket_queue<ICHSearchNode> q;
void Graph::EnqueAndUpdate_d(ContractionHierarchy &ch, ContractionIndex &ci, NodeID v, uint16_t i, distance_t dist, uint16_t path_count) {

    //store original values in queue
    FlatCutIndex cv = ci.get_contraction_label(v).cut_index;
    if((cv.paths()[i] & (1 << 15)) == 0) {
        q.push(ICHSearchNode(v, i, cv.distances()[i], cv.paths()[i]), ch.nodes[v].dist_index);
        // setting the highest bit
        cv.paths()[i] = cv.paths()[i] | (1 << 15);
    }

    // update values in index
    if(cv.distances()[i] > dist) {
        cv.distances()[i] = dist;
        cv.paths()[i] = path_count | (1 << 15);
    } else {
        cv.paths()[i] = cv.paths()[i] + path_count;
    }
}

void Graph::EnqueAndUpdate_i(ContractionHierarchy &ch, ContractionIndex &ci, NodeID v, uint16_t i, uint16_t path_count) {

    //store original values in queue
    FlatCutIndex cv = ci.get_contraction_label(v).cut_index;
    if((cv.paths()[i] & (1 << 15)) == 0) {
        q.push(ICHSearchNode(v, i, cv.distances()[i], cv.paths()[i]), ch.nodes[v].dist_index);
        // setting the highest bit
        cv.paths()[i] = cv.paths()[i] | (1 << 15);
    }

    // update count in index
    cv.paths()[i] = cv.paths()[i] - path_count;
}

void Graph::DCL_Dec_Opt(ContractionHierarchy &ch, ContractionIndex &ci, vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > >& updates) {

    vector<pair<edge_t, edata_t> > C;
    GS_Dec(ch, updates, C);

    //update distances involving ancestors
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first <= a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];

                if(a.distances()[i] >= dist) {
                    uint16_t path_count = iter.second.second * b.paths()[i];
                    EnqueAndUpdate_d(ch, ci, iter.first.first, i, dist, path_count);
                }
            }
        }
    }

    // update distances involving descendants
    while(!q.empty()) {
        ICHSearchNode next = q.pop();

        uint16_t convex_path_count = 0;
        FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
        // resetting the highest bit
        cv.paths()[next.i] &= ~(1 << 15);
        if(cv.distances()[next.i] == next.distance) {
            convex_path_count = cv.paths()[next.i] - next.path_count;
        } else if(cv.distances()[next.i] < next.distance) {
            convex_path_count = cv.paths()[next.i];
        } else
            continue;

        // queue updates for descendants
        for(NodeID u: ch.nodes[next.v].down_neighbors) {
            Neighbor &x = UpNeighbor(ch, u, next.v);
            distance_t dist = x.distance + cv.distances()[next.i];

            FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
            if(cu.distances()[next.i] >= dist) {
                uint16_t path_count = x.path_count * convex_path_count;
                EnqueAndUpdate_d(ch, ci, u, next.i, dist, path_count);
            }
        }
    }
}

void Graph::DCL_Inc_Opt(ContractionHierarchy &ch, ContractionIndex &ci, std::vector<std::pair<std::pair<distance_t, distance_t>, std::pair<NodeID, NodeID> > >& updates) {

    vector<pair<edge_t, edata_t> > C;
    GS_Inc(ch, updates, C);

    //update distances involving ancestors
    for(pair<edge_t, edata_t> iter: C) {
        FlatCutIndex a = ci.get_contraction_label(iter.first.first).cut_index;
        if(iter.second.first == a.distances()[ch.nodes[iter.first.second].dist_index]) {

            FlatCutIndex b = ci.get_contraction_label(iter.first.second).cut_index;
            for(size_t i = 0; i <= ch.nodes[iter.first.second].dist_index; i++) {
                distance_t dist = iter.second.first + b.distances()[i];

                if(dist == a.distances()[i]) {
                    uint16_t path_count = iter.second.second * b.paths()[i];
                    EnqueAndUpdate_i(ch, ci, iter.first.first, i, path_count);
                }
            }
        }
    }

    // update distances involving descendants
    while(!q.empty()) {
        ICHSearchNode next = q.pop();

        FlatCutIndex cv = ci.get_contraction_label(next.v).cut_index;
        // resetting the highest bit
        cv.paths()[next.i] &= ~(1 << 15);
        uint16_t convex_path_count = next.path_count - cv.paths()[next.i];

        // update descendants
        for(NodeID u: ch.nodes[next.v].down_neighbors) {
            Neighbor &x = UpNeighbor(ch, u, next.v);
            FlatCutIndex cu = ci.get_contraction_label(u).cut_index;
            distance_t dist = x.distance + cv.distances()[next.i];

            if(dist == cu.distances()[next.i]) {
                uint16_t path_count = x.path_count * convex_path_count;
                EnqueAndUpdate_i(ch, ci, u, next.i, path_count);
            }
        }

        if(cv.paths()[next.i] == 0) {
            cv.distances()[next.i] = infinity;
            for(Neighbor &w: ch.nodes[next.v].up_neighbors) {
                if(ch.nodes[w.node].dist_index >= next.i) {
                    Neighbor &x = UpNeighbor(ch, next.v, w.node);
                    FlatCutIndex cw = ci.get_contraction_label(w.node).cut_index;
                    distance_t dist = x.distance + cw.distances()[next.i];
                    uint16_t path_count = x.path_count * cw.paths()[next.i];

                    if(dist < cv.distances()[next.i]) {
                        cv.distances()[next.i] = dist;
                        cv.paths()[next.i] = path_count;
                    } else if(dist == cv.distances()[next.i]) {
                        cv.paths()[next.i] = cv.paths()[next.i] + path_count;
                    }
                }
            }
        }
    }
}

void Graph::contract_seq(ContractionIndex &ci, vector<pair<pair<distance_t,distance_t>, NodeID> >& contracted_updates) {

    // we start searches in order of original distance and cancel searches for nodes already updated
    sort(contracted_updates.begin(), contracted_updates.end());
    // search from each update using DFS
    vector<SearchNode> stack;
    for (auto& it: contracted_updates) {
        // if already processed, distance offset has (probably) changed
        if (it.first.first != ci.get_contraction_label(it.second).distance_offset)
            continue;
        stack.push_back(SearchNode(it.first.second, it.second));
        while (!stack.empty()) {
            SearchNode next = stack.back(); stack.pop_back();

            // update label
            ci.update_distance_offset(next.node, next.distance);
            // enqueue neighbors
            for (Neighbor n : node_data[next.node].neighbors) {
                if (ci.get_contraction_label(n.node).parent == next.node)
                    stack.push_back(SearchNode(next.distance + n.distance, n.node));
            }
        }
    }
}

//--------------------------- Graph debug ---------------------------

bool Graph::is_consistent() const
{
    // all nodes in subgraph have correct subgraph ID assigned
    for (NodeID node : nodes)
        if (node_data[node].subgraph_id != subgraph_id)
        {
            DEBUG("wrong subgraph ID for " << node << " in " << *this);
            return false;
        }
    // number of nodes with subgraph_id of subgraph equals number of nodes in subgraph
    size_t count = 0;
    for (NodeID node = 0; node < node_data.size(); node++)
        if (contains(node))
            count++;
    if (count != nodes.size())
    {
        DEBUG(count << "/" << nodes.size() << " nodes contained in " << *this);
        return false;
    }
    return true;
}

bool Graph::is_undirected() const
{
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].neighbors)
        {
            bool found = false;
            for (Neighbor nn : node_data[n.node].neighbors)
                if (nn.node == node && nn.distance == n.distance)
                {
                    found = true;
                    break;
                }
            if (!found)
                return false;
        }
    return true;
}

vector<pair<distance_t,distance_t>> Graph::distances() const
{
    vector<pair<distance_t,distance_t>> d;
    for (const Node &n : node_data)
        d.push_back(pair(n.distance, n.outcopy_distance));
    return d;
}

vector<pair<NodeID,NodeID>> Graph::flow() const
{
    vector<pair<NodeID,NodeID>> f;
    for (const Node &n : node_data)
        f.push_back(pair(n.inflow, n.outflow));
    return f;
}

NodeID Graph::random_node() const
{
    return nodes[rand() % nodes.size()];
}

pair<NodeID,NodeID> Graph::random_pair(size_t steps) const
{
    if (steps < 1)
        return make_pair(random_node(), random_node());
    NodeID start = random_node();
    NodeID stop = start;
    for (size_t i = 0; i < steps; i++)
    {
        NodeID n = NO_NODE;
        do
        {
            n = util::random(node_data[stop].neighbors).node;
        } while (!contains(n));
        stop = n;
    }
    return make_pair(start, stop);
}

// generate batch of random node pairs, filtered into buckets by distance (as for H2H/P2H)
void Graph::random_pairs(vector<vector<pair<NodeID,NodeID>>> &buckets, distance_t min_dist, size_t bucket_size, const ContractionIndex &ci)
{
    assert(buckets.size() > 0);
    const distance_t max_dist = diameter(true);
    const double x = pow(static_cast<double>(max_dist) / min_dist, 1.0 / buckets.size());
    vector<distance_t> bucket_caps;
    // don't push last cap - implied and works nicely with std::upper_bound
    for (size_t i = 1; i < buckets.size(); i++)
        bucket_caps.push_back(min_dist * pow(x, i));
    size_t todo = buckets.size();
    cout << "|";
    size_t counter = 0;
    while (todo)
    {
        // generate some queries using random walks for speedup
        pair<NodeID, NodeID> q = ++counter % 5 ? make_pair(random_node(), random_node()) : random_pair(1 + rand() % 100);
        distance_t d = ci.get_distance(q.first, q.second);
        if (d >= min_dist)
        {
            size_t bucket = upper_bound(bucket_caps.begin(), bucket_caps.end(), d) - bucket_caps.begin();
            if (buckets[bucket].size() < bucket_size)
            {
                buckets[bucket].push_back(q);
                if (buckets[bucket].size() == bucket_size)
                {
                    todo--;
                    cout << bucket << "|" << flush;
                }
            }
        }
    }
}

void Graph::randomize()
{
    shuffle(nodes.begin(), nodes.end(), default_random_engine());
    for (NodeID node : nodes)
        shuffle(node_data[node].neighbors.begin(), node_data[node].neighbors.end(), default_random_engine());
}

void print_graph(const Graph &g, ostream &os)
{
    vector<Edge> edges;
    g.get_edges(edges);
    sort(edges.begin(), edges.end());
    os << "p sp " << Graph::super_node_count() << " " << edges.size() << endl;
    for (Edge e : edges)
        os << "a " << e.a << ' ' << e.b << ' ' << e.d << endl;
}

void read_graph(Graph &g, istream &in)
{
    char line_id;
    uint32_t v, w, d;

    while (in >> line_id) {
        switch (line_id)
        {
        case 'p':
            in.ignore(3);
            in >> v;
            in.ignore(1000, '\n');
            g.resize(v);
            break;
        case 'a':
            in >> v >> w >> d;
            g.add_edge(v, w, d, true);
            break;
        default:
            in.ignore(1000, '\n');
        }
    }
    g.remove_isolated();
}

//--------------------------- ostream -------------------------------

// for easy distance printing
struct Dist
{
    distance_t d;
    Dist(distance_t d) : d(d) {}
};

static ostream& operator<<(ostream& os, Dist distance)
{
    if (distance.d == infinity)
        return os << "inf";
    else
        return os << distance.d;
}

// for easy bit string printing
struct BitString
{
    uint64_t bs;
    BitString(uint64_t bs) : bs(bs) {}
};

static ostream& operator<<(ostream& os, BitString bs)
{
    size_t len = bs.bs & 63ul;
    uint64_t bits = bs.bs >> 6;
    for (size_t i = 0; i < len; i++)
    {
        os << (bits & 1);
        bits >>= 1;
    }
    return os;
}

ostream& operator<<(ostream& os, const CutIndex &ci)
{
    return os << "CI(p=" << bitset<8>(ci.partition) << ",c=" << (int)ci.cut_level
        << ",di=" << ci.dist_index << ",d=" << ci.distances << ")";
}

ostream& operator<<(ostream& os, const FlatCutIndex &ci)
{
    uint64_t partition_bitvector = *ci.partition_bitvector();
    vector<uint16_t> dist_index(ci.dist_index(), ci.dist_index() + ci.cut_level() + 1);
    vector<distance_t> distances(ci.distances(), ci.distances() + ci.label_count());
    return os << "FCI(pb=" << BitString(partition_bitvector) << ",di=" << dist_index << ",d=" << distances << ")";
}

ostream& operator<<(ostream& os, const ContractionLabel &cl)
{
    return os << "CL(" << cl.cut_index << ",d=" << cl.distance_offset << ",p=" << cl.parent << ")";
}

ostream& operator<<(ostream& os, const Neighbor &n)
{
    if (n.distance == 1)
        return os << n.node;
    else
        return os << n.node << "@" << Dist(n.distance);
}

ostream& operator<<(ostream& os, const Node &n)
{
    return os << "N(" << n.subgraph_id << "#" << n.neighbors << ")";
}

ostream& operator<<(ostream& os, const Partition &p)
{
    return os << "P(" << p.left << "|" << p.cut << "|" << p.right << ")";
}

ostream& operator<<(ostream& os, const DiffData &dd)
{
    return os << "D(" << dd.node << "@" << dd.dist_a << "-" << dd.dist_b << "=" << dd.diff() << ")";
}

ostream& operator<<(ostream& os, const Graph &g)
{
#ifdef MULTI_THREAD
    g.node_data.normalize();
#endif
    return os << "G(" << g.subgraph_id << "#" << g.nodes << " over " << g.node_data << ")";
}

} // road_network
