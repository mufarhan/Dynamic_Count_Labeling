#include "road_network.h"
#include "util.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace road_network;

int main(int argc, char** argv)
{
    ifstream ifs(argv[1]);
    Graph g;
    read_graph(g, ifs);
    ifs.close();

    ifs.open(string(argv[2]) + string("_cl"));
    ContractionIndex con_index(ifs);
    ifs.close();
    ifs.open(string(argv[2]) + string("_gs"));
    ContractionHierarchy ch(ifs);
    ifs.close();

    vector<pair<pair<distance_t,distance_t>, NodeID> > contracted_updates;
    vector<pair<pair<distance_t, distance_t>, pair<NodeID, NodeID> > > updates;
    NodeID a, b; distance_t weight;

    ifs.open(argv[3]);
    while(ifs >> a >> b >> weight) {

        distance_t new_weight;
        if(argv[4][0] == 'd')
            new_weight = weight * 0.5;
        else if(argv[4][0] == 'i')
            new_weight = weight * 1.5;

        g.update_edge(a, b, new_weight);
        g.update_edge(b, a, new_weight);

        ContractionLabel x = con_index.get_contraction_label(a), y = con_index.get_contraction_label(b);
        if (con_index.is_contracted(a) || con_index.is_contracted(b)) {
            if (x.distance_offset > y.distance_offset)
                contracted_updates.push_back(make_pair(make_pair(x.distance_offset, y.distance_offset + new_weight), a));
            else if(x.distance_offset < y.distance_offset)
                contracted_updates.push_back(make_pair(make_pair(y.distance_offset, x.distance_offset + new_weight), b));
            continue;
        }
        updates.push_back(make_pair(make_pair(weight, new_weight), make_pair(a, b)));
    }
    ifs.close();

    util::start_timer();
    if(argv[4][0] == 'd') {
      //sequential
      g.DCL_Dec(ch, con_index, updates);
      // optimized
      //g.DCL_Dec_Opt(ch, con_index, updates);
      //parallel
      //g.DCL_Dec_Par(ch, con_index, updates);
    } else if(argv[4][0] == 'i') {
      //sequential
      g.DCL_Inc(ch, con_index, updates);
      // optimized
      //g.DCL_Inc_Opt(ch, con_index, updates);
      //parallel
      //g.DCL_Inc_Par(ch, con_index, updates);
    }
    g.contract_seq(con_index, contracted_updates);

    double random_update_time = util::stop_timer();
    cout << "ran " << updates.size() << " random updates in " << random_update_time << endl;
    
    return 0;
}
