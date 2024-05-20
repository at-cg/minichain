#include "graphUtils.h"
#include <assert.h>

graphUtils::graphUtils(gfa_t *g)
{
    this->g = g;
}

// Read and store the Graph in Adjacency list
void graphUtils::read_graph()
{
    uint v;
    n_vtx = gfa_n_vtx(g);
    // Resize node_len
    node_len.resize(n_vtx, 0);
    adj_.resize(n_vtx); 
    /* Node_len */
    for (int v = 0; v < n_vtx/2; v++)
    {
        gfa_seg_t segment = (g)->seg[v];
        int len =  segment.len;
        node_len[2*v] = len;
        node_len[2*v + 1] = len;
    }
    // look for all the edges , if the sum of all the 
    // edges are zero then, that's a linear reference
    u_int num_edges = 0;
    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        int v_ = av->v_lv >> 32;
        // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::endl; 
        for (int i = 0; i < n_edges; i++)
        {
            num_edges++;
        }
    }


    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        if (num_edges == 0) // Linear Reference
        {
            lin_ref = 1; // Mark as a linear reference
        }else
        {
            int v_ = av->v_lv >> 32;
            // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::endl; 
            for (int i = 0; i < n_edges; i++)
            {
                uint w = av[i].w;
                adj_[v_].push_back(w);
            }
        }
    }


    // Read walks
    if (g->n_walk > 0) is_hprc_gfa = true;

    // flip gfa walks

    for (int n = 0; n < g->walk[0].n_v; n++)
    {
        int v = g->walk[0].v[n];
    }
    
    num_walks = g->n_walk; // count number of walks
    for (size_t w = 0; w < g->n_walk; w++)
    {
        // store the walk in vector
        std::vector<int> walk;
        std::vector<int> rev_walk;
        for (size_t n = 0; n < g->walk[w].n_v; n++)
        {
            int v = g->walk[w].v[n];
            // walk.push_back((v)); // Forward Strand
            haps[v].push_back(g->walk[w].sample); // Map walk to sample
            // rev_walk.push_back((v|1)); // Reverse Strand
            fwd_walk_map[g->walk[w].sample].push_back(v); // Map reverse walk to sample
            rev_walk_map[g->walk[w].sample].push_back((v|1)); // Map reverse walk to sample
        }

        walk_ids.push_back(g->walk[w].sample);

        walk_map[w] = g->walk[w].sample;
        // walk_map[w + g->n_walk] = g->walk[w].sample;
        // fwd_walk_map[g->walk[w].sample].push_back(w);
        // fwd_walk_map[g->walk[w].sample].push_back(w + g->n_walk);

        // std::cerr << " haplotype : " << g->walk[w].sample << " walk id : " << w << std::endl; 

        if (param_z) std::cerr << "Walk id : " << walk_map[w]<< std::endl;
        
        // std::reverse(walk.begin(), walk.end()); // Just checking whether reverse is needed or not (We don't need)
        // std::reverse(rev_walk.begin(), rev_walk.end()); // Just checking whether reverse is needed or not (We don't need)
        // walks.push_back(walk);
        // rev_walks.push_back(rev_walk);
        walk.clear();
        rev_walk.clear();
    }

    // reverse walks and rev_walks
    std::reverse(walks.begin(), walks.end());
    std::reverse(rev_walks.begin(), rev_walks.end());

}

void graphUtils::print_graph()
{
    std::cerr << " This is Graph " << std::endl;
    for (size_t i = 0; i < n_vtx; i++)
    {
        std::cerr << i << "->";
        for (int &x : adj_[i])
        {
            std::cerr << x << "->";
        }
        std::cerr << std::endl;
    }
}

void DFS(std::vector<std::vector<int>> un_adj, std::vector<bool> &visited, int u, std::vector<int> &component, int num_components)
{
    std::stack<int> stack;

    // Push the source node
    stack.push(u);

    while (!stack.empty())
    {
        int s = stack.top();
        stack.pop();
        component[s] = num_components;

        if (!visited[s])
        {
            visited[s] = true;
        }
        
        for (auto v:un_adj[s])
        {
            if (!visited[v])
            {
                stack.push(v);
            }
        }
    } 

}

void graphUtils::Connected_components()
{
    size_t num_components;
    if (lin_ref == 1)
    {
        num_components = n_vtx;
        component.resize(n_vtx);
        for (int u = 0; u < n_vtx; u++) // Put unique nodes in unique cids
        {
            component[u] = u;
        }
    }else
    {

        // Create Adjacency list
        std::vector<std::vector<int>> un_adj;
        un_adj.resize(n_vtx); // u -> v
        for (int u = 0; u < n_vtx; u++)
        {
            for (auto v:adj_[u])
            {
                un_adj[u].push_back(v);
                un_adj[v].push_back(u);
            }
        }

        num_components = 0;
        component.resize(n_vtx);  
        // run DFS
        std::vector<bool> visited(n_vtx,false);
        for (int u = 0; u < n_vtx; u++)
        {
            if (visited[u] == false)
            {
                DFS(un_adj,visited,u,component,num_components);
                num_components++;
            }
        } 
        un_adj.clear(); 
        visited.clear(); 
    }

    num_cid = num_components;
    // Storing Connected Components
    conn_comp.resize(num_components);
    for (int i = 0; i < n_vtx; i++)
    {
        conn_comp[component[i]].push_back(i); // Add Vertex to it's component
    }

    
    component_idx.resize(num_cid);
    idx_component.resize(num_cid);
    // Map Components
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        component_idx[cid].resize(conn_comp[cid].size());
        idx_component[cid].resize(n_vtx);
        for (size_t j = 0; j < conn_comp[cid].size(); j++)
        {
            component_idx[cid][j] = conn_comp[cid][j]; // map global to local
            idx_component[cid][conn_comp[cid][j]] = j; // map local to global
        }
    }

    // Create Adjacency list for all connected components
    adj_cc = new std::vector<std::vector<int>>[num_components];
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        adj_cc[cid].resize(conn_comp[cid].size());
        for (auto const &j : conn_comp[cid])
        {
            for (int k = 0; k < adj_[j].size(); k++)
            {
                adj_cc[cid][idx_component[cid][j]].push_back(idx_component[cid][adj_[j][k]]); // Map local index
            }
        }
    }

    // HPRC GFA
    if (is_hprc_gfa)
    {
        local_walks.resize(num_cid, std::vector<std::vector<int>>(num_walks));
        for (auto walk : walk_map)
        {
            std::reverse(rev_walk_map[walk.second].begin(), rev_walk_map[walk.second].end());
        }

        for (auto walk : walk_map)
        {
            int cid = component[fwd_walk_map[walk.second][0]];
            std::vector<int> temp;
            for (auto u_w : fwd_walk_map[walk.second])
            {
                temp.push_back(idx_component[cid][u_w]);
            }
            local_walks[cid][walk.first] = temp;
            temp.clear();

            cid = component[rev_walk_map[walk.second][0]];
            for (auto u_w : rev_walk_map[walk.second])
            {
                temp.push_back(idx_component[cid][u_w]);
            }
            local_walks[cid][walk.first] = temp;
            temp.clear();
        }
    }
}


int graphUtils::is_cyclic() // Check cyclicity of component.

{
    // construct
    std::vector<std::vector<int>> in_degree;
    std::vector<std::vector<int>> out_degree;
    in_degree.resize(num_cid);
    out_degree.resize(num_cid);
    // queue
    std::vector<std::queue<int>> q;
    q.resize(num_cid);
    std::vector<bool> cid_cycle;
    cid_cycle.resize(num_cid);
    int cycle_count = 0;
    top_order.resize(num_cid);
    omp_set_dynamic(1);
    omp_set_num_threads(0);
    #pragma omp parallel for
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        // Intialize
        in_degree[cid].resize(adj_cc[cid].size(), 0);
        out_degree[cid].resize(adj_cc[cid].size(), 0);
        // Compute
        for (int v = 0; v < adj_cc[cid].size(); v++)
        {
            out_degree[cid][v] = adj_cc[cid][v].size(); // [0]-> 1,3 (in adj_cc[cid][0] so out_degree of 0 is 2)
        }
        for (int u = 0; u < adj_cc[cid].size(); u++)
        {
            for (auto const &v : adj_cc[cid][u])
            {
                in_degree[cid][v]++; // [0]-> 1,3 (in adj_cc[cid][0] so in_degree of 1 and 3 is 1 for u = 0 and compute for all u's in ad_cc[cid])
            }
        }
        // kahn's Topological sort
        for (int v = 0; v < adj_cc[cid].size(); v++)
        {
            if (in_degree[cid][v] == 0)
            {
                q[cid].push(v);
            }
        }
        int count = 0;
        while (!q[cid].empty())
        {
            int u = q[cid].front();
            q[cid].pop();
            for (auto const &v : adj_cc[cid][u])
            {
                if (--in_degree[cid][v] == 0)
                {
                    q[cid].push(v); // front -->[0,1,2]--> back , (when you push, queue adds v in incrementing order)
                }
            }
            count++; // make count of vertices in top_order[cid]
        }
        // verify component is acyclic
        if (count != adj_cc[cid].size())
        {
            cycle_count++;
        }
    }
    q.clear();
    in_degree.clear();
    out_degree.clear();
    // std::cerr << "[Connected components : " << num_cid << ", components with cycle : " << cycle_count << "]"<< std::endl;
    return cycle_count;
}

void graphUtils::topologicat_sort()
{

    // construct
    std::vector<std::vector<int>> in_degree;
    std::vector<std::vector<int>> out_degree;
    in_degree.resize(num_cid);
    out_degree.resize(num_cid);
    map_top_sort.resize(num_cid);

    // In node and out node
    in_node.resize(num_cid);
    out_node.resize(num_cid);
    // queue
    std::vector<std::queue<int>> q;
    q.resize(num_cid);
    omp_set_dynamic(1);
    omp_set_num_threads(0);
    #pragma omp parallel for
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        //    /*
        //     ###########################
        //     # kahn's Toplogical sort  #
        //     ###########################
        //     */

        // In degree and out degree computation

        // Intialize
        in_degree[cid].resize(adj_cc[cid].size(), 0);
        out_degree[cid].resize(adj_cc[cid].size(), 0);
        // Compute
        for (int v = 0; v < adj_cc[cid].size(); v++)
        {
            out_degree[cid][v] = adj_cc[cid][v].size(); // [0]-> 1,3 (in adj_cc[cid][0] so out_degree of 0 is 2)
        }
        for (int u = 0; u < adj_cc[cid].size(); u++)
        {
            for (auto const &v : adj_cc[cid][u])
            {
                in_degree[cid][v]++; // [0]-> 1,3 (in adj_cc[cid][0] so in_degree of 1 and 3 is 1 for u = 0 and compute for all u's in ad_cc[cid])
            }
        }
        // kahn's Topological sort
        for (int v = 0; v < adj_cc[cid].size(); v++)
        {
            if (in_degree[cid][v] == 0)
            {
                q[cid].push(v);
            }
        }
        int count = 0;
        while (!q[cid].empty())
        {
            int u = q[cid].front();
            q[cid].pop();
            top_order[cid].push_back(u);
            for (auto const &v : adj_cc[cid][u])
            {
                if (--in_degree[cid][v] == 0)
                {
                    q[cid].push(v); // front -->[0,1,2]--> back , (when you push, queue add v in incrementing order)
                }
            }
            count++; // make count of vertices in top_order[cid]
        }
        // verify component is acyclic
        if (count != adj_cc[cid].size())
        {
            std::cout << " Can't do Topological ordering : cycle exist " << std::endl; 
        }
        // std::cerr << " Top sort for cid : " << cid << std::endl;
        // for (size_t i = 0; i < top_order[cid].size(); i++)
        // {
        //     std::cerr << component_idx[cid][top_order[cid][i]] << std::endl;
        // }
        // std::cerr << std::endl;
        
    }

    // in_node and out_node computation
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        // Compute In-node and Out-node
        in_node[cid].resize(adj_cc[cid].size());  // all nodes
        out_node[cid].resize(adj_cc[cid].size()); // all nodes
        // computing in_nodes and out_nodes
        for (int u = 0; u < adj_cc[cid].size(); u++)
        {
            for (auto const &v : adj_cc[cid][u])
            {
            in_node[cid][v].push_back(u);
            out_node[cid][u].push_back(v);
            }
        }
    }

    /* Mapping for Top_Sort */
    for (int cid = 0; cid < num_cid; cid++)
    {
        map_top_sort[cid].resize(top_order[cid].size());
        
        for (int v = 0; v < top_order[cid].size(); v++)
        {
            map_top_sort[cid][top_order[cid][v]] = v;
        }
    }
    in_degree.clear();
    out_degree.clear();
}

std::vector<std::vector<int>> graphUtils::shrink(int cid)
{
    /*
    ######################
    #     shrinking      #
    ######################
    */
    // Compute Shrinking    
    size_t N = adj_cc[cid].size();
    std::vector<int> cids;
    // fill with vertex id
    for (size_t i = 0; i < N; i++)
    {
        cids.push_back(i);
    }

    std::vector<int> covered(N, 0);
    std::vector<std::vector<int>> ret;
    int K = path_cover[cid].size(), inf = path_cover[cid].size();
    std::vector<int> starts(N, 0), ends(N, 0);
    std::map<std::pair<int, int>, int> edge_covered;
    for (auto path : path_cover[cid]) {
        for (int i = 0; i < path.size(); i++) {
            covered[path[i]]++;
        if (i > 0)
            edge_covered[{ path[i - 1], path[i] }]++;
        }
        starts[path[0]]++;
        ends[path.back()]++;
    }
    flowGraph fg(N * 2);
    // i_in = i, i_out = i + N
    // add r(i, j) = c(j,i) + f(i,j) - l(i,j)
    auto add = [&](int i, int j, int cap, int l, int ff) {
    // std::cerr << "add edge " << i << " " << j << "  " << cap << "  " << l << " " <<ff << std::endl;
    fg.add_edge(i, j, 0 + ff - l);
    fg.add_edge(j, i, cap - ff);
    };
    for (int i = 0; i < N; i++)
    for (size_t jid : out_node[cid][i]) {
            size_t j = jid;
            int ff = edge_covered.count({i, j}) ? edge_covered[{i, j}] : 0;
            add(i + N, j, inf, 0, ff);
        }
    for (int i = 0; i < N; i++) {
        add(i, i + N, inf, 1, covered[i]);
        add(fg.S, i, inf, 0, starts[i]);
        add(i + N, fg.T, inf, 0, ends[i]);
    }

    int total = inf;
    std::vector<int> Q(fg.N, 0), pre(fg.N, -1), d(fg.N, 0);
    while (1) {
        int Qsize = 0;
        Q[Qsize++] = fg.S;
        for (int i = 0; i < fg.N; i++) {
            pre[i] = -1;
            d[i] = 0;
        }
        d[fg.S] = 1;
        for (int idx = 0; idx < Qsize && d[fg.T] == 0;) {
            int i = Q[idx++];
            for (int e = fg.f[i]; e; e = fg.t[e]) {
                int j = fg.p[e]; 
                if (fg.c[e] > 0 && d[j] == 0) {
                    d[j] = 1;
                    pre[j] = e;
                    Q[Qsize++] = j;
                }
            }
        }
        if (d[fg.T] == 0) break;
        std::vector<int> tmp;
        int flow = fg.c[pre[fg.T]];
        for (int i = fg.T; ;) {
            tmp.push_back(i);
            int e = pre[i];
            if (e == -1) break;
            flow = std::min(flow, fg.c[e]);
            i = fg.p[e ^ 1];
        }
        for (int i = fg.T; ;) {
            int e = pre[i];
            if (e == -1) break;
            fg.c[e] -= flow;
            fg.c[e^1] += flow;
            i = fg.p[e ^ 1];
        }
        if (flow == 0) exit(1);
        total -= flow;
    // std::cerr << " Now shrink by " << flow << " to " << total << std::endl;
    }

    // std::cerr << " Minimum Flow :  " << total << std::endl;

    // convert flow back to path cover
    // ret = pc;
    // ret.resize(total);

    for (int itr = 0; itr < total; itr++) {
        std::vector<int> tmp;
        for (int i = fg.S; i != fg.T; ) {
            if (0 <= i && i < N)
                tmp.push_back(cids[i]);
            int nxt = -1;
            for (int e = fg.f[i]; e; e = fg.t[e]) {
                int j = fg.p[e];
                int ff = fg.c[e] + ((i < N && i + N == j) ? 1 : 0);
                if ((e & 1) == 0 && ff > 0) {
                    nxt = j;
                    fg.c[e]--;
                    break;
                }
            }
            if (nxt == -1) {
                std::cerr << i << " not found nxt " << std::endl;
                return ret; // return ret
            }
            i = nxt;
        }
        ret.push_back(tmp);
    }
    return ret;
}

bool check_MPC(std::vector<std::vector<int>> adj, std::vector<int> path_verify, int k, int cid)
{
    int count = 0;
    for (int i = 0; i < path_verify.size() - 1; i++)
    {
        int u = path_verify[i];
        int v = path_verify[i+1];
        for (auto x:adj[u])
        {
            if (x == v)
            {
                count++; // count the edges which will be "vertex - 1"
            }
            
        }
        
    }
    if (count + 1 ==  path_verify.size())
    {
        std::cerr << "cid = " << cid << " path #" << k+1<< " MPC is OK! " << std::endl;
        return true;
    }
}

void graphUtils::MPC()
{
    if(is_hprc_gfa)
    {
        // use GFA v1.1 walks as path cover
        path_cover = local_walks;

    }else {
        /*
        ######################
        # MPC (greedy cover) #
        ######################
        */
        path_cover.resize(num_cid);
        omp_set_dynamic(1);
        omp_set_num_threads(0);
        #pragma omp parallel for
        for (size_t cid = 0; cid < num_cid; cid++)
        {
            /* Greedy MPC! */
            int T = 0;
            int covered_count = 0;
            std::vector<int> covered;
            std::vector<int> max_cover;
            std::vector<int> pre;
            std::vector<int> path;
            int V = adj_cc[cid].size(); 
            covered.resize(V,0);
            max_cover.resize(V,0);
            pre.resize(V,-1); // Initialise as "None"
            while (covered_count < V)
            {
                for (size_t i = 0; i < max_cover.size(); i++)
                {
                    max_cover[i] = 0;
                    pre[i] = -1;
                }
                
                for (auto const & v: top_order[cid])
                {
                    if (covered[v] == 0)
                    {
                        max_cover[v]++;
                    }
                    
                    for (auto const& u : adj_cc[cid][v])
                    {
                        if (max_cover[u] < max_cover[v])
                        {
                            max_cover[u] = max_cover[v];
                            pre[u] = v;
                        }   
                    }
                }
                auto max = std::max_element(max_cover.begin(),max_cover.end());
                T = std::distance(max_cover.begin(),max); // argmax(max_cover[v])
                
                int new_covered = 0;
                while (covered_count < V)
                {
                    if (T == -1)
                    {
                        break;
                    }else
                    {
                        if (covered[T] == 0)
                        {
                            covered_count++;
                            covered[T] = 1;
                            new_covered++;
                        }
                        path.push_back(T);
                        T = pre[T];
                    }
                    
                }
                std::reverse(path.begin(),path.end());
                path_cover[cid].push_back(path);
                // std::cerr << "cid = " << cid << " path #" << path_cover[cid].size() << " : " << path.size() << " " << new_covered << " " << (V - covered_count) << std::endl;
                path.clear();
            }

            // Shrink Path Cover
            path_cover[cid] = shrink(cid);
            std::reverse(path_cover[cid].begin(),path_cover[cid].end());
        }
    }

    if (param_z)
    {
        /* Check MPC */
            for (size_t cid = 0; cid < num_cid; cid++)
        {
            bool check_path;
            for (int k = 0; k < path_cover[cid].size(); k++)
            {
                check_path = check_MPC(adj_cc[cid], path_cover[cid][k], k, cid);
                /* If MPC is correct then only proceed */
                if (check_path == false)
                {
                    std::cerr << " MPC Computation is wrong! " << std::endl;
                    exit(0);
                }
            }
        }
    }
    // Compute [min-max] paths for cids
    int min = path_cover[0].size();
    int max = min;
    for(int cid = 1; cid < num_cid; cid++)
    {
        min = min<path_cover[cid].size()?min:path_cover[cid].size();
        max = max>path_cover[cid].size()?max:path_cover[cid].size();
    }
    fprintf(stderr, "[M::%s] range Haplotypes [%d-%d] \n", __func__,min,max);
}

void graphUtils::MPC_index()
{
    index.resize(num_cid);
    rev_index.resize(num_cid);
    dist2begin.resize(num_cid);
    last2reach.resize(num_cid);
    Distance.resize(num_cid);
    paths.resize(num_cid);
    omp_set_dynamic(1);
    omp_set_num_threads(0);
    #pragma omp parallel for
    for (size_t cid = 0; cid < num_cid; cid++)
    {
        /*
        ######################
        # last2reach & index #
        ######################
        */
        int K = path_cover[cid].size();
        int N = adj_cc[cid].size();

        // Compute index ( arranged in topological order by MPC , no need to sort)
        index[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),-1));
        paths[cid].resize(N);
        rev_index[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),-1));
        dist2begin[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),0));
        for (size_t k = 0; k < path_cover[cid].size(); k++)
        {
            int i = 0;
            int temp = 0;
            for (auto idx : path_cover[cid][k])
            {
                index[cid][k][idx] = i; // assuming topologically sorted in origional component
                rev_index[cid][k][i++] = idx;
                dist2begin[cid][k][idx] = temp;
                temp += (int)node_len[component_idx[cid][idx]];
            }
        }

        // last2reach computation
        last2reach[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),-1)); // Initialise last2reach
        Distance[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),std::numeric_limits<int>::max()/2)); // Initilaise Distance

        for (int k = 0; k < K; k++)
        {
            int i = 0;
            for (auto v:path_cover[cid][k])
            {
                last2reach[cid][k][v] = i++;
            }
        }

        // last2reach computation
        for (int k = 0; k < K; k++)
        {
            for (int v : top_order[cid]) {
                for (size_t u : in_node[cid][v]) {
                    last2reach[cid][k][v] = std::max(last2reach[cid][k][v], last2reach[cid][k][u]);
                }
            }
        }

        // Set all values in path to zero (for DP to be correct)
        for (int k = 0; k < K; k++)
        {
            for (auto v : path_cover[cid][k])
            {
                Distance[cid][k][v] = (int)0;
            }
            
        }

        for (int k = 0; k < K; k++)
        {
            for (int v : top_order[cid]) {
                for (size_t u : in_node[cid][v]) {
                    if (last2reach[cid][k][v] == last2reach[cid][k][u] && last2reach[cid][k][v] != -1)
                    {
                        Distance[cid][k][v] = ((int)node_len[component_idx[cid][u]] + std::min(Distance[cid][k][v], Distance[cid][k][u]));
                    }
                }
            }
        }

        // Correct Distance for which las2reach = -1
        for (int k = 0; k < K; k++)
        {
            for (size_t v = 0; v < N; v++)
            {
                if (last2reach[cid][k][v] == -1)
                {
                    Distance[cid][k][v] = (int)0;
                }
            }
        }
        

        // // Print last2reach 
        // std::cerr << " last2reach " << std::endl;
        // for (size_t i = 0; i < K; i++)
        // {
        //     for (size_t j = 0; j < N; j++)
        //     {
        //         std::cerr << last2reach[cid][i][j] << " " ;
        //     }
        //     std::cerr << std::endl;
        // }
        // std::cerr << " index " << std::endl;
        // for (size_t i = 0; i < K; i++)
        // {
        //     for (size_t j = 0; j < N; j++)
        //     {
        //         std::cerr <<  index[cid][i][j] << " " ;
        //     }
        //     std::cerr << std::endl;
        // }
        // std::cerr << " dist2begin " << std::endl;
        // for (size_t i = 0; i < K; i++)
        // {
        //     for (size_t j = 0; j < N; j++)
        //     {
        //         std::cerr << dist2begin[cid][i][j] << "\t" ;
        //     }
        //     std::cerr << std::endl;
        // }

        // std::cerr << " Distance " << std::endl;
        // for (size_t i = 0; i < 1; i++)
        // {
        //     for (size_t j = 0; j < N; j++)
        //     {
        //         std::cerr << Distance[cid][i][j] << " " ;
        //     }
        //     std::cerr << std::endl;
        // }
    }
}

bool compare_T(const Tuples &a, const Tuples &b){
    return std::tie( a.top_v , a.pos, a.task) < std::tie( b.top_v , b.pos, b.task);
};

std::vector<mg128_t> graphUtils::Chaining(std::vector<mg128_t> anchors, std::string query_name)
{
    if (param_z)
    {
        std::cerr << "Number of Anchors : " << anchors.size() << "\n";
    }

    std::vector<std::set<std::string>> chain_haps(num_cid);
    std::vector<mg128_t> best; // Best Anchors
    std::vector<std::vector<Anchors>> M; // Anchors
    M.resize(num_cid);
    std::vector<std::vector<mg128_t>> idx_Anchor;
    idx_Anchor.resize(num_cid); 
    for (int j = 0; j < anchors.size(); j++) // For each anchor Divide Anchors corresponding to their cids
    {
        Anchors M_; // Anchor
        int minimizer_len = (int)(anchors[j].y>>32&0xff);
        int node = (int)(anchors[j].x>>32);
        M_.v = node;
        M_.c = (int)(anchors[j].y) - minimizer_len + 1;
        M_.d =  (int)(anchors[j].y);
        M_.x =  (int)(anchors[j].x) - minimizer_len + 1;
        M_.y =  (int)(anchors[j].x);
        M_.c_ = (int)(M_.c - G);
        M[component[node]].push_back(M_);
        idx_Anchor[component[node]].push_back(anchors[j]);
    }
    std::vector<std::pair<std::vector<mg128_t>,int>> best_chains; // Best Chains
    best_chains.resize(num_cid);
    std::vector<int> chain_count(num_cid,0);

    if (is_hprc_gfa)
    {
        if (benchmark) tau_1 = 0.99999f;
        // chain vtxs
        std::vector<std::vector<std::vector<int32_t>>> chains(num_cid);

        // Recombinations count
        std::vector<int> min_loc(num_cid, std::numeric_limits<int>::max()), max_loc(num_cid, std::numeric_limits<int>::min());
        count++; // Read count
        std::vector<std::string> haps_seq_cid;
        haps_seq_cid.resize(num_cid);
        std::vector<float> acc_cid(num_cid);
        std::vector<float> precision_vec(num_cid);
        std::vector<float> recall_vec(num_cid);

        bool is_haplotype = false;
        for (int cid = 0; cid < num_cid; cid++)
        {
            int N = M[cid].size(); // #Anchors
            int K = path_cover[cid].size(); // #Paths
            
            // minichain
            if (param_z) std::cerr << "Query length : " << seq_len << std::endl;
            int sf = scale_factor;    
            if (N == 0) continue;

            /* Initialise T */
            std::vector<Tuples> T; // Tuples of Anchors
            int cost = (M[cid][0].d - M[cid][0].c + 1) * scale_factor;

            // RMQ (AVL Tree)
            std::pair<int, int> defaul_value = {std::numeric_limits<int>::min(), -1};
            std::vector<IndexX> I(K, IndexX());

            // Initialise C as dynamic array
            Score** C = new Score*[N];
            Score* C_m = new Score[N];
            Score* C_p = new Score[N];
            for (int i = 0; i < N; i++)
            {
                C[i] = new Score[K];
            }

            // Initialise C, C_m and C_p
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    C[i][j].s = (int64_t)cost;
                    C[i][j].i = -1;
                    C[i][j].j = -1;
                }
                C_m[i].s = (int64_t)cost;
                C_m[i].i = -1;
                C_m[i].j = -1;
                C_p[i].s = (int64_t)cost;
                C_p[i].i = -1;
                C_p[i].j = -1;
            }
            
            
            for (int j = 0; j < N; j++) // Add Tuple of Anchors
            {
                int node = M[cid][j].v;
                int v = idx_component[cid][node];
                for (int k = 0; k < K; k++)
                {
                    if (index[cid][k][v] != -1) // v is on the path "k"
                    {
                        Tuples t;
                        // for task 0
                        t.anchor = j; // anchor id
                        t.path = k; // path id
                        t.pos = M[cid][j].x; // position of the anchor on a graph
                        t.task = 0; // task id
                        t.d = M[cid][j].d; // position of the anchor on a query
                        t.v = v; // sorting purpose
                        t.w = v; // vertex in which the anchor lies.
                        t.top_v = map_top_sort[cid][v]; // Topological sorting of the vertex
                        T.push_back(t); // Add Tuple

                        // for task 1
                        t.anchor = j; // anchor id
                        t.path = k; // path id
                        t.pos = M[cid][j].x; // position of the anchor on a graph
                        t.task = 1; // task id
                        t.d = M[cid][j].d; // position of the anchor on a query
                        t.v = v; // sorting purpose
                        t.w = v; // vertex in which the anchor lies.
                        t.top_v = map_top_sort[cid][v]; // Topological sorting of the vertex
                        T.push_back(t); // Add Tuple

                        // for task 2
                        t.anchor = j;
                        t.path = k;
                        t.pos = M[cid][j].y;
                        t.task = 2;
                        t.d = M[cid][j].d;
                        t.v = v; 
                        t.w = v;
                        t.top_v = map_top_sort[cid][v];
                        T.push_back(t);

                    }else if (last2reach[cid][k][v] != -1) // v is not on the path "k" and if last2reach exist
                    {
                        int w = last2reach[cid][k][v]; // vertex -> Index
                        w = rev_index[cid][k][w]; // index -> vertex
                        Tuples t;
                        // for task 0
                        t.anchor = j;
                        t.path = k;
                        t.pos = std::numeric_limits<int>::max(); // Anchor lies at infinity distance on the graph
                        t.task = 0;
                        t.d = M[cid][j].d;
                        t.v = w;
                        t.w = v;
                        t.top_v = map_top_sort[cid][w]; // Topological Sort
                        T.push_back(t);
                    }   
                }
            }

            // Sort the tuples in ascending order of (rank(v), pos, task)
            std::sort(T.begin(), T.end(), compare_T); // Sort Tuples
            std::pair<int64_t, int> rmq = {std::numeric_limits<int64_t>::min(), -1};
            int64_t wt = (int64_t)cost;

            int64_t infty= std::numeric_limits<int64_t>::min();

            // Backward Pointers
            int64_t max_score = -1;
            int idx = 0;
            int path = 0;

            // Optimized Backtracking
            std::vector<std::pair<std::pair<int64_t, int>, std::pair<int, int>>> D(N, std::pair<std::pair<int64_t, int>, \ 
                        std::pair<int, int>>({std::numeric_limits<int>::min(), -1},{-1, -1})); // (score, index), index

            // Initialize a pointer and a array.
            std::vector<int> x(K,0); // pointer for the path
            std::vector<int> rmq_coor(K,0); // current index of the anchor which lies outside the window of G
            std::vector<std::vector<std::pair<int, std::pair<int, int>>>> path_distance(K); // path_distance[path][l] = (M[cid][t.anchor].y + dist2begin[cid][t.path][t.v] , anchor)
            std::pair<int, int> key;

            if(param_z) std::cerr << "Chaining started for cid : " << cid << "\n";

            int i = 0, j = 0, v = 0, w = 0;
            for (auto t:T) //  Process Tuples in the Lexicographic order of (rank(v), pos, task)
            {
                i = t.anchor; j = t.path; v = t.v; w = t.w;
                if (t.task == 0)
                {
                    int64_t p = M[cid][i].x + dist2begin[cid][j][t.v] + M[cid][i].c + Distance[cid][j][w]- 2;
                    if (index[cid][j][w] != -1) // j \in paths(M[cid][i].v)
                    {
                        // Query after anchor deletion
                        rmq = I[j].RMQ(0, M[cid][i].c - 1);

                        if (rmq.first > infty)
                        {
                            // Dynamic Array
                            int64_t loc_score = rmq.first + wt - p;
                            if (C_p[i].s > C[i][j].s)
                            {
                                C[i][j].s = C_p[i].s;
                                C[i][j].i = C_p[i].i;
                                C[i][j].j = C_p[i].j;
                            }

                            if (loc_score > C[i][j].s)
                            {
                                C[i][j].s = loc_score;
                                C[i][j].i = rmq.second;
                                C[i][j].j = j;
                            }

                            if (C[i][j].s > C_m[i].s)
                            {
                                C_m[i].s = C[i][j].s;
                                C_m[i].i = C[i][j].i;
                                C_m[i].j = C[i][j].j;
                            }
                            
                        }
                    }else // last2reach
                    {
                        rmq = I[j].RMQ(0, M[cid][i].c - 1);
                        if (rmq.first > infty){
                            // Dynamic Array
                            int64_t loc_score = rmq.first + wt - p - (int64_t)recomb;
                            if (loc_score > C_p[i].s)
                            {
                                C_p[i].s = loc_score;
                                C_p[i].i = rmq.second;
                                C_p[i].j = j;
                            }
                        }

                        if (param_z) if (Distance[cid][j][w] > 0) std::cerr << "Distance : " << Distance[cid][j][w] << std::endl;
                    }

                    // Debug
                    if(param_z) std::cerr << "CID: " << cid << " i: " << t.anchor << " j: " << t.path << " v : " << v << " top_v : "<<  t.top_v << " t.pos : " << t.pos << " t.task : " << t.task << " C[i][j]: " \ 
                    << C[t.anchor][t.path].s << " p: " << p << " wt: " << wt << " recomb: " << recomb << " Distance: " << Distance[cid][j][w] << " dist2begin :" << dist2begin[cid][j][t.v] <<  "\n";

                } else if (t.task == 1) // Update scores with maximum score with recombination
                {
                    // Dynamic Array
                    int64_t loc_score = C_m[i].s - (int64_t)recomb;
                    if (loc_score > C[i][j].s)
                    {
                        C[i][j].s = loc_score;
                        C[i][j].i = C_m[i].i;
                        C[i][j].j = C_m[i].j;
                    }

                } else
                {
                    int64_t q = C[i][j].s + (int)M[cid][i].y + dist2begin[cid][j][v] + (int)M[cid][i].d; // krsna
                    I[j].add(M[cid][i].d, {q, i});
                }
            }

            if(param_z) std::cerr << "Chaining done for cid : " << cid << "\n";

            // delete dynamic array
            delete [] C_m;
            delete [] C_p;
            // Backtracking for read mapping and graph generation
            if (N!=0)
            {
                // Fill D Array
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < K; j++)
                    {
                        D.push_back({{C[i][j].s, C[i][j].i}, {i, j}}); // (score, index), index
                    }
                }
                std::vector<int> temp_chain; // temporary chain
                std::vector<int> chain; // final chain
                std::vector<bool> visited(N, false); // visited array

                // // Sort D by a pair (score, index)
                std::sort(D.begin(), D.end(), [](const std::pair<std::pair<int64_t, int>, std::pair<int, int>> &a, const std::pair<std::pair<int64_t, int>, std::pair<int, int>> &b) -> bool
                {
                    return std::tie(a.first.first, a.first.second, a.second.first, a.second.second) > std::tie(b.first.first, b.first.second, a.second.first, a.second.second); // sort by score and index
                });

                max_score = D[0].first.first;
                idx = D[0].second.first;
                path = D[0].second.second;

                int prev_idx = 0;

                int64_t score = max_score;
                int64_t min_score = scale_factor * (float)(M[cid][0].d - M[cid][0].c + 1); // minimum score for a contig
                int64_t threshold_score = (float)tau_1 * score;

                // print score and threshold_score
                if(param_z) std::cerr << "CID: " << cid << " score: " << score << " threshold_score: " << threshold_score << "\n";

                int64_t best_cid_score = score; // best score for a contig
                std::pair<std::vector<mg128_t>, int64_t> chain_pair; // (chain, score)

                bool flag; // chain is disjoint

                int count_chains = 0; // count of anchors in chain
                int count_anchors = 0; // count of anchors in chain
                std::pair<std::pair<int64_t, int>, int> anchor_id;

                if (param_z) std::cerr << "Backtracking started for cid : " << cid << "\n";
                accuracy = -1.0f;
                float sum_acc = 0.0f;
                float precision_ = 0.0f;
                float recall_ = 0.0f;
                int count_acc = 0;
                while (max_score >= threshold_score && count_anchors < N && max_score >  min_score)
                {
                    int loc_min = std::numeric_limits<int>::max();
                    int loc_max = std::numeric_limits<int>::min();
                    count_chains++;
                    int recombination = 0;
                    // std::string prev_hap = walk_map[path];
                    int prev_hap = path;
                    int prev_i = idx;
                    std::vector<std::string> chained_haps;
                    for (anchor_id = {{C[idx][path].s, idx}, path}; anchor_id.first.second != -1 && \
                        anchor_id.second != -1; anchor_id = {{C[anchor_id.first.second][anchor_id.second].s, \
                        C[anchor_id.first.second][anchor_id.second].i}, C[anchor_id.first.second][anchor_id.second].j}) // backtracking
                    {
                        flag = true; // chain is disjoint
                        count_anchors++; // count anchors
                        // count Recombination
                        if (prev_hap != anchor_id.second)
                        {
                            recombination++;
                            prev_hap = anchor_id.second;
                        }
                        chained_haps.push_back(walk_map[anchor_id.second]); // push haplotype to chained_haps
                        if (visited[anchor_id.first.second] == true)
                        {
                            recombination = -1; // don't use recombinations for which there is no chain 
                            flag = false; // chain is not disjoint
                            temp_chain.clear(); // clear the chains which are not disjoint
                            break;
                        }else
                        {
                            temp_chain.push_back(anchor_id.first.second); // push anchor to temp_chain
                            visited[anchor_id.first.second] = true; // mark anchor as visited
                            chain_haps[cid].insert(walk_map[anchor_id.second]); // insert haplotype to chain_haps
                        }

                        if (param_z) // debug
                        {
                            std::cerr << " cid  : " << cid << " idx : " << anchor_id.first.second << " Hap : " << walk_map[anchor_id.second] <<  " parent : " << C[anchor_id.first.second][anchor_id.second].i << "," << C[anchor_id.first.second][anchor_id.second].j <<  " node : " \ 
                            << M[cid][anchor_id.first.second].v << " C[i][j] : " << C[anchor_id.first.second][anchor_id.second].s << " M.y : " << M[cid][anchor_id.first.second].y  <<  \
                            " M.d : " << M[cid][anchor_id.first.second].d << " size of chain :" << temp_chain.size() << " score : " << max_score << " threshold : " \
                            << threshold_score << " Recombination : " << recombination <<  " Chain Count : " << count_chains << "\n";
                        }
                    }

                    if (benchmark) chains[cid].push_back(temp_chain); // push anchor to chains

                    // Store the anchors if chain is disjoint and clear the keys
                    if (flag == true && temp_chain.size() > 0 ) // minimap2 min_cnt
                    {
                        for (int i = 0; i < temp_chain.size(); i++)
                        {
                            chain.push_back(temp_chain[i]); // push anchor to chain
                        }
                        temp_chain.clear(); // clear the temp_chain

                        // Here we will compute the chained haps

                        // compute precision and recall for recombinations
                        int true_rec = 0;
                        int chain_rec = 0;
                        std::vector<std::string> true_haplotypes;
                        chain_rec = max_sum;
                        // Split the string into vector of strings by !
                        std::vector<std::string> qname_;
                        std::string delimiter = "!";
                        size_t pos = 0;
                        std::string token;
                        while ((pos = query_name.find(delimiter)) != std::string::npos) {
                            token = query_name.substr(0, pos);
                            qname_.push_back(token);
                            query_name.erase(0, pos + delimiter.length());
                        }
                        qname_.push_back(query_name);

                        // // print qname_
                        // for (auto q:qname_)
                        // {
                        //     std::cerr << "Qname : " << q << "\n";
                        // }
                        // exit(0);
                        if (qname_.size() == 3)
                        {
                            is_haplotype = true;
                            float f1_score = 0.0f;
                            // Split the string into vector of strings by >
                            std::string delimiter2 = ">";
                            size_t pos2 = 0;
                            std::string token2;
                            while ((pos2 = qname_[2].find(delimiter2)) != std::string::npos) {
                                token2 = qname_[2].substr(0, pos2);
                                true_haplotypes.push_back(token2);
                                qname_[2].erase(0, pos2 + delimiter2.length());
                            }
                            true_haplotypes.push_back(qname_[2]);

                            // remove NULL from true_haplotypes
                            true_haplotypes.erase(std::remove(true_haplotypes.begin(), true_haplotypes.end(), ""), true_haplotypes.end());

                            true_rec = true_haplotypes.size() - 1;
                            assert(true_rec == atoi(qname_[1].c_str()));

                            if (recomb == 0 || true)
                            {
                                chained_haps.clear();
                                std::vector<mg128_t> chained_anchors;
                                for (auto idx:chain)
                                {
                                    chained_anchors.push_back(idx_Anchor[cid][idx]);
                                }
                                // reverse the chain
                                std::reverse(chained_anchors.begin(), chained_anchors.end());

                                std::set<int32_t> set_chained_vtxs;
                                std::vector<int32_t> chained_vtx;
                                for (auto idx:chained_anchors)
                                {
                                    int32_t v = idx.x>>32;
                                    set_chained_vtxs.insert(v);
                                }
                                // fill the chained_vtx
                                for (auto v:set_chained_vtxs){ chained_vtx.push_back(v); }
                                set_chained_vtxs.clear();

                                // sort chained vtx by topological order
                                std::sort(chained_vtx.begin(), chained_vtx.end(), [&](int32_t a, int32_t b) -> bool \
                                        { return map_top_sort[cid][idx_component[cid][a]] < map_top_sort[cid][idx_component[cid][b]]; });
                                
                                // Do DP to find minimum recombination
                                int min_recomb = std::numeric_limits<int>::max();
                                std::map<std::string, std::vector<std::pair<int, std::pair<std::string, int>>>> C_;
                                // Inialize C as \inf
                                for (auto hap:walk_ids)
                                {
                                    for (int i = 0; i < chained_vtx.size(); i++)
                                    {
                                        C_[hap].push_back({std::numeric_limits<int>::max(), {"-1", -1}});
                                    }
                                }

                                // Inialize C
                                int node = chained_vtx[0];
                                for (auto hap:haps[node])
                                {
                                    C_[hap][0] = {0, {"-1", -1}};
                                }

                                // Do DP
                                for (int i = 1; i < chained_vtx.size(); i++)
                                {
                                    int node_1 = chained_vtx[i];
                                    for (auto hap:haps[node_1])
                                    {
                                        int node_2 = chained_vtx[i-1];
                                        for (auto hap2:haps[node_2])
                                        {
                                            if (hap == hap2)
                                            {
                                                C_[hap][i] = std::min(C_[hap][i], {C_[hap2][i-1].first, {hap2, i-1}});
                                            }else
                                            {
                                                C_[hap][i] = std::min(C_[hap][i], {C_[hap2][i-1].first + 1, {hap2, i-1}});
                                            }
                                        }
                                    }
                                }

                                std::pair<int, std::pair<std::string, int>> min_dp_recomb = {std::numeric_limits<int>::max(), {"-1", -1}};
                                // Find the mimum recombination
                                for (auto hap:walk_ids)
                                {
                                    if (min_dp_recomb > C_[hap][chained_vtx.size() - 1])
                                    {
                                        min_dp_recomb = C_[hap][chained_vtx.size() - 1];
                                    }
                                }

                                // recombination = min_dp_recomb.first;
                                if (min_loc[cid] == std::numeric_limits<int>::max())
                                {
                                    min_loc[cid] = min_dp_recomb.first;
                                }else
                                {
                                    min_loc[cid] += min_dp_recomb.first;
                                }
                                

                                // std::vector<std::string> chained_haps;  
                                // Backtrack to find the haplotypes
                                for (auto parent = min_dp_recomb.second; parent.first != "-1" && parent.second != -1; parent = C_[parent.first][parent.second].second)
                                {
                                    chained_haps.push_back(parent.first);
                                }


                                // std::vector<std::string> chained_haps;
                                std::reverse(chained_haps.begin(), chained_haps.end());
                                // find unique haplotypes preseving the order
                                std::vector<std::string> unique_haps;
                                for (auto hap:chained_haps)
                                {
                                    if (std::find(unique_haps.begin(), unique_haps.end(), hap) == unique_haps.end())
                                    {
                                        unique_haps.push_back(hap);
                                    }
                                }

                                // push unique haplotypes to chained_haps
                                chained_haps = unique_haps; 
                                unique_haps.clear();

                                // for (auto hap:chained_haps)
                                // {
                                //     hap_seqs += ">" + hap;
                                // }
                            }else
                            {
                                // std::vector<std::string> chained_haps;
                                std::reverse(chained_haps.begin(), chained_haps.end());
                                // find unique haplotypes preseving the order
                                std::vector<std::string> unique_haps;
                                for (auto hap:chained_haps)
                                {
                                    if (std::find(unique_haps.begin(), unique_haps.end(), hap) == unique_haps.end())
                                    {
                                        unique_haps.push_back(hap);
                                    }
                                }

                                if (min_loc[cid] == std::numeric_limits<int>::max())
                                {
                                    min_loc[cid] = recombination;
                                }else
                                {
                                    min_loc[cid] += recombination;
                                }

                                // push unique haplotypes to chained_haps
                                chained_haps = unique_haps; 
                                unique_haps.clear();

                                // for (auto hap:chained_haps)
                                // {
                                //     hap_seqs += ">" + hap;
                                // }
                                // exit(0);
                            }

                            // make_pair
                            std::set<std::pair<std::string, std::string>> chain_hap_pair;
                            if (chained_haps.size() == 1) // When there is only one haplotype
                            {
                                chain_hap_pair.insert({"S", chained_haps[0]}); 
                                chain_hap_pair.insert({chained_haps[0], "E"});
                            }

                            for (int i = 0; i < chained_haps.size() - 1; i++)
                            {
                                if (i == 0)
                                {
                                    chain_hap_pair.insert({"S", chained_haps[i]});
                                    chain_hap_pair.insert({chained_haps[chained_haps.size() - 1], "E"});
                                }
                                chain_hap_pair.insert({chained_haps[i], chained_haps[i + 1]});
                            }

                            std::set<std::pair<std::string, std::string>> true_hap_pair;
                            if (true_haplotypes.size() == 1) // When there is only one haplotype
                            {
                                true_hap_pair.insert({"S", true_haplotypes[0]});
                                true_hap_pair.insert({true_haplotypes[0], "E"});
                            }
                            
                            for (int i = 0; i < true_haplotypes.size() - 1; i++)
                            {
                                if (i == 0)
                                {
                                    true_hap_pair.insert({"S", true_haplotypes[i]});
                                    true_hap_pair.insert({true_haplotypes[true_haplotypes.size() - 1], "E"});
                                }
                                true_hap_pair.insert({true_haplotypes[i], true_haplotypes[i + 1]});
                            }

                            // Compute precision and recall
                            int true_pos = 0;
                            int false_pos = 0;
                            int false_neg = 0;
                            
                            for (auto hap_pair:chain_hap_pair)
                            {
                                if (true_hap_pair.find(hap_pair) != true_hap_pair.end())
                                {
                                    true_pos++; // chain haplotype is in the true haplotype
                                }else
                                {
                                    false_pos++; // chain haplotype is not in the true haplotype
                                }
                            }

                            for (auto hap_pair:true_hap_pair)
                            {
                                if (chain_hap_pair.find(hap_pair) == chain_hap_pair.end())
                                {
                                    false_neg++; // true haplotype is not in the chain
                                }
                            }


                            precision_ = (float)true_pos/(float)(true_pos + false_pos);
                            recall_ = (float)true_pos/(float)(true_pos + false_neg);

                            // compute F1 
                            float loc_acc = 0.0f;
                            f1_score = 2.0f * ((precision_ * recall_)/(precision_ + recall_));
                            if (precision_ != 0.0f && recall_ != 0.0f)
                            {
                                sum_acc = f1_score;
                                loc_acc = f1_score;
                            }else
                            {
                                sum_acc = 0.0f;
                                loc_acc = 0.0f;
                            }
                            count_acc++;

                            // add chained_haps to hap_seqs with S and E and > sign 
                            haps_seq_cid[cid] += ">S";
                            for (auto hap:true_haplotypes)
                            {
                                haps_seq_cid[cid] += ">" + hap;
                            }

                            haps_seq_cid[cid] += ">E";

                            haps_seq_cid[cid] += "!";

                            haps_seq_cid[cid] += ">S";

                            for (auto hap:chained_haps)
                            {
                                haps_seq_cid[cid] += ">" + hap;
                            }

                            haps_seq_cid[cid] += ">E";
                            haps_seq_cid[cid] += "!";
                            haps_seq_cid[cid] += std::to_string(loc_acc);

                            
                        }
                    }

                    max_score = std::numeric_limits<int>::min(); //{{std::numeric_limits<int>::min(), -1}, -1};
                    for (int i = prev_idx; i < D.size(); i++) // O(N) time maximum score search
                    {
                        if (visited[D[i].second.first] == false) // only look on univisited anchors
                        {
                            max_score = D[i].first.first; // get the next maximum score
                            prev_idx = i; // update the prev_idx
                            idx = D[i].second.first; // update the idx
                            path = D[i].second.second; // update the path
                            break;
                        }
                    }
                    if(param_z) std::cerr << " Second pass for max_score computation finished ... " << std::endl;
                }


                acc_cid[cid] = sum_acc/(float)count_acc;
                precision_vec[cid] = precision_;
                recall_vec[cid] = recall_;
                
                for (int i = 0; i < chain.size(); i++)
                {
                    chain_pair.first.push_back(idx_Anchor[cid][chain[i]]); // push anchor to chain_pair
                }
                // add max score
                chain_pair.second = best_cid_score; // push max score to chain_pair
                best_chains[cid] = chain_pair; //   push chain_pair to best_chains
            }

            for (int i = 0; i < N; i++)
            {
                delete [] C[i];
            }
        }
        int max_f1_cid_ = 0;
        float max_precision = 0.0f;
        float max_recall = 0.0f;
        float max_f1 = std::numeric_limits<float>::min();
        for (int cid = 0; cid < num_cid; cid++)
        {
            if (max_f1 < acc_cid[cid])
            {
                max_f1 = acc_cid[cid];
                max_f1_cid_ = cid;
            }
        }
        accuracy = max_f1;
        precision = precision_vec[max_f1_cid_];
        recall = recall_vec[max_f1_cid_];
        hap_seqs = haps_seq_cid[max_f1_cid_];
        max = std::max(max, min_loc[max_f1_cid_]);
        min = std::min(min, min_loc[max_f1_cid_]);
        max_sum += max;

        if (benchmark && !is_haplotype)
        {
            accuracy = -1.0f;
            max_f1_cid_ = 0;
            int64_t max_score = std::numeric_limits<int64_t>::min();
            for (int cid = 0; cid < num_cid; cid++)
            {
                if (max_score < best_chains[cid].second)
                {
                    max_score = best_chains[cid].second;
                    max_f1_cid_ = cid;
                }
            }
            // find count correct over all cids
            std::set<int32_t> walk_vtx;
            for (auto walk:fwd_walk_map[walk_map[0]])
            {
                walk_vtx.insert(walk);
            }

            for (auto walk:rev_walk_map[walk_map[0]])
            {
                walk_vtx.insert(walk);
            }

            int32_t loc_correct = 0;
            int32_t loc_not_correct = 0;
            float frac_correct_loc = 0.0f;
            std::set<int32_t> chain_vtx;
            for (auto chain_num:chains[max_f1_cid_]) // If any chain is part of CHM13#0 then do count_correct++
            {
                for (auto idx:chain_num)
                {
                    chain_vtx.insert(idx_Anchor[max_f1_cid_][idx].x>>32);
                }

                // find the intersection of chain_vtx and walk_vtx
                std::vector<int32_t> common_vtx;
                // do set_minus
                std::set_difference(chain_vtx.begin(), chain_vtx.end(), walk_vtx.begin(), walk_vtx.end(), std::inserter(common_vtx, common_vtx.begin()));
                // If the intersection is empty then do count_correct++;
                float frac_correct_ = ((float)chain_vtx.size() - (float)common_vtx.size())/(float)chain_vtx.size();
                if (common_vtx.size() == 0)
                {
                    loc_correct = 1;
                }else
                {
                    frac_correct_loc = frac_correct_;
                    loc_not_correct = 1;
                }
            }
            count_correct += loc_correct;
            count_not_correct += loc_not_correct;
            frac_correct += frac_correct_loc;
        }
        
    } else {
        for (int cid = 0; cid < num_cid; cid++)
        {
            int N = M[cid].size(); // #Anchors
            int K = path_cover[cid].size(); // #Paths
            /* Initialise Search Trees */
            /* Initialise T */
            std::vector<Tuples> T; // Tuples of Anchors
            std::vector<std::pair<int, int>> C; // Array of Cost of each Anchor
            C.resize(N);
            int cost =(int) (M[cid][0].d - M[cid][0].c + 1); // Cost of each Anchor
            for (int j = 0; j < N; j++) // Add Tuple of Anchors
            {
                int node = M[cid][j].v;
                int v = idx_component[cid][node];
                for (int k = 0; k < K; k++)
                {
                    if (index[cid][k][v] != -1) // v is on the path "k"
                    {
                        Tuples t;
                        // for task 0
                        t.anchor = j; // anchor id
                        t.path = k; // path id
                        t.pos = M[cid][j].x; // position of the anchor on a graph
                        t.task = 0; // task id
                        t.d = M[cid][j].d; // position of the anchor on a query
                        t.v = v; // sorting purpose
                        t.w = v; // vertex in which the anchor lies.
                        t.top_v = map_top_sort[cid][v]; // Topological sorting of the vertex
                        T.push_back(t); // Add Tuple
                        // for task 1
                        t.anchor = j;
                        t.path = k;
                        t.pos = M[cid][j].y;
                        t.task = 1;
                        t.d = M[cid][j].d;
                        t.v = v; 
                        t.w = v;
                        t.top_v = map_top_sort[cid][v];
                        T.push_back(t);
                    }else if (last2reach[cid][k][v] != -1) // v is not on the path "k" and if last2reach exist
                    {
                        int w = last2reach[cid][k][v]; // vertex -> Index
                        w = rev_index[cid][k][w]; // index -> vertex
                        int len = node_len[component_idx[cid][w]];
                        // if first anchor exceeds distance then remove it
                        int gap = Distance[cid][k][v] + M[cid][j].x - len + 1; // gap from last2reach node to the anchor on different node
                        if (gap <= G)
                        {
                            Tuples t;
                            // for task 0
                            t.anchor = j;
                            t.path = k;
                            t.pos = std::numeric_limits<int>::max(); // Anchor lies at infinity distance on the graph
                            t.task = 0;
                            t.d = M[cid][j].d;
                            t.v = w;
                            t.w = v;
                            t.top_v = map_top_sort[cid][w]; // Topological Sort
                            T.push_back(t);
                        }
                    }   
                }
                /* Initialise C */
                C[j] = {cost , -1};
            }

            // Initialize a pointer and a array.
            std::vector<int> x(K,0); // pointer for the path
            std::vector<int> rmq_coor(K,0); // current index of the anchor which lies outside the window
            /* Erase redundant Tuples and Sort the Tuples by T.v, T.pos, T.task */ 
            std::sort(T.begin(),T.end(),compare_T); // Sort Tuples
            std::vector<std::vector<std::pair<int, std::pair<int, int>>>> path_distance(K); // path_distance[path][l] = (M[cid][t.anchor].y + dist2begin[cid][t.path][t.v] , anchor)
            int infty_int = std::numeric_limits<int>::max();
            int _infty_int64 = std::numeric_limits<int>::min();
            std::pair<int, int> rmq; // rmq
            std::pair<int, int> key;

            float c_1 = expf(-div*(float)kmer_len), c_2 = 0.05*c_1; // penalities
            if(!is_ggen) // minimap2 gap cost
            {
                c_1 = 0.01f * kmer_len; // minimap2
                c_2 = 0.0f;
            }

            std::vector<std::map<std::pair<int, int>, std::pair<int, int>>> MAP(K);
            std::map<std::pair<int, int>, std::pair<int, int>>::iterator start_it, end_it;

            for (auto t:T) //  Process Tuples in the Lexicographic order of (rank(v), pos, task)
            {
                if(t.task == 0) // update the score
                {
                    int weight = (M[cid][t.anchor].d - M[cid][t.anchor].c + 1);
                    int range = (dist2begin[cid][t.path][t.v] + Distance[cid][t.path][t.w] + M[cid][t.anchor].x - G - 1);
                    int Rd_1 = (dist2begin[cid][t.path][t.v] + Distance[cid][t.path][t.w] + M[cid][t.anchor].x);
                    int Qd_1 = (M[cid][t.anchor].c);
                    if (index[cid][t.path][t.w] != -1) // same path
                    {
                        for (int l = x[t.path]; l < path_distance[t.path].size(); l++)
                        {
                            if (path_distance[t.path][l].first > range)
                            {
                                rmq_coor[t.path] = l;
                                break;
                            }
                            
                        }
                        // delete anchors from search tree which are not in a range of G
                        for (int l = x[t.path]; l < rmq_coor[t.path]; l++)
                        {
                            key = path_distance[t.path][l].second;
                            MAP[t.path].erase(key);
                        }

                        // Update the pointer
                        x[t.path] = rmq_coor[t.path];   
                    }
                    // Minigraph gap cost
                    start_it = MAP[t.path].lower_bound({M[cid][t.anchor].c_, infty_int});
                    end_it = MAP[t.path].upper_bound({M[cid][t.anchor].d - 1, infty_int}); // overlap 
                    // // Maximize the score function from window query
                    std::pair<int, int> max_cost = std::make_pair(_infty_int64, -1);
                    int h = 0; 
                    for (auto i = std::reverse_iterator(end_it); i != std::reverse_iterator(start_it) && h < max_itr; i++)
                    {
                        int idx = i->second.first; // anchor_idx
                        int v = i->second.second; // vertex
                        int Rd_2 = dist2begin[cid][t.path][v] + M[cid][idx].y + 1; // count of bases hence +1
                        int Qd_2 = M[cid][idx].d + 1;
                        if (Rd_2 - 1 > range)
                        {
                            int g = abs((Rd_1 - Rd_2) - (Qd_1 - Qd_2));
                            int d = std::min((Rd_1 - Rd_2), (Qd_1 - Qd_2));
                            float log_pen = g >= 2 ? mg_log2(g) : 0.0f; // mg_log2() only works for dd>=2
                            int gap = (c_1*(float)g + c_2*(float)d + log_pen);
                            weight = std::min(d, kmer_len);
                            max_cost = max_cost.first > (C[idx].first - gap) ? max_cost : std::make_pair(C[idx].first - gap, idx); 
                            C[t.anchor] = std::max(C[t.anchor], { weight + max_cost.first, max_cost.second});
                            h++;
                        }
                    }

                }else // update the priority of the anchor
                {
                    MAP[t.path][{M[cid][t.anchor].d, t.anchor}] = {t.anchor, t.v};
                    path_distance[t.path].push_back({(int)(M[cid][t.anchor].y + dist2begin[cid][t.path][t.v]), {M[cid][t.anchor].d, t.anchor}}); // (value, key) pair
                }
            }

            // Backtracking for read mapping and graph generation
            if (N!=0)
            {
                int max_score = 0; // maximum score
                std::vector<mg128_t> temp_chain; // temporary chain
                std::vector<mg128_t> chain; // final chain
                std::vector<bool> visited(N, false); // visited array
                // Create an array D and copy C to D
                std::vector<std::pair<std::pair<int, int>, int>> D; // (score, index), index
                for (size_t i = 0; i < C.size(); i++)
                {
                    D.push_back({C[i], i}); // (score, index), index
                }
                // Sort D by a pair (score, index)
                std::sort(D.begin(), D.end(), [](const std::pair<std::pair<int, int>, int> &a, const std::pair<std::pair<int, int>, int> &b) -> bool
                {
                    return std::tie(a.first.first, a.first.second) > std::tie(b.first.first, b.first.second); // sort by score and index
                });
                // Find the minimum score
                int min_score =  (M[cid][0].d - M[cid][0].c + 1); // minimum score for a contig
                // Find the maximum score
                std::pair<int, int> chain_rmq = D[0].first; // Max value of (score, index)
                int prev_idx = 0; // index of the maximum score
                max_score = chain_rmq.first; // maximum score
                int threshold_score = is_ggen > 0 ? 1000 : tau_1*(float)max_score; // threshold score
                int best_cid_score = max_score; // best score for a contig
                std::pair<std::vector<mg128_t>, int> chain_pair; // (chain, score)
                bool flag; // flag for disjoint chain
                int count_anchors = 0; // count anchors
                while (max_score >= threshold_score && count_anchors < N && max_score >  min_score) // at most N times
                {
                    for (int anchor_id = D[prev_idx].second; anchor_id != -1 ; anchor_id = C[anchor_id].second) // backtracking
                    {
                        flag = true; // chain is disjoint
                        count_anchors++; // count anchors keeps count of visted anchors
                        if (visited[anchor_id] == true)
                        {
                            flag = false; // chain is not disjoint
                            temp_chain.clear(); // clear the chains which are not disjoint
                            break;
                        }else
                        {
                            temp_chain.push_back(idx_Anchor[cid][anchor_id]); // push anchor to temp_chain
                            visited[anchor_id] = true; // mark anchor as visited
                        }

                        if (param_z) // debug
                        {
                            std::cerr << " cid  : " << cid << " idx : " << anchor_id <<  " parent : " << C[anchor_id].second  <<  " node : " << M[cid][anchor_id].v << " C[j] : " << C[anchor_id].first << " M.y : " << M[cid][anchor_id].y  <<  " M.d : " << M[cid][anchor_id].d << " size of chain :" << temp_chain.size() << " score : " << max_score << " threshold : " << threshold_score << "\n";
                        }
                    }
                    // Store the anchors if chain is disjoint and clear the keys
                    if (flag == true && temp_chain.size()) // minimap2 min_cnt
                    {
                        for (int i = 0; i < temp_chain.size(); i++)
                        {
                            chain.push_back(temp_chain[i]); // push anchor to chain
                        }
                        chain_count[cid]++;
                        temp_chain.clear(); // clear the temp_chain
                    }

                    for (int i = prev_idx; i < D.size(); i++) // O(N) time maximum score search
                    {
                        if (visited[D[i].second] == false)
                        {
                            chain_rmq = D[i].first; // get the next maximum score
                            prev_idx = i; // update the prev_idx
                            break;
                        }
                    }
                    max_score = chain_rmq.first; // update the max_score
                }
                // push chain to chain_pair
                for (int i = 0; i < chain.size(); i++)
                {
                    chain_pair.first.push_back(chain[i]); // push anchor to chain_pair
                }
                // add max score
                chain_pair.second = best_cid_score; // push max score to chain_pair
                best_chains[cid] = chain_pair; //   push chain_pair to best_chains
            }
        }
    }
    // Out of all "cids" pick one which has maximum score
    int max_score_ = best_chains[0].second; // Initializing max_score as the score of the first cid
    int best_cid = 0; // Initializing best_cid as the first cid
    for (int cid = 0; cid < num_cid; cid++)
    {
        // Here we are checking if the score of the current cid is greater than the previously stored max_score
        if (max_score_ <= best_chains[cid].second)
        {
            // If the score of the current cid is greater than the previously stored max_score, update max_score
            max_score_ = best_chains[cid].second;
            best_cid = cid; // Update best_cid
        }
    }
    if(is_ggen){ // For graph generation
        // std::cerr << " Graph Generation! " << std::endl;
        // exit(0);
        best = best_chains[best_cid].first; // best is the chain of the best_cid
        std::reverse(best.begin(),best.end());
    }else
    {
        // Here we are using tau_2 as a threshold value and multiplying it with max_score_
        int threshold = tau_2 *(float)max_score_;

        // Now pick all the chains from all the cid which satisfies threshold
        std::vector<int> best_cids;
        for (int cid = 0; cid < num_cid; cid++)
        {
            //Here we are checking if the score of the chains is greater than or equal to threshold
            if (best_chains[cid].second >= threshold)
            {
                // If the score of the chains is greater than or equal to threshold then push the cid in best_cids
                best_cids.push_back(cid);
                if (param_z)
                {
                    std::cerr << " Best cid : " << cid << " Chain count : " << chain_count[cid] << "\n"; 
                }
            }
        }
        // Now pick all the chains from best_cids
        for (auto cid:best_cids)
        {
            // Here we are iterating through all the chains of the cid
            for (int i = 0; i < best_chains[cid].first.size(); i++)
            {
                // Now we are pushing the chains in best vector
                best.push_back(best_chains[cid].first[i]);
            }
        }
        // Reverse the order of the elements in best
        std::reverse(best.begin(), best.end());
        std::vector<int> red_idx;
        std::vector<int> count(node_len.size(), 0); // Initialize the count array with all 0s
        /* Count the frequency of each node */
        for (int i = 0; i < best.size(); i++) {
            int node = (int)(best[i].x >> 32);
            count[node]++;
        }
        /* Find the indices of anchors with frequency <= 5 */
        if(!is_hprc_gfa && !is_hap)
        {
            for (int i = 0; i < best.size(); i++) {
                int node = best[i].x >> 32;
                if (count[node] <= 5) {
                    red_idx.push_back(i);
                }
            }
            if(param_z) std::cerr << "Not HPRC GFA!" << std::endl;
        }else
        {
            for (int i = 0; i < best.size(); i++) {
                int node = (int)(best[i].x >> 32);
                if (count[node] < 1) { // pass everything
                    red_idx.push_back(i);
                }
            }
            if(param_z) std::cerr << "HPRC GFA!" << std::endl;
        }
        
        /* Remove anchors from collected indices */
        for (int i = red_idx.size() - 1; i >= 0; i--) {
            best.erase(best.begin()+red_idx[i]);
        }
    }
    if (param_z)
    {
        std::cerr << "Number of Best Anchors : " << best.size() << "\n";
    }
   return best; // return union of all the disjoint set of chains

}
