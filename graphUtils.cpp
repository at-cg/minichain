#include "graphUtils.h"
#include <assert.h>

graphUtils::graphUtils(gfa_t *g)
{
    this->g = g;
}

// Read and store the Graph in Adjacency list
void graphUtils::read_graph()
{
    uint32_t v;
    n_vtx = gfa_n_vtx(g);
    // Resize node_len
    node_len.resize(n_vtx, 0);
    // Array of vectors
    adj_ = new std::vector<int>[n_vtx];
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
    u_int32_t num_edges = 0;
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
                uint32_t w = av[i].w;
                adj_[v_].push_back(w);
            }
        }
    }
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
    // In degree and Outdegree computation
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
    // std::cout << " Computed Topological Order " << std::endl;
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
    /*
    ######################
    # MPC (greedy cover) #
    ######################
    */
    path_cover.resize(num_cid);
    in_node.resize(num_cid);
    out_node.resize(num_cid);
    omp_set_dynamic(1);
    omp_set_num_threads(0);
    #pragma omp parallel for
    for (size_t cid = 0; cid < num_cid; cid++)
   {

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
        if (param_z)
        {
            /* Check MPC */
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
  fprintf(stderr, "[M::%s] range wcc [%d-%d] \n", __func__,min,max);
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
        dist2begin[cid].resize(path_cover[cid].size(),std::vector<int64_t>(adj_cc[cid].size(),0));
        for (size_t k = 0; k < path_cover[cid].size(); k++)
        {
            int i = 0;
            int64_t temp = 0;
            for (auto idx : path_cover[cid][k])
            {
                index[cid][k][idx] = i; // assuming topologically sorted in origional component
                rev_index[cid][k][i++] = idx;
                dist2begin[cid][k][idx] = temp;
                temp += (int64_t)node_len[component_idx[cid][idx]];
            }
        }

        // last2reach computation
        last2reach[cid].resize(path_cover[cid].size(),std::vector<int>(adj_cc[cid].size(),-1)); // Initialise last2reach
        Distance[cid].resize(path_cover[cid].size(),std::vector<int64_t>(adj_cc[cid].size(),std::numeric_limits<int64_t>::max()/2)); // Initilaise Distance

        for (int k = 0; k < K; k++)
        {
            int i = 0;
            for (auto v:path_cover[cid][k])
            {
                last2reach[cid][k][v] = i++;
            }
        }

        // Topological order (kahn's Algorithm)
        std::vector<int> incd(N, 0), Q;
        for (int i = 0; i < N; i++)
        {
            incd[i] = in_node[cid][i].size();
            if (incd[i] == 0)
            {
                Q.push_back(i);
            }
        }
        for (int i = 0; i < Q.size(); ) {
        int s = Q[i++];
        for (size_t t : out_node[cid][s]) {
            incd[t]--;
            if (incd[t] == 0)
                Q.push_back(t);
            }
        }

        // last2reach computation
        for (int k = 0; k < K; k++)
        {
            for (int v : Q) {
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
                Distance[cid][k][v] = (int64_t)0;
            }
            
        }

        for (int k = 0; k < K; k++)
        {
            for (int v : Q) {
                for (size_t u : in_node[cid][v]) {
                    if (last2reach[cid][k][v] == last2reach[cid][k][u] && last2reach[cid][k][v] != -1)
                    {
                        Distance[cid][k][v] = ((int64_t)node_len[component_idx[cid][u]] + std::min(Distance[cid][k][v], Distance[cid][k][u]));
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
                    Distance[cid][k][v] = (int64_t)0;
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
        // for (size_t i = 0; i < K; i++)
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

std::vector<mg128_t> graphUtils::Chaining(std::vector<mg128_t> anchors)
{
    if (param_z)
    {
        std::cerr << " Number of Anchors : " << anchors.size() << "\n";
    }
    std::vector<mg128_t> best; // Best Anchors
    std::vector<std::vector<Anchors>> M; // Anchors
    M.resize(num_cid);
    std::vector<std::vector<mg128_t>> idx_Anchor;
    idx_Anchor.resize(num_cid); 
    for (int j = 0; j < anchors.size(); j++) // For each anchor Divide Anchors corresponding to their cids
    {
        Anchors M_; // Anchor
        int minimizer_len = (int32_t)(anchors[j].y>>32&0xff);
        int node = (int)(anchors[j].x>>32);
        M_.v = node;
        M_.c = (int32_t)(anchors[j].y) - minimizer_len + 1;
        M_.d =  (int32_t)(anchors[j].y);
        M_.x =  (int32_t)(anchors[j].x) - minimizer_len + 1;
        M_.y =  (int32_t)(anchors[j].x);
        M_.c_ = (int)(M_.c - G);
        M[component[node]].push_back(M_);
        idx_Anchor[component[node]].push_back(anchors[j]);
    }
    std::vector<std::pair<std::vector<mg128_t>,int64_t>> best_chains; // Best Chains
    best_chains.resize(num_cid);
    std::vector<int> chain_count(num_cid,0);
    for (int cid = 0; cid < num_cid; cid++)
    {
        int N = M[cid].size(); // #Anchors
        int K = path_cover[cid].size(); // #Paths
        /* Initialise Search Trees */
        std::pair<std::pair<int64_t, int> , int64_t> default_value = {{std::numeric_limits<int64_t>::min(), -1}, (int64_t)-1};
        std::vector<SearchTree> I(K, SearchTree(default_value)); // Search Tree
        /* Initialise T */
        std::vector<Tuples> T; // Tuples of Anchors
        std::vector<std::pair<int64_t,int>> C; // Array of Cost of each Anchor
        C.resize(N);
	    int64_t sf = scale_factor; // Scale Factor
	    int64_t cost =(int64_t) (M[cid][0].d - M[cid][0].c + 1)*sf; // Cost of each Anchor
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
                    int64_t gap = Distance[cid][k][v] + M[cid][j].x - len; // gap from last2reach node to the anchor on different node
                    if (gap < G)
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
        std::vector<std::vector<std::pair<int64_t, std::pair<int, int>>>> y_array(K); // y_array[l] = (M[cid][t.anchor].y + dist2begin[cid][t.path][t.v] , anchor)
        int infty_int = std::numeric_limits<int>::max();
        int64_t _infty_int64 = std::numeric_limits<int64_t>::min();
        // std::pair<int64_t,int> rmq; // rmq
        std::pair<std::pair<int64_t, int> , int64_t> rmq; // rmq
        std::pair<int, int> key;
        for (auto t:T) //  Process Tuples in the Lexicographic order of (rank(v), pos, task)
        {
            if(t.task == 0) // update the score
            {
                int64_t val_1 = ( M[cid][t.anchor].c - 1 + M[cid][t.anchor].x - 1 + dist2begin[cid][t.path][t.v] + Distance[cid][t.path][t.w]);
                int64_t val_2 = sf*(M[cid][t.anchor].d - M[cid][t.anchor].c + 1);
                int64_t sum_d_D = (int64_t)(dist2begin[cid][t.path][t.v] + M[cid][t.anchor].x - G - 1);
                for (int l = x[t.path]; l < y_array[t.path].size(); l++)
                {
                    if (y_array[t.path][l].first >= sum_d_D)
                    {
                        rmq_coor[t.path] = l;
                        break;
                    }
                    
                }
                // delete anchors from search tree which are not in a range of G
                for (int l = x[t.path]; l < rmq_coor[t.path]; l++)
                {
                    key = y_array[t.path][l].second;
                    I[t.path].remove(key);
                }
                
                // Extend the chain
                // rmq = I[t.path].RMQ({M[cid][t.anchor].c_, infty_int}, {M[cid][t.anchor].c - 1, infty_int});
                rmq = I[t.path].RMQ_3({M[cid][t.anchor].c_, infty_int}, {M[cid][t.anchor].c - 1, infty_int}, sum_d_D);
                if (rmq.first.first > _infty_int64)
                {
                    C[t.anchor] = std::max(C[t.anchor], { rmq.first.first - val_1 + val_2, rmq.first.second});
                }

                // Update the pointer
                x[t.path] = rmq_coor[t.path];

                if (param_z) // debug
                {
                    std::cerr << " cid  : " << cid << " idx : " << t.anchor << " top_v :" << t.top_v << " pos : " << t.pos << " task : " << t.task << " path : " << t.path <<  " parent : " << C[t.anchor].second  <<  " node : " << M[cid][t.anchor].v << " index : " << index[cid][t.path][t.w] << " C[j] : " << C[t.anchor].first << " update_C : " << (rmq.first.first - val_1 + val_2) << " rmq.first : " << rmq.first.first  << " val_1 : " << val_1 << " dist2begin : " << dist2begin[cid][t.path][t.v] << " Distance : " <<  Distance[cid][t.path][t.w] <<  " M[i].d : " << t.d << "\n"; 
                }

            }else // update the priority of the anchor
            {
                int64_t val_3 = (M[cid][t.anchor].d + M[cid][t.anchor].y + dist2begin[cid][t.path][t.v]); // M[i].d + M[i].y + dist2begin[v]
                I[t.path].add({M[cid][t.anchor].d, t.anchor}, std::make_pair(std::make_pair(C[t.anchor].first + val_3 , t.anchor), (int64_t)(M[cid][t.anchor].y + dist2begin[cid][t.path][t.v]))); // C[j].first + M[i].d + M[i].y + dist2begin[v] (priority, anchor)
                y_array[t.path].push_back({(int64_t)(M[cid][t.anchor].y + dist2begin[cid][t.path][t.v]), {M[cid][t.anchor].d, t.anchor}}); // (value, key) pair
                if (param_z) // debug
                {
                    std::cerr << " cid  : " << cid << " idx : " << t.anchor << " top_v :" << t.top_v << " pos : " << t.pos << " task : " << t.task << " path : " << t.path << " val_3 : " << val_3  << " M.y : " << M[cid][t.anchor].y <<  " dist2begin : "  << dist2begin[cid][t.path][t.v] << " M[i].d : " << t.d << "\n"; 
                }
            }
        }

        // Backtracking for read mapping and graph generation
        if (N!=0)
        {
            int64_t max_score = 0; // maximum score
            std::vector<mg128_t> temp_chain; // temporary chain
            std::vector<mg128_t> chain; // final chain
            std::vector<bool> visited(N, false); // visited array
            // Create an array D and copy C to D
            std::vector<std::pair<std::pair<int64_t, int>, int>> D; // (score, index), index
            for (size_t i = 0; i < C.size(); i++)
            {
                D.push_back({C[i], i}); // (score, index), index
            }
            
            // Sort D by a pair (score, index)
            std::sort(D.begin(), D.end(), [](const std::pair<std::pair<int64_t, int>, int> &a, const std::pair<std::pair<int64_t, int>, int> &b) -> bool
            {
                return std::tie(a.first.first, a.first.second) > std::tie(b.first.first, b.first.second); // sort by score and index
            });

            // // Create one-to-one mapping
            // std::vector<int> mapping(D.size()); // In case if someone wants to experiment with the original index of anchors
            // for (size_t i = 0; i < D.size(); i++)
            // {
            //     mapping[D[i].second] = i;
            // }

            int min_score = scale_factor * (float)(M[cid][0].d - M[cid][0].c + 1); // minimum score for a contig
            // Find the maximum score
            std::pair<int64_t, int> chain_rmq = D[0].first; // Max valur of (score, index)
            int prev_idx = 0; // index of the maximum score
            max_score = chain_rmq.first; // maximum score
            int64_t threshold_score;
            if (is_ggen) // For graph generation
            {
                threshold_score =  min_score; // threshold score for graph generation (at least 2 anchors with gap-cost=0)
            }else // For read mapping
            {
                threshold_score = tau_1 * (float)max_score; // threshold score
            }
            int64_t best_cid_score = max_score; // best score for a contig
            std::pair<std::vector<mg128_t>, int64_t> chain_pair; // (chain, score)
            bool flag; // flag for disjoint chain
            int count_anchors = 0; // count anchors
            // while (max_score >= threshold_score && count_anchors < N && max_score >  min_score ) // at most N times
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
                if (flag == true)
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


    // Out of all "cids" pick one which has maximum score
    int64_t max_score_ = best_chains[0].second; // Initializing max_score as the score of the first cid
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
        best = best_chains[best_cid].first; // best is the chain of the best_cid
        // Reverse the order of the elements in best
        std::reverse(best.begin(),best.end());
    }else
    {
        // Here we are using tau_2 as a threshold value and multiplying it with max_score_
        int64_t threshold = tau_2 *(float)max_score_;

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
        std::reverse(best.begin(),best.end());


        std::vector<int> red_idx;
        std::vector<int> count(node_len.size(), 0); // Initialize the count array with all 0s

        /* Count the frequency of each node */
        for (int i = 0; i < best.size(); i++) {
            int node = (int)(best[i].x >> 32);
            count[node]++;
        }

        /* Find the indices of anchors with frequency <= 5 */
        for (int i = 0; i < best.size(); i++) {
            int node = (int)(best[i].x >> 32);
            if (count[node] <= 5) {
                red_idx.push_back(i);
            }
        }

        /* Remove anchors from collected indices */
        for (int i = red_idx.size() - 1; i >= 0; i--) {
            best.erase(best.begin()+red_idx[i]);
        }
    }

    if (param_z)
    {
        std::cerr << " Number of Best Anchors : " << best.size() << "\n";
    }

   return best; // return union of all the disjoint set of chains

}
