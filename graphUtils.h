#ifndef __graphUtils__
#define __graphUtils__
#include "gfa.h"
#include "gfa-priv.h"
#include "kalloc.h"
#include "ksort.h"
#include "kvec.h"
#include <iostream>
#include <vector>
#include <stack>
#include <map>
#include <queue>
#include <algorithm>
#include <limits>
#include "minigraph.h"
#include "mgpriv.h"
#include <omp.h>
#include <math.h>


// Anchors
struct Anchors {
	int v;
	int x;
	int y;
	int c;
	int d;
    int c_;
};

struct Tuples {
	int v;
	int w;
	int pos;
	int path;
	int anchor;
	int task;
	int top_v;
	int d;
};


// for flow graph
struct flowGraph {
	int N, M, S, T;
	std::vector<int> f, p, t, c;

	flowGraph(int NN) : N(NN+2) {
		init(N);
		S = NN;
		T = NN + 1;
	}

	void init(int N) {
		f.clear();
		f.resize(N, 0);
		t.clear();
		t.resize(2);
		p = t;
		c = t;
	}

	void add_edge(int i, int j, int cap) {
		// std::cerr << "add edge " << i << " " << j << "  " << cap << std::endl;
		p.push_back(j);
		t.push_back(f[i]);
		c.push_back(cap);
		f[i] = t.size() - 1;
	}
};


class graphUtils
{
	public:
		gfa_t * g;	// This is Graph 
        std::vector<std::vector<int>> adj_;	// This is the adjacency list
        std::vector<std::vector < int>> conn_comp;	// connected components
        std::vector<int> component;	// component id
        int num_comp;	// number of connected components
        int n_vtx;	// number of vertices
        std::vector<std::vector < int>> *adj_cc;	// This is the graph of connected components
        int num_cid;	// number of connected components
        std::vector<std::vector < int>> top_order;	// topological order
        //count of characters in a node 
        std::vector<int> node_len;	// node length
        // Chaining (Revised Algorithm)
        std::vector<std::vector<std::vector< int>>> index;	// index
        std::vector<std::vector<std::vector< int>>> rev_index;	// rev_index
        std::vector<std::vector<std::vector< int>>> last2reach;	// last2reach
        std::vector<std::vector<std::vector< int64_t>>> dist2begin;	// dist2begin
        std::vector<std::vector<std::vector< int64_t>>> Distance;	// Distance
        std::vector<std::vector < int>> component_idx;	// mapping between origional index and local index
        std::vector<std::vector < int>> idx_component;	// mapping between local index and origional index
        std::vector<std::vector<std::vector< int>>> path_cover;	// Path Cover
        std::vector<std::vector<std::vector< int>>> paths;	// Path Cover
        // in_node and out_node computation 
        std::vector<std::vector<std::vector< int>>> in_node;	// in_node
        std::vector<std::vector<std::vector< int>>> out_node;	// out_node
        /*Map Top_Sort */
        std::vector<std::vector < int>> map_top_sort;
        float scale_factor;
        bool param_z;
        int G;
        int lin_ref = 0;

        // tau_1 : intra cid threshold, and tau_2 : inter cid threshold.
        float tau_1;
        float tau_2;
        bool is_ggen;
        int64_t kmer_len;
        float div;
        int max_itr;

        graphUtils(gfa_t *g);	// This is constructor

        void read_graph(); // This is to read the graph

        void print_graph(); // This is to print the graph

        int is_cyclic(); // This is to check if the graph is cyclic

        void Connected_components(); // This is to find connected components

        void topologicat_sort(); // This is to find topological order

        void MPC(); // This is to find the minimum path cover

        std::vector<std::vector < int>> shrink(int cid); // This is to shrink the graph

        std::vector<std::vector < mg128_t>> get_anchors(); // This is to get anchors from the Minigraph

        void MPC_index(); // This is to compute indexing with the minimum path cover

        std::vector<mg128_t> Chaining(std::vector<mg128_t> anchors); // This is to chain anchors

};


// map-algo.c
void get_Op(graphUtils *graphOp);

#endif