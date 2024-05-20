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
#include <sstream>
#include <fstream>
#include <set>

// Score
struct Score {
	int64_t s;
    int i;
    int j;
};

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
        std::vector<std::vector<std::vector< int>>> dist2begin;	// dist2begin
        std::vector<std::vector<std::vector< int>>> Distance;	// Distance
        std::vector<std::vector < int>> component_idx;	// mapping between origional index and local index
        std::vector<std::vector < int>> idx_component;	// mapping between local index and origional index
        std::vector<std::vector<std::vector< int>>> path_cover;	// Path Cover
        std::vector<std::vector<std::vector< int>>> paths;	// Path Cover
        // in_node and out_node computation 
        std::vector<std::vector<std::vector< int>>> in_node;	// in_node
        std::vector<std::vector<std::vector< int>>> out_node;	// out_node
        /*Map Top_Sort */
        std::vector<std::vector < int>> map_top_sort;
        int scale_factor;
        bool param_z;
        int G;
        int recomb;
        int lin_ref = 0;
        int seq_len = 0;

        // tau_1 : intra cid threshold, and tau_2 : inter cid threshold.
        float tau_1;
        float tau_2;
        bool is_ggen;
        bool is_hap;
        int kmer_len;
        float div;
        int max_itr;
        std::string hap_seqs;
        float precision = 0.0f, recall = 0.0f;

        // for recombinations count
        int min = std::numeric_limits<int>::max(), max = std::numeric_limits<int>::min(), max_sum = 0, count = 0;

        // Read walks from graph
        std::string graph_name;
        bool is_hprc_gfa = false;
        std::vector<std::vector<int>> walks;
        std::vector<std::vector<int>> rev_walks;
        std::vector<std::vector<std::vector<int>>> local_walks;
        std::vector<std::vector<int>> hprc_adj;
        std::vector<std::vector<int>> Linear_Order;
        std::vector<std::vector<int>> Cycles;
        std::map<int, std::string> walk_map;
        std::vector<std::string> walk_ids;
        std::map<int, std::vector<std::string>> haps;
        std::map<std::string, std::vector<int>> rev_walk_map;
        std::map<std::string, std::vector<int>> fwd_walk_map;
        std::vector<int> ref_query;
        float accuracy = 0.0f;
        int num_walks = 0;
        bool benchmark;
        int32_t count_correct = 0;
        int32_t count_not_correct = 0;
        float frac_correct = 0.0f;

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

        std::vector<mg128_t> Chaining(std::vector<mg128_t> anchors, std::string query_name); // This is to chain anchors

};


// map-algo.c
void get_Op(graphUtils *graphOp);

// AVL Tree
template<typename T, typename V>
class AVLTree {
   public:
    AVLTree(const V &default_value = V()) : default_value(default_value) , root_(-1) {}

    void add(T key, V value) { root_ = add(root_, key, value); }

    void update(T key, V value) {
        auto node = find(root_, key);
        if (node != -1) {
            nodes_[node].value = value;
        }
    }

    V RMQ(T key1, T key2) {
        return RMQ(root_, key1, key2);
    }

    V RMQ_1(T key1, T key2) {
        return RMQ_1(root_, key1, key2);
    }

    V RMQ_2(T key1, T key2, int  range) {
        return RMQ_2(root_, key1, key2, range);
    }

    void remove(T key) { root_ = remove(root_, key); }

    V get(T key) {
        int node = find(root_, key);
        if (node != -1) {
            return nodes_[node].value;
        }
        return default_value;
    }

   private:
    struct Node {
        T key;
        V value;
        int height;
        int left;
        int right;

        Node(T k, V v) : key(k), value(v), height(1), left(-1), right(-1) {}
    };
    int root_;
    std::vector<Node> nodes_;
    V default_value;

    int height(int node) { return node != -1 ? nodes_[node].height : 0; }

    int balanceFactor(int node) { return height(nodes_[node].left) - height(nodes_[node].right); }

    void updateHeight(int node) { nodes_[node].height = std::max(height(nodes_[node].left), height(nodes_[node].right)) + 1; }

    int rotateRight(int node) {
        int left = nodes_[node].left;
        nodes_[node].left = nodes_[left].right;
        nodes_[left].right = node;
        updateHeight(node);
        updateHeight(left);
        return left;
    }

    int rotateLeft(int node) {
        int right = nodes_[node].right;
        nodes_[node].right = nodes_[right].left;
        nodes_[right].left = node;
        updateHeight(node);
        updateHeight(right);
        return right;
    }

    int balance(int node) {
        updateHeight(node);
        if (balanceFactor(node) == 2) {
            if (balanceFactor(nodes_[node].left) < 0) {
                nodes_[node].left = rotateLeft(nodes_[node].left);
            }
            return rotateRight(node);
            } else if (balanceFactor(node) == -2) {
            if (balanceFactor(nodes_[node].right) > 0) {
                nodes_[node].right = rotateRight(nodes_[node].right);
            }
            return rotateLeft(node);
        }
        return node;
    }

    int add(int node, T key, V value) {
        if (node == -1) {
            nodes_.emplace_back(key, value);
            return (int)nodes_.size() - 1;
        }
        if (key < nodes_[node].key) {
            nodes_[node].left = add(nodes_[node].left, key, value);
        } else if (key > nodes_[node].key) {
            nodes_[node].right = add(nodes_[node].right, key, value);
        } else {
            nodes_[node].value = std::max(nodes_[node].value, value);
        }
        return balance(node);
    }

    int find(int node, T key) {
        if (node == -1) {
            return -1;
        }
        if (key < nodes_[node].key) {
            return find(nodes_[node].left, key);
        } else if (key > nodes_[node].key) {
            return find(nodes_[node].right, key);
        } else {
            return node;
        }
    }

    int findMin(int node) {
      if (nodes_[node].left == -1) {
          return node;
      }
      return findMin(nodes_[node].left);
    }


    int remove(int node, T key) {
        if (node == -1) {
            return -1;
        }
        if (key < nodes_[node].key) {
            nodes_[node].left = remove(nodes_[node].left, key);
        } else if (key > nodes_[node].key) {
            nodes_[node].right = remove(nodes_[node].right, key);
        } else {
            if (nodes_[node].left == -1 && nodes_[node].right == -1) {
                return -1;
            }
            if (nodes_[node].left == -1) {
                return nodes_[node].right;
            }
            if (nodes_[node].right == -1) {
                return nodes_[node].left;
            }
            int next = findMin(nodes_[node].right);
            nodes_[node].key = nodes_[next].key;
            nodes_[node].value = nodes_[next].value;
            nodes_[node].right = remove(nodes_[node].right, nodes_[next].key);
        }
        return balance(node);
    }

    V RMQ(int node, T key1, T key2) {
        if (node == -1) {
            return default_value;
        }
        if (key1 <= nodes_[node].key && key2 >= nodes_[node].key) {
            V leftMax = RMQ(nodes_[node].left, key1, key2);
            V rightMax = RMQ(nodes_[node].right, key1, key2);
            return std::max(nodes_[node].value, std::max(leftMax, rightMax));
          } else if (key1 > nodes_[node].key) {
                return RMQ(nodes_[node].right, key1, key2);
          } else {
                return RMQ(nodes_[node].left, key1, key2);
          }
    }

    V RMQ_1(int node, T key1, T key2) {
        if (node == -1) {
            return default_value;
        }
        V maxValue = default_value;
        std::stack<int> s;
        s.push(node);
        while (!s.empty()) {
            node = s.top();
            s.pop();
            if (key1 <= nodes_[node].key && key2 >= nodes_[node].key) {
                maxValue = std::max(maxValue, nodes_[node].value);
            }
            if (nodes_[node].left != -1 && key1 <= nodes_[node].key) {
                s.push(nodes_[node].left);
            }
            if (nodes_[node].right != -1 && key2 > nodes_[node].key) {
                s.push(nodes_[node].right);
            }
        }
        return maxValue;
    }

    V RMQ_2(int node, T key1, T key2, int range) {
        if (node == -1) {
            return default_value;
        }
        V maxValue = default_value;
        std::stack<int> s;
        s.push(node);
        while (!s.empty()) {
            node = s.top();
            std::pair<std::pair<int, int>, int> value = nodes_[node].value;
            s.pop();
            if (key1 <= nodes_[node].key && key2 >= nodes_[node].key && value.second > range) {
                maxValue = std::max(maxValue, nodes_[node].value);
            }
            if (nodes_[node].left != -1 && key1 <= nodes_[node].key) {
                s.push(nodes_[node].left);
            }
            if (nodes_[node].right != -1 && key2 > nodes_[node].key) {
                s.push(nodes_[node].right);
            }
        }
        return maxValue;
    }

};

// typedef AVLTree<std::pair<int, int>, std::pair<int, int>> SearchTree_H;
typedef AVLTree<std::pair<int, int>, std::pair<std::pair<int, int>, int>> SearchTree;
typedef AVLTree<std::pair<int, int>, std::pair<int, int>> SearchTree_H;


template<typename T, typename V>
struct Treap {
	struct Node {
		int ls, rs, size, pri;
		T key;
		V value, max;
	};
	std::vector<Node> t;
	int root;
	V default_value;
	Treap(const V &default_value = V()) : default_value(default_value) {
		root = 0;
		t.resize(1);
	}
	inline int randomm() {
		static int seed = 703; 
		return seed = int(seed * 48271LL % 2147483647);
	}
	inline int update(int now) {
		t[now].size = 1;
		t[now].max = t[now].value;
		if (t[now].ls) {
			t[now].size += t[t[now].ls].size;
			t[now].max = max(t[now].max, t[t[now].ls].max);
		}
		if (t[now].rs) {
			t[now].size += t[t[now].rs].size;
			t[now].max = max(t[now].max, t[t[now].rs].max);
		}
		return now;
	}
	inline int new_node (T key, V value) {
		t.push_back(Node({ 0, 0, 1, randomm(), key, value, value }));
		return t.size() - 1;
	}
	int merge(int x, int y) {
		if (!x || !y) return x + y;
		if (t[x].pri > t[y].pri) {
			t[x].rs = merge(t[x].rs, y);
			return update(x);
		}
		else {
			t[y].ls = merge(x, t[y].ls);
			return update(y);
		}
	}
	void split(int now, T key, int &x, int &y) {
		if (!now) {
			x = y = 0; 
			return;
		}
		if (t[now].key <= key) {
			x = now;
			split(t[now].rs, key, t[now].rs, y);
			update(x);
		}
		else {
			y = now;
			split(t[now].ls, key, x, t[now].ls);
			update(y);
		}
	}
	// void Del(int &root, int key) {
	// 	int x = 0, y = 0, z = 0;
	// 	split(root, key, x, z);
	// 	split(x, key - 1, x, y);
	// 	y = merge(t[y].ls, t[y].rs);
	// 	root = merge(merge(x, y), z);
	// }
	void add(T key, V value) {
		int x = 0, y = 0, z = 0;
		split(root, key, x, y);
		root = merge(merge(x, new_node(key, value)), y);
	}	
	V RMQ(T l, T r) {
		int now = root;
		while (now != 0 && (t[now].key < l || t[now].key > r)) {
			if (t[now].key < l)
				now = t[now].rs;
			else
				now = t[now].ls;
		}
		if (now == 0) {
			return default_value;
		}
		V ret = t[now].value;
		int x = t[now].ls;
		while (x != 0) {
			if (t[x].key >= l) {
				ret = max(ret, t[x].value);
				if (t[x].rs != 0)
					ret = max(ret, t[t[x].rs].max);
				x = t[x].ls;
			}
			else
				x = t[x].rs;
		}
		int y = t[now].rs;
		while (y != 0) {
			if (t[y].key <= r) {
				ret = max(ret, t[y].value);
				if (t[y].ls != 0)
					ret = max(ret, t[t[y].ls].max);
				y = t[y].rs;
			}
			else
				y = t[y].ls;
		}
		return ret;
	}
};

typedef Treap<std::pair<int, int>, std::pair<std::pair<int, int>, int>> IndexT;

typedef Treap<int, std::pair<int64_t, int>> IndexX;
#endif