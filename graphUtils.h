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
#include <immintrin.h>

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

    V RMQ_2(T key1, T key2, int64_t  sum_d_D) {
        return RMQ_2(root_, key1, key2, sum_d_D);
    }

    V RMQ_3(T key1, T key2, int64_t  sum_d_D) {
        return RMQ_3(root_, key1, key2, sum_d_D);
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

    // int add(int node, T key, V value) {
    //     if (node == -1) {
    //         nodes_.emplace_back(key, value);
    //         return (int)nodes_.size() - 1;
    //     }
    //     int prev = -1, curr = node;
    //     while (curr != -1) {
    //         prev = curr;
    //         if (key < nodes_[curr].key) {
    //             curr = nodes_[curr].left;
    //         } else if (key > nodes_[curr].key) {
    //             curr = nodes_[curr].right;
    //         } else {
    //             nodes_[curr].value = std::max(nodes_[curr].value, value);
    //             return node;
    //         }
    //     }
    //     nodes_.emplace_back(key, value);
    //     int new_node = (int)nodes_.size() - 1;
    //     if (key < nodes_[prev].key) {
    //         nodes_[prev].left = new_node;
    //     } else {
    //         nodes_[prev].right = new_node;
    //     }
    //     node = balance(node);
    //     return node;
    // }



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

    V RMQ_2(int node, T key1, T key2, int64_t sum_d_D ) {
        if (node == -1) {
            return default_value;
        }
        std::pair<std::pair<int64_t, int>, int64_t> value = nodes_[node].value;
        if (key1 <= nodes_[node].key && key2 >= nodes_[node].key && value.second >= sum_d_D) {
            V leftMax = RMQ_2(nodes_[node].left, key1, key2, sum_d_D);
            V rightMax = RMQ_2(nodes_[node].right, key1, key2, sum_d_D);
            return std::max(nodes_[node].value, std::max(leftMax, rightMax));
          } else if (key1 > nodes_[node].key) {
                return RMQ_2(nodes_[node].right, key1, key2, sum_d_D);
          } else {
                return RMQ_2(nodes_[node].left, key1, key2, sum_d_D);
          }
    }

    V RMQ_3(int node, T key1, T key2, int64_t sum_d_D) {
        if (node == -1) {
            return default_value;
        }
        V maxValue = default_value;
        std::stack<int> s;
        s.push(node);
        while (!s.empty()) {
            node = s.top();
            std::pair<std::pair<int64_t, int>, int64_t> value = nodes_[node].value;
            s.pop();
            if (key1 <= nodes_[node].key && key2 >= nodes_[node].key && value.second >= sum_d_D) {
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

typedef AVLTree<std::pair<int, int>, std::pair<std::pair<int64_t, int>, int64_t>> SearchTree;


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
        std::vector<int> *adj_;	// This is the adjacency list
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