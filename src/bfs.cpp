#include <queue>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "graph_utils.hpp"

using namespace std;

vector<int> reconstruct_path(int meeting_node, int id_start, int id_end, 
    const vector<int>& parents_start, const vector<int>& parents_end);

/* Implements a bidirectional BFS algorithm to find the shortest path between the starting node and the ending node
'id_start' and 'id_end' are the one-dimensional IDs for the starting and ending lattice coordinates 
'adj' is the adjacency list used to reconstruct the random lattice
'n' is the dimension (or horizontal length + 1) of the lattice
*/
vector<int> bfs(int id_start, int id_end, const vector<vector<int>>& adj, int n) { 
    
    int N = n * n; // number of nodes after 2D lattice has been flattened, i.e. number of one-dimensional ID's
    queue<int> to_visit_start; // queue of nodes to visit from the starting node
    vector<int> parents_start(N, -1); 
    vector<bool> visited_start(N, false); // vector indicating truth value of whether each node has been visited from the starting node

    queue<int> to_visit_end; // queue of nodes to visit from the ending node
    vector<int> parents_end(N, -1);
    vector<bool> visited_end(N, false); // vector indicating truth value of whether each node has been visited from the ending node

    int meeting_node = -1;

    bool found = false; // boolean variable indicating whether the 

    to_visit_start.push(id_start);
    to_visit_end.push(id_end);
    visited_start[id_start] = true;
    visited_end[id_end] = true;
    while (!to_visit_start.empty() and !to_visit_end.empty()) {
        int x = to_visit_start.front();
        to_visit_start.pop();
        for (int i : adj[x]) {
            if (!visited_start[i]) {
                to_visit_start.push(i);
                visited_start[i] = true;
                parents_start[i] = x;
                if (visited_end[i]) {
                    meeting_node = i;
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
        int y = to_visit_end.front();
        to_visit_end.pop();
        for (int i : adj[y]) {
            if (!visited_end[i]) {
                to_visit_end.push(i);
                visited_end[i] = true;
                parents_end[i] = y;
                if (visited_start[i]) {
                    meeting_node = i;
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
    }
    return reconstruct_path(meeting_node, id_start, id_end, parents_start, parents_end);
}

vector<int> reconstruct_path(int meeting_node, int id_start, int id_end, 
    const vector<int>& parents_start, const vector<int>& parents_end) {
    vector<int> path;

    int node = meeting_node;
    while (node != id_start) {
        path.push_back(node);
        node = parents_start[node];
    }
    path.push_back(id_start);
    reverse(path.begin(), path.end());

    node = meeting_node;
    while (node != id_end) {
        node = parents_end[node];
        path.push_back(node);
    }

    return path;
}
