#pragma once
#include <vector>

using namespace std;

/** 
struct AdaptedPath {
    vector<int> path;
    int actual_start;
    int actual_end;
};

AdaptedPath bfs_adaptive(int id_start, int id_end, const vector<vector<int>>& adj, int n);


vector<int> neighbors(int id, int n);
 **/

vector<int> bfs(int id_start, int id_end, const vector<vector<int>>& adj, int n);
