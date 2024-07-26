#pragma once

#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using namespace std;

/// @brief Simple directed graph using an adjacency list.
/// @tparam VertexT vertex type
/// @tparam WeightT edge weight type
template <typename VertexT, typename WeightT>
class graph {
   private:
    // TODO_STUDENT
    size_t numV; // number of vertices in graph
    size_t numE; // number of edges in graph 
    map<VertexT, map<VertexT, WeightT>>adjacencyList; // adjacency list declared as a map within a map, an element would look like {vertex a , {vertex b, distance from a to b}}
    set<VertexT>vertices; // a set of vertices (we use a set to prevent getting duplicate vertices)

   public:
    /// Default constructor
    graph() {
        // TODO_STUDENT
        numV = 0;
        numE = 0;
        adjacencyList = {};
        vertices = {};
    }

    /// @brief Add the vertex `v` to the graph, must run in at most O(log |V|).
    /// @param v
    /// @return true if successfully added; false if it existed already
    bool addVertex(VertexT v) {
        // TODO_STUDENT
        if(vertices.count(v) == 1) { return false; } // if v is already in set of vertices, return false (count() big-O: logN)
        else {
            // else add the vertex
            vertices.emplace(v); // insert() big-O: logN
            numV += 1; // increase count of vertices
            return true; // return true as we added a vertex
        }
        // tot runtime = 3logN = logN
    }

    /// @brief Add or overwrite directed edge in the graph, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight edge weight / label
    /// @return true if successfully added or overwritten;
    ///         false if either vertices isn't in graph
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        // TODO_STUDENT
        if(vertices.count(from) == 0 || vertices.count(to) == 0) { return false; } // if either vertice is not present in the graph, return false (count() big-O: logN)
        else {
            if( adjacencyList[from].count(to) == 0 ) { // count() runtime : logN;
                // if theres no edge present, increase count  as adding a new edge
                numE += 1;
            }
            adjacencyList[from][to] = weight; // if edge didnt exist, add it, if edge existed, overwrite it
            return true; // return true as we set a new weight
        }
        // tot runtime = 2logN = logN
    }

    /// @brief Maybe get the weight associated with a given edge, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight output parameter
    /// @return true if the edge exists, and `weight` is set;
    ///         false if the edge does not exist
    bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
        // TODO_STUDENT
        if(vertices.count(from) == 0 || vertices.count(to) == 0) { return false; } // if either vertice is not present in the graph, return false (count() big-O: logN)
        else {
            if ( adjacencyList.at(from).count(to) == 0 ) { return false; } // if vertex from does not have vertex to in its adjaceny list, return false (count() big-O: logN)
            else {
               weight = adjacencyList.at(from).at(to); // else both vertices exist and we obtain the weight between the edges
               return true; // return true as we got the weight
            }
        }
        // tot runtime = 2logN = logN
    }

    /// @brief Get the out-neighbors of `v`. Mu st run in at most O(|V|).
    /// @param v
    /// @return vertices that v has an edge to
    set<VertexT> neighbors(VertexT v) const {
        set<VertexT> S;
        // TODO_STUDENT
        if(vertices.count(v) == 0) {
            S = {}; // if the vertex v isnt in the set of vertices, it has no neighbors, set the set to an empty set
        }
        else { // else the vertex exists in the vertices set, this doesnt necessarily mean that the vertex is connected to the graph, we could very well have a disconnected graph so we have another chec below
            if (adjacencyList.find(v) != adjacencyList.end()){ // this line checks if vertex is in the adjacency list BECAUSE if the vertex doesnt exist in the adjacency list then it has no neighbors
                for (auto currPair : adjacencyList.at(v)) { // iterate through list of neighbors of vertex v
                    // currPair is a pair of an adjacent vertex and the weighted edge between this vertex and v (currPair.first = neighbor vertex, currPair.second = distance between vertex v and neighbor)
                    S.emplace(currPair.first); // currPair.first is the first piece of the pair which is the neighbor vertex, so we add that to our set 
                }
            }
            else {
                S = {}; // if vertex v wasnt in the adjacency list, it has no neighbors so set the set to an empty set
            }
        }
        return S; // return set
    }

    /// @brief Return a vector containing all vertices in the graph
    vector<VertexT> getVertices() const {
        // TODO_STUDENT
        // return vertices;
        // return vector<VertexT>{}; rechec
        vector<VertexT>out;
        for(auto curr : vertices) { // iterate through the vertices set adding each element to a vector
            out.push_back(curr);
        }
        return out; // return the vector
    }

    /// @brief Get the number of vertices in the graph. Runs in O(1).
    size_t NumVertices() const {
        // TODO_STUDENT
        return numV; // return number of vertices
    }

    /// @brief Get the number of directed edges in the graph. Runs in at most O(|V|).
    size_t NumEdges() const {
        // TODO_STUDENT
        return numE; // return number of edges
    }
};
