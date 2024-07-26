#include "application.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <algorithm>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

double INF = numeric_limits<double>::max();

class prioritize { // this class is declared to help our worklist map order the values in a min queue, where the first element in the queue is the one with the least distance from start
    public:
    bool operator() (const pair<long long, double>& p1, const pair<long long, double>& p2) const {
        return p1.second > p2.second;
    }
};

graph<long long, double> buildGraph(
    const map<long long, Coordinates>& Nodes,
    const vector<FootwayInfo>& Footways,
    const vector<BuildingInfo>& Buildings) {
    
    graph<long long, double> G;

    // TODO_STUDENT
    for(auto n : Nodes) {
        G.addVertex(n.first); // add all vertices to the vertices set , we get 23459 (46 buildings missing, we add them later)
    }

    // now add edges between (non building) vertices using the footways
    size_t i = 0;
    for(i = 0; i < Footways.size(); ++i) { // for each footway in the Footways vector
        size_t j = 0;
        for(j = 0; j < (Footways.at(i)).Nodes.size() - 1; ++j) { // for each node in the current footway
            long long from = (Footways.at(i)).Nodes.at(j); // get the i'th node in the current footway (node from)
            long long to = (Footways.at(i)).Nodes.at(j+1); // get the i+1'th node in the current footway, this node is the immediately next node to i'th node (node to)
            double weight = distBetween2Points((Nodes.at(from)).Lat, (Nodes.at(from)).Lon, (Nodes.at(to)).Lat, (Nodes.at(to)).Lon); // calculate the distance between these two nodes using the given function and passing in their latitudes and longitudes
            G.addEdge(from, to, weight); // add a directed edge from -> to
            G.addEdge(to, from, weight); // add a directed edge to -> from
        }
    }

    for(auto b : Buildings) {
        G.addVertex(b.Coords.ID); // gets the building vertices // 23459 + 46 = 23505
        for (auto n : Nodes) { // now add edges between building vertices & footways
            double dist = distBetween2Points(b.Coords.Lat, b.Coords.Lon, n.second.Lat, n.second.Lon); // calculate the distance between the building vertex and the current node (from our nodes map)
            if (dist <= 0.041 && n.second.OnFootway == true) { // if the distance calculated is within 0.041 AND the node we used to connect is on a footway, then this code bloc executes
                G.addEdge(b.Coords.ID, n.second.ID, dist); // add a directed edge building -> node
                G.addEdge(n.second.ID, b.Coords.ID, dist); // add a directed edge node -> building
            }
        }
    }

    return G; // THEN WE RETURN GRAPH YAYAYAYAY :)
}

vector<long long> dijkstra(
    const graph<long long, double>& G,
    long long start,
    long long target,
    const set<long long>& ignoreNodes) {
    
    vector<long long> path;

    // TODO_STUDENT

    if (start == target) {
        path.push_back(start); // if the start and target are the same then we just add start to our path and return the path
        return path;
    }

    priority_queue<pair<long long, double> , vector<pair<long long, double>> , prioritize> worklist = {}; // worklist to store nodes we want to process
    map<long long, double>distances = {};  // map to store node and its distance from the start
    map<long long, long long>nodeAndprevNode = {}; // {node, nodes predecessor}

    pair<long long, double>currPair = {}; // {node, distance from start to node}

    for (auto vertex : G.getVertices()) { // for vertices in the graph, add them into the worklist along with their distances from the start
        if (vertex == start) {
            distances[vertex] = 0.00; // if vertex we are adding is the start, distance is start - start = 0
        }
        else {
            distances[vertex] = INF; // else distance is INF for right now, later we will actually calculate
        }
    }

    worklist.push({start, distances[start]}); // add the start vertex and its distance from strt (which is 0) as a pair to the worklist

    while (!worklist.empty()) { // as long as list isnt empty

        currPair = worklist.top(); // get first element in queue, top returns a pair
        // currPair.first will be the current vertex
        // currPair.second will be distance from the start vertex to the current vertex

        worklist.pop(); // remove first element in queue

        // if vertex removed right now is the target, yay we found our way, leave loop
        if (currPair.first == target) {
           break;
        }

        // if vertex removed is in ignoreNodes then we continue BUT if the vertex is start then we dont want to skip it as we need to proccess the vertex
        if (ignoreNodes.count(currPair.first) == 1 && currPair.first != start) {
            continue; 
        }

        // else its not in ignore nodes, and not target then find the shortest distance from vertex to one of its neighbors // find closest neighbor
        for (auto currNeighbor : G.neighbors(currPair.first)) {
            
            if(currNeighbor != target && ignoreNodes.count(currNeighbor) == 1) {
                continue; // if neighbor is target, we dont want to skip it & if neighbor is in ignore nodes, dont process it
            }

            double tempDist; // declare temp dist variable to hold new distance (if set)
            bool weightSet = false; // declare bool variable to chec if weight was set
            weightSet = G.getWeight(currPair.first, currNeighbor, tempDist); // call getWeight() function and pass temp dist to store weight of edge between current vertex and current neighbor
            tempDist = tempDist + currPair.second; // add temp dist (distance from current vertex to neighbor) with distance from start to current vertex

            if (weightSet == true && tempDist < distances[currNeighbor]) { // if the weight was set and the new distance (temp dist) is greater than previously set distance then this bloc executes
                distances[currNeighbor] = tempDist; // reset distance of current neighbor from start to the new distance we found (temp dist)
                nodeAndprevNode[currNeighbor] = currPair.first; // set predecessor of current neighbor to our current vertex 
                worklist.push({currNeighbor, tempDist}); // add the current neighbor and its distance from start as a pair to the worklist to be processed later on
            }
        }    
    }

    // at this point we may have found the shortest path from start to target AND if we did we didnt store it in path, so we retrace it
    
    if(distances[target] == INF) {
        path.clear();
        return path; // if distance to target remains INF then we didnt find a path to target, return an empty path
    }
    else {
        // go back in the nodeAndPrevNode map starting from target, tracking predecessors all the way to start
        long long currNode = target;
        while (currNode != start) {
            path.push_back(currNode); // add currNode to our path
            currNode = nodeAndprevNode[currNode]; // set currNode to the predecessor of currNode
        }
        path.push_back(start); // lastly add start to the path, at this point the path is [target -> targets predecessor -> other nodes -> start]
 
        // after while loop finishes we have a path stored from start to target but the path vector has it stored in the reverse order (i.e : from target to start) so we reverse it to get start -> target
        reverse(path.begin(), path.end());
    }

    return path; // THEN WE RETURN PATH YAYAYAYAYAY :)
}

double pathLength(const graph<long long, double>& G, const vector<long long>& path) {
    double length = 0.0;
    double weight;
    for (size_t i = 0; i + 1 < path.size(); i++) {
        bool res = G.getWeight(path.at(i), path.at(i + 1), weight);
        assert(res);
        length += weight;
    }
    return length;
}

void outputPath(const vector<long long>& path) {
    for (size_t i = 0; i < path.size(); i++) {
        cout << path.at(i);
        if (i != path.size() - 1) {
            cout << "->";
        }
    }
    cout << endl;
}

void application(
    const vector<BuildingInfo>& Buildings,
    const graph<long long, double>& G) {
    string person1Building, person2Building;

    set<long long> buildingNodes;
    for (const auto& building : Buildings) {
        buildingNodes.insert(building.Coords.ID);
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
        cout << "Enter person 2's building (partial name or abbreviation)> ";
        getline(cin, person2Building);

        //
        // find the building coordinates
        //
        bool foundP1 = false;
        bool foundP2 = false;
        Coordinates P1Coords, P2Coords;
        string P1Name, P2Name;

        for (const BuildingInfo& building : Buildings) {
            if (building.Abbrev == person1Building) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (building.Abbrev == person2Building) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        for (const BuildingInfo& building : Buildings) {
            if (!foundP1 &&
                building.Fullname.find(person1Building) != string::npos) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (!foundP2 && building.Fullname.find(person2Building) != string::npos) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        if (!foundP1) {
            cout << "Person 1's building not found" << endl;
        } else if (!foundP2) {
            cout << "Person 2's building not found" << endl;
        } else {
            cout << endl;
            cout << "Person 1's point:" << endl;
            cout << " " << P1Name << endl;
            cout << " (" << P1Coords.Lat << ", " << P1Coords.Lon << ")" << endl;
            cout << "Person 2's point:" << endl;
            cout << " " << P2Name << endl;
            cout << " (" << P2Coords.Lat << ", " << P2Coords.Lon << ")" << endl;

            string destName;
            Coordinates destCoords;

            Coordinates centerCoords = centerBetween2Points(
                P1Coords.Lat, P1Coords.Lon, P2Coords.Lat, P2Coords.Lon);

            double minDestDist = numeric_limits<double>::max();

            for (const BuildingInfo& building : Buildings) {
                double dist = distBetween2Points(
                    centerCoords.Lat, centerCoords.Lon,
                    building.Coords.Lat, building.Coords.Lon);
                if (dist < minDestDist) {
                    minDestDist = dist;
                    destCoords = building.Coords;
                    destName = building.Fullname;
                }
            }

            cout << "Destination Building:" << endl;
            cout << " " << destName << endl;
            cout << " (" << destCoords.Lat << ", " << destCoords.Lon << ")" << endl;

            vector<long long> P1Path = dijkstra(G, P1Coords.ID, destCoords.ID, buildingNodes);
            vector<long long> P2Path = dijkstra(G, P2Coords.ID, destCoords.ID, buildingNodes);

            // This should NEVER happen with how the graph is built
            if (P1Path.empty() || P2Path.empty()) {
                cout << endl;
                cout << "At least one person was unable to reach the destination building. Is an edge missing?" << endl;
                cout << endl;
            } else {
                cout << endl;
                cout << "Person 1's distance to dest: " << pathLength(G, P1Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P1Path);
                cout << endl;
                cout << "Person 2's distance to dest: " << pathLength(G, P2Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P2Path);
            }
        }

        //
        // another navigation?
        //
        cout << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or #> ";
        getline(cin, person1Building);
    }
}
