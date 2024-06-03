// Graph.h
// header file for a graph class
// max glover
// it 279

#ifndef GRAPH_H
#define GRAPH_H


#include <iostream>
#include <vector>

using namespace std;

class Graph {
public:
    // Constructor and Destructor
    Graph();
    ~Graph();


    // Public methods
    bool readGraph(const string& filename);
    void printGraph();
    void computeTopologicalSort();
    void computeShortestPaths(const string& nodeName);
    void computeMinimumSpanningTree();


private:
    int numVertices;
    int numEdges;

// struct to represent an edge
    struct Edge {
        int from;
        int to;
        int weight;
    };


    // Private members
    vector<vector<int>> adjacencyMatrix;
    vector<string> nodeNames;
    // takes more space to have extra vector of edges but saves time complexity in printgraph and kruskals algorithm
    vector<Edge> edges;



    // Helper functions
    int getNodeIndex(const string& nodeName);
    void printPathAndCost(int startNode, int targetNode, const vector<int>& parent, const vector<int>& distance);
    void sortEdges(vector<Edge>& edges);


};


#endif // GRAPH_H





