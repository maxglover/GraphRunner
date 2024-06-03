// Graph.cpp
// by max glover
// it 279
// implementation of a graph with corresponding methods

#include "Graph.h"
#include <fstream>
#include <sstream>
#include <queue>
#include <algorithm>
#include "DisjointSet.h"



using namespace std;

//constructor, not sure if implementation is needed but starts with inintalizing size to 0
Graph::Graph() {
    numVertices=0;
    numEdges=0;
}

// do not have dynamically allocated memory so no cleanup necessary
Graph::~Graph() {
}

// reads in graph from input file
bool Graph::readGraph(const string& filename) {
   
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return false;
    }


    // Read the number of vertices from the file
    inputFile >> numVertices;


    // clears and replaces adjacency matrix with new matrix
    adjacencyMatrix.clear();
    adjacencyMatrix.resize(numVertices, vector<int>(numVertices, 0));
    nodeNames.resize(numVertices);


    // Read the node names
    for (int i = 0; i < numVertices; ++i) {
        inputFile >> nodeNames[i];
    }


    // Read the number of edges
    inputFile >> numEdges;
    edges.resize(0);


    // Read the edge data and fill the adjacency matrix
    for (int i = 0; i < numEdges; ++i) {
        string fromString, toString;
        int weight;
        inputFile >> fromString >> toString >> weight;
        int from = getNodeIndex(fromString);
        int to = getNodeIndex(toString);

        edges.push_back({from, to, weight});

        // make sure that from and to indices are within bounds- fixed an error for me
        if (from >= 0 && from < numVertices && to >= 0 && to < numVertices) {
            adjacencyMatrix[from][to] = weight;
        } 
        else {
            inputFile.close();
            return false;
        }
    }


    inputFile.close();
    return true;
}



//doesnt display the graph like the sample output order but it displays the graph correctly in the same order as read in
void Graph::printGraph() {
    cout << numVertices<<endl;
    for (int vertex =0; vertex<nodeNames.size();vertex++)
    {
       cout << nodeNames[vertex]<<endl;
    }
    cout<< numEdges<<endl;
    for (int edge = 0; edge<edges.size(); edge++)
    {
       cout << nodeNames[edges[edge].from]<< " " << nodeNames[edges[edge].to]<< " "<< edges[edge].weight<<endl;
    }   
}

// topological sort of graph if possible
void Graph::computeTopologicalSort() {
    vector<int> inDegree(adjacencyMatrix.size(), 0);

    // fill in degrees
    for (int row = 0; row < adjacencyMatrix.size(); row++) {
        for (int col = 0; col < adjacencyMatrix[row].size(); col++) {
            if (adjacencyMatrix[row][col]!=0)
            {
                inDegree[col]++;
            }
            
        }
       
    }

    queue<int> inDegreeQueue;
    //puts vertex on queue if at 0
    for (int i = 0; i < inDegree.size(); i++)
    {
        if(inDegree[i]==0){
            inDegreeQueue.push(i);
        }
    }

    vector<int> sorted;

    //loop, pops off first from queue and adds to our sort list
    //decrements in degrees of nodes the current points to when popped off
    while(!inDegreeQueue.empty()){
        int current = inDegreeQueue.front();
        inDegreeQueue.pop();
        sorted.push_back(current);

        for(int i =0; i<adjacencyMatrix.size(); i++){
            if (adjacencyMatrix[current][i]!=0)
            {
                inDegree[i]--;
            
                if (inDegree[i]==0)
                {
                    inDegreeQueue.push(i);
                }
            }
            
        }

    }

    // not a DAG
    if (sorted.size() != adjacencyMatrix.size()) {
        cout << "This graph cannot be topologically sorted." << endl;
        return;
    }
    else{
        cout << "Topological Sort:" << endl;
        for(int vertex : sorted){
            cout << nodeNames[vertex];
            if (vertex!=sorted.back()){
                cout << " --> ";
            }
        }
        cout << endl;
    }
}


void Graph::computeShortestPaths(const string& nodeName) {
    int startNode = getNodeIndex(nodeName);
    // usually we set costs to unknown nodes to infinity but I used int max(32 bit i think is right?)
    vector<int> cost(adjacencyMatrix.size(), INT32_MAX);
    vector<int> parent(adjacencyMatrix.size(), -1);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;


    cost[startNode] = 0;
    pq.push({0, startNode});

    //djikstras algorithm
    while (!pq.empty()) {
        int current = pq.top().second;
        pq.pop();


        for (int v = 0; v < adjacencyMatrix.size(); v++) {
            if (adjacencyMatrix[current][v] && cost[current] + adjacencyMatrix[current][v] < cost[v]) {
                cost[v] = cost[current] + adjacencyMatrix[current][v];
                parent[v] = current;
                pq.push({cost[v], v});
            }
        }
    }
    
    


    cout << "Shortest paths from " << nodeName << ":" << endl;
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        if (i != startNode) {
            printPathAndCost(startNode, i, parent, cost);
        }
    }
}

// kruskals algorithm
// did not check if graph is connected as precondition, unclear on what i am supposed to do if graph is not connected
void Graph::computeMinimumSpanningTree() {

    vector<Edge> mstEdges = edges;
    // Sort edges by weight
    sortEdges(mstEdges);
    DisjointSet mySet(adjacencyMatrix.size());


    vector<Edge> mst;
    int count =0;

    // for each edge, check if the edges are connected and if not add them to the minimum spanning tree while updating disjoint set
    while (mstEdges.size()!=count)
    {
        if (mySet.find(mstEdges[count].from) != mySet.find(mstEdges[count].to))
        {
            mst.push_back(mstEdges[count]);
            mySet.doUnion(mstEdges[count].from, mstEdges[count].to);
        }

        // if all vertices are connected break the loop
        if (mySet.find(mstEdges[count].from)==adjacencyMatrix.size())
        {
            break;
        }
        
        count++;
    }

    int totalCost=0;
    cout << "Minimum Spanning Tree:" << endl;
    for (int i = 0; i < mst.size(); i++)
    {
        cout << nodeNames[mst[i].from] << " -- " << nodeNames[mst[i].to] << " || Weight: " << mst[i].weight << endl;
        totalCost+=mst[i].weight;
    }
    cout << "Total Cost: "<< totalCost << endl;
    
}

//used in read graph
int Graph::getNodeIndex(const string& nodeName) {
    for (int i = 0; i < nodeNames.size(); i++)
    {
        if (nodeNames[i] == nodeName){
            return i;
        }
    }
    return -1;
}


//used in djikstras algorithm for shortest paths
void Graph::printPathAndCost(int startNode, int targetNode, const vector<int>& parent, const vector<int>& cost) {
    // If there is no path to the target node
    if (cost[targetNode] == INT32_MAX) {
        cout << "No path from " << nodeNames[startNode] << " to " << nodeNames[targetNode] << " found." << endl;
        return;
    }


    vector<int> path;


    // Traverse the parent array to make the path
    for (int v = targetNode; v != -1; v = parent[v]) {
        path.push_back(v);
    }


    cout << nodeNames[startNode];
    // Print the path in reverse order
    for (int i = path.size() - 2; i >= 0; i--) {
        cout << " --> " << nodeNames[path[i]];
    }


    cout << " || Weight: " << cost[targetNode] << endl;

    
}

// selection sort to sort edges
void Graph::sortEdges(vector<Edge>& edges) {
    int n = edges.size();


    for (int i = 0; i < n - 1; ++i) {
        // Find the minimum element in unsorted part of the vector
        int minIndex = i;
        for (int j = i + 1; j < n; ++j) {
            if (edges[j].weight < edges[minIndex].weight) {
                minIndex = j;
            }
            // if weights equal, sort by from and to index
            else if (edges[j].weight == edges[minIndex].weight)
            {
                if (edges[j].from<edges[minIndex].from)
                {
                    minIndex =j;
                }
                else if (edges[j].from==edges[minIndex].from && edges[j].to < edges[minIndex].to)
                {
                    minIndex = j;
                }
                
                
            }
            
        }


        // Swap the found minimum element with the first element
        swap(edges[minIndex], edges[i]);
    }
}




