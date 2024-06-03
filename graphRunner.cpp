//IT 279 F2020 - Nick Dunne and __max glover____
//Short program to handle the menu interface to the Graph class
#include "Graph.h"
#include "TextMenu.h"
#include <iostream>
#include <string>

//Prints a short welcome message to the user
void printWelcomeMessage();

//Executes the program with the given menu and graph
void handleMenu(TextMenu &menu, Graph &graph);

//Input the file into the graph for the user
void readFileIntoGraph(Graph &graph);

//Asks the user to input a file name for a graph, and returns that input to the user
string promptFileName();

//Asks the user for the node to get shortest paths from, and returns the input to the user
string promptShortestPathsNode();

const int QUIT_CHOICE = 6;

int main()
{
    //Create menu object to use as program input later
    TextMenu menu = TextMenu();

    //Create options to put in the menu
    string menuChoices[6] = {
        /*1*/ "Read in a new graph",
        /*2*/ "Print out the graph",
        /*3*/ "Topological Sort",
        /*4*/ "Shortest Path (from a node of your choice to all other nodes)",
        /*5*/ "Minimum Spanning Tree",
        /*6*/ "Exit the program"};

    //Add the menu choices to the TextMenu
    menu.setUp(6, menuChoices);

    //Create graph object to use later
    Graph graph = Graph();

    //Print welcome message
    printWelcomeMessage();

    //Start the menu
    handleMenu(menu, graph);

    return 0;
}

void printWelcomeMessage()
{
    //Short welcome message to display to users
    cout << "Welcome to Program 5 - Working with Graph Algorithms." << endl;
}

void handleMenu(TextMenu &menu, Graph &graph)
{

    //Read a user-provided file into the Graph class
    readFileIntoGraph(graph);

    //String value to keep the node for the shortest path
    string node;
    // get first menu choice from user
    int menuChoice = menu.getChoice();

    while (menuChoice != QUIT_CHOICE)
    {
        switch (menuChoice)
        {
        case 1:
            //Read in another graph here
            readFileIntoGraph(graph);
            break;
        case 2:
            //Print out the graph that is in the file
            // << your call here >>
            graph.printGraph();
            break;
        case 3:
            //Call topological sort
            // << your call here >>
            graph.computeTopologicalSort();
            break;
        case 4:
            //Retrieve input from user to ask for a node
            node = promptShortestPathsNode();

            //Call shortest paths
            // << your call here >>
            graph.computeShortestPaths(node);
            break;
        case 5:
            //Call minimum spanning tree
            // << your call here >>
            graph.computeMinimumSpanningTree();
            break;
        }
        menuChoice = menu.getChoice();
    }
}

void readFileIntoGraph(Graph &graph)
{
    //Get file name from user
    string fileName = promptFileName();

    //Attempt to read in the file
    bool fileReadSuccessfully;
    
    // << your call here >>
    if (graph.readGraph(fileName))
    {
        fileReadSuccessfully = true;
    }


    //If the read doesn't work, keep asking for a file until it does.
    while (!fileReadSuccessfully)
    {
        //Provide error message and new prompt to the user
        cout << "That input file does not exist, please try again with an existing input file." << endl;
        fileName = promptFileName();

        //Attempt to read the file in again
        // << your call here >>
        if (graph.readGraph(fileName))
    {
        fileReadSuccessfully = true;
    }
    }
}

/*
    ########## You shouldn't need to change anything below this line ##########
*/

string promptFileName()
{
    //String to store the file name
    string fileName;

    //Generate the prompt for the user
    cout << "Please enter the name of a text file to retrieve your graph from: " << flush;

    //Input retrieved from user
    getline(cin, fileName);

    //Return the file name to the calling method (the menu method)
    return fileName;
}

string promptShortestPathsNode()
{
    //String to store the node name given by the user
    string node;

    //Generate the prompt for the user
    cout << "Please type the name of a node to get shortest paths for: " << flush;

    //Input retrieved from the user
    getline(cin, node);

    //Return the typed in selection made by the user
    return node;
}
