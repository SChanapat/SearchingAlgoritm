
#include <iostream>
#include <time.h>
#include <limits>
#include <vector>

using namespace std;

const int MaxDistance = std::numeric_limits<int>::max();
const int graphSize = 5;
const float graphDensity = 0.7;
const int maxEdgeDistance = 20;
const int minEdgeDistance = 1;

enum edgeState
{
    visited,    // edge is in the tree
    unvisited,  // edge is out of the tree
    halfvisited // edge has one vertext in the tree and one vertex out of the tree
};
bool vertecesvisited[graphSize];

void printVerteces()
{
    for (int i = 0; i < graphSize; i++)
    {
        cout << "Vertex " << i << ": " << vertecesvisited[i] << endl;
    }
}

class edge
{

public:
    int _vertex1, _vertex2, _distance;
    bool trueEdge;
    edge(int vertex1, int vertex2, int distance)
    {
        _vertex1 = vertex1;
        _vertex2 = vertex2;
        _distance = distance;
        trueEdge = true;
    }

    bool hasVertex(int vertex)
    {
        if (_vertex1 == vertex || _vertex2 == vertex)
            return true;
    }

    edgeState getState()
    {
        edgeState result;
        if (vertecesvisited[_vertex1] && vertecesvisited[_vertex2])
            result = visited;
        else if (!vertecesvisited[_vertex1] && !vertecesvisited[_vertex2])
            result = unvisited;
        else
            result = halfvisited;
        return result;
    }

    bool isLooping()
    {
        if (vertecesvisited[_vertex1] && vertecesvisited[_vertex2])
            return true;
        else
            return false;
    }

    bool equals(edge otherEdge)
    {
        if ((_vertex1 == otherEdge._vertex1 && _vertex2 == otherEdge._vertex2) || (_vertex1 == otherEdge._vertex2 && _vertex2 == otherEdge._vertex1))
            return true;
        else
            return false;
    }

    void visit()
    {
        cout << "visiting: " << _vertex1 << " , " << _vertex2 << endl;
        vertecesvisited[_vertex1] = true;
        vertecesvisited[_vertex2] = true;
    }

    void print()
    {
        cout << _vertex1 << " | " << _vertex2 << " | " << _distance << " | " << getState() << endl;
    }
};

class edgeList
{
public:
    vector<edge> edges;

    edgeList()
    {

        for (int i = 0; i < graphSize; i++)
        {
            vertecesvisited[i] = false;
        }
    }

    void addEdge(edge next)
    {
        edges.push_back(next);
    }

    edge shortestEdge()
    {
        edge lowest = edges.front();

        for (edge n : edges)
        {
            if (n._distance < lowest._distance)
            {
                lowest = n;
            }
        }

        return lowest;
    }

    void visitNextEdge()
    {
        edge lowest = edge(0, 0, maxEdgeDistance + 1);
        lowest.trueEdge = false;

        for (edge n : edges)
        {
            if (n.getState() == halfvisited && (lowest.trueEdge == false || n._distance < lowest._distance))
            {
                lowest = n;
                cout << "----"; // TODO: this marks current shortest edge but also appears on first edge when *lowest isnt null anymore
            }
        }

        if (lowest.trueEdge == false)
        {
            cout << "Tree broken!" << endl;
            abort();
        }
        else
        {
            lowest.visit();
        }
    }

    edge *shortestEdgeHalfvisited()
    {
        edge *lowest = new edge(0, 0, 1);
        lowest->trueEdge = false;

        for (edge n : edges)
        {
            if (n.getState() == halfvisited && (lowest->trueEdge == false || n._distance < lowest->_distance))
            {
                lowest = &n;
                cout << "----";
            }
            n.print();
        }

        cout << "Testing for breakage. 0" << endl;
        if (lowest->trueEdge == false){
            cout << "Tree broken!" << endl;
            abort();
        }
        cout << "Testing for breakage. 1" << endl;
        return lowest;
    }

    void printEdges(){
        cout << "===============1" << endl;
        for (edge n : edges){
            // if (n._distance < lowest._distance)
            n.print();
        }
        cout << "===============2" << endl;
    }

    // Used to detect if the tree is complete and all verteces connected
    bool allVertecesvisited(){
        cout << "Testing verteces visited." << endl;
        for (int i = 0; i < graphSize; i++){
            if (!vertecesvisited[i])
                return false;
        }

        return true;
    }

    void PrimMST()
    {

        int firstVertex = rand() % graphSize;
        vertecesvisited[firstVertex] = true;

        bool keepRunning = true;

        while (keepRunning)
        {
            printVerteces();
            visitNextEdge();

            if (allVertecesvisited()){
                cout << "Tree Complete!" << endl;
                keepRunning = false;
            }

        }
    }
};
class pathFinder
{

private:
    int vertexCount;
    int vertecesValues[graphSize]; // Keeps values for the Dijksta algorithm
    bool exploredVerteces[graphSize];
    bool vertexvisited[graphSize]; // Prim algorithm data - which vertexes are in the tree (visited)

    int generateEdge()
    {
        float chanceForEdge = (double)(rand() % 11) / 10;
        int edgeDistance = -1;

        if (chanceForEdge < graphDensity) // If chance for edge falls within the density range, generate an edge.
        {
            edgeDistance = minEdgeDistance + rand() % (maxEdgeDistance - minEdgeDistance);
        }

        return edgeDistance;
    }

    // Used in shortestPath() - Finds  the lowest value unexplored vertex for the Dijkstra algorithm
    int getLowestValueUnexploredVertex()
    {
        int lowestIndex = -1;
        int lowestValue = MaxDistance;

        for (int i = 0; i < graphSize; i++)
        {
            // If vertex is unexplored and of lowest value, keep its index.
            if (!exploredVerteces[i] && vertecesValues[i] < lowestValue)
            {
                lowestIndex = i;
                lowestValue = vertecesValues[i];
            }
        }
        return lowestIndex;
    }

    void clearPathValues()
    {
        for (int i = 0; i < graphSize; i++)
        {
            vertecesValues[i] = MaxDistance;
        }
    }

public:
    int graphMatrixDistances[graphSize][graphSize]; // Keeps original graph
    //  pathFinder Constructor
    pathFinder() : vertexCount(graphSize)
    {
        std::fill_n(exploredVerteces, graphSize, false);
        std::fill_n(vertexvisited, graphSize, false);

        // graphDensity = 0.7;
        // maxEdgeDistance = 20;
        // minEdgeDistance = 1;
    }

    void generateGraph()
    {
        // For each vertex
        for (int i = 0; i < vertexCount; i++)
        {
            // In relation to each other vertex
            for (int j = 0; j < vertexCount; j++)
            {
                int currentEdgeDistance = -1; // Default value is -1, meaning no connection between the two verteces.

                if (i == j) // All verteces have 0 distance to themselves.
                {
                    currentEdgeDistance = 0;
                }
                else if (i < j) // If i < j this is the first time this edge distance is set. I'm using generateEdge() to set the distance.
                {
                    currentEdgeDistance = generateEdge();
                }
                else // in this case, j is bigger than i, and this edge already exists (this being an undirected graph).
                {
                    currentEdgeDistance = graphMatrixDistances[j][i];
                }

                cout << currentEdgeDistance << " [" << i << ", " << j << "]"
                     << "\t";
                graphMatrixDistances[i][j] = currentEdgeDistance;
            }

            cout << endl;
        }
    }

    void printGraph()
    {

        cout << endl;

        for (int i = 0; i < vertexCount; i++)
        {
            // In relation to each other vertex
            for (int j = 0; j < vertexCount; j++)
            {
                cout << graphMatrixDistances[i][j] << "\t"
                     << "[" << i << ", " << j << "]"
                     << "\t";
            }

            cout << endl;
        }

        cout << endl;
    }

    // shortestPath finds shortest path from origin vertex to target vertex using Dijkstra's Algorithm.

    void shortestPath(int origin, int target)
    {
        clearPathValues();

        int cursor = origin;
        vertecesValues[origin] = 0;

        bool keepSearching = true;

        while (keepSearching)
        {

            exploredVerteces[cursor] = true;
            // update values on linked verteces
            for (int i = 0; i < graphSize; i++)
            {

                // This is meaningless target isn't linked to cursor
                int totalPathToTarget = vertecesValues[cursor] + graphMatrixDistances[cursor][i];

                if (i != cursor && graphMatrixDistances[origin][i] != -1 && totalPathToTarget < vertecesValues[i])
                {
                    // update taget value
                    vertecesValues[i] = totalPathToTarget;

                    if (i == target)
                    {
                        cout << "Shortest Path Found !!!" << endl;
                        cout << "Destination shortest distance: " << vertecesValues[i] << endl;
                        keepSearching = false;
                        break;
                    }
                }
            }

            if (cursor == getLowestValueUnexploredVertex())
            {
                cout << "No Path between origin and destination." << endl;
            }
            else
            {
                cursor = getLowestValueUnexploredVertex();
                cout << "Next Cursor: " << cursor << endl;
                break; // keepSearching = false;
            }
        }
    }
};

int main()
{
    srand(time(0));

    pathFinder p = pathFinder();
    p.generateGraph();
    p.printGraph();
    // p.shortestPath(0, 4);

    edgeList prim = edgeList();

    for (int i = 0; i < graphSize; i++)
    {
        for (int j = 0; j < graphSize; j++)
        {
            if (j < i) // Take edges only once.
            {
                if (p.graphMatrixDistances[i][j] != -1)
                {
                    edge newEdge = edge(i, j, p.graphMatrixDistances[i][j]);
                    prim.addEdge(newEdge);
                }
            }
        }
    }
    cout << "Printing edges before MST:" << endl;
    prim.printEdges();
    cout << "===============" << endl;

    prim.PrimMST();

    prim.printEdges();

    return 0;
}