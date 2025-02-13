#pragma once

#include <vector>

class UndirectedGraph
{
public:
    UndirectedGraph(int numVertices, std::vector<int> verticesA, std::vector<int> verticesB, std::vector<double> edgeWeights);

    void QuicksortByEdgeWeight();
    int GetNumVertices();
    int GetNumEdges();
    int GetFirstVertexAtIndex(int index);
    int GetSecondVertexAtIndex(int index);
    double GetEdgeWeightAtIndex(int index);
    std::vector<int>& GetEdgeListForVertex(int vertex);

private:
    int SelectPivotIndex(int startIndex, int endIndex);
    int Partition(int startIndex, int endIndex, int pivotIndex);
    void SwapEdges(int indexOne, int indexTwo);
\
    int m_numVertices;
    std::vector<int> m_verticesA;
    std::vector<int> m_verticesB;
    std::vector<double> m_edgeWeights;
    std::vector<std::vector<int> > m_edges;
};

