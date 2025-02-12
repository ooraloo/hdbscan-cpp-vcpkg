#include "UndirectedGraph.h"

UndirectedGraph::UndirectedGraph(int numVertices, std::vector<int> verticesA, std::vector<int> verticesB, std::vector<double> edgeWeights) :
    m_numVertices(numVertices),
    m_verticesA(verticesA),
    m_verticesB(verticesB),
    m_edgeWeights(edgeWeights)
{
    m_edges.resize(m_numVertices);
    for (int i = 0; i < m_edgeWeights.size(); i++)
    {
        m_edges[m_verticesA[i]].push_back(m_verticesB[i]);

        if (m_verticesA[i] != m_verticesB[i])
        {
            m_edges[m_verticesB[i]].push_back(m_verticesA[i]);
        }
    }

}

void UndirectedGraph::QuicksortByEdgeWeight()
{
    int edgeWeightsLength = static_cast<int>(m_edgeWeights.size());
    if (edgeWeightsLength <= 1)
    {
        return;
    }

    std::vector<int> startIndexStack(edgeWeightsLength / 2);
    std::vector<int> endIndexStack(edgeWeightsLength / 2);

    startIndexStack[0] = 0;
    endIndexStack[0] = edgeWeightsLength - 1;

    int stackTop = 0;
    while (stackTop >= 0)
    {
        int startIndex = startIndexStack[stackTop];
        int endIndex = endIndexStack[stackTop];
        stackTop--;
        int pivotIndex = SelectPivotIndex(startIndex, endIndex);
        pivotIndex = Partition(startIndex, endIndex, pivotIndex);
        if (pivotIndex > startIndex + 1)
        {
            startIndexStack[stackTop + 1] = startIndex;
            endIndexStack[stackTop + 1] = pivotIndex - 1;
            stackTop++;
        }
        if (pivotIndex < endIndex - 1)
        {
            startIndexStack[stackTop + 1] = pivotIndex + 1;
            endIndexStack[stackTop + 1] = endIndex;
            stackTop++;
        }
    }
}
int UndirectedGraph::SelectPivotIndex(int startIndex, int endIndex)
{
    if (startIndex - endIndex <= 1)
        return startIndex;

    double first = m_edgeWeights[startIndex];
    double middle = m_edgeWeights[startIndex + (endIndex - startIndex) / 2];
    double last = m_edgeWeights[endIndex];

    if (first <= middle)
    {
        if (middle <= last)
            return startIndex + (endIndex - startIndex) / 2;

        if (last >= first)
            return endIndex;

        return startIndex;
    }

    if (first <= last)
        return startIndex;

    if (last >= middle)
        return endIndex;

    return startIndex + (endIndex - startIndex) / 2;
}

int UndirectedGraph::Partition(int startIndex, int endIndex, int pivotIndex)
{
    double pivotValue = m_edgeWeights[pivotIndex];
    SwapEdges(pivotIndex, endIndex);
    int lowIndex = startIndex;
    for (int i = startIndex; i < endIndex; i++)
    {
        if (m_edgeWeights[i] < pivotValue)
        {
            SwapEdges(i, lowIndex);
            lowIndex++;
        }
    }
    SwapEdges(lowIndex, endIndex);
    return lowIndex;
}

void UndirectedGraph::SwapEdges(int indexOne, int indexTwo)
{
    if (indexOne == indexTwo)
        return;

    int tempVertexA = m_verticesA[indexOne];
    int tempVertexB = m_verticesB[indexOne];
    double tempEdgeDistance = m_edgeWeights[indexOne];
    m_verticesA[indexOne] = m_verticesA[indexTwo];
    m_verticesB[indexOne] = m_verticesB[indexTwo];
    m_edgeWeights[indexOne] = m_edgeWeights[indexTwo];
    m_verticesA[indexTwo] = tempVertexA;
    m_verticesB[indexTwo] = tempVertexB;
    m_edgeWeights[indexTwo] = tempEdgeDistance;
}

int UndirectedGraph::GetNumVertices()
{
    return m_numVertices;
}

int UndirectedGraph::GetNumEdges()
{
    return static_cast<int>(m_edgeWeights.size());
}

int UndirectedGraph::GetFirstVertexAtIndex(int index)
{
    return m_verticesA[index];
}

int UndirectedGraph::GetSecondVertexAtIndex(int index)
{
    return m_verticesB[index];
}

double UndirectedGraph::GetEdgeWeightAtIndex(int index)
{
    return m_edgeWeights[index];
}

std::vector<int>& UndirectedGraph::GetEdgeListForVertex(int vertex)
{
    return m_edges[vertex];
}
