#include <limits>
#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <list>

#include "UndirectedGraph.h"
#include "HdbscanCluster.h"
#include "HdbscanAlgorithm.h"

template<class InputIt, class T = typename std::iterator_traits<InputIt>::value_type>
constexpr InputIt std_find(InputIt first, InputIt last, const T& value)
{
    for (; first != last; ++first)
        if (*first == value)
            return first;

    return last;
}

std::vector<double> HdbscanAlgorithm::CalculateCoreDistances(std::vector<std::vector<double> > distances, int minClusterSize)
{
    auto length = distances.size();
    std::vector<double> coreDistances(length);
    int numNeighbors = minClusterSize - 1;

    if (minClusterSize == 1)
    {
        for (int point = 0; point < length; point++)
        {
            coreDistances[point] = 0;
        }
        return coreDistances;
    }

    for (int point = 0; point < length; point++)
    {
        std::vector<double> pointDistances;
        pointDistances.reserve(length - 1);
        for (int neighbor = 0; neighbor < length; neighbor++)
        {
            if (point == neighbor)
            {
                continue;
            }

            pointDistances.push_back(distances[point][neighbor]);
        }
        std::sort(pointDistances.begin(), pointDistances.end());
        coreDistances[point] = pointDistances[numNeighbors - 1];
    }

    return coreDistances;
}

UndirectedGraph HdbscanAlgorithm::ConstructMst(std::vector<std::vector<double> > distances, std::vector<double> coreDistances)
{
    auto length = static_cast<int>(distances.size());
    auto numMRDNeighbors = 2 * length - 1;
    std::vector<bool> attachedPoints(length);

    std::vector<int> nearestMRDNeighbors(numMRDNeighbors);
    std::vector<double> nearestMRDDistances(numMRDNeighbors);

    for (int i = 0; i < length - 1; i++)
    {
        nearestMRDDistances[i] = std::numeric_limits<double>::max();
    }

    auto currentPoint = length - 1;
    int numAttachedPoints = 1;
    attachedPoints[length - 1] = true;

    while (numAttachedPoints < length)
    {
        int nearestMRDPoint = -1;
        double nearestMRDDistance = std::numeric_limits<double>::max();
        for (int neighbor = 0; neighbor < length; neighbor++)
        {
            if (currentPoint == neighbor)
                continue;
            if (attachedPoints[neighbor])
                continue;

            double distance = distances[currentPoint][neighbor];
            double mutualReachabiltiyDistance = distance;
            if (coreDistances[currentPoint] > mutualReachabiltiyDistance)
                mutualReachabiltiyDistance = coreDistances[currentPoint];

            if (coreDistances[neighbor] > mutualReachabiltiyDistance)
                mutualReachabiltiyDistance = coreDistances[neighbor];

            if (mutualReachabiltiyDistance < nearestMRDDistances[neighbor])
            {
                nearestMRDDistances[neighbor] = mutualReachabiltiyDistance;
                nearestMRDNeighbors[neighbor] = currentPoint;
            }

            if (nearestMRDDistances[neighbor] <= nearestMRDDistance)
            {
                nearestMRDDistance = nearestMRDDistances[neighbor];
                nearestMRDPoint = neighbor;
            }

        }
        attachedPoints[nearestMRDPoint] = true;
        numAttachedPoints++;
        currentPoint = nearestMRDPoint;
    }

    std::vector<int> otherVertexIndices(numMRDNeighbors);
    for (int i = 0; i < length - 1; i++)
    {
        otherVertexIndices[i] = i;
    }

    // self edges
    for (auto i = length - 1; i < numMRDNeighbors; i++)
    {
        auto vertex = i - (length - 1);
        nearestMRDNeighbors[i] = vertex;
        otherVertexIndices[i] = vertex;
        nearestMRDDistances[i] = coreDistances[vertex];
    }

    UndirectedGraph undirectedGraphObject(length, nearestMRDNeighbors, otherVertexIndices, nearestMRDDistances);
    return undirectedGraphObject;
}

void HdbscanAlgorithm::ComputeHierarchyAndClusterTree(
    UndirectedGraph* mst,
    int minClusterSize,
    std::vector<std::vector<int> >& hierarchy,
    std::vector<double>& pointNoiseLevels,
    std::vector<int>& pointLastClusters,
    std::vector<HdbscanCluster*>& clusters)
{
    int hierarchyPosition = 0;

    //The current edge being removed from the MST:
    int currentEdgeIndex = mst->GetNumEdges() - 1;
    int nextClusterLabel = 2;
    bool nextLevelSignificant = true;

    //The previous and current cluster numbers of each point in the data set:
    std::vector<int> previousClusterLabels(mst->GetNumVertices());
    std::vector<int> currentClusterLabels(mst->GetNumVertices());

    for (int i = 0; i < currentClusterLabels.size(); i++)
    {
        currentClusterLabels[i] = 1;
        previousClusterLabels[i] = 1;
    }

    clusters.push_back(NULL);
    clusters.push_back(new HdbscanCluster(1, NULL, std::numeric_limits<double>::quiet_NaN(), mst->GetNumVertices()));

    std::set<int> clusterOne;
    clusterOne.insert(1);

    std::set<int> affectedClusterLabels;
    std::set<int> affectedVertices;
    while (currentEdgeIndex >= 0)
    {
        double currentEdgeWeight = mst->GetEdgeWeightAtIndex(currentEdgeIndex);
        std::vector<HdbscanCluster*> newClusters;
        while (currentEdgeIndex >= 0 && mst->GetEdgeWeightAtIndex(currentEdgeIndex) == currentEdgeWeight)
        {
            int firstVertex = mst->GetFirstVertexAtIndex(currentEdgeIndex);
            int secondVertex = mst->GetSecondVertexAtIndex(currentEdgeIndex);
            std::vector<int>& firstVertexEdgeList = mst->GetEdgeListForVertex(firstVertex);
            std::vector<int>::iterator secondVertexInFirstEdgeList = std_find(firstVertexEdgeList.begin(), firstVertexEdgeList.end(), secondVertex);
            if (secondVertexInFirstEdgeList != mst->GetEdgeListForVertex(firstVertex).end())
                mst->GetEdgeListForVertex(firstVertex).erase(secondVertexInFirstEdgeList);
            std::vector<int>& secondVertexEdgeList = mst->GetEdgeListForVertex(secondVertex);
            std::vector<int>::iterator firstVertexInSecondEdgeList = std_find(secondVertexEdgeList.begin(), secondVertexEdgeList.end(), firstVertex);
            if (firstVertexInSecondEdgeList != mst->GetEdgeListForVertex(secondVertex).end())
                mst->GetEdgeListForVertex(secondVertex).erase(firstVertexInSecondEdgeList);

            if (currentClusterLabels[firstVertex] == 0)
            {
                currentEdgeIndex--;
                continue;
            }

            affectedVertices.insert(firstVertex);
            affectedVertices.insert(secondVertex);
            affectedClusterLabels.insert(currentClusterLabels[firstVertex]);
            currentEdgeIndex--;
        }
        if (!affectedClusterLabels.size())
            continue;
        while (affectedClusterLabels.size())
        {
            int examinedClusterLabel = *prev(affectedClusterLabels.end());
            affectedClusterLabels.erase(prev(affectedClusterLabels.end()));
            std::set<int> examinedVertices;
            //std::set<int>::iterator affectedIt;
            for (auto affectedIt = affectedVertices.begin(); affectedIt != affectedVertices.end();)
            {
                int vertex = *affectedIt;
                if (currentClusterLabels[vertex] == examinedClusterLabel)
                {
                    examinedVertices.insert(vertex);
                    affectedIt = affectedVertices.erase(affectedIt);

                }
                else
                {
                    ++affectedIt;
                }
            }
            std::set<int> firstChildCluster;
            std::list<int> unexploredFirstChildClusterPoints;
            int numChildClusters = 0;
            while (examinedVertices.size())
            {
                std::set<int> constructingSubCluster;
                std::list<int> unexploredSubClusterPoints;
                bool anyEdges = false;
                bool incrementedChildCount = false;
                int rootVertex = *prev(examinedVertices.end());
                constructingSubCluster.insert(rootVertex);
                unexploredSubClusterPoints.push_back(rootVertex);
                examinedVertices.erase(prev(examinedVertices.end()));
                while (unexploredSubClusterPoints.size())
                {
                    int vertexToExplore = *unexploredSubClusterPoints.begin();
                    unexploredSubClusterPoints.erase(unexploredSubClusterPoints.begin());
                    std::vector<int>& vertexToExploreEdgeList = mst->GetEdgeListForVertex(vertexToExplore);
                    for (std::vector<int>::iterator it = vertexToExploreEdgeList.begin(); it != vertexToExploreEdgeList.end();)
                    {
                        int neighbor = *it;
                        anyEdges = true;
                        if (std_find(constructingSubCluster.begin(), constructingSubCluster.end(), neighbor) == constructingSubCluster.end())
                        {
                            constructingSubCluster.insert(neighbor);
                            unexploredSubClusterPoints.push_back(neighbor);
                            if (std_find(examinedVertices.begin(), examinedVertices.end(), neighbor) != examinedVertices.end())
                                examinedVertices.erase(std_find(examinedVertices.begin(), examinedVertices.end(), neighbor));

                        }
                        else
                        {
                            ++it;
                        }
                    }
                    if (!incrementedChildCount && constructingSubCluster.size() >= minClusterSize && anyEdges)
                    {
                        incrementedChildCount = true;
                        numChildClusters++;

                        //If this is the first valid child cluster, stop exploring it:
                        if (firstChildCluster.size() == 0)
                        {
                            firstChildCluster = constructingSubCluster;
                            unexploredFirstChildClusterPoints = unexploredSubClusterPoints;
                            break;
                        }
                    }

                }
                //If there could be a split, and this child cluster is valid:
                if (numChildClusters >= 2 && constructingSubCluster.size() >= minClusterSize && anyEdges)
                {
                    //Check this child cluster is not equal to the unexplored first child cluster:
                    int firstChildClusterMember = *prev(firstChildCluster.end());
                    if (std_find(constructingSubCluster.begin(), constructingSubCluster.end(), firstChildClusterMember) != constructingSubCluster.end())
                        numChildClusters--;
                    //Otherwise, c a new cluster:
                    else
                    {
                        HdbscanCluster* newCluster = CreateNewCluster(constructingSubCluster, currentClusterLabels,
                            clusters[examinedClusterLabel], nextClusterLabel, currentEdgeWeight);
                        newClusters.push_back(newCluster);
                        clusters.push_back(newCluster);
                        nextClusterLabel++;
                    }
                }
                else if (constructingSubCluster.size() < minClusterSize || !anyEdges)
                {
                    CreateNewCluster(constructingSubCluster, currentClusterLabels,
                        clusters[examinedClusterLabel], 0, currentEdgeWeight);

                    for (std::set<int>::iterator it = constructingSubCluster.begin(); it != constructingSubCluster.end(); it++)
                    {
                        int point = *it;
                        pointNoiseLevels[point] = currentEdgeWeight;
                        pointLastClusters[point] = examinedClusterLabel;
                    }
                }
            }
            if (numChildClusters >= 2 && currentClusterLabels[*firstChildCluster.begin()] == examinedClusterLabel)
            {
                while (unexploredFirstChildClusterPoints.size())
                {
                    int vertexToExplore = *unexploredFirstChildClusterPoints.begin();
                    unexploredFirstChildClusterPoints.pop_front();
                    for (std::vector<int>::iterator it = mst->GetEdgeListForVertex(vertexToExplore).begin(); it != mst->GetEdgeListForVertex(vertexToExplore).end(); it++)
                    {
                        int neighbor = *it;
                        if (std_find(firstChildCluster.begin(), firstChildCluster.end(), neighbor) == firstChildCluster.end())
                        {
                            firstChildCluster.insert(neighbor);
                            unexploredFirstChildClusterPoints.push_back(neighbor);
                        }
                    }
                }
                HdbscanCluster* newCluster = CreateNewCluster(firstChildCluster, currentClusterLabels,
                    clusters[examinedClusterLabel], nextClusterLabel, currentEdgeWeight);
                newClusters.push_back(newCluster);
                clusters.push_back(newCluster);
                nextClusterLabel++;
            }
        }
        if (nextLevelSignificant || newClusters.size())
        {
            std::vector<int> lineContents(previousClusterLabels.size());
            for (int i = 0; i < previousClusterLabels.size(); i++)
                lineContents[i] = previousClusterLabels[i];
            hierarchy.push_back(lineContents);
            hierarchyPosition++;
        }
        std::set<int> newClusterLabels;
        for (std::vector<HdbscanCluster*>::iterator it = newClusters.begin(); it != newClusters.end(); it++)
        {
            HdbscanCluster* newCluster = *it;
            newCluster->HierarchyPosition = hierarchyPosition;
            newClusterLabels.insert(newCluster->Label);
        }

        for (int i = 0; i < previousClusterLabels.size(); i++)
        {
            previousClusterLabels[i] = currentClusterLabels[i];
        }
        if (!newClusters.size())
            nextLevelSignificant = false;
        else
            nextLevelSignificant = true;
    }

    {
        std::vector<int> lineContents(previousClusterLabels.size() + 1);
        for (int i = 0; i < previousClusterLabels.size(); i++)
            lineContents[i] = 0;
        hierarchy.push_back(lineContents);
    }
}

void HdbscanAlgorithm::PropagateTree(std::vector<HdbscanCluster*>& clusters)
{
    std::map<int, HdbscanCluster*> clustersToExamine;
    std::unordered_set<int> markedClusters;

    //Find all leaf clusters in the cluster tree:
    for (HdbscanCluster* cluster : clusters)
    {
        if (cluster != NULL && !cluster->HasChildren)
        {
            int label = cluster->Label;
            clustersToExamine.erase(label);
            clustersToExamine.insert({ label, cluster });
            markedClusters.insert(label);
        }
    }
    //Iterate through every cluster, propagating stability from children to parents:
    while (clustersToExamine.size())
    {
        std::map<int, HdbscanCluster*>::iterator currentKeyValue = prev(clustersToExamine.end());
        HdbscanCluster* currentCluster = currentKeyValue->second;
        clustersToExamine.erase(currentKeyValue->first);
        currentCluster->Propagate();

        if (currentCluster->Parent != NULL)
        {
            HdbscanCluster* parent = currentCluster->Parent;
            int label = parent->Label;

            if (markedClusters.find(label) == markedClusters.end())
            {
                clustersToExamine.erase(label);
                clustersToExamine.insert({ label, parent });
                markedClusters.insert(label);
            }
        }
    }
}

std::vector<int> HdbscanAlgorithm::FindProminentClusters(std::vector<HdbscanCluster*>& clusters, std::vector<std::vector<int> >& hierarchy, int numPoints)
{
    //Take the list of propagated clusters from the root cluster:
    std::vector<HdbscanCluster*> solution = clusters[1]->PropagatedDescendants;
    std::vector<int> flatPartitioning(numPoints);

    //Store all the hierarchy positions at which to find the birth points for the flat clustering:
    std::map<int, std::vector<int> > significantHierarchyPositions;

    std::vector<HdbscanCluster*>::iterator it = solution.begin();
    while (it != solution.end())
    {
        int hierarchyPosition = (*it)->HierarchyPosition;
        if (significantHierarchyPositions.count(hierarchyPosition) > 0)
            significantHierarchyPositions[hierarchyPosition].push_back((*it)->Label);
        else
            significantHierarchyPositions[hierarchyPosition].push_back((*it)->Label);
        it++;
    }

    //Go through the hierarchy file, setting labels for the flat clustering:
    while (significantHierarchyPositions.size())
    {
        std::map<int, std::vector<int> >::iterator entry = significantHierarchyPositions.begin();
        std::vector<int> clusterList = entry->second;
        int hierarchyPosition = entry->first;
        significantHierarchyPositions.erase(entry->first);

        std::vector<int> lineContents = hierarchy[hierarchyPosition];

        for (int i = 0; i < lineContents.size(); i++)
        {
            int label = lineContents[i];
            if (std_find(clusterList.begin(), clusterList.end(), label) != clusterList.end())
                flatPartitioning[i] = label;
        }
    }
    return flatPartitioning;
}
std::vector<double> HdbscanAlgorithm::FindMembershipScore(std::vector<int> clusterids, std::vector<double> coreDistances)
{
    auto length = clusterids.size();
    std::vector<double> prob(length, std::numeric_limits<double>::max());
    int i = 0;

    while (i < length)
    {
        if (prob[i] == std::numeric_limits<double>::max())
        {
            int clusterno = clusterids[i];
            std::vector<int>::iterator iter = clusterids.begin() + i;
            std::vector<int> indices;
            while ((iter = std_find(iter, clusterids.end(), clusterno)) != clusterids.end())
            {
                indices.push_back(static_cast<int>(distance(clusterids.begin(), iter)));
                iter++;
                if (iter == clusterids.end())
                    break;
            }
            if (clusterno == 0)
            {
                for (int j = 0; j < indices.size(); j++)
                {
                    prob[indices[j]] = 0;
                }
                i++;
                continue;
            }
            std::vector<double> tempCoreDistances(indices.size());
            for (int j = 0; j < indices.size(); j++)
            {
                tempCoreDistances[j] = coreDistances[j];
            }
            double maxCoreDistance = *max_element(tempCoreDistances.begin(), tempCoreDistances.end());
            for (int j = 0; j < tempCoreDistances.size(); j++)
            {
                prob[indices[j]] = (maxCoreDistance - tempCoreDistances[j]) / maxCoreDistance;
            }
        }

        i++;
    }
    return prob;

}

/// <summary>
/// Removes the set of points from their parent Cluster, and creates a new Cluster, provided the
/// clusterId is not 0 (noise).
/// </summary>
/// <param name="points">The set of points to be in the new Cluster</param>
/// <param name="clusterLabels">An array of cluster labels, which will be modified</param>
/// <param name="parentCluster">The parent Cluster of the new Cluster being created</param>
/// <param name="clusterLabel">The label of the new Cluster </param>
/// <param name="edgeWeight">The edge weight at which to remove the points from their previous Cluster</param>
/// <returns>The new Cluster, or null if the clusterId was 0</returns>
HdbscanCluster* HdbscanAlgorithm::CreateNewCluster(
    std::set<int>& points,
    std::vector<int>& clusterLabels,
    HdbscanCluster* parentCluster,
    int clusterLabel,
    double edgeWeight)
{
    std::set<int>::iterator it = points.begin();
    while (it != points.end())
    {
        clusterLabels[*it] = clusterLabel;
        ++it;
    }
    parentCluster->DetachPoints(static_cast<int>(points.size()), edgeWeight);

    if (clusterLabel != 0)
    {
        return new HdbscanCluster(clusterLabel, parentCluster, edgeWeight, static_cast<int>(points.size()));
    }

    parentCluster->AddPointsToVirtualChildCluster(points);
    return NULL;
}
