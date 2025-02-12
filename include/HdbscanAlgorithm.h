#pragma once

#include <limits>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <list>

#include "UndirectedGraph.h"
#include "HdbscanCluster.h"

class HdbscanAlgorithm
{
public:
    /// <summary>
    /// Calculates the core distances for each point in the data set, given some value for k.
    /// </summary>
    /// <param name="distances">A vector of vectors where index [i][j] indicates the jth attribute of data point i</param>
    /// <param name="k">Each point's core distance will be it's distance to the kth nearest neighbor</param>
    /// <returns> An array of core distances</returns>
    static std::vector<double> CalculateCoreDistances(std::vector<std::vector<double>> distances, int minClusterSize);

    static UndirectedGraph ConstructMst(std::vector<std::vector<double>> distances, std::vector<double> coreDistances);

    /// <summary>
    /// Propagates constraint satisfaction, stability, and lowest child death level from each child
    /// cluster to each parent cluster in the tree.  This method must be called before calling
    /// findProminentClusters() or calculateOutlierScores().
    /// </summary>
    /// <param name="clusters">A list of Clusters forming a cluster tree</param>
    static void ComputeHierarchyAndClusterTree(
        UndirectedGraph* mst,
        int minClusterSize,
        std::vector<std::vector<int>>& hierarchy,
        std::vector<double>& pointNoiseLevels,
        std::vector<int>& pointLastClusters,
        std::vector<HdbscanCluster*>& clusters);

    static void PropagateTree(std::vector<HdbscanCluster*>& clusters);

    // Label of 0 represents noise
    static std::vector<int> FindProminentClusters(std::vector<HdbscanCluster*>& clusters, std::vector<std::vector<int>>& hierarchy, int numPoints);

    static std::vector<double> FindMembershipScore(std::vector<int> clusterids, std::vector<double> coreDistances);

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
    static HdbscanCluster* CreateNewCluster(
        std::set<int>& points,
        std::vector<int>& clusterLabels,
        HdbscanCluster* parentCluster,
        int clusterLabel,
        double edgeWeight);
};
