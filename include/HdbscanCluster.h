#pragma once

#include <set>
#include <vector>
#include <limits>
#include <stdexcept>

class HdbscanCluster
{
public:
    HdbscanCluster(int label, HdbscanCluster* parent, double birthLevel, int numPoints);

    bool operator==(const HdbscanCluster& other) const;

    void DetachPoints(int numPoints, double level);
    void Propagate();
    void AddPointsToVirtualChildCluster(std::set<int> points);
    bool VirtualChildClusterConstraintsPoint(int point);
    void ReleaseVirtualChildCluster();
    int GetClusterId();

    int Label;
    HdbscanCluster* Parent;
    std::vector<HdbscanCluster*> PropagatedDescendants;
    double PropagatedLowestChildDeathLevel = std::numeric_limits<double>::max();
    double Stability = 0;
    bool HasChildren = 0;
    int HierarchyPosition = 0;// First level where points with this cluster's label appear

private:
    void PropagateSelf();
    void PropagateDescendants();

    int m_id;
    int m_numPoints;

    double m_birthLevel;
    double m_deathLevel = 0;
    double m_propagatedStability = 0;
    std::set<int> m_virtualChildCluster;

    static int m_counter;
};
