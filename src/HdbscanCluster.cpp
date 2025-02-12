#include "HdbscanCluster.h"

int HdbscanCluster::m_counter = 0;

HdbscanCluster::HdbscanCluster(int label, HdbscanCluster* parent, double birthLevel, int numPoints) :
    Label(label),
    Parent(parent),
    m_birthLevel(birthLevel),
    m_numPoints(numPoints)
{
    m_id = ++m_counter;

    if (Parent != NULL)
    {
        Parent->HasChildren = true;
    }

    PropagatedDescendants.resize(0);
}
bool HdbscanCluster ::operator==(const HdbscanCluster& other) const
{
    return (this->m_id == other.m_id);
}
void HdbscanCluster::DetachPoints(int numPoints, double level)
{
    m_numPoints -= numPoints;
    Stability += (numPoints * (1 / level - 1 / m_birthLevel));

    if (m_numPoints == 0)
        m_deathLevel = level;
    else if (m_numPoints < 0)
        throw std::invalid_argument("Cluster cannot have less than 0 points.");
}

void HdbscanCluster::Propagate()
{
    if (Parent != NULL)
    {
        if (PropagatedLowestChildDeathLevel == std::numeric_limits<double>::max())
        {
            PropagatedLowestChildDeathLevel = m_deathLevel;
        }

        if (PropagatedLowestChildDeathLevel < Parent->PropagatedLowestChildDeathLevel)
        {
            Parent->PropagatedLowestChildDeathLevel = PropagatedLowestChildDeathLevel;
        }

        if (!HasChildren)
        {
            PropagateSelf();
        }
        else
        {
            // Choose the parent over descendants if there is a tie in stability:
            if (Stability >= m_propagatedStability)
            {
                PropagateSelf();
            }
            else
            {
                PropagateDescendants();
            }
        }
    }
}

void HdbscanCluster::PropagateSelf()
{
    Parent->m_propagatedStability += Stability;
    Parent->PropagatedDescendants.push_back(this);
}

void HdbscanCluster::PropagateDescendants()
{
    Parent->m_propagatedStability += m_propagatedStability;
    Parent->PropagatedDescendants.insert(Parent->PropagatedDescendants.end(), PropagatedDescendants.begin(), PropagatedDescendants.end());
}

void HdbscanCluster::AddPointsToVirtualChildCluster(std::set<int> points)
{
    for (std::set<int>::iterator it = points.begin(); it != points.end(); ++it)
    {
        m_virtualChildCluster.insert(*it);
    }
}
bool HdbscanCluster::VirtualChildClusterConstraintsPoint(int point)
{
    return (m_virtualChildCluster.find(point) != m_virtualChildCluster.end());
}

void HdbscanCluster::ReleaseVirtualChildCluster()
{
    m_virtualChildCluster.clear();
}

int HdbscanCluster::GetClusterId()
{
    return m_id;
}
