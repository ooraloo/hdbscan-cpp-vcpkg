#include<iostream>
#include<cstdio>
#include<cstdint>
#include <cmath>
#include <string>

#include "HdbscanAlgorithm.h"

const int NUM_INPUTS = 300;
const int MIN_POINTS = 5;
const int MIN_CLUSTER_SIZE = 5;
const std::vector<double> m_testDataset({
	0.837, 2.136,
	-1.414, 7.410,
	1.155, 5.100,
	-1.019, 7.815,
	1.271, 1.893,
	3.438, 0.262,
	-1.808, 1.597,
	1.414, 4.381,
	-0.205, 8.432,
	-0.711, 8.660,
	-1.712, 2.778,
	-2.670, 8.354,
	1.243, 4.504,
	-2.228, 6.895,
	1.455, -0.029,
	0.454, 3.956,
	1.069, 4.531,
	2.569, 0.507,
	-1.067, 3.132,
	-1.079, 2.205,
	2.715, 1.291,
	1.777, 1.187,
	0.734, 5.037,
	-1.996, 2.852,
	-1.918, 2.605,
	-0.556, 4.696,
	1.697, 0.866,
	0.595, 4.710,
	-2.880, 2.304,
	0.187, 4.027,
	-0.513, 7.874,
	-2.057, 7.379,
	1.873, 4.181,
	-1.131, 6.767,
	-1.644, 7.941,
	-2.419, 7.440,
	-2.016, 7.484,
	-2.621, 7.986,
	-2.203, 2.479,
	1.664, 0.663,
	0.670, 3.595,
	-1.985, 2.055,
	-0.047, 5.474,
	1.545, 4.211,
	-1.702, 2.461,
	-1.022, 2.768,
	-1.378, 8.103,
	-1.552, 2.746,
	-1.479, 7.569,
	1.989, 1.510,
	-1.950, 3.484,
	2.455, 0.621,
	-0.894, 7.617,
	1.697, 0.755,
	1.756, 2.055,
	-1.111, 2.822,
	-0.042, 7.809,
	-1.141, 1.976,
	-1.806, 7.728,
	1.393, 0.929,
	-2.257, 7.302,
	0.572, 4.323,
	-1.550, 9.283,
	-1.038, 2.953,
	-2.110, 3.107,
	-1.187, 2.784,
	-2.458, 7.512,
	2.370, 0.951,
	-2.667, 7.848,
	-1.497, 3.214,
	1.322, 4.179,
	-0.487, 3.329,
	-1.037, 8.063,
	-1.605, 2.974,
	-1.504, 1.924,
	-0.785, 8.453,
	-1.758, 2.974,
	1.190, 4.728,
	2.140, 0.706,
	-1.035, 8.206,
	1.255, 0.090,
	0.596, 4.086,
	1.280, 1.058,
	1.730, 1.147,
	-0.949, 8.464,
	0.935, 5.332,
	2.369, 0.795,
	0.429, 4.974,
	-2.048, 6.654,
	-1.457, 7.487,
	0.529, 3.808,
	1.782, 0.908,
	-1.956, 8.616,
	-1.746, 3.012,
	-1.180, 3.128,
	1.164, 3.791,
	1.362, 1.366,
	2.601, 1.088,
	0.272, 5.470,
	-3.122, 3.282,
	-0.588, 8.614,
	1.669, -0.436,
	-0.683, 7.675,
	2.368, 0.552,
	1.052, 4.545,
	2.227, 1.263,
	2.439, -0.073,
	1.345, 4.857,
	-1.315, 6.839,
	0.983, 5.375,
	-1.063, 2.208,
	-1.607, 3.565,
	1.573, 0.484,
	-2.179, 8.086,
	1.834, 0.754,
	2.106, 3.495,
	-1.643, 7.527,
	1.106, 1.264,
	1.612, 1.823,
	0.460, 5.450,
	-0.538, 3.016,
	1.678, 0.609,
	-1.012, 3.603,
	1.342, 0.594,
	1.428, 1.624,
	2.045, 1.125,
	1.673, 0.659,
	-1.359, 2.322,
	1.131, 0.936,
	-1.739, 1.948,
	-0.340, 8.167,
	-1.638, 2.433,
	-1.688, 2.241,
	2.430, -0.064,
	-1.380, 7.185,
	-1.252, 2.339,
	-2.395, 3.398,
	-2.092, 7.481,
	0.488, 3.268,
	-0.539, 7.456,
	-2.592, 8.076,
	-1.047, 2.965,
	1.256, 3.382,
	-1.622, 4.272,
	1.869, 5.441,
	-1.764, 2.222,
	-1.382, 7.288,
	0.008, 4.176,
	-1.103, 7.302,
	-1.794, 7.581,
	-1.512, 7.944,
	0.959, 4.561,
	-0.601, 6.300,
	0.225, 4.770,
	1.567, 0.018,
	-1.034, 2.921,
	-0.922, 8.099,
	-1.886, 2.248,
	1.869, 0.956,
	1.101, 4.890,
	-1.932, 8.306,
	0.670, 4.041,
	0.744, 4.122,
	1.640, 1.819,
	0.815, 4.785,
	-2.633, 2.631,
	-0.961, 1.274,
	0.214, 4.885,
	1.435, 1.307,
	1.214, 3.648,
	1.083, 4.063,
	-1.226, 8.296,
	1.482, 0.690,
	1.896, 5.185,
	-1.324, 4.131,
	-1.150, 7.893,
	2.469, 1.679,
	2.311, 1.304,
	0.573, 4.088,
	-0.968, 3.122,
	2.625, 0.950,
	1.684, 4.196,
	-2.221, 2.731,
	-1.578, 3.034,
	0.082, 4.567,
	1.433, 4.377,
	1.063, 5.176,
	0.768, 4.398,
	2.470, 1.315,
	-1.732, 7.164,
	0.347, 3.452,
	-1.001, 2.849,
	1.016, 4.485,
	0.560, 4.214,
	-2.118, 2.035,
	-1.362, 2.383,
	-2.784, 2.992,
	1.652, 3.656,
	-1.940, 2.189,
	-1.815, 7.978,
	1.202, 3.644,
	-0.969, 3.267,
	1.870, -0.108,
	-1.807, 2.068,
	1.218, 3.893,
	-1.484, 6.008,
	-1.564, 2.853,
	-0.686, 8.683,
	1.076, 4.685,
	-0.976, 6.738,
	1.380, 4.548,
	-1.641, 2.681,
	-0.002, 4.581,
	1.714, 5.025,
	-1.405, 7.726,
	-0.708, 2.504,
	-0.886, 2.646,
	1.984, 0.490,
	2.952, -0.344,
	0.432, 4.335,
	-1.866, 7.625,
	2.527, 0.618,
	2.041, 0.455,
	-2.580, 3.188,
	1.620, 0.068,
	-2.588, 3.131,
	0.444, 3.115,
	-0.457, 7.306,
	-1.129, 7.805,
	2.130, 5.192,
	1.004, 4.191,
	-1.393, 8.746,
	0.728, 3.855,
	0.893, 1.011,
	-1.108, 2.920,
	0.789, 4.337,
	1.976, 0.719,
	-1.249, 3.085,
	-1.078, 8.881,
	-1.868, 3.080,
	2.768, 1.088,
	0.277, 4.844,
	3.411, 0.872,
	-1.581, 7.553,
	-1.530, 7.705,
	-1.825, 7.360,
	-1.686, 7.953,
	-1.651, 3.446,
	-1.304, 3.003,
	-0.731, 6.242,
	2.406, 4.870,
	-1.536, 3.014,
	1.489, 0.652,
	0.514, 4.627,
	-1.815, 3.290,
	-1.937, 3.914,
	-0.615, 3.950,
	2.032, 0.197,
	2.149, 1.037,
	-1.370, 7.770,
	0.914, 4.550,
	0.334, 4.936,
	-2.160, 3.410,
	1.367, 0.635,
	-0.571, 8.133,
	-1.006, 3.084,
	1.495, 3.858,
	-0.590, 7.695,
	0.715, 5.413,
	2.114, 1.247,
	1.201, 0.602,
	-2.546, 3.150,
	-1.959, 2.430,
	2.338, 3.431,
	3.353, 1.700,
	1.843, 0.073,
	1.320, 1.404,
	2.097, 4.847,
	-1.243, 8.152,
	-1.859, 7.789,
	2.747, 1.545,
	2.608, 1.089,
	1.660, 3.563,
	2.352, 0.828,
	2.223, 0.839,
	3.229, 1.132,
	-1.559, 7.248,
	-0.647, 3.429,
	-1.327, 8.515,
	0.917, 3.906,
	2.295, -0.766,
	1.816, 1.120,
	-1.120, 7.110,
	-1.655, 8.614,
	-1.276, 7.968,
	1.974, 1.580,
	2.518, 1.392,
	0.439, 4.536,
	0.369, 7.791,
	-1.791, 2.750
});

std::vector<int> m_correctClusters({
	2, 5, 6, 5, 2, 2, 7, 6, 5, 5, 7, 5, 6, 5, 2, 6, 6, 2, 7, 7, 2, 2, 6, 7, 7, 0, 2, 6, 0, 6, 5, 5, 6, 5, 
	5, 5, 5, 5, 7, 2, 6, 7, 6, 6, 7, 7, 5, 7, 5, 2, 7, 2, 5, 2, 2, 7, 5, 7, 5, 2, 5, 6, 5, 7, 7, 7, 5, 2, 
	5, 7, 6, 7, 5, 7, 7, 5, 7, 6, 2, 5, 2, 6, 2, 2, 5, 6, 2, 6, 5, 5, 6, 2, 5, 7, 7, 6, 2, 2, 6, 7, 5, 2, 
	5, 2, 6, 2, 2, 6, 5, 6, 7, 7, 2, 5, 2, 0, 5, 2, 2, 6, 7, 2, 7, 2, 2, 2, 2, 7, 2, 7, 5, 7, 7, 2, 5, 7, 
	7, 5, 6, 5, 5, 7, 6, 0, 6, 7, 5, 6, 5, 5, 5, 6, 5, 6, 2, 7, 5, 7, 2, 6, 5, 6, 6, 2, 6, 7, 0, 6, 2, 6, 
	6, 5, 2, 6, 7, 5, 2, 2, 6, 7, 2, 6, 7, 7, 6, 6, 6, 6, 2, 5, 6, 7, 6, 6, 7, 7, 7, 6, 7, 5, 6, 7, 2, 7, 
	6, 5, 7, 5, 6, 5, 6, 7, 6, 6, 5, 7, 7, 2, 2, 6, 5, 2, 2, 7, 2, 7, 0, 5, 5, 6, 6, 5, 6, 2, 7, 6, 2, 7, 
	5, 7, 2, 6, 2, 5, 5, 5, 5, 7, 7, 5, 0, 7, 2, 6, 7, 7, 0, 2, 2, 5, 6, 6, 7, 2, 5, 7, 6, 5, 6, 2, 2, 7, 
	7, 0, 2, 2, 2, 6, 5, 5, 2, 2, 6, 2, 2, 2, 5, 7, 5, 6, 2, 2, 5, 5, 5, 2, 2, 6, 5, 7
});

double euclidean_distance(const double* x, const double* y, int n)
{
	double distance = 0;
	for (int i = 0; i < n; ++i, ++x, ++y)
	{
		double difference = (*x) - (*y);
		distance += difference * difference;
	}

	return std::sqrt(distance);
}

int main() {
	std::cout << "Testing HDBSCAN library...\n";

	std::vector<std::vector<double> > distanceMatrix(NUM_INPUTS);
    auto ndim = m_testDataset.size() / NUM_INPUTS;

    for (int i = 0; i < NUM_INPUTS; i++)
    {
        distanceMatrix[i].resize(NUM_INPUTS);
        for (int j = 0; j < i; j++)
        {
            auto distance = euclidean_distance(&m_testDataset[ndim * i], &m_testDataset[ndim * j], ndim);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
    }
	
	auto coreDistances = HdbscanAlgorithm::CalculateCoreDistances(distanceMatrix, MIN_POINTS);
    auto minimumSpanningTree = HdbscanAlgorithm::ConstructMst(distanceMatrix, coreDistances);
    minimumSpanningTree.QuicksortByEdgeWeight();

    std::vector<double> pointNoiseLevels(NUM_INPUTS);
    std::vector<int> pointLastClusters(NUM_INPUTS);
    std::vector<std::vector<int> > hierarchy;
    std::vector<HdbscanCluster*> clusters;

    HdbscanAlgorithm::ComputeHierarchyAndClusterTree(&minimumSpanningTree, MIN_CLUSTER_SIZE, hierarchy, pointNoiseLevels, pointLastClusters, clusters);
    HdbscanAlgorithm::PropagateTree(clusters);
    std::vector<int> clusterLabels = HdbscanAlgorithm::FindProminentClusters(clusters, hierarchy, NUM_INPUTS);

	std::cout << "Cluster assignments: \n";

	for (int i = 0; i < NUM_INPUTS; i++)
	{
		if (clusterLabels[i] != m_correctClusters[i])
		{
			std::cout << "\nIncorrect cluster found.\n";
			return 1;
		}
		
		std::cout << std::to_string(clusterLabels[i]) << " ";
	}

	return 0;
}
