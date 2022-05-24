#include <iostream>
#include <vector>
#include "misc_utils.hpp"
#include "algo_utils.hpp"
#include <chrono>

using namespace std;

class conv_kmeans{
    // T1 would be double or float
    // T2 would be an integer
    template <typename Tdouble, typename Tint>
    Tint kmeans(vector<vector <Tdouble> > &, 
    Tint, Tdouble, Tint, Tint , Tint, 
    vector<vector <Tdouble> > &, vector<Tint> &);

};

template <typename Tdouble, typename Tint>
Tint kmeans(vector<vector <Tdouble> > &dataset, 
Tint num_clusters, Tdouble threshold, Tint num_iterations, 
Tint numRows, Tint numCols, vector<vector <Tdouble> > &new_centroids, 
vector<Tint> &assigned_clusters){

    Tint loop_counter = 0;
    vector<vector<Tdouble> > centroids(num_clusters, vector<Tdouble>(numCols));
    vector<vector <Tdouble> > dist_matrix(dataset.size(), vector<Tdouble>(num_clusters));
    vector<Tint> cluster_size(num_clusters);   

    // double centroids[num_clusters][numCols];
    // double dist_matrix[dataset.size()][num_clusters];
 
    // Create objects
    algorithm_utils alg_utils;
    print_utils pu;

    // auto t1 = std::chrono::high_resolution_clock::now();
    
    // Initialize centroids
    alg_utils.init_centroids(centroids, dataset, num_clusters);

    // auto t2 = std::chrono::high_resolution_clock::now();
    // auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    // cout << "Time for centroid init: " << ms.count() << "\n";

    // Print centroids
    // print_2d_vector(centroids, centroids.size(), "Initial Centroids");

    // Assign data to nearest center
    auto t3 = std::chrono::high_resolution_clock::now();

    alg_utils.calculate_distances(dataset, centroids, dist_matrix, 
    num_clusters, assigned_clusters, cluster_size);

    auto t4 = std::chrono::high_resolution_clock::now();
    auto temp1 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
    cout << "Time for first center-point distance calc: " << temp1.count() << "\n";

    // Print dist_matrix
    // print_2d_vector(dist_matrix, 5, "Distance matrix");
    // print_vector(assigned_clusters, 5, "Assigned clusters");
    // print_2d_vector(cluster_size, 5, "Cluster Assignments");

    int cent_int = 0;
    int dist_int = 0;

    while (loop_counter < num_iterations){
        
        loop_counter++;
        // print_map(assign_dict, 5, "Before Erasing");

        auto t5 = std::chrono::high_resolution_clock::now();

        // Calculate new centroids
        alg_utils.update_centroids(dataset, new_centroids, assigned_clusters, 
        cluster_size, numCols);
        
        auto t6 = std::chrono::high_resolution_clock::now();
        auto temp2 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5);
        cent_int = cent_int + temp2.count();

        // print_2d_vector(new_centroids, centroids.size(), "Updated Centroids");
        // print_2d_vector(dist_matrix, 3, "Before: Distance matrix");

        // Check Convergence
        if (alg_utils.check_convergence(new_centroids, centroids, threshold)){
                cout << "Convergence at iteration: " << loop_counter << "\n";
                break;
        }

        // Re-calculate distances
        auto t11 = std::chrono::high_resolution_clock::now();
        alg_utils.calculate_distances(dataset, new_centroids, dist_matrix,
                                    num_clusters, assigned_clusters, cluster_size);

        auto t12 = std::chrono::high_resolution_clock::now();
        auto temp5 = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11);
        dist_int = dist_int + temp5.count();

        // Print for testing
        // print_2d_vector(cluster_size, 5, "After: Cluster Size");
        // print_2d_vector(new_centroids, new_centroids.size(), "After: Updated centroids");

        // Move the new centroids to older
        centroids = new_centroids;
        
        // reset centroids
        alg_utils.reinit(new_centroids);
    
        // break;
    }

    cout << "Total time for centroid updation: " << cent_int << " milliseconds. \n";
    cout << "Total time for distance calc: " << dist_int << " milliseconds. \n";

    // print_2d_vector(new_centroids, new_centroids.size(), "Final centroids");

    return loop_counter;

}
