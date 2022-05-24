#include <iostream>
#include <vector>
#include "math.h"
#include <map>
#include <algorithm>
#include <unordered_set>
#include "algo_utils.hpp"
#include "misc_utils.hpp"
#pragma once

using namespace std;

class dckm_utils{

    template <typename TD>
    void find_midpoints(vector<TD> &center1, vector<TD> &center2, 
    vector<TD> &midpoint);

    template <typename TD>
    void find_affine_vector(vector<TD> &midpoint, vector<TD> &center, 
    vector<TD> &affine);

    template <typename TD>
    bool find_context_direction(vector<TD> &centroid_vector, vector<TD> &midpoint,
    vector<TD> &actual_point, string chk_type);

    template <typename TD, typename TI>
    void restore_radius(vector<vector <TD> > &dist_matrix,
    vector<unordered_set<TI> > &assign_dict, 
    vector<TD> &cluster_radius);

    template <typename TD, typename TI>
    void find_neighbors(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<TD> &cluster_radius, 
    vector<vector<TI> > &neighbors, map<string, vector<TD> > &mid_points, 
    map<string, vector<TD> > &affine_vectors);

    template <typename TD, typename TI>
    void determine_data_expression(vector<vector<TD> > &dataset, 
    vector<vector <TD> > &curr_centroids, 
    vector<vector<TI> > &neighbors, vector<unordered_set<TI>> &assign_dict,
    map<string, vector<TD> > &affine_vectors, 
    map<string, vector<TD> > &mid_points, 
    vector<vector<TI> > &he_data, vector<vector <TD> > &center_dist_mat);

    template <typename TDouble, typename Tint>
    void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
    vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
    vector<unordered_set<Tint>> &assign_dict, 
    vector<TDouble> &cluster_radius,  
    vector<vector <Tint> > &he_data);

    template <typename T1, typename T2>
    void calculate_dckm_distances(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat,
    T2 num_clusters, vector<unordered_set<T2>> &assign_dict, 
    vector<T1> &cluster_radius);

    template <typename TDouble, typename TInt>
    void update_dckm_centroids(vector<vector <TDouble> > &dataset, 
    vector<vector<TDouble> > &new_centroids, 
    vector<unordered_set<TInt> > &assign_dict);

};


template <typename TD>
void find_midpoints(vector<TD> &center1, vector<TD> &center2, 
vector<TD> &midpoint){

    for (int i=0; i<center1.size(); i++){
        midpoint[i] = (center1[i] + center2[i])/2;
    }
}


template <typename TD>
void find_affine_vector(vector<TD> &midpoint, vector<TD> &ot_point, vector<TD> &affine){
    for (int i=0; i<ot_point.size(); i++)
        affine[i] = ot_point[i] - midpoint[i];
}


template <typename TD>
bool find_context_direction(vector<TD> &centroid_vector, vector<TD> &midpoint,
vector<TD> &actual_point, string chk_type){

    int mysize = midpoint.size(); 
    TD vec_sum = 0.0;
    TD temp_holder = 0;
    
    if (chk_type == "validity"){
            for (int i=0; i<mysize; i++){
                temp_holder = actual_point[i] - midpoint[i];
                vec_sum =  vec_sum + (temp_holder * (-centroid_vector[i]));
            }
    }
    
    else{
        for (int i=0; i<mysize; i++){
            temp_holder = actual_point[i] - midpoint[i];
            vec_sum =  vec_sum + (temp_holder * centroid_vector[i]);
        }
    }

    if (vec_sum>0)
        return true;
    
    return false;
}


template <typename TD, typename TI>
void restore_radius(vector<vector <TD> > &dist_matrix,
vector<unordered_set<TI> > &assign_dict, 
vector<TD> &cluster_radius){

    for (int i=0; i<assign_dict.size(); i++){
        for (unordered_set<int>::iterator it = assign_dict[i].begin(); 
                it!=assign_dict[i].end() ; it++){
                if ((dist_matrix[(*it)][i] > cluster_radius[i]))
                    cluster_radius[i] = dist_matrix[(*it)][i];
        }
    }
}


template <typename TD, typename TI>
void find_neighbors(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector <TD> &cluster_radius, 
vector<vector<TI> > &neighbors, map<string, vector<TD> > &mid_points, 
map<string, vector<TD> > &affine_vectors){

    TD dist = 0;
    TD radius = 0;

    algorithm_utils alg_utils;

    vector<TD> temp_midpoint(centroids[0].size());
    vector<TD> temp_affine(centroids[0].size());
    string key = "";
    vector<vector<TD> > temp_master;
    
    // Calculate inter-centroid distances
    for(int curr_center=0; curr_center<centroids.size(); curr_center++){
        
        radius = cluster_radius[curr_center];
        vector<TD> temp;
        vector<TI> temp1;
        
        for (int ot_center=0; ot_center<centroids.size(); ot_center++){    
            
            // Do only k calculation :) save some 
            if (curr_center == ot_center){
                center_dist_mat[curr_center][ot_center] = 0;
            }
            else if (center_dist_mat[curr_center][ot_center] == 0){
                dist = alg_utils.calc_euclidean(centroids[curr_center], centroids[ot_center]);
                center_dist_mat[curr_center][ot_center] = (dist/2);
                center_dist_mat[ot_center][curr_center] = center_dist_mat[curr_center][ot_center];
            }

            if ((curr_center != ot_center) && (center_dist_mat[curr_center][ot_center] < radius)){   
                // The following is the neighbor for the current center
                // neighbors[curr_center].push_back(ot_center); 
                temp.push_back(center_dist_mat[curr_center][ot_center]);
                temp.push_back(ot_center);
                temp_master.push_back(temp);
                temp.clear();
                    
                // Get the mid-point coordinates for this pair of centroids
                find_midpoints(centroids[curr_center], centroids[ot_center], 
                                    temp_midpoint);
                // Determine the affine vector                
                find_affine_vector(temp_midpoint, centroids[ot_center], temp_affine);
                
                // Update the primary containers
                key = std::to_string(ot_center) + std::to_string(curr_center);
                mid_points.insert_or_assign(key, temp_midpoint);
                affine_vectors.insert_or_assign(key, temp_affine);
            }
        }

        if (temp_master.size()>1){
            sort(temp_master.begin(), temp_master.end(), [](const std::vector<TD>& a, const std::vector<TD>& b) {
                return a[0] < b[0];});
            
            for(int i = 0; i<temp_master.size();i++)
                    temp1.push_back(trunc(temp_master[i][1]));
            
            neighbors[curr_center] = temp1;
        }

        else if (temp_master.size() == 1){
            temp1.push_back(temp_master[0][1]);
            neighbors[curr_center] = temp1;
        }
        
        temp1.clear();
        temp_master.clear();
    }
}


// template <typename TD, typename TI>
// void determine_data_expression(vector<vector<TD> > &dataset, 
// vector<vector <TD> > &curr_centroids,
// vector<vector<TI> > &neighbors, vector<unordered_set<TI>> &assign_dict,
// map<string, vector<TD> > &affine_vectors, 
// map<string, vector<TD> > &mid_points, 
// vector<vector<TI> > &he_data){

//     TI data_point = 0;
//     TI curr_cluster = 0;
//     vector<TI> * curr_neighbors;
//     TI closest_nei;
//     string key = "";
//     vector<TI> temp;

//     // For each cluster
//     for (int i=0; i<assign_dict.size(); i++){

//         curr_cluster = i;
//         curr_neighbors = &neighbors[curr_cluster];

//         // For each point in the current cluster
//         for(unordered_set<int>::iterator ref = assign_dict[i].begin();
//         ref != assign_dict[i].end(); ++ref){

//             // Detemine HE-ness
//             data_point = (*ref);
            
//             // If the current point is within the 
//             //shortest radius then it can be ignored.
//             // if ((*curr_neighbors).size()>=1){
                
//             // closest_nei = (*curr_neighbors)[0];
//             // key = std::to_string(closest_nei) + std::to_string(i);

//             // cout << "current neighbor: " << i << " closest neighbor: " << closest_nei << "\n";
                
//             // if(find_context_direction(affine_vectors[key], 
//             //     mid_points[key], dataset[data_point], "validity")){
                    
//                     // cout  << "Not detected: " << data_point << "\t :" << curr_cluster << "\n";
//                     // if (data_point == 67){
//                     //     cout << "Key" << key << "\n";
//                     //     print_vector(dataset[data_point], 2, "data");
//                     //     print_vector(mid_points[key], 2, "mid points");
//                     //     vector<TD> tem1 = curr_centroids[curr_cluster];
//                     //     print_vector(tem1, 2, "center-1");
//                     //     tem1 = curr_centroids[closest_nei];
//                     //     print_vector(tem1, 2, "center-2");
                        
//                     // }
                    
//                 //     continue;
//                 // }
            
//             // else {
                
//                 for (int j=2; j<(*curr_neighbors).size(); j++){            
                    
//                     key = std::to_string((*curr_neighbors)[j]) + std::to_string(curr_cluster);
                    
//                     if(find_context_direction(affine_vectors[key], 
//                     mid_points[key], dataset[data_point], "validity")){
//                         continue;
//                         // break;
//                     }

//                     else {
                    
//                     if(find_context_direction(affine_vectors[key], 
//                     mid_points[key], dataset[data_point], "redundant")){
                        
//                         if (temp.size() == 0){
//                             temp.push_back(data_point);
//                             temp.push_back(curr_cluster);
//                         }
                        
//                         temp.push_back((*curr_neighbors)[j]);
//                         }
//                     }
//                 }

//                 if (temp.size() > 0)
//                     he_data.push_back(temp);
//                 temp.clear();
//             }    

//            // }  
//         }
//     }

template <typename TD, typename TI>
void determine_data_expression(vector<vector<TD> > &dataset, 
vector<vector <TD> > &curr_centroids, vector<vector<TI> > &neighbors, 
vector<unordered_set<TI> > &assign_dict,
map<string, vector<TD> > &affine_vectors, 
map<string, vector<TD> > &mid_points, 
vector<vector<TI> > &he_data, vector<vector <TD> > &center_dist_mat){

    TI data_point = 0;
    TI curr_cluster = 0;
    vector<TI> * curr_neighbors;
    TI closest_nei;
    // bool status;
    string key = "";
    vector<TI> temp;
    TD tempe =0;

    algorithm_utils alg;

    // For each cluster
    for (int i=0; i<assign_dict.size(); i++){

        curr_cluster = i;
        curr_neighbors = &neighbors[curr_cluster];
        // cout << "current cluster: " << curr_cluster << "\n";
        // print_2d_vector(neighbors, neighbors.size(), "Neighbors");

        // For each point in the current cluster
        for(unordered_set<int>::iterator ref = assign_dict[i].begin();
        ref != assign_dict[i].end(); ++ref){

            // Detemine HE-ness
            data_point = (*ref);
            
            // If the current point is within the 
            //shortest radius then it can be ignored.
            // if ((*curr_neighbors).size()>=1){
                
            // closest_nei = (*curr_neighbors)[0];
            
            // key = std::to_string(closest_nei) + std::to_string(i);

            // cout << "current neighbor: " << i << " closest neighbor: " << closest_nei << "\n";
                
            // if(find_context_direction(affine_vectors[key], 
            //     mid_points[key], dataset[data_point], "validity")){
            //         // cout  << "Not detected: " << data_point << "\t :" << curr_cluster << "\n";
                    
            //         if (data_point == 67){
            //             cout << "Key" << key << "\n";
            //             print_vector(dataset[data_point], 2, "data");
            //             print_vector(mid_points[key], 2, "mid points");
            //             vector<TD> tem1 = curr_centroids[curr_cluster];
            //             print_vector(tem1, 2, "center-1");
            //             tem1 = curr_centroids[closest_nei];
            //             print_vector(tem1, 2, "center-2");
                        
            //         }
                    
            //         continue;
            //     }
            
            // else {
                
                for (int j=0; j<(*curr_neighbors).size(); j++){            
                    
                    key = std::to_string((*curr_neighbors)[j]) + std::to_string(curr_cluster);
                    // tempe = alg.calc_euclidean(dataset[data_point], curr_centroids[(*curr_neighbors)[j]]);
                    // if (tempe < center_dist_mat[curr_cluster][(*curr_neighbors)[j]]){
                    //     continue;
                    // }

                    if(find_context_direction(affine_vectors[key], 
                    mid_points[key], dataset[data_point], "validity")){
                        continue;
                        // break;
                    }

                    else {
                    
                    // if(find_context_direction(affine_vectors[key], 
                    // mid_points[key], dataset[data_point], "redundant")){
                        
                        if (temp.size() == 0){
                            temp.push_back(data_point);
                            temp.push_back(curr_cluster);
                        }
                        
                        temp.push_back((*curr_neighbors)[j]);
                        //}
                    }
                }

                if (temp.size() > 0)
                    he_data.push_back(temp);
                temp.clear();
            //} 

            }  
        }
    }




template <typename TDouble, typename Tint>
void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
vector<unordered_set<Tint> > &assign_dict, 
vector<TDouble> &cluster_radius,  
vector<vector <Tint> > &he_data){

    Tint curr_center = 0;
    Tint new_center = 0;
    Tint data_point = 0;
    TDouble temp = 0.0;

    Tint refer_size = 0;
    algorithm_utils alg_utils;

    for (int i=0; i < he_data.size(); i++){

        refer_size = he_data[i].size();
        
        if (refer_size>0){

            data_point = he_data[i][0];
            curr_center = he_data[i][1];

        for (int j = 1; j<refer_size; j++){
            
            new_center = he_data[i][j];
            temp = alg_utils.calc_euclidean(dataset[data_point], centroids[new_center]);
            dist_mat[data_point][new_center] = temp;

            if(temp <= dist_mat[data_point][curr_center]){
                
                dist_mat[data_point][curr_center] = temp;

                cluster_radius[curr_center] = 0.0;
                
                // Move the data in the dict
                assign_dict[curr_center].erase(data_point);
                assign_dict[new_center].insert(data_point);
                curr_center = new_center;
            } 
        }

    }

    }

}


template <typename T1, typename T2>
void calculate_dckm_distances(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat,
T2 num_clusters, vector<unordered_set<T2>> &assign_dict, 
vector<T1> &cluster_radius){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters);
    double temp = 0.0;
    algorithm_utils alg;

    // Calculate the distance of points to nearest center
    for (int i=0; i < dataset.size(); i++){
        
        double shortestDist2 = std::numeric_limits<double>::max();
        
        for (int j=0; j < centroids.size(); j++){ 
            temp = alg.calc_euclidean(dataset[i], centroids[j]);
            temp_dist[j] = temp;
            
            if (temp < shortestDist2){
                shortestDist2 = temp;
                current_center = j;
            }
        }

        dist_mat[i] = temp_dist;
        assign_dict[current_center].insert(i);

        // Store the max so far
        if (shortestDist2 > cluster_radius[current_center])
            cluster_radius[current_center] = shortestDist2;
    }
}


template <typename TDouble, typename TInt>
void update_dckm_centroids(vector<vector <TDouble> > &dataset, 
vector<vector<TDouble> > &new_centroids, 
vector<unordered_set<TInt> > &assign_dict){

    unordered_set<int>::iterator ref ;

    for (TInt i=0; i<assign_dict.size(); i++){
        
        if (assign_dict[i].size() > 0){
            for(ref = assign_dict[i].begin(); ref!=assign_dict[i].end(); ref++){
                for (TInt j=0; j<new_centroids[0].size();j++)
                    new_centroids[i][j] = new_centroids[i][j] + dataset[(*ref)][j];
            }
        }

        for (TInt k=0; k<new_centroids[0].size(); ++k){
            new_centroids[i][k] = new_centroids[i][k]/assign_dict[i].size();
        }
    }
}