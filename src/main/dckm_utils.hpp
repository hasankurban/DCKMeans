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
    vector<TI> &assigned_clusters, 
    vector<vector <TD> > &cluster_size);

    template <typename TD, typename TI>
    void find_neighbors(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<vector <TD> > &cluster_size, 
    vector<vector<TI> > &neighbors, map<string, vector<TD> > &mid_points, 
    map<string, vector<TD> > &affine_vectors);
 
    // template <typename TD, typename TI>
    // void determine_data_expression(vector<vector<TD> > &dataset, 
    // vector<vector <TD> > &centroids, vector<TI> &assigned_clusters, 
    // vector<vector<TI> > &neighbors,
    // map<string, vector<TD> > &affine_vectors, 
    // map<string, vector<TD> > &mid_points, 
    // vector<vector<TI> > &he_data);

    template <typename TD, typename TI>
    void determine_data_expression(vector<vector<TD> > &dataset, 
    vector<vector <TD> > &curr_centroids, vector<TI> &assigned_clusters, 
    vector<vector<TI> > &neighbors, vector<unordered_set<TI>> &assign_dict,
    map<string, vector<TD> > &affine_vectors, 
    map<string, vector<TD> > &mid_points, 
    vector<vector<TI> > &he_data);

    template <typename TD, typename TI>
    void find_HE_data(vector<vector<TD> > &dataset, vector<vector <TD> > &centroids, 
    vector<vector <TD> > &distances, vector<TI> &assigned_clusters,
    vector<vector<TD> > &cluster_size, vector<vector <TD> > &center_dist_mat, 
    map<string, vector<TD> > &mid_points, 
    map<string, vector<TD> > &affine_vectors,
    vector<vector<TI> > &neighbors, 
    vector<vector<TI> > &he_data);

    // template <typename TDouble, typename Tint>
    // void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
    // vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
    // Tint num_clusters, vector<Tint> &assigned_clusters, 
    // vector<vector<TDouble> > &cluster_size, 
    // vector<vector <Tint> > &he_data);

    template <typename TDouble, typename Tint>
    void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
    vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
    Tint num_clusters, vector<Tint> &assigned_clusters,
    vector<unordered_set<Tint>> &assign_dict, 
    vector<vector<TDouble> > &cluster_size,  
    vector<vector <Tint> > &he_data);

    int chk_sign(double &val1);

    template <typename Tint>
    void get_assign_dict(vector<unordered_set<Tint> > &assign_dict, 
    vector<Tint> &assigned_clusters);

    bool chk_validity(vector<double> &mid_points, vector<double> &centroid, 
    vector<double> &actual_point, vector<double> &affine_vec);

};


template <typename TD>
void find_midpoints(vector<TD> &center1, vector<TD> &center2, 
vector<TD> &midpoint){

    // print_vector(center1, center1.size(), "Test1");
    // print_vector(center2, center2.size(), "Test2");

    for (int i=0; i<center1.size(); i++){
        midpoint[i] = (center1[i] + center2[i])/2;
        // midpoint[i] = midpoint[i]/2;
    }

    // print_vector(midpoint, midpoint.size(), "Test3");
}


template <typename TD>
void find_affine_vector(vector<TD> &midpoint, vector<TD> &ot_point, vector<TD> &affine){
    for (int i=0; i<ot_point.size(); i++)
        affine[i] = ot_point[i] - midpoint[i];
}


int chk_sign(double &val1){

    int sign = 0;
    if (val1>0)
        sign = 1;
    else if (val1<0)
        sign = -1;
    return sign;
}


template <typename TD>
bool find_context_direction(vector<TD> &centroid_vector, vector<TD> &midpoint,
vector<TD> &actual_point, string chk_type){

    int mysize = midpoint.size(); 
    vector<TD> cent_point_vec(mysize);
    TD vec_sum = 0.0;
    
    find_affine_vector(midpoint, actual_point, cent_point_vec);

    if (chk_type == "validity"){
            for (int i=0; i<mysize; i++)
                vec_sum =  vec_sum + (cent_point_vec[i] * (-centroid_vector[i]));
    }
    
    else{
        for (int i=0; i<mysize; i++)
            vec_sum =  vec_sum + (cent_point_vec[i] * centroid_vector[i]);
    }

    if (vec_sum>0)
        return true;

    // cout << "Inner product: " << vec_sum << "\n";
    // print_vector(cent_point_vec, cent_point_vec.size(), "point vector");
    // print_vector(centroid_vector, cent_point_vec.size(), "Centroid vector");
    
    return false;
    // return status;
}


template <typename TD, typename TI>
void restore_radius(vector<vector <TD> > &dist_matrix,
vector<TI> &assigned_clusters, 
vector<vector <TD> > &cluster_size){

    for (int i=0; i<cluster_size.size(); i++){
        for (int j=0; j< assigned_clusters.size();j ++){
                if ((assigned_clusters[j] == i) && (dist_matrix[j][i] > cluster_size[i][1]))
                    cluster_size[i][1] = dist_matrix[j][i];
        }
    }
}


template <typename TD, typename TI>
void find_neighbors(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector<vector <TD> > &cluster_size, 
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
        
        radius = cluster_size[curr_center][1];
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

            // Start neighbor finding
            // if (center_dist_mat[curr_center][ot_center] < radius){

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
                // mid_points[std::to_string(curr_center) + std::to_string(ot_center)] = temp_midpoint;
                // affine_vectors[std::to_string(curr_center) + std::to_string(ot_center)] = temp_affine;
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
// vector<vector <TD> > &centroids, vector<TI> &assigned_clusters, 
// map<TI, vector<vector<TI> > > &neighbors,
// map<string, vector<TD> > &affine_vectors, 
// map<string, vector<TD> > &mid_points, 
// vector<vector<TI> > &he_data){
    
//     TI my_cluster = 0;
//     string key = "";
//     bool status = false;
//     map<TI, vector<vector<TI> > >::iterator it;
//     vector<vector<TI> > &mapit;
//     TI my_nei;

//     for (int i = 0; i < dataset.size(); i++){

//         vector<TI> temp;
//         my_cluster = assigned_clusters[i];
//         status = false;

//         it = neighbors[my_cluster].begin();
//         mapit = (*it).second;
        
//         if (mapit.size()!=0){
            
//             my_nei = mapit[0][1];
//             key = std::to_string(my_nei) + std::to_string(my_cluster);
            
//             if(find_context_direction(affine_vectors[key], 
//                 mid_points[key], dataset[i], "validity")){
//                 continue;
//             }
//         }

//         // for (int j=0; j<neighbors[my_cluster].size(); j++){

//         for (it = neighbors[my_cluster].begin(); it != neighbors[my_cluster].end(); it++){
            
//             key = std::to_string(neighbors[my_cluster][j]) + std::to_string(my_cluster);
            
//             if ( (my_cluster != neighbors[my_cluster][j]) && find_context_direction(affine_vectors[key], 
//             mid_points[key], dataset[i], "redundant")){
                
//                 // cout << "Same direction" << "\n" ;
//                 if (temp.size() == 0)
//                     temp.push_back(i);
//                     temp.push_back(my_cluster);
                
//                 temp.push_back(neighbors[my_cluster][j]);
//                 status = true;
//             } 
//         }

//         if (status){
//             he_data.push_back(temp);
//         }
//     }
// }


// template <typename TD, typename TI>
// void determine_data_expression(vector<vector<TD> > &dataset, 
// vector<vector <TD> > &centroids, vector<TI> &assigned_clusters, 
// vector<vector<TI> >  &neighbors,
// map<string, vector<TD> > &affine_vectors, 
// map<string, vector<TD> > &mid_points, 
// vector<vector<TI> > &he_data){
    
//     TI my_cluster = 0;
//     string key = "";
//     bool status = false;
//     vector<TI> * curr_neighbors;
//     TI my_nei;
//     vector<TI> temp;

//     for (int i = 0; i < dataset.size(); i++){

//         my_cluster = assigned_clusters[i];
//         status = false;

//         curr_neighbors = &neighbors[my_cluster];
        
//         if ((*curr_neighbors).size()>1){
            
//             my_nei = (*curr_neighbors)[0];
//             key = std::to_string(my_nei) + std::to_string(my_cluster);
            
//             if(find_context_direction(affine_vectors[key], 
//                 mid_points[key], dataset[i], "validity")){
//                 continue;
//             }

//         for (int j=0; j<(*curr_neighbors).size(); j++){            
            
//             key = std::to_string((*curr_neighbors)[j]) + std::to_string(my_cluster);

//             // if ( (my_cluster != neighbors[my_cluster][j]) && find_context_direction(affine_vectors[key], 
//             // mid_points[key], dataset[i], "redundant")){
            
//             if (find_context_direction(affine_vectors[key], 
//             mid_points[key], dataset[i], "redundant")){
                
//                 // cout << "Same direction" << "\n" ;
//                 if (temp.size() == 0)
//                     temp.push_back(i);
//                     temp.push_back(my_cluster);
                
//                 temp.push_back((*curr_neighbors)[j]);
//                 status = true;
//             } 
//         }

//         if (status)
//             he_data.push_back(temp);
//             temp.clear();
//         }
    
//     } 
// }


template <typename TD, typename TI>
void determine_data_expression(vector<vector<TD> > &dataset, 
vector<vector <TD> > &curr_centroids, vector<TI> &assigned_clusters, 
vector<vector<TI> > &neighbors, vector<unordered_set<TI>> &assign_dict,
map<string, vector<TD> > &affine_vectors, 
map<string, vector<TD> > &mid_points, 
vector<vector<TI> > &he_data){

    TI data_point = 0;
    TI curr_cluster = 0;
    vector<TI> * curr_neighbors;
    TI closest_nei;
    // bool status;
    string key = "";
    vector<TI> temp;

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
                    
                    if(find_context_direction(affine_vectors[key], 
                    mid_points[key], dataset[data_point], "validity")){
                        continue;
                    }

                    else {
                    
                    if(find_context_direction(affine_vectors[key], 
                    mid_points[key], dataset[data_point], "redundant")){
                        
                        if (temp.size() == 0){
                            temp.push_back(data_point);
                            temp.push_back(curr_cluster);
                        }
                        
                        temp.push_back((*curr_neighbors)[j]);
                        }
                    }
                }

                if (temp.size() > 0)
                    he_data.push_back(temp);
                temp.clear();
            //} 

            }  
        }
    }


template <typename TD, typename TI>
void find_HE_data(vector<vector<TD> > &dataset, vector<vector <TD> > &centroids, 
vector<vector <TD> > &distances, vector<TI> &assigned_clusters,
vector<vector<TD> > &cluster_size, vector<vector <TD> > &center_dist_mat, 
map<string, vector<TD> > &mid_points, 
map<string, vector<TD> > &affine_vectors,
vector<vector<TD> > &neighbors, 
vector<vector<TI> > &he_data){

        // Find neighbors
        find_neighbors(centroids, center_dist_mat, cluster_size, 
        neighbors, mid_points, affine_vectors);

        // Calculate data expression
        determine_data_expression(dataset, centroids, assigned_clusters, 
        neighbors, affine_vectors, mid_points, he_data);

        // cout << "Number of HE points: " << he_data.size() << "\n";
}


// template <typename TDouble, typename Tint>
// void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
// vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
// Tint num_clusters, vector<Tint> &assigned_clusters, 
// vector<vector<TDouble> > &cluster_size,  
// vector<vector <Tint> > &he_data){

//     Tint current_center = 0;
//     vector<TDouble> temp_dist (num_clusters);
//     TDouble temp = 0.0;

//     // vector<Tint> * refer;
//     Tint refer_size = 0;
//     algorithm_utils alg_utils;

//     // print_2d_vector(neighbors, neighbors.size(), "Neighbors");

//     for (int i=0; i < he_data.size(); i++){

//         // refer = &he_data[i];
//         refer_size = he_data[i].size();

//         for (int j = 1; j<refer_size; j++){
//             temp = alg_utils.calc_euclidean(dataset[he_data[i][0]], centroids[he_data[i][j]]);
//             dist_mat[he_data[i][0]][he_data[i][j]] = temp;

//             // cout << "Point: " << he_data[i][0] << " previous center" << assigned_clusters[he_data[i][0]] << " curr calcu: " << he_data[i][j] << "\n";
//             // cout << "previous dist: " << dist_mat[he_data[i][0]][assigned_clusters[he_data[i][0]]] << 
//             // "curr calc: " << temp << "\n";
            
//             // cout << he_data[i][0] << "<-->" << he_data[i][j] << ": " << temp << "\n";    

//             if(temp < dist_mat[he_data[i][0]][assigned_clusters[he_data[i][0]]]){
//             // if(temp < dist_mat[he_data[i][0]][he_data[i][j]]){
//                 dist_mat[he_data[i][0]][he_data[i][j]] = temp;
//                 cluster_size[assigned_clusters[he_data[i][0]]][0] = cluster_size[assigned_clusters[he_data[i][0]]][0] - 1;
                
//                 // cout << "Point: " << he_data[i][0] << "\t" << assigned_clusters[he_data[i][0]] << "-->" << he_data[i][j] << "\n";
                
//                 cluster_size[assigned_clusters[he_data[i][0]]][1] = 0.0;

//                 // Update the dict
//                 // assign_dict[assigned_clusters[he_data[i][0]]].de

//                 assigned_clusters[he_data[i][0]] = he_data[i][j];
//                 cluster_size[he_data[i][j]][0] = cluster_size[he_data[i][j]][0] + 1; 
//             } 
//         }
//     }

// }


template <typename TDouble, typename Tint>
void calculate_HE_distances(const vector<vector<TDouble> > &dataset, 
vector<vector<TDouble> > &centroids, vector<vector<TDouble> > &dist_mat,
Tint num_clusters, vector<Tint> &assigned_clusters,
vector<unordered_set<Tint>> &assign_dict, 
vector<vector<TDouble> > &cluster_size,  
vector<vector <Tint> > &he_data){

    Tint new_center = 0;
    Tint curr_center = 0;
    Tint data_point = 0;
    vector<TDouble> temp_dist (num_clusters);
    TDouble temp = 0.0;

    Tint refer_size = 0;
    algorithm_utils alg_utils;

    for (int i=0; i < he_data.size(); i++){

        refer_size = he_data[i].size();
        
        if (refer_size>0){

            data_point = he_data[i][0];
            curr_center = assigned_clusters[data_point];

        for (int j = 1; j<refer_size; j++){
            
            new_center = he_data[i][j];

            temp = alg_utils.calc_euclidean(dataset[data_point], centroids[new_center]);
            dist_mat[data_point][new_center] = temp;

            if(temp < dist_mat[data_point][assigned_clusters[data_point]]){
                
                dist_mat[data_point][new_center] = temp;

                cluster_size[assigned_clusters[data_point]][0] = cluster_size[assigned_clusters[data_point]][0] - 1; 
                cluster_size[assigned_clusters[data_point]][1] = 0.0;

                assigned_clusters[data_point] = new_center;
                cluster_size[new_center][0] = cluster_size[new_center][0] + 1;

                // Move the data in the dict
                assign_dict[curr_center].erase(data_point);
                assign_dict[new_center].insert(data_point);
            } 
        }

    }

    }

}


template <typename Tint>
void get_assign_dict(vector<unordered_set<Tint> > &assign_dict, 
vector<Tint> &assigned_clusters){

    for (int i=0; i<assigned_clusters.size(); i++){
        assign_dict[assigned_clusters[i]].insert(i);
    }
}