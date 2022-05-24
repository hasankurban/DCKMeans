from utils.kmeans_utils import *
from utils.vis_utils import *


def DCKMeans(data, num_clusters, threshold, num_iterations, seed):

    loop_counter = 0
    assign_dict = {}
    dist_mat = np.zeros((num_clusters, num_clusters))

    centroids = init_centroids(data, num_clusters, seed)
    # new_centroids = np.zeros(centroids.shape)

    # Calculate the cluster assignments for data points
    assigned_clusters, distances = calculate_distances(data, centroids)
    dckm_calc = num_clusters * data.shape[0]

    print("Initial centroids: ", centroids, "\n")
    # print(assigned_clusters[0:5])

    while loop_counter < num_iterations:

        # Re-calculate the centroids
        new_centroids = calculate_centroids(data, assigned_clusters)

        if check_convergence(new_centroids, centroids, threshold):
            print("Kmeans: Convergence at iteration: ", loop_counter)
            break

        # print(loop_counter, " Updated centroids: ", new_centroids)

        assign_dict = get_membership(assigned_clusters, assign_dict, num_clusters)

        # for i in assign_dict:
        #     print(i, ": ", len(assign_dict[i]))

        neighbors, he_indices_dict = find_all_he_indices_neighbor(data, new_centroids, distances,
                                                assign_dict, dist_mat)

        # temp = []
        # for i in he_indices_dict.keys():
        #     print("iter: ", loop_counter, "\t center: ", i, "HE datas: ", len(he_indices_dict[i]),
        #           he_indices_dict[i])

        for center in neighbors:
            if len(he_indices_dict[center]) > 0:
                # print("Neighbors for ", center, ": ", neighbors[center])
                # print(loop_counter, ": HE points for cluster: ", center, len(he_indices_dict[center]))
                temp, dist123 = calculate_distances_specific(data[he_indices_dict[center]],
                                             new_centroids[neighbors[center]], neighbors[center])

                # for index in range(len(he_indices_dict[center])):
                #     point = he_indices_dict[center][index]
                #     print("Point: ", point, "\t", assigned_clusters[point], "-->", temp[index],
                #           " old dist: ", distances[point], " new distance: ", dist123[index])

                distances[he_indices_dict[center]] = dist123
                assigned_clusters[he_indices_dict[center]] = temp
                dckm_calc += len(he_indices_dict[center]) * len(neighbors[center])

        # Calculate the cluster assignments for data points
        centroids[:] = new_centroids[:]

        loop_counter += 1

    print("KMeans exiting at: ", loop_counter, " iterations")
    return new_centroids, loop_counter, dckm_calc






