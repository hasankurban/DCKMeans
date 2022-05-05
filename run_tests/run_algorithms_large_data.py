from utils.assign_clusters import calc_raw_accuracy
from base.expr_DEKmeans import *
from base.kmeans import *
import time
from sklearn.cluster import KMeans as km


def run_algorithms(data, labels, dataset_name, result_dictionary, num_iterations,
                   epsilon, num_clusters, i_seed):

    algorithms = ['KMeans-SciPy', 'DEKMeans-LB', 'DEKMeans-PROB', 'DEKMeans-STO']


    ###########
    # Conventional KMeans
    ###########

    km_start_time = time.time()
    # km_centroids, km_iter = Kmeans(data, num_clusters, epsilon, num_iterations, [], False, i_seed)
    km_centroids, km_iter = km_clustering(data, num_clusters, num_iterations, epsilon, i_seed)
    centroids2 = init_centroids(data, num_clusters, i_seed)
    res = km(n_clusters=num_clusters, init=centroids2, tol=epsilon).fit(data)
    km_centroids = res.cluster_centers_
    km_iter = res.n_iter_
    km_TraningTime = round(time.time() - km_start_time, 2)

    km_acc, _ = calc_raw_accuracy(labels, pred_membership(data, km_centroids), data)
    km_ari = check_ARI(pred_membership(data, km_centroids), labels)
    km_amis = check_amis(pred_membership(data, km_centroids), labels)

    result_dictionary['Accuracy'].append(km_acc)
    result_dictionary['ARI'].append(km_ari)
    result_dictionary['AMI'].append(km_amis)
    result_dictionary['Deviation'].append(0)

    result_dictionary['Runtime'].append(km_TraningTime)
    result_dictionary['Num_iterations'].append(km_iter)
    result_dictionary['Dataset'].append(dataset_name)
    result_dictionary['Num_clusters'].append(num_clusters)


    ###########
    # DEKMeans Tree
    ###########

    # kmtree_start_time = time.time()
    # kmtree_centroids, kmtree_iter = DEKmeans_tree(data, num_clusters, data_threshold,
    #                                               num_iterations, i_seed)
    # kmtree_TraningTime = round(time.time() - kmtree_start_time, 2)
    #
    # kmtree_acc, _ = calc_raw_accuracy(labels, pred_membership(data, kmtree_centroids), data)
    # kmtree_ari = check_ARI(pred_membership(data, kmtree_centroids), labels)
    # kmtree_amis = check_amis(pred_membership(data, kmtree_centroids), labels)
    # kmtree_dev = round(np.sqrt(np.mean(np.square(km_centroids - kmtree_centroids))), 3)
    #
    # result_dictionary['Accuracy'].append(kmtree_acc)
    # result_dictionary['ARI'].append(kmtree_ari)
    # result_dictionary['AMI'].append(kmtree_amis)
    # result_dictionary['Deviation'].append(kmtree_dev)
    #
    # result_dictionary['Runtime'].append(kmtree_TraningTime)
    # result_dictionary['Num_iterations'].append(kmtree_iter)
    # result_dictionary['Dataset'].append(dataset_name)
    # result_dictionary['Num_clusters'].append(num_clusters)


    #$#########
    # DEKM Look back
    ###########

    kmlb_start_time = time.time()
    kmlb_centroids, kmlb_iter = DEKmeans_lb(data, num_clusters, num_iterations, i_seed)
    kmlb_TraningTime = round(time.time() - kmlb_start_time, 2)

    result_dictionary['Runtime'].append(kmlb_TraningTime)
    result_dictionary['Num_iterations'].append(kmlb_iter)
    result_dictionary['Dataset'].append(dataset_name)
    result_dictionary['Num_clusters'].append(num_clusters)

    kmlb_acc, _ = calc_raw_accuracy(labels, pred_membership(data, kmlb_centroids), data)
    kmlb_ari = check_ARI(pred_membership(data, kmlb_centroids), labels)
    kmlb_amis = check_amis(pred_membership(data, kmlb_centroids), labels)
    kmlb_dev = round(np.sqrt(np.mean(np.square(km_centroids - kmlb_centroids))), 3)

    result_dictionary['Accuracy'].append(kmlb_acc)
    result_dictionary['ARI'].append(kmlb_ari)
    result_dictionary['AMI'].append(kmlb_amis)
    result_dictionary['Deviation'].append(kmlb_dev)


    ###########
    # DEKM look back probabilistic
    ###########

    kmprob_start_time = time.time()
    kmprob_centroids, kmprob_iter = DEKmeans_prob(data, num_clusters, num_iterations, epsilon, i_seed)
    kmprob_TraningTime = round(time.time() - kmprob_start_time, 2)

    result_dictionary['Runtime'].append(kmprob_TraningTime)
    result_dictionary['Num_iterations'].append(kmprob_iter)
    result_dictionary['Dataset'].append(dataset_name)
    result_dictionary['Num_clusters'].append(num_clusters)

    kmprob_acc, _ = calc_raw_accuracy(labels, pred_membership(data, kmprob_centroids), data)
    kmprob_ari = check_ARI(pred_membership(data, kmprob_centroids), labels)
    kmprob_amis = check_amis(pred_membership(data, kmprob_centroids), labels)
    kmprob_dev = round(np.sqrt(np.mean(np.square(km_centroids - kmprob_centroids))), 3)

    result_dictionary['Accuracy'].append(kmprob_acc)
    result_dictionary['ARI'].append(kmprob_ari)
    result_dictionary['AMI'].append(kmprob_amis)
    result_dictionary['Deviation'].append(kmprob_dev)


    ###########
    # DEKM look back stochastic
    ###########

    kmsto_start_time = time.time()
    kmsto_centroids, kmsto_iter = DEKmeans_sto(data, num_clusters, num_iterations, i_seed)
    kmsto_TraningTime = round(time.time() - kmsto_start_time, 2)

    result_dictionary['Runtime'].append(kmsto_TraningTime)
    result_dictionary['Num_iterations'].append(kmsto_iter)
    result_dictionary['Dataset'].append(dataset_name)
    result_dictionary['Num_clusters'].append(num_clusters)

    kmsto_acc, _ = calc_raw_accuracy(labels, pred_membership(data, kmsto_centroids), data)
    kmsto_ari = check_ARI(pred_membership(data, kmsto_centroids), labels)
    kmsto_amis = check_amis(pred_membership(data, kmsto_centroids), labels)
    kmsto_dev = round(np.sqrt(np.mean(np.square(km_centroids - kmsto_centroids))), 3)

    result_dictionary['Accuracy'].append(kmsto_acc)
    result_dictionary['ARI'].append(kmsto_ari)
    result_dictionary['AMI'].append(kmsto_amis)
    result_dictionary['Deviation'].append(kmsto_dev)

    result_dictionary['Algorithm'] += algorithms

    return result_dictionary
