import sys
sys.path.append("../")

from utils.dataIO import *
from base.KMeans import *
from base.DCKmeans import *
from pathlib import Path
import time
from utils.assign_clusters import *
import pandas as pd

'''
Algo1: KMeans
Algo2: DCKMeans
'''

# Set parameters
threshold = 0.001
num_iterations = 1

# file_list = ['test_data_case1.csv']
file_list = ['test_100_2_3.csv']
# file_list = ['crop.csv']
# file_list = ['magic.csv']
# file_list = ['user_knowledge_train.csv']
# file_list = ['hapt_train.csv']
# file_list = ['covertype.csv']
# file_list = ['spambase.csv']
# file_list = ['10_clusters.csv']

DATA_PATH = "/Users/schmuck/Documents/Box Sync/Ph.D./DATASETS"

# Make changes for adjusting the current directory here
file_path = os.path.join(DATA_PATH, "benchmark", "clustering_data")
# file_path = os.path.join(Path(__file__).parents[1], "benchmark", "scal_data")
file_path = os.path.join(DATA_PATH, "real_data")
file_path = os.path.join(DATA_PATH, "sample_data")

num_clusters = 3
seed = 12

for data_file in file_list:

    data, labels = read_simulated_data(os.path.join(file_path, data_file))
    # data = np.array(pd.read_csv(os.path.join(file_path, data_file), sep=" "))
    # data, labels = read_real_data(os.path.join(file_path, data_file), ["user_knowledge_train.csv"])
    # vis_PCA(data, labels)
    # exit(0)

    print(data.shape)

    # km_start_time = time.time()
    # km_centroids, km_iter = Kmeans(data, num_clusters, threshold, num_iterations, seed)
    # km_TraningTime = round(time.time() - km_start_time, 5)

    kmlb_start_time = time.time()
    # _, _, _ = DCKMeans(data, num_clusters, threshold, num_iterations, seed)
    kmlb_centroids, kmlb_iter, dckm_calc = DCKMeans(data, num_clusters, threshold, num_iterations, seed)
    kmlb_TraningTime = round(time.time() - kmlb_start_time, 2)

    # print(kmlb_centroids)
    # print("Distance calculations by KMeans: ", num_clusters*data.shape[0]*km_iter)
    # print("Distance calculations by DCKMeans: ", dckm_calc)

    # acc, new_labels1 = calc_raw_accuracy(labels, assign1, data)
    # print("KM: ", acc, check_ARI(new_labels1, labels))
    #
    # acc, new_labels2 = calc_raw_accuracy(labels, assign2, data)
    # print("KMLB: ", acc, check_ARI(new_labels2, labels))
    #
    # print("Diff in clustering: ", len(np.where(new_labels1 != new_labels2)[0]))
    # print("Diff in clustering: ", len(np.where(assign1 != assign2)[0]))
    # print(assign1[np.where(assign1 != assign2)[0]], assign2[np.where(assign1 != assign2)[0]])

    # print(km_TraningTime, kmlb_TraningTime)
    # print(km_TraningTime)
    # print(kmlb_TraningTime)
    # print("Dev: ", round(np.sqrt(np.mean(np.square(km_centroids - kmlb_centroids))), 3))
    # print(km_centroids)