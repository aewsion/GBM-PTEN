import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

folder = ["p_ten_combine_null"] # or "p_ten_combine_wild"
ll = len(folder)
num = 20
num_start = 0
savedir = './results/index'

if not os.path.exists(savedir):
    os.makedirs(savedir)

def expected_nn_distance(n):
    # A = 10,000
    return 100 / (2 * np.sqrt(n))

def plot_index(index_data, index_type, folder, case_list, n_min_list, sim_num):

    plt.figure(figsize=(10, 6))
    cmap = plt.cm.plasma 
    num_colors = len(folder)
    num_case = len(case_list)
    for i in range(num_case):
        for k in range(num_colors):
            color = cmap(k / num_colors)
            column_mean = np.mean(index_data[k][n_min_list[i]:n_min_list[i]+sim_num], axis=0) 
            plt.plot(column_mean, label=f"{folder[k]} Mean", linestyle='--', color=color)
        
        plt.title(f'{index_type} with {case_list[i]}')
        plt.xlabel('Index')
        plt.ylabel(index_type)
        plt.legend()
        plt.xlim([0, 89])
        savedir_fig = os.path.join(savedir, f'{index_type}')
        if not os.path.exists(savedir_fig):
            os.makedirs(savedir_fig)
        plot_filename = os.path.join(savedir_fig, f'{index_type}_with_{case_list[i]}_plot.pdf')
        plt.savefig(plot_filename, format='pdf')
        plt.close()
        
    
    for k in range(num_colors):
        for i in range(num_case):
            color = cmap(i / num_case)
            column_mean = np.mean(index_data[k][n_min_list[i]:n_min_list[i]+sim_num], axis=0) 
            plt.plot(column_mean, label=f"{case_list[i]}", linestyle='--', color=color)
        
        plt.title(f'{folder[k]}')
        plt.xlabel('Index')
        plt.ylabel(index_type)
        plt.legend()
        plt.xlim([0, 89])
        plot_filename = os.path.join(savedir_fig, f'{index_type}_with_{folder[k]}_plot.pdf')
        plt.savefig(plot_filename, format='pdf')
        plt.close()
        
def calculate_index():    
    for k in range(ll):
        savedir_k = os.path.join(savedir, folder[k].replace(" ", ""))
        
        if not os.path.exists(savedir_k):
            os.makedirs(savedir_k)
        
        nnd_index = np.zeros((num, 200)) 
        nni_index = np.zeros((num, 200)) 
        nnd_index2 = np.zeros((num, 200)) 
    
        for i in range(num):
            fileroad = os.path.join(folder[k], f'set_{i+num_start}', 'tc_num.csv')
            tmp = pd.read_csv(fileroad, header=None)
            
            l = tmp.shape[1]
            print(f'l = {l}')
            if l > 48*90:
                l = 48*90
            
            tc_density = tmp / 2154
            
            for j in range(l):
                if j % 48 == 0:
                    address = os.path.join(folder[k], f'set_{i+num_start}', 'spatial', f'{j}', 'cg.csv')
                    sparse_matrix = np.array(pd.read_csv(address, header=None))
                    
                    # ANND
                    m_matrix = sparse_matrix.copy()
                    tc_matrix = sparse_matrix.copy()
                    tc_matrix[(tc_matrix != 1)] = 0
                    nonzero_index = np.argwhere(tc_matrix != 0)
                    m_matrix[m_matrix != 3] = 0
                    m_index = np.argwhere(m_matrix != 0)
    
                    nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree')
                    nn.fit(m_index)
                    distances, _ = nn.kneighbors(nonzero_index)
                    nnd_dist = distances.flatten()
                    Cell_num = len(nonzero_index)
                    sum_p = np.sum(nnd_dist) / Cell_num
                    nnd_index[i, j // 48] = sum_p
                    
                    # NNI
                    tc_matrix = sparse_matrix.copy()
                    tc_matrix[(tc_matrix != 1) & (tc_matrix != 7)] = 0
                    nonzero_index = np.argwhere(tc_matrix != 0)
                    
                    nn = NearestNeighbors(n_neighbors=2, algorithm='ball_tree')
                    nn.fit(nonzero_index)
                    distances, _ = nn.kneighbors(nonzero_index)
                    nnd_dist = distances[:,1].flatten()
                    #print("nnd_dist2: ", nnd_dist)
                    Cell_num = len(nonzero_index)
                    sum_p = np.sum(nnd_dist) / Cell_num
                    R_e = expected_nn_distance(Cell_num)
                    NNI = sum_p / R_e
                    nni_index[i, j // 48] = NNI 
                    
                    # ANND
                    m_matrix = sparse_matrix.copy()
                    tc_matrix = sparse_matrix.copy()
                    tc_matrix[(tc_matrix != 1)] = 0
                    nonzero_index = np.argwhere(tc_matrix != 0)
                    m_matrix[m_matrix != 4] = 0
                    m_index = np.argwhere(m_matrix != 0)
    
                    nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree')
                    nn.fit(m_index)
                    distances, _ = nn.kneighbors(nonzero_index)
                    nnd_dist = distances.flatten()
                    Cell_num = len(nonzero_index)
                    sum_p = np.sum(nnd_dist) / Cell_num
                    nnd_index2[i, j // 48] = sum_p
                    
        nnd_indexwrite = os.path.join(savedir_k, 'ANND_ATCM1_all_48.csv')
        nni_indexwrite = os.path.join(savedir_k, 'NNI_ATCM1_all_48.csv')
        nnd_index2write = os.path.join(savedir_k, 'ANND_ATCM2_all_48.csv')
        
        pd.DataFrame(nnd_index).to_csv(nnd_indexwrite, header=False, index=False)
        pd.DataFrame(nni_index).to_csv(nni_indexwrite, header=False, index=False)
        pd.DataFrame(nnd_index2).to_csv(nnd_index2write, header=False, index=False)
        
def read_index(folder, file):
    index_list = []
    for k in range(len(folder)):
        readdir_k = os.path.join(savedir, folder[k].replace(" ", ""), file)
        index = pd.read_csv(readdir_k, header=None).to_numpy()
        index_list.append(index)
    return index_list

def main():
    calculate_index()
    nn_index_list = read_index(folder, 'ANND_ATCM1_all_48.csv')
    nni_index_list = read_index(folder, 'NNI_ATCM1_all_48.csv')
    nn_index2_list = read_index(folder, 'ANND_ATCM2_all_48.csv')
    
    case_list = ["CSF1R_I_CSF1R_I", "CSF1R_I_IGF1R_I", 
                    "IGF1R_I_AKT_I", "CSF1R_I_AKT_I"]
                    
    n_min_list = [0, 5, 10, 15]
    sim_num = 5
    plot_index(nn_index_list, 'ANND', folder, case_list, n_min_list, sim_num)
    plot_index(nni_index_list, 'NNI', folder, case_list, n_min_list, sim_num)
    plot_index(nn_index2_list, 'ANND2', folder, case_list, n_min_list, sim_num)
    print("done")
    
if __name__ == "__main__":
    main()
