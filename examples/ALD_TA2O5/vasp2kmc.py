import numpy as np
from copy import deepcopy
from tqdm import tqdm
import pickle

def run(filename,outfile,strict_distance, stric_z, element_dict):
    data_dic = read(filename)
    data_dic['neighbor_matrix'] = generate_neighbor_matrix(data_dic,strict_distance)
    write(data_dic,outfile,stric_z,element_dict)


def read(filename):
    with open(filename) as f:
        content = f.readlines()
        nrow = 1
        vectors = [] #晶胞参数
        coord = []
        # element_type = []
        for line in content:
            line = line.strip()
            if nrow >= 3 and nrow <= 5:
                vector = [float(num_str) for num_str in line.split()]
                vectors.append(vector)
            if nrow == 6:
                element_type = line.split()
            if nrow == 7:
                element_num = [int(num_str) for num_str in line.split()]
            ## 读取所有的坐标
            if nrow >=9 :
                coord.append([float(num_str) for num_str in line.split()])
            nrow += 1
        fraction_coord = []
        # Z_max = np.max(coord, axis=0)[2] #去掉真空层
        # vectors[2][2] = Z_max
        # 周期性性边界条件
        for point in coord:
            for i in range(0,2):
                if point[i] < 0:
                    point[i] += 1
            fraction_coord.append(deepcopy(point))
    # print(coord)
    return {'vectors':vectors, 'element_type':element_type, 'element_num':element_num, 
            'coord':coord, 'fraction_coord':fraction_coord}

def generate_neighbor_matrix(data_dic,strict_distance):
    vectors = deepcopy(data_dic['vectors'])
    coord = deepcopy(data_dic['fraction_coord'])
    # 计算两个格点之间的距离
    element_num = data_dic['element_num']
    all_neighbor_list = []
    ## 构建两个距离的矩阵
    distance_matrix = np.zeros((element_num[1]+element_num[0],element_num[1]+element_num[0]))
    for i in tqdm(range(0, element_num[0])):
        neighbor_list = []
        for j in range(element_num[0], element_num[1]+element_num[0]):
            distance = cal_distance(deepcopy(coord[i]),deepcopy(coord[j]),deepcopy(vectors))
            if distance < strict_distance:
                distance_matrix[i][j] = 1
                distance_matrix[j][i] = 1
    return distance_matrix


def cal_distance(coord1, coord2, vectors):
    if coord1[0]-coord2[0] > 0.5:
        coord2[0] = coord2[0] + 1
    elif coord2[0]-coord1[0] > 0.5:
        coord1[0] = coord1[0] + 1
    if coord1[1]-coord2[1] > 0.5:
        coord2[1] = coord2[1] + 1
    elif coord2[1]-coord1[1] > 0.5:
        coord1[1] = coord1[1] + 1
    # if coord1[2]-coord2[2] > 0.5:
    #     coord2[2] = coord2[2] + 1
    # elif coord2[2]-coord1[2] > 0.5:
    #     coord1[2] = coord1[2] + 1

    vector_matrix = np.array([vectors[0],vectors[1],vectors[2]])
    coord1_ =  np.dot(vector_matrix,coord1)
    coord2_ = np.dot(vector_matrix,coord2)
    return np.linalg.norm(coord1_ - coord2_)

def write(data_dic, filename,  stric_z, element_dict):
    fraction_coord = data_dic['fraction_coord']
    coord = data_dic['coord']
    data_write = ''
    # Sites
    data_write += 'Sites' + '\n'
    data_write += '\n'
    for i in range(len(fraction_coord)):
        data_write += str(i+1).ljust(3, ' ')  + ' ' + str(fraction_coord[i][0]).ljust(5, ' ')  + ' ' + str(fraction_coord[i][1]).ljust(5, ' ')  + ' ' + str(fraction_coord[i][2]).ljust(5, ' ')  + '\n'
    # Neighbors
    data_write += '\n'
    data_write += 'Neighbors' + '\n'
    data_write += '\n'
    neighbor_matrix = data_dic['neighbor_matrix']
    neighbor_list = []
    for i,row in enumerate(neighbor_matrix):
        col_indices = np.where(row == 1)[0]
        neighbor_list.append(deepcopy(col_indices))
        formatted_elements = [str(i+1).ljust(3, ' ') for i in col_indices]
        data_write += str(i+1) + ' ' + ' '.join(formatted_elements) +  '\n'
    data_dic['neighbor_list'] = neighbor_list
    # Values
    ## 0:Vac, 1:O, 2:Ta
    element_num = data_dic['element_num']
    data_write += '\n'
    data_write += 'Values' + '\n'
    data_write += '\n'
    element_cumulative = [element_num[0]]
    # 使用循环来计算累加和
    for num in element_num[1:]:
        element_cumulative.append(element_cumulative[-1] + num)
    # 元素对应的数字列表
    ele_list = [ element_dict[key] for key in data_dic['element_type'] if key in element_dict]
    for i in range(len(fraction_coord)):
        element_index = -1
        if fraction_coord[i][2] < stric_z:
            for j in range(len(element_cumulative)):
                if i > element_cumulative[j]:
                    element_index = j
                    break
            data_write += str(i+1).ljust(3, ' ') + ' ' + str(ele_list[element_index+1]) + ' ' + str(cal_coordination(data_dic,stric_z,i)).ljust(3, ' ')  + '\n'
        else:
            data_write += str(i+1).ljust(3, ' ') + ' ' + '0' + ' ' + str(cal_coordination(data_dic,stric_z,i)).ljust(3, ' ')  + '\n'

    with open (filename, 'w') as f:
            f.write(data_write)

def cal_coordination(data_dic, stric_z, index):
    coordnation = 0
    print(index)
    print(len(data_dic['neighbor_list']))
    neighbors = data_dic['neighbor_list'][index]
    fraction_coord = data_dic['fraction_coord']
    for i in range(len(neighbors)):
        neighbor = neighbors[i]
        if fraction_coord[neighbor][2] < stric_z:
            coordnation += 1
    return coordnation

if __name__ == '__main__':
    filename = './Ta2o5.vasp'
    ctritical_z = 0.1 #小于ctritical_z填充原子
    outfile = 'data.ald'
    element_dict = {'Ta':2,'O':1} #元素对应的数字表示
    ctritical_distance = 2.5 #近邻列表距离
    run(filename,outfile,ctritical_distance, ctritical_z,element_dict)