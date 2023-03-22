import os, sys, shutil, subprocess
import random
from train_GA import train_GA
from sklearn.model_selection import KFold
from distutils.dir_util import copy_tree

def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False


def random_hyper_param(input_file):
    random_param = [random.choice(i) for i in grid.values()]

    os.rename(input_file, 'input_write')

    with open('input_write','r') as fin:
        with open('input_train.ini', 'w') as fout:
            lines = [line.split() for line in fin]
            for i in lines:
                if i[0] == 'num_poses':
                    i[2] = str(random_param[0])
                if i[0] == 'ga_num_mating':
                    i[2] = str(random_param[1])
                if i[0] == 'par_type':
                    i[2] = random_param[2]
                if i[0] == 'keep_par':
                    i[2] = str(random_param[3])
                if i[0] == 'crossover_type':
                    i[2] = random_param[4]
                if i[0] == 'cross_prob':
                    i[2] = str(random_param[5])
                if i[0] == 'mut_prob':
                    i[2] = str(random_param[6])

                
                fout.write(' '.join(i))
                fout.write('\n')

    os.remove('input_write')
    return

def train_to_test(input_file, run_script):
    with open(input_file, 'r') as fin:
        with open('input_test.ini', 'w') as fout:
            lines = [line.split() for line in fin]
            for i in lines:
                if i[0] == 'method':
                    i[2] = 'test'
                
                fout.write(' '.join(i))
                fout.write('\n')

    with open(run_script, 'r') as fin:
        with open('test_run_script', 'w') as fout:
            lines = [line.split() for line in fin]
            for i in lines:
                if i[0] == '/home/matthak/.conda/envs/metal_dock/bin/python3':
                    i[4] = 'input_test.ini'
                
                fout.write(' '.join(i))
                fout.write('\n')
    
    os.remove(input_file)
    os.remove(run_script)
    return

def protein_numbers():
    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    data_set = sorted(dir_list)

    return data_set


if __name__ == '__main__':

    grid = {'num_poses': [5,10,15,20,25,30,35,40,45,50],
            'ga_num_mating': [10,20,30,40,50],
            'par_type': ['sss','rws','sus','rank','random'],
            'keep_par': [-1, 0, 10, 20, 30, 40],
            'crossover_type': ['single_point','two_points','uniform','scattered'],
            'cross_prob': [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
            'mut_prob': [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    }


    n_cvs = 30

    input_dir = os.getcwd()

    os.chdir(f'{input_dir}/data_set')

    data_set = protein_numbers()

    os.chdir(f'..')
    ##### K FOLD SPLIT #####
    kf = KFold(n_splits=4, shuffle=True)

    for i in range(kf.n_splits):
        os.mkdir(f'{input_dir}/fold_{i+1}')
        os.mkdir(f'{input_dir}/fold_{i+1}/train')
        os.mkdir(f'{input_dir}/fold_{i+1}/train/data_set')

        os.mkdir(f'{input_dir}/fold_{i+1}/test')
        os.mkdir(f'{input_dir}/fold_{i+1}/test/data_set')

    count=1

    for train_index, test_index in  kf.split(data_set):
        print("TRAIN:", train_index,"TEST:", test_index)
        for i in train_index:
            copy_tree(f'{input_dir}/data_set/protein_{data_set[i]}', f'{input_dir}/fold_{count}/train/data_set/protein_{data_set[i]}')

        for i in test_index:
            copy_tree(f'{input_dir}/data_set/protein_{data_set[i]}', f'{input_dir}/fold_{count}/test/data_set/protein_{data_set[i]}')

        count+=1

    ##### HYPER PARAMETER COMPARISON ###### 
    os.mkdir('CV_search')
    os.chdir('CV_search')

    for i in range(n_cvs):
        os.chdir(f'{input_dir}/CV_search')
        os.mkdir(f'CV_param_{i+1}')
        os.chdir(f'CV_param_{i+1}')
        shutil.copyfile(f'{input_dir}/input_train.ini', 'input_train.ini')
        random_hyper_param('input_train.ini')
        #train_to_test('input_train.ini', 'train_run_script')

        for j in range(kf.n_splits):
            os.chdir(f'{input_dir}/CV_search/CV_param_{i+1}')
            copy_tree(f'{input_dir}/fold_{j+1}',f'fold_{j+1}')
            os.chdir(f'{input_dir}/CV_search/CV_param_{i+1}/fold_{j+1}/train')
            shutil.copyfile(f'{input_dir}/CV_search/CV_param_{i+1}/input_train.ini', 'input_train.ini')
            shutil.copyfile(f'{input_dir}/train_run_script', 'train_run_script')

            os.chdir(f'{input_dir}/CV_search/CV_param_{i+1}/fold_{j+1}/test')

            os.chdir(f'{input_dir}/CV_search/CV_param_{i+1}/fold_{j+1}')
