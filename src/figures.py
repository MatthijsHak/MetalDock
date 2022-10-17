import statistics

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def plot_each_protein(protein_numbers, output_dir):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    first_gen = []
    last_gen = []

    protein_number = protein_numbers

    first_gen_dict = {}
    last_gen_dict = {}

    first_protein_rmsd_dict = {}
    last_protein_rmsd_dict = {}

    for x in range(0,len(protein_number)):
        first_gen_dict[x] = protein_number[x]
        first_protein_rmsd_dict[protein_number[x]] = []
        
        first_gen_dict[x] = f'{output_dir}/first_gen/protein_'+str(first_gen_dict[x])
        
        protein = open(first_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            first_gen.append(float(protein_rmsd[i][5]))
            first_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))
            
        last_gen_dict[x] = protein_number[x]
        last_protein_rmsd_dict[protein_number[x]] = []
        
        last_gen_dict[x] = f'{output_dir}/last_gen/protein_'+str(last_gen_dict[x])
        
        protein = open(last_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            last_gen.append(float(protein_rmsd[i][5]))
            last_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))

        mu_first = statistics.mean(first_protein_rmsd_dict[protein_number[x]])
        sigma_first = statistics.stdev(first_protein_rmsd_dict[protein_number[x]])
        
        mu_last = statistics.mean(last_protein_rmsd_dict[protein_number[x]])
        sigma_last = statistics.stdev(last_protein_rmsd_dict[protein_number[x]])
        
        label_first = [f'First Generation: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
        label_last = [f'Last Generation: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

        bins = np.linspace(0, 15, 100)

        fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)

        ax.hist(first_protein_rmsd_dict[protein_number[x]], density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
        ax.hist(last_protein_rmsd_dict[protein_number[x]], density=False, color='r', alpha=0.5, bins =bins, label=label_last)

        ax.set_title('Protein {}'.format(protein_number[x]))

        ax.set_xticks(range(0,15,2))
        
        ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

        ax.set_xlim(0,15)
        ax.set_ylabel('Conformations (N)')
        ax.set_xlabel('RMSD (Å)')

        plt.legend(loc='upper right')
    
        plt.tight_layout()
        plt.savefig(f"protein_"+str(protein_number[x])+".png", bbox_inches='tight')
        plt.close()

def plot_total_conformations(output_dir):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    first_gen_all = []
    last_gen_all = []

    first_gen = open(f'{output_dir}/first_gen/all_conf_first_gen')
    first_rmsd = [line.split() for line in first_gen]
        
    for i in range(0,len(first_rmsd)):
        first_gen_all.append(float(first_rmsd[i][5]))

    last_gen = open(f'{output_dir}/last_gen/all_conf_last_gen')
    last_rmsd = [line.split() for line in last_gen]
        
    for i in range(0,len(last_rmsd)):
        last_gen_all.append(float(last_rmsd[i][5]))

    mu_first = statistics.mean(first_gen_all)
    sigma_first = statistics.stdev(first_gen_all)
    
    mu_last = statistics.mean(last_gen_all)
    sigma_last = statistics.stdev(last_gen_all)
    
    label_first = [f'First Generation: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
    label_last = [f'Last Generation: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

    bins = np.linspace(0, 15, 100)

    fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)
    
    ax.hist(first_gen_all, density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
    ax.hist(last_gen_all, density=False, color='r', alpha=0.5, bins =bins, label=label_last)
    
    ax.legend(loc="upper right")
    ax.set_title('All Conformations')
    
    ax.set_xlabel('RMSD (Å)')
    ax.set_ylabel('Conformations(N)')
    
    ax.set_xticks(range(0,15,2))

    ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"all_conformations.png", bbox_inches='tight')
    plt.close()

def plot_parameters(metal, output_dir):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    data_set = open(f'{output_dir}/parameter_history','r')
    data_set_line = [ lines.split() for lines in data_set]

    fitness = []
    rmsd_min = []

    r_OA = []
    e_OA = []

    r_SA = []
    e_SA = []

    r_HD = []
    e_HD = []

    r_NA = []
    e_NA = []

    r_N = []
    e_N = []

    r_M = []
    e_M = []
    

    for i in range(1,len(data_set_line)):
        r_OA.append(float(data_set_line[i][4]))
        e_OA.append(float(data_set_line[i][5]))
        
        r_SA.append(float(data_set_line[i][6]))
        e_SA.append(float(data_set_line[i][7]))
        
        r_HD.append(float(data_set_line[i][8]))
        e_HD.append(float(data_set_line[i][9]))
        
        r_NA.append(float(data_set_line[i][8]))
        e_NA.append(float(data_set_line[i][9]))
        
        r_N.append(float(data_set_line[i][10]))
        e_N.append(float(data_set_line[i][11]))
        
        r_M.append(float(data_set_line[i][12]))
        e_M.append(float(data_set_line[i][13]))
        
        fitness.append(float(data_set_line[i][16]))
        rmsd_min.append(float(data_set_line[i][18]))

    x_axis = range(0,len(r_OA))

    fig, ax = plt.subplots(6,2, figsize=(15,15), sharex=True)

    ax[0][0].plot(x_axis,r_OA, c='b')
    ax[0][1].plot(x_axis,e_OA, c='b')

    ax[1][0].plot(x_axis,r_SA, c='b')
    ax[1][1].plot(x_axis,e_SA, c='b')

    ax[2][0].plot(x_axis,r_HD, c='b')
    ax[2][1].plot(x_axis,e_HD, c='b')

    ax[2][0].plot(x_axis,r_NA, c='b')
    ax[2][1].plot(x_axis,e_NA, c='b')

    ax[3][0].plot(x_axis,r_N, c='b')
    ax[3][1].plot(x_axis,e_N, c='b')

    ax[4][0].plot(x_axis,r_N, c='b')
    ax[4][1].plot(x_axis,e_N, c='b')

    ax[5][0].plot(x_axis,r_M, c='b')
    ax[5][1].plot(x_axis,e_M, c='b')

    #########################################################
    ax[0][0].set_ylabel('r {}-OA (Å)'.format(metal))
    ax[0][1].set_ylabel('$\epsilon$ {}-OA (kcal/mol)'.format(metal))

    ax[1][0].set_ylabel('r {}-SA (Å)'.format(metal))
    ax[1][1].set_ylabel('$\epsilon$ {}-SA (kcal/mol)'.format(metal))

    ax[2][0].set_ylabel('r {}-HD (Å)'.format(metal))
    ax[2][1].set_ylabel('$\epsilon$ {}-HD (kcal/mol)'.format(metal))

    ax[3][0].set_ylabel('r {}-NA (Å)'.format(metal))
    ax[3][1].set_ylabel('$\epsilon$ {}-NA (kcal/mol)'.format(metal))

    ax[4][0].set_ylabel('r {}-N (Å)'.format(metal))
    ax[4][1].set_ylabel('$\epsilon$ {}-N (kcal/mol)'.format(metal))
        
    ax[5][0].set_ylabel('r {}-H-bond (Å)'.format(metal,metal))
    ax[5][1].set_ylabel('$\epsilon$ {}-H-bond (kcal/mol)'.format(metal,metal))

    ax[5][0].set_xlabel('N parents')
    ax[5][1].set_xlabel('N parents')

    # for i in range(0,6):
    #     for j in range(0,2):
    #         ax[i][j].set_xlim(0,500)

    plt.tight_layout()
    plt.savefig(f"parameters.png", bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(1,2, figsize=(15,5), sharex=True)

    ax[0].plot(x_axis,fitness, c='b', label='Data set')

    ax[1].plot(x_axis,rmsd_min, c='b', label='{}'.format(metal))

    ax[0].set_ylabel('Fitness Function')
    ax[1].set_ylabel('RMSD (Å)')

    ax[0].set_xlabel('N parentes')
    ax[1].set_xlabel('N parents')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"RMSD_fitness.png", bbox_inches='tight')
    plt.close()
