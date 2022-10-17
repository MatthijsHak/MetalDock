
import os, glob, subprocess
import shutil
import numpy as np
import uuid
import prepare_dock as d


import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool
import time
from time import perf_counter

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


def numerical_derviative(tmp_dir, dir_list, par, parameter_set, index, x, dx, learning_rate):
    args = [(tmp_dir, dir_list, par, parameter_set, index, x+i) for i in [dx,-dx]]

    pool = multiprocessing.Pool(processes=2)
    output = pool.starmap(func_dock, args)
    
    pool.close()
    pool.join()

    df = output[0] - output[1] / (2*dx)
    # df = ( func_dock(tmp_dir, dir_list, par, parameter_set, index, x+dx) -  func_dock(tmp_dir, dir_list, par, parameter_set, index, x-dx)) / (2 * dx)
    
    parameter_set[index] = parameter_set[index] - learning_rate * df

    return parameter_set[index]


def gradient_descent(tmp_dir, dir_list, par, parameter_set, learning_rate=0.5, n_iter=50):
    iteration = 0
    while iteration < n_iter:
        args = [(tmp_dir, dir_list, par, parameter_set, idx, i, 0.5, learning_rate) for idx, i in enumerate(parameter_set)]
        with MyPool(processes=12) as pool:
            start = perf_counter()

            pool.starmap(numerical_derviative, args)
            pool.close()
            pool.join()
            end = perf_counter()
            execution_time = (end - start)
            # for idx, i in enumerate(parameter_set):
            #      pool.starmap(numerical_derviative, zip(tmp_dir, dir_list, par, parameter_set, idx, i, 0.5))
        iteration+=1 
        print('{}\n'.format(parameter_set))
        print('{}'.format(execution_time))
    return parameter_set



def func_dock(tmp_dir, dir_list, par, parameter_set, index, par_value):
    rmsd_avg = []
    rmsd_list = []
    avg_list = []
    min_list = []
    rmsd_print_list = []
    output = []

    parameter_set[index] = par_value

    input_dir ='/home/matthak/AutoDock_GA_par/autodock_metal_GA/src/test/gradient_descent'
    dir_name = str(uuid.uuid4())
    os.mkdir(f'{tmp_dir}/'+dir_name)
    os.chdir(f'{tmp_dir}/'+dir_name)

    tmp_dir_GA=os.getcwd()


    for n_prot in dir_list:
        os.chdir(f'{tmp_dir}/'+dir_name)
        os.mkdir(f'protein_{n_prot}')
        os.chdir(f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}')

        os.system(f"cp -r  ../../../../data_set/protein_{n_prot}/output/docking .")
        os.chdir(f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}/docking')
        output_temp_dir=f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}/docking'

        # Obtain ligand and protein names
        for files in glob.glob("*_1.pdbqt"):
            file_list = files.split('_1.pdbqt')
            name_ligand = file_list[0]

        for files in glob.glob("clean_*.pdb"):
            file_list = files.split('clean_')
            file_list = file_list[1].split('.')
            name_protein = file_list[0]

        dock = d.get_coordinates(f'ref.xyz',par.metal_symbol)
        e = open(f'{output_temp_dir}/energy','r')
        lines = [line.split() for line in e]
        energy = lines[0][4]

        ##### Fitness function ######
        subprocess.call([f'cp {input_dir}/'+par.parameter_file+' .'], shell=True)
        # insert parameters for R and epsilon for H-bond
        subprocess.call([r'''awk '{ if ($2 == "'''+par.metal_symbol.upper()+'''" || $2 == "'''+par.metal_symbol+'''") ($7 = '''+str(parameter_set[10])+''') && ($8 = '''+str(parameter_set[11])+'''); print $0}' '''+par.parameter_file+''' > file_1'''], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call([r'''awk '{ if ($2 == "'''+par.metal_symbol.upper()+'''" || $2 == "'''+par.metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+par.parameter_file], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call(['rm file_1'], shell=True)

        #create_gpf():
        subprocess.call([os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+par.parameter_file+" -p npts='{},{},{}'".format(par.box_size,par.box_size,par.box_size)+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2])], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        gpf = open('clean_'+name_protein+'.gpf', 'a')
        gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[0],parameter_set[1])+'    12 6 OA '+par.metal_symbol+'\n')
        gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[2],parameter_set[3])+'    12 6 SA '+par.metal_symbol+'\n')
        gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[4],parameter_set[5])+'    12 6 HD '+par.metal_symbol+'\n')
        gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[6],parameter_set[7])+'    12 6 NA '+par.metal_symbol+'\n')
        gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[8],parameter_set[9])+'    12 6  N '+par.metal_symbol+'\n')
        gpf.close()

        #autogrid()
        subprocess.call([os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        #create_dpf()
        d.write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, par.parameter_file, energy, par.num_poses, par.dock_algorithm, random_pos=par.random_pos, SA=par.sa_dock, GA=par.ga_dock)

        #autodock()
        subprocess.call([os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        #write_all_conformations()
        subprocess.call([os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg"], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


        i = 1
        while os.path.exists(name_ligand+"_%i.pdbqt" % i):
            subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
            rmsd_list.append(rmsd)
            i += 1

    avg_output = np.mean(rmsd_list)
    avg_list.append(avg_output)
    output.append(f"Average RMSD: {avg_output:.4}\n")

    minimum_rmsd = min(rmsd_list)
    min_list.append(minimum_rmsd)
    output.append(f"Lowest RMSD: {minimum_rmsd:.4}\n")

    stdv_rmsd = np.std(rmsd_list)
    output.append(f"Standard Deviation RMSD: {stdv_rmsd:.4}\n")

    var_rmsd = np.var(rmsd_list)
    output.append(f"Variance RMSD: {var_rmsd:.4}\n")
    output.append("-----------------------------------------------------------------------------------------------------------\n")

    print(''.join(output))
    shutil.rmtree(f'{tmp_dir_GA}', ignore_errors=True)
    return avg_output