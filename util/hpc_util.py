import os, subprocess,sys
import math

def execute_code(exec_str, log=None, make_call=True, bg=False):
    if log is not None:
        exec_str = "{0:s} >& {1:s}".format(exec_str, log)
    if bg:
        exec_str = "{0:s} &".format(exec_str)

    # print('... exec_str: ',exec_str)
    # if make_call: subprocess.call(exec_str, shell=True,timeout=30,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if make_call: subprocess.call(exec_str, shell=True)
    # print('... code exectued')

def execute_str(bin, infile, bg=False):

    exec_str = "{0:s} {1:s}".format(bin, infile)
    if bg: exec_str = "{0:s} &".format(exec_str)
    return exec_str

def mpi_execute_str(bin, infile, nproc=1, mpi=None, bg=False, machinefile=None):

    if mpi is None: mpi = os.path.expandvars('$MPIEXEC')

    # Allow bin, infile, ncpu to be lists of length 2 for case where MPI
    # communicates between two codes
    if type(bin)    is str: bin    = [bin]
    if type(infile) is str: infile = [infile]
    if type(nproc)  is int: nproc  = [nproc]

    # Execution string
    exec_str = "{0:s} -np {1:d} {2:s} {3:s}".format(mpi, nproc[0], bin[0],
                                                   infile[0])
    if len(bin) == 2 and len(infile) == 2 and len(nproc) == 2:
        exec_str = "{0:s} : -np {1:d} {2:s} {3:s}".format(exec_str, nproc[1],
                                                         bin[1], infile[1])
    if machinefile is not None:
        exec_str = "{0:s} -machinefile {1:s}".format(exec_str, machinefile)
    if bg: exec_str = "{0:s} &".format(exec_str)

    return exec_str

def batch_maui_pbs():
    pass

def batch_slurm_pbs(name, hpc, exec_dir, bin, infile, log):
    f = open("batch.sh", "w")
    if hpc.partition:
      partition = "#SBATCH -p {0:s}\n".format(hpc.partition)
    else:
      partition = "" 
    f.write("#!/bin/bash\n"\
    +"#SBATCH --job-name={0:s}\n".format(name)\
    +"#SBATCH --output={0:s}.output\n".format(name)\
    +"#SBATCH --error={0:s}.error\n".format(name)\
    +"#SBATCH --time=48:00:00\n"\
    +"#SBATCH --nodes={0:d}\n".format(int(math.ceil(float(hpc.nproc)/hpc.ppn)))\
    +"#SBATCH --ntasks-per-node={0:d}\n".format(hpc.ppn)\
    +partition
    +"cd $SLURM_SUBMIT_DIR\n"\
    +"cd {0:s}\n".format(exec_dir)\
    +"module load openmpi/1.8.3/gcc\n"\
    +"mpirun -n {0:d} {1:s} {2:s} >& {3:s}\n".format(hpc.nproc, bin, infile, log))
    f.close()
    os.system("sbatch batch.sh")
    os.remove("batch.sh")

def batch_machine_specific():
    pass
