from mpi4py import MPI
import subprocess
import numpy as np
import os
import sys

class TaskManager():
    def __init__(self,tasks_list_fname,outdir):
        self.tasks_list_fname = tasks_list_fname
        self.outdir = outdir

    def generate_tasks(self):
        tasks_list = []

        num = 1
        with open(self.tasks_list_fname,'r') as file:
            for line in file:
                command = line
                out     = os.path.join(self.outdir,"progress_%s.out" % num)
                num += 1

                tasks_list.append([command,out])

        return tasks_list
    
    def execute_tasks(self,tasks_list):
        #loop over tasks in list
        for task in tasks_list:         
            command = task[0]
            out_fname = task[1]
        
            out = open(out_fname,"w")
            
            #initiate subprocess for calculations
            process = subprocess.Popen(command.split(),
                                       stdin=subprocess.PIPE,
                                       stdout=out,
                                       bufsize=1,universal_newlines=True)
            
            #wait till the end of calculations
            process.wait()
        
            out.close()
        
    
if __name__ == '__main__':
    tasks_list_fname = sys.argv[1]
    outdir = sys.argv[2]

    #create output directory
    os.makedirs(outdir,exist_ok=True)

    #initialize mpi variables
    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()  #number of processes
    mpi_rank = mpi_comm.Get_rank()  #rank of current process
    mpi_master = 0                  #master process

    task_manager = TaskManager(tasks_list_fname,outdir)
    
    job_list = []
    
    #prepare task list
    if mpi_rank == mpi_master:  #preparation is performed on master node
        print(
        "------------------------\n"
        "--- MPI task manager ---\n"
        "------------------------\n")
        
        #here we form total list of tasks
        tasks_list_total = task_manager.generate_tasks()

        #calculation of numbers of jobs per node
        n_tasks_all = len(tasks_list_total)                     #total number of tasks
        n_remain = n_tasks_all % mpi_size                       #number of unnisigned tasks
        n_main = int((n_tasks_all - n_remain) / mpi_size)       #number of tasks per process (without remaining)
        
        print("Total number of jobs: %s" % n_tasks_all)
        print("Number of processes requested: %s" % mpi_size)
        print("Job list size: %s, Remaining jobs: %s" % (n_main,n_remain))
        
        n_local_jobs = np.zeros(mpi_size,int)
        i_job = 0

        #loop over processes
        for i_process in range(mpi_size):
            job_local = []

            #loop over block size
            for i_el in range(n_main):
                job_local.append(tasks_list_total[i_job])
                n_local_jobs[i_process] += 1
                i_job += 1

            #append remaining
            if i_process < n_remain:
                job_local.append(tasks_list_total[i_job])
                n_local_jobs[i_process] += 1
                i_job += 1

            print("Length of task list: %s for process number: %s" % (n_local_jobs[i_process],i_process+1))

            #append task to task list
            job_list.append(job_local)
    else:
        job_list = None

    #distribute tasks between processes
    local_job_list = mpi_comm.scatter(job_list,root=mpi_master)
    
    #perform calculations with local tasks list
    local_results = task_manager.execute_tasks(local_job_list)

    if mpi_rank == mpi_master:
        print("Task accomplished!")
