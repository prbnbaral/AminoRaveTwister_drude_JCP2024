from __future__ import print_function
import sys
import os
sys.executable


import numpy
numpy.set_printoptions(threshold=sys.maxsize)
import pandas as pd


## MD related libraries (these should be in a folder called files)
from input.openmm.omm_readinputs import *
from input.openmm.omm_readparams import *
from input.openmm.omm_vfswitch import *
from input.openmm.omm_barostat import *
from input.openmm.omm_restraints import *
#from input.openmm.omm_rewrap import *

import MDAnalysis as mda
import MDAnalysis.analysis.nuclinfo


# Openn MM libraries
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

# ParmEd libraries
from parmed import load_file
from parmed.openmm.reporters import NetCDFReporter
from parmed import unit as u

# Plumed plugin libraries
#Platform.loadPluginsFromDirectory('/usr/local/openmm/lib')
Platform.loadPluginsFromDirectory('/home/username/Modules/openmm/lib/plugins')
from openmmplumed import PlumedForce
print(Platform.getPluginLoadFailures())

# Own libraries
## path to the python script for data generation
Path2script = "/home/username/python_codes"
sys.path.append(Path2script)

from folding_package.traj_gen import *
from folding_package.distance_cal import *
from folding_package.col_ids import *
from folding_package.com import *
from folding_package.com_distance import *
from folding_package.traj_gen_com import *
from folding_package.res_ids import *
from folding_package.plumed_script import *
from folding_package.wc_pair_finder import *
from folding_package.plumed_script import *

# ML libraries

from keras import backend as K
from keras.models import Model
from keras.layers import Input, Dense, Lambda
from keras.initializers import RandomUniform, Constant
#from keras.optimizers import RMSprop
from keras.constraints import unit_norm
from keras import regularizers
from keras.callbacks import Callback
from keras.losses import mean_squared_error

import tensorflow as tf
tf.compat.v1.disable_eager_execution()

from folding_package.COLVAR2npy import *
from folding_package.Analyze_prave import *
from folding_package.rave import *
from folding_package.amino import *

import argparse

def arpm(iteration):
    psffile='input/md/'+'system.psf'
    crdfile='input/md/'+system+'.crd'
    topparfile='input/md/'+'toppar_drude.str'
    
    if iteration == 0:
        dcdfile=system+'.dcd' #This would be the trajectory file obtained from unbiased MD, to be used for training the OPs
        print("DCD FILE NAME: "+dcdfile) 
    if iteration == 1:
        dcdfile='output/md/'+'equilibration_0.dcd'
    
    if iteration > 1:
        dcdfile='output/md/'+'output_'+str(iteration-1)+'.dcd'


    # Load parameters
    print("Loading parameters")
    #inputs = read_inputs(inpfile)
    params = read_params(topparfile)
    top = read_top(psffile)
    gro = read_crd(crdfile)
    print('Parameters loaded')

    #Define universe for MDAnalysis
    rna = mda.Universe(psffile, dcdfile)
    traj_len = len(rna.trajectory)
    print("Length of trajectory is: "+str(traj_len))

    ## Atoms selection subroutine
    ura_wc_n3=rna.select_atoms('resname URA and name N3')
    ura_wc_o4=rna.select_atoms('resname URA and name O4')
    ade_wc_n1=rna.select_atoms('resname ADE and name N1')
    ade_wc_n6=rna.select_atoms('resname ADE and name N6')
    gua_wc_o6=rna.select_atoms('resname GUA and name O6')
    gua_wc_n1=rna.select_atoms('resname GUA and name N1')
    gua_wc_n2=rna.select_atoms('resname GUA and name N2')
    cyt_wc_n4=rna.select_atoms('resname CYT and name N4')
    cyt_wc_n3=rna.select_atoms('resname CYT and name N3')
    cyt_wc_o2=rna.select_atoms('resname CYT and name O2')

    num_ura_sel = ura_wc_n3.positions.shape[0]
    num_ade_sel = ade_wc_n1.positions.shape[0]
    num_gua_sel = gua_wc_o6.positions.shape[0]
    num_cyt_sel = cyt_wc_n4.positions.shape[0]

    print("Number of URACIL: "+str(num_ura_sel))
    print("Number of ADENINE: "+str(num_ade_sel))
    print("Number of GUANINE: "+str(num_gua_sel))
    print("Number of CYTOSINE: "+str(num_cyt_sel))

    # Calculates the number of bases, atom id range of bases
    number_of_bases = int((num_ade_sel+num_cyt_sel+num_gua_sel+num_ura_sel))
    max_id = max(rna.select_atoms('nucleic').ids)
    min_id =min(rna.select_atoms('nucleic').ids)

    # COM related OPs generation

    if os.path.isfile('output/ml/'+'data_rc1.csv'):
        print('Previosly generated order parameters are imported')
        OPs_1 = pd.read_csv('output/ml/'+'data_rc1.csv', sep=' ')
        OPs_1 = OPs_1.drop(OPs_1.columns[0], axis=1)
    else:
        print('New order parameters are generated')
        # WC pairs related OP generation      
        pairs = wc_pair_finder(number_of_bases, rna)
        OPs_1 = traj_gen_com_n(number_of_bases, max_id, min_id, 1, rna, 4)
        OPs_1.to_csv('output/ml/'+'data_rc1.csv', sep=' ')
    

    #initialize AMINO package
    if not os.path.isfile('input/data_1_0'):
        print('data_1_0 file does not exist.')    
        if iteration==0:
        
             all_ops = []
             for i in OPs_1.columns:
                 all_ops.append(OrderParameter(i, OPs_1[i]))
             final_ops = find_ops(all_ops, 2, 50)
        
             print("\nAMINO order parameters:")
             cols = []
             for i in final_ops:
                 print(i)
                 cols.append(str(i))
        

        OPs_1_reduced=OPs_1.filter(cols)
        OPs_1_reduced.to_csv('input/data_1_0', sep=' ')
    
        print("OPs trajectory length: "+str(OPs_1_reduced.shape[0]))
        with open('input/data_1_0', "r+") as myfile:
            myfile.seek(0)
            myfile.write("# ")
            myfile.close()
    else:
        OPs_1_reduced=pd.read_csv('input/data_1_0', sep=' ')
        OPs_1_reduced=OPs_1_reduced.drop(OPs_1_reduced.columns[0], axis=1)
        print('data_2_0 file exists.')
    # COM related OPs generation

    if os.path.isfile('output/ml/'+'data_rc2.csv'):
        print('Previosly generated order parameters are imported')
        OPs_2 = pd.read_csv('output/ml/'+'data_rc2.csv', sep=' ')
        OPs_2 = OPs_2.drop(OPs_2.columns[0], axis=1)
    else:
        print('New order parameters are generated')
        OPs_2 =traj_gen_com_n(number_of_bases, max_id, min_id, 1, rna, 4)
        OPs_2.to_csv('output/ml/'+'data_rc2.csv', sep=' ')
    #OPs =traj_gen_com(number_of_bases, max_id, min_id, 10, rna)


    #initialize AMINO package
    if not os.path.isfile('input/data_2_0'):
        print('data_2_0 file does not exist.')    
        if iteration==0:
    
            all_ops = []
            for i in OPs_2.columns:
                all_ops.append(OrderParameter(i, OPs_2[i]))
            final_ops = find_ops(all_ops, 2, 50)
    
            print("\nAMINO order parameters:")
            cols = []
            for i in final_ops:
                print(i)
                cols.append(str(i))

        OPs_2_reduced=OPs_2.filter(cols)
        OPs_2_reduced.to_csv('input/data_2_0', sep=' ')
        print("OPs trajectory length: "+str(OPs_2_reduced.shape[0]))
        with open('input/data_2_0', "r+") as myfile:
            myfile.seek(0)
            myfile.write("# ")
            myfile.close()
    else:
        OPs_2_reduced=pd.read_csv('input/data_2_0', sep=' ')
        OPs_2_reduced=OPs_2_reduced.drop(OPs_2_reduced.columns[0], axis=1)
        print('data_2_0 file exists.')
    

    # Defining inputs for RAVE

    RCs = {}
    for i in range(0, 2):
        nRC = i+1
        # Defining inputs for RAVE
        system_name = 'data_'+str(nRC)
        n_trajs = 1   
        save_path = 'output/rave/'  
        input_path = 'input/'
        T = 300     
        bias = True  


        #predictive time delay #
        time_delay= list(range(0, 2, 1)) #predictive time delay
        trials = range(4) 

        #network variables
        batch_size = 1
        training_size = 4
        #training_size = int(OPs_reduced.shape[0]/batch_size)*batch_size-batch_size
        print("training size: "+str(training_size))

        #module = sys.modules[__name__]
        #func = getattr(module, "OPs_{}_reduced".format(nRC))
        #op_dim = func.shape[1] 
        op_dim= locals()["OPs_{}_reduced".format(nRC)].shape[1]
        rc_dim = 1  
        int_dim = 128  
        s_vari = 0.005 
        learning_rate = 0.0002 
        decay = 0.0       
        epochs = 20 

        # End of input definitions for RAVE

        rave(system_name, n_trajs, save_path, input_path, T, bias, time_delay, training_size, batch_size, op_dim, rc_dim, int_dim, s_vari, learning_rate, decay, epochs, trials, iteration)

        network_info = '_int_dim'+str(int_dim)+'_lr'+str(learning_rate)+'_decay'+str(decay)+'_batch_size'+str(batch_size)
        save_result(system_name, op_dim, time_delay, trials, s_vari, training_size, network_info, save_path)  

        # Extracting reactions coordinate weights from RAVE output
        RCs[nRC] = pd.read_csv('output/rave/final_result_'+system_name+'_svar'+str(s_vari)+'_train_size'+str(training_size)+network_info+'.txt', delimiter= '\s+')
    
    # MD part
    if iteration == 0:
        print("EQUILIBRATION STARTS")
        print('Building the system')
        
        inpfile='input/md/'+'step4_equilibration.inp'
        inputs = read_inputs(inpfile)
        gen_box(top, gro)
        #gen_box(top, gro).setBox(6.0635, 6.0503, 6.0604)
        
        
        # Defining the metadynamics inputs
#        from folding_package.plumed_script import *
#       bias_factor=15
#       height=0
#       pace=100
#       RCs_1 = RCs[1]
#       RCs_2 = RCs[2]
#       OPs_1_r = OPs_1_reduced
#       OPs_2_r = OPs_2_reduced
#       
#       sigma_1= abs((OPs_1_reduced.std()[0:RCs_1.shape[1]-1].values*RCs_1.iloc[len(time_delay)-1,1:RCs_1.shape[1]].values).sum())*0.1
#       print("Metadynamics sigma value for RC 1:" + str(sigma_1))
#       
#       sigma_2= abs((OPs_2_reduced.std()[0:RCs_2.shape[1]-1].values*RCs_2.iloc[len(time_delay)-1,1:RCs_2.shape[1]].values).sum())*0.1
#       print("Metadynamics sigma value for RC 2:" + str(sigma_2))
#       
#       plumed_script_double_same(OPs_1_r, OPs_2_r, number_of_bases, min_id, max_id, time_delay, rna, RCs_1, RCs_2, iteration, T, bias_factor, sigma_1, sigma_2, height, pace, 'output/meta/')
#       
#       file = open('output/meta/'+'plumed_s_'+str(iteration)+'.dat',"r+")
#       script = file.read()
#       file.close()
        
        system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                  nonbondedCutoff=inputs.r_off*nanometers,
                                  constraints=inputs.cons,
                                  ewaldErrorTolerance=inputs.ewald_Tol,
                                 )
        print("SIMULATION BOX VECTORS: "+str(system.getDefaultPeriodicBoxVectors()))
        if inputs.vdw == 'Force-switch':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
            system = vfswitch(system, top, inputs)
        if inputs.vdw == 'LJPME':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
        if inputs.vdw == 'Switch':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
                                      
        #if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs)
        if inputs.pcouple == 'yes':      system = barostat(system, inputs)
        if inputs.rest == 'yes':         system = restraints(system, top, gro, inputs)
        print('System has been built')
        
        
        #system.addForce(PlumedForce(script))
        
        integrator = DrudeLangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.drude_temp*kelvin, inputs.drude_fric_coeff/picosecond, inputs.dt*picoseconds)
        integrator.setMaxDrudeDistance(inputs.drude_hardwall) # Drude Hardwall
        
        print("Setting platform")
        platform = Platform.getPlatformByName('CUDA')
        
        #prop = dict(CudaPrecision='mixed', DeviceIndex="0")
        prop = dict(CudaPrecision='mixed')
        print("Platform has been set")
        
        print("Building simulation context")
        simulation = Simulation(top.topology, system, integrator, platform, prop)
        simulation.context.setPositions(gro.positions)
                
        
        # Drude VirtualSites
        simulation.context.computeVirtualSites()
        
        
        # Energy minimization
        
        if inputs.mini_nstep > 0:
            print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
            simulation.minimizeEnergy(tolerance=inputs.mini_Tol*kilojoule/mole, maxIterations=inputs.mini_nstep)
            print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
        
        
        # Generate initial velocities
        if inputs.gen_vel == 'yes':
            print("\nGenerate initial velocities")
            if inputs.gen_seed:
                simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
            else:
                simulation.context.setVelocitiesToTemperature(inputs.gen_temp) 
      
        # Production
        if inputs.nstep > 0:
            print("\nMD run: %s steps" % inputs.nstep)
            if inputs.nstdcd > 0:
                odcd = 'output/md/'+'equilibration_'+str(iteration)+'.dcd'
                simulation.reporters.append(DCDReporter(odcd, inputs.nstdcd))
            simulation.reporters.append(
                StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, progress=True,
                                  remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
            )
            simulation.reporters.append(PDBReporter('output/md/'+'equilibration_'+str(iteration)+'.pdb', inputs.nstep))
            simulation.step(inputs.nstep) 
        # Write restart file
        orst = 'output/md/'+'equilibration_'+str(iteration)+'.rst'
        
        if 'orst' in locals() or 'orst' in globals():
            state = simulation.context.getState( getPositions=True, getVelocities=True )
            with open(orst, 'w') as f:
                f.write(XmlSerializer.serialize(state))
        if 'ochk' in locals() or 'ochk' in globals():
            with open(ochk, 'wb') as f:
                f.write(simulation.context.createCheckpoint())
                
        print("EQUILIBRATION ENDS")  
    
    if iteration > 0:
        print("PRODUCTION STARTS")
        print('Building the system')
        inpfile='input/md/'+'step5_production.inp'
        inputs = read_inputs(inpfile)
        
        gen_box(top, gro)
        system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                  nonbondedCutoff=inputs.r_off*nanometers,
                                  constraints=inputs.cons,
                                  ewaldErrorTolerance=inputs.ewald_Tol,
                                 )
        if inputs.vdw == 'Force-switch':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
            system = vfswitch(system, top, inputs)
        if inputs.vdw == 'LJPME':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
        if inputs.vdw == 'Switch':
            system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                                      nonbondedCutoff=inputs.r_off*nanometers,
                                      constraints=inputs.cons,
                                      ewaldErrorTolerance=inputs.ewald_Tol)
                                      
        #if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs)
        if inputs.pcouple == 'yes':      system = barostat(system, inputs)
        if inputs.rest == 'yes':         system = restraints(system, top, gro, inputs)
        print('System has been built')
                
        # Defining the metadynamics inputs
        #from folding_package.plumed_script import *
        bias_factor=15
        height=1.5
        pace=500
        RCs_1 = RCs[1]
        RCs_2 = RCs[2]
        OPs_1_r = OPs_1_reduced
        OPs_2_r = OPs_2_reduced
        
        sigma_1= abs((OPs_1_reduced.std()[0:RCs_1.shape[1]-1].values*RCs_1.iloc[len(time_delay)-1,1:RCs_1.shape[1]].values).sum())*0.1
        print("Metadynamics sigma value for RC 1:" + str(sigma_1))
        
        sigma_2= abs((OPs_2_reduced.std()[0:RCs_2.shape[1]-1].values*RCs_2.iloc[len(time_delay)-1,1:RCs_2.shape[1]].values).sum())*0.1
        print("Metadynamics sigma value for RC 2:" + str(sigma_2))
        
        if os.path.isfile('output/meta/'+'plumed_s_'+str(iteration)+'.dat'):
            
            file = open('output/meta/'+'plumed_s_'+str(iteration)+'.dat',"r+")
            script = file.read()
            file.close()
            
        else:
            
            plumed_script_double_same(OPs_1_r, OPs_2_r, number_of_bases, min_id, max_id, time_delay, rna, RCs_1, RCs_2, iteration, 300, bias_factor, sigma_1, sigma_2, height, pace, 'output/meta/')
            file = open('output/meta/'+'plumed_s_'+str(iteration)+'.dat',"r+")
            script = file.read()
            file.close()
        
        system.addForce(PlumedForce(script))
        
        integrator = DrudeLangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.drude_temp*kelvin, inputs.drude_fric_coeff/picosecond, inputs.dt*picoseconds)
        integrator.setMaxDrudeDistance(inputs.drude_hardwall) # Drude Hardwall
        
        print("Setting platform")
        platform = Platform.getPlatformByName('CUDA')
        prop = dict(CudaPrecision='mixed')
        print("Platform has been set")
        
        print("Building simulation context")
        simulation = Simulation(top.topology, system, integrator, platform, prop)
        simulation.context.setPositions(gro.positions)
        
        # Drude VirtualSites
        simulation.context.computeVirtualSites()
        
        if iteration == 1:
            irst = 'output/md/'+'equilibration_0.rst'
        if iteration > 1:
            irst = 'output/md/'+'output_'+str(iteration-1)+'.rst'
        
        if 'irst' in locals() or 'irst' in globals():
            with open(irst, 'r') as f:
                simulation.context.setState(XmlSerializer.deserialize(f.read()))
        if 'ichk' in locals() or 'ichk' in globals():
            with open(ichk, 'rb') as f:
                simulation.context.loadCheckpoint(f.read())
        
        
        # Production
        if inputs.nstep > 0:
            print("\nMD run: %s steps" % inputs.nstep)
            if inputs.nstdcd > 0:
                odcd = 'output/md/'+'output_'+str(iteration)+'.dcd'
                simulation.reporters.append(DCDReporter(odcd, inputs.nstdcd))
            simulation.reporters.append(
                StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, progress=True,
                                  remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
            )
            simulation.reporters.append(PDBReporter('output/md/'+'output_'+str(iteration)+'.pdb', inputs.nstep))
            simulation.step(inputs.nstep) 
        # Write restart file
        orst = 'output/md/'+'output_'+str(iteration)+'.rst'
        
        if 'orst' in locals() or 'orst' in globals():
            state = simulation.context.getState( getPositions=True, getVelocities=True )
            with open(orst, 'w') as f:
                f.write(XmlSerializer.serialize(state))
        if 'ochk' in locals() or 'ochk' in globals():
            with open(ochk, 'wb') as f:
                f.write(simulation.context.createCheckpoint())
                
        print("PRODUCTION ENDS")  
          
        print("End of iteration: "+str(iteration))


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='iteration', help='iteration', required=True)
args = parser.parse_args()

arpm(int(args.iteration))


#if __name__ == "__main__":
#    num_iter=50
#    for iteration in range(0, num_iter):
#        arpm(iteration)
