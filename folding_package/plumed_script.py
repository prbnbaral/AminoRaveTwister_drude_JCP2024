import os
import numpy
import pandas as pd
import MDAnalysis as mda
from folding_package.com_distance import  com_distance
from folding_package.res_ids import *

## THIS SUBROUTINE HAS BEEN LAST MODIFIED ON 24TH OF MAY


def plumed_script(OPs_r, RCs, time_delay, iteration, T, bias_factor, label, sigma, height, pace, output_path):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            #file.write(section.rjust(2) + "  " + type.rjust(2) + "  " + parameter.rjust(2) + "  " + value.rjust(8) + "                    \n")
            #file.write("This script is automatically created. \n\n\n\n")
            file.write("""""")
            # Defines the OPs to be combined
            for i in range(0, len(OPs_r.columns)):
                           file.write(str(OPs_r.columns[i])+": DISTANCE ATOMS="
                                      +str(int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('-')])+1)+','
                                      +str(int(OPs_r.columns[i][OPs_r.columns[i].find('-')+len('-'):OPs_r.columns[i].rfind('')])+1)+'\n' 
                           )
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_r.columns)):
                file.write(str(OPs_r.columns[i]))
                if i!=len(OPs_r.columns)-1:
                    file.write(',')
                if i==len(OPs_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs.shape[1]):
                file.write(str(RCs.iloc[len(time_delay)-1,i]))
                if i!=RCs.shape[1]-1:
                    file.write(',')
                if i==RCs.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # Define the settings of metadynamics run

            file.write("METAD ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG="+label+" ")
            file.write("SIGMA="+str(sigma)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("LABEL=metadynamics ")
            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            for i in range(0, len(OPs_r.columns)):
                file.write(str(OPs_r.columns[i]))
                file.write(',')
            file.write(",metadynamics.bias ")
            file.write("STRIDE=1 ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass


def plumed_script_com(OPs_r, number_of_bases, min_atom_id, max_atom_id, time_delay, universe_name, RCs, iteration, T, bias_factor, label, sigma, height, pace, output_path):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            #file.write(section.rjust(2) + "  " + type.rjust(2) + "  " + parameter.rjust(2) + "  " + value.rjust(8) + "                    \n")
            #file.write("This script is automatically created. \n\n\n\n")
            file.write("""""")
            # Defines the OPs to be combined
            #for i in range(0, len(OPs_r.columns)):
            #    n_1 = int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('_')])
            #    n_2 = int(OPs_r.columns[i][OPs_r.columns[i].find('_')+len('_'):OPs_r.columns[i].rfind('')])
            #    file.write('C'+str(n_1)+": COM ATOMS=")
            #    file.write(','.join([str(x+1) for x in res_ids(n_1, universe_name)]))
            #    file.write('\n')
            #    file.write('C'+str(n_2)+": COM ATOMS=")
            #    file.write(','.join([str(x+1) for x in res_ids(n_2, universe_name)]))
            #    file.write('\n')
            # Combines the defined OPs
            
            for i in range(0, number_of_bases):                
                file.write('C'+str(i+1)+": COM ATOMS=")
                file.write(','.join([str(x+1) for x in res_ids(i+1, min_atom_id, max_atom_id, universe_name)]))
                file.write('\n')
            
            
            # Defines the OPs to be combined
            for i in range(0, len(OPs_r.columns)):
                n_1 = int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('_')])
                n_2 = int(OPs_r.columns[i][OPs_r.columns[i].find('_')+len('_'):OPs_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2))+": DISTANCE ATOMS="+'C'+str(n_1)+','+'C'+str(n_2)+'\n')
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_r.columns)):
                n_1 = int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('_')])
                n_2 = int(OPs_r.columns[i][OPs_r.columns[i].find('_')+len('_'):OPs_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_r.columns)-1:
                    file.write(',')
                if i==len(OPs_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs.shape[1]):
                file.write(str(RCs.iloc[len(time_delay)-1,i]))
                if i!=RCs.shape[1]-1:
                    file.write(',')
                if i==RCs.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # Define the settings of metadynamics run

            file.write("METAD ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG="+label+" ")
            file.write("SIGMA="+str(sigma)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("LABEL=metadynamics ")
            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            for i in range(0, len(OPs_r.columns)):
                n_1 = int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('_')])
                n_2 = int(OPs_r.columns[i][OPs_r.columns[i].find('_')+len('_'):OPs_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            file.write(",metadynamics.bias ")
            file.write("STRIDE=1 ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass
        
        
        
        
        
        
def plumed_script_double(OPs_1_r, OPs_2_r, number_of_bases, min_atom_id, max_atom_id, time_delay, universe_name, RCs_1, RCs_2, iteration, T, bias_factor, sigma_1, sigma_2, height, pace, output_path):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            file.write("""""")

            # RC 1 OPs to be combined
            for i in range(0, len(OPs_1_r.columns)):
                           file.write(str(OPs_1_r.columns[i])+": DISTANCE ATOMS="
                                      +str(int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('-')])+1)+','
                                      +str(int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('-')+len('-'):OPs_1_r.columns[i].rfind('')])+1)+'\n' 
                           )
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_1_r.columns)):
                file.write(str(OPs_1_r.columns[i]))
                if i!=len(OPs_1_r.columns)-1:
                    file.write(',')
                if i==len(OPs_1_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_1.shape[1]):
                file.write(str(RCs_1.iloc[len(time_delay)-1,i]))
                if i!=RCs_1.shape[1]-1:
                    file.write(',')
                if i==RCs_1.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # RC 2 definitions
            for i in range(0, number_of_bases):                
                file.write('C'+str(i+1)+": COM ATOMS=")
                file.write(','.join([str(x+1) for x in res_ids(i+1, min_atom_id, max_atom_id, universe_name)]))
                file.write('\n')
            
            
            # Defines the OPs to be combined
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2))+": DISTANCE ATOMS="+'C'+str(n_1)+','+'C'+str(n_2)+'\n')
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_2"+" ARG=")
            
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_2_r.columns)-1:
                    file.write(',')
                if i==len(OPs_2_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_2.shape[1]):
                file.write(str(RCs_2.iloc[len(time_delay)-1,i]))
                if i!=RCs_2.shape[1]-1:
                    file.write(',')
                if i==RCs_2.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # Define the settings of metadynamics run

            file.write("METAD ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG=RC_1,RC_2 ")
            file.write("SIGMA="+str(sigma_1)+","+str(sigma_2)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("LABEL=metadynamics ")
            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            for i in range(0, len(OPs_1_r.columns)):
                file.write(str(OPs_1_r.columns[i]))
                file.write(',')
            file.write("metadynamics.bias ")
            file.write("STRIDE=1 ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass



def plumed_script_double_n(OPs_1_r, OPs_2_r, number_of_bases, min_atom_id, max_atom_id, time_delay, universe_name, RCs_1, RCs_2, iteration, T, bias_factor, sigma_1, sigma_2, height, pace, output_path):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            file.write("""""")
            
            file.write("WHOLEMOLECULES ENTITY="+str(min_atom_id+1)+"-"+str(max_atom_id+1)+'\n')

            # RC 1 OPs to be combined
            for i in range(0, len(OPs_1_r.columns)):
                           file.write(str(OPs_1_r.columns[i])+": DISTANCE ATOMS="
                                      +str(int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('-')])+1)+','
                                      +str(int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('-')+len('-'):OPs_1_r.columns[i].rfind('')])+1)+'\n' 
                           )
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_1_r.columns)):
                file.write(str(OPs_1_r.columns[i]))
                if i!=len(OPs_1_r.columns)-1:
                    file.write(',')
                if i==len(OPs_1_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_1.shape[1]):
                file.write(str(RCs_1.iloc[len(time_delay)-1,i]))
                if i!=RCs_1.shape[1]-1:
                    file.write(',')
                if i==RCs_1.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # RC 2 definitions
            rc_name_l = []
            for i in range(0, len(OPs_2_r.columns)):    
                rc_name = OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')]
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')]+": COM ATOMS=")
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')
                    
                
                rc_name = OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')]
                
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')]+": COM ATOMS=")
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')                
            
            
            # Defines the OPs to be combined
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2))+": DISTANCE ATOMS="+'C'+str(n_1)+','+'C'+str(n_2)+'\n')
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_2"+" ARG=")
            
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_2_r.columns)-1:
                    file.write(',')
                if i==len(OPs_2_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_2.shape[1]):
                file.write(str(RCs_2.iloc[len(time_delay)-1,i]))
                if i!=RCs_2.shape[1]-1:
                    file.write(',')
                if i==RCs_2.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # Define the settings of metadynamics run

            file.write("METAD ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG=RC_1,RC_2 ")
            file.write("SIGMA="+str(sigma_1)+","+str(sigma_2)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("LABEL=metadynamics ")
            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            for i in range(0, len(OPs_1_r.columns)):
                file.write(str(OPs_1_r.columns[i]))
                file.write(',')
            file.write("metadynamics.bias ")
            file.write("STRIDE=1 ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass


def plumed_script_double_same(OPs_1_r, OPs_2_r, number_of_bases, min_atom_id, max_atom_id, time_delay, universe_name, RCs_1, RCs_2, iteration, T, bias_factor, sigma_1, sigma_2, height, pace, output_path):
    ## this part of the script has been updated on 22 September 2020
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            file.write("""""")
            file.write("RESTART"+"\n")
            
            file.write("WHOLEMOLECULES ENTITY0="+str(min_atom_id+1)+"-"+str(max_atom_id+1)+'\n')
            
            
            rc_name_l = []
            for i in range(0, len(OPs_1_r.columns)):    
                
                rc_name = OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]
                
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]+": COM ATOMS=")
                    
                    
                    ##subroutine to seperate whole string to numbers
                    
                    str1=OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]
                    
                    j=0
                    dec=0
                    resid_nums=[]
                    #print(len(list(str3))-2)
                    for k in range(0, len(str1)):
                        if j<len(list(str1))-2:
                            #print('JJ: '+str(j))
                            #print("The values is "+str(list(str3)[j]))
                            if int(list(str1)[j+1])%10==9:
                                if int(list(str1)[j]+list(str1)[j+1])!=89:
                                    print(list(str1)[j]+list(str1)[j+1])
                                    resid_nums.append(list(str1)[j]+list(str1)[j+1])
                                    j=j+2
                                    #print('J1: '+str(j))
                                else:
                                    print(list(str1)[j])
                                    resid_nums.append(list(str1)[j])
                                    j=j+1
                                    #print('J1: '+str(j))
                            elif list(str1)[j]==list(str1)[j+2]:
                                print(list(str1)[j]+list(str1)[j+1])
                                resid_nums.append(list(str1)[j]+list(str1)[j+1])
                                j=j+2
                                dec=1
                                print('two decimal starts')
                                #print('J3: '+str(j))
            
                            else:
                                #list(str3)[j]!=list(str3)[j+2]:
                                print(list(str1)[j])
                                resid_nums.append(list(str1)[j])
                                j=j+1
                                #print('J2: '+str(j))
            
                        if j>=len(list(str1))-2:
                            if dec==1:
                                print(list(str1)[j]+list(str1)[j+1])
                                resid_nums.append(list(str1)[j]+list(str1)[j+1])
                                break
                            else:
                                print(list(str1)[j])
                                resid_nums.append(list(str1)[j])
                                j=j+1
                                print(list(str1)[j])
                                resid_nums.append(list(str1)[j])
                                break

                    
                    #file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]), min_atom_id, max_atom_id, universe_name)]))
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join([i for i in resid_nums]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')               
                      
                resid_nums.clear()
                rc_name = OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]+": COM ATOMS=")
                    
                    ##subroutine to seperate whole string to numbers
                    str2=OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]
                    j=0
                    dec=0
                    resid_nums=[]
                    #print(len(list(str3))-2)
                    for k in range(0, len(str2)):
                        if j<len(list(str2))-2:
                            #print('JJ: '+str(j))
                            #print("The values is "+str(list(str3)[j]))
                            if int(list(str2)[j+1])%10==9:
                                if int(list(str2)[j]+list(str2)[j+1])!=89:
                                    print(list(str2)[j]+list(str2)[j+1])
                                    resid_nums.append(list(str2)[j]+list(str2)[j+1])
                                    j=j+2
                                    #print('J1: '+str(j))
                                else:
                                    print(list(str2)[j])
                                    resid_nums.append(list(str2)[j])
                                    j=j+1
                                    #print('J1: '+str(j))
                            elif list(str2)[j]==list(str2)[j+2]:
                                print(list(str2)[j]+list(str2)[j+1])
                                resid_nums.append(list(str2)[j]+list(str2)[j+1])
                                j=j+2
                                dec=1
                                print('two decimal starts')
                                #print('J3: '+str(j))
            
                            else:
                                #list(str3)[j]!=list(str3)[j+2]:
                                print(list(str2)[j])
                                resid_nums.append(list(str2)[j])
                                j=j+1
                                #print('J2: '+str(j))
                        if j>=len(list(str2))-2:
                            if dec==1:
                                print(list(str2)[j]+list(str2)[j+1])
                                resid_nums.append(list(str2)[j]+list(str2)[j+1])
                                break
                            else:
                                print(list(str2)[j])
                                resid_nums.append(list(str2)[j])
                                j=j+1
                                print(list(str2)[j])
                                resid_nums.append(list(str2)[j])
                                break
                                        
                    #file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]), min_atom_id, max_atom_id, universe_name)]))
                    print('RESID NUMS')
                    print(" ".join([i for i in resid_nums]))
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join([i for i in resid_nums]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')
            
            resid_nums.clear()

            # RC 1 OPs to be combined
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2))+": DISTANCE ATOMS="+'C'+str(n_1)+','+'C'+str(n_2)+' NOPBC'+'\n')
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_1_r.columns)-1:
                    file.write(',')
                if i==len(OPs_1_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_1.shape[1]):
                file.write(str(RCs_1.iloc[len(time_delay)-1,i]))
                if i!=RCs_1.shape[1]-1:
                    file.write(',')
                if i==RCs_1.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # RC 2 definitions 
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_2"+" ARG=")
            
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_2_r.columns)-1:
                    file.write(',')
                if i==len(OPs_2_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_2.shape[1]):
                file.write(str(RCs_2.iloc[len(time_delay)-1,i]))
                if i!=RCs_2.shape[1]-1:
                    file.write(',')
                if i==RCs_2.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            #Define RADIUS OF GYRATION
                        
            file.write("\nGYRATION TYPE=RADIUS ATOMS=")
            file.write(str(min_atom_id+1)+"-"+str(max_atom_id+1)+' '+'LABEL=rg '+'NOPBC'+'\n')
            
            # Define the settings of metadynamics run
            file.write("METAD ")
            file.write("LABEL=metadynamics ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG=RC_1,RC_2 ")
            file.write("SIGMA="+str(sigma_1)+","+str(sigma_2)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("GRID_MIN=-15,-15 ")
            file.write("GRID_MAX=15,15 ")
            file.write("GRID_BIN=256,256 ")
            #file.write("REWEIGHTING_NGRID=256,256 ")
            #file.write("REWEIGHTING_NHILLS=10 ")
            
            file.write("CALC_RCT ")
            file.write("RCT_USTRIDE=1 ")

            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            file.write("RC_1,RC_2,")
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            file.write("rg,")
            file.write("metadynamics.bias,")
            file.write("metadynamics.rct ")            
            file.write("STRIDE="+str(pace)+" ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass


# PLUMED_SCRIPT_STACK HAS BEEN ADDED ON 24TH OF MAY
    
def plumed_script_stack(OPs_1_r, OPs_2_r, number_of_bases, min_atom_id, max_atom_id, time_delay, universe_name, RCs_1, RCs_2, iteration, T, bias_factor, sigma_1, sigma_2, height, pace, output_path):
    ## this part of the script has been updated on 22 September 2020
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            file.write("""""")
            file.write("RESTART"+"\n")
            
            file.write("WHOLEMOLECULES ENTITY0="+str(min_atom_id+1)+"-"+str(max_atom_id+1)+'\n')
            
            
            rc_name_l = []
            for i in range(0, len(OPs_1_r.columns)):    
                
                rc_name = OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]
                
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]+": COM ATOMS=")
                    
                    
                    ##subroutine to seperate whole string to numbers
                    
                    str1=OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]
                    print(str1)
                    j=0
                    dec=0
                    resid_nums=[]
                    #print(len(list(str3))-2)
                    for k in range(0, len(str1)):
                        print(list(str1)[k])
                        resid_nums.append(list(str1)[k])

                    
                    #file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')]), min_atom_id, max_atom_id, universe_name)]))
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join([i for i in resid_nums]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')               
                      
                resid_nums.clear()
                rc_name = OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]
                
                if rc_name not in rc_name_l:
                    
                    rc_name_l.append(rc_name)
                    
                    file.write('C'+OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]+": COM ATOMS=")
                    
                    ##subroutine to seperate whole string to numbers
                    str2=OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]
                    j=0
                    dec=0
                    resid_nums=[]
                    #print(len(list(str3))-2)
                    for k in range(0, len(str2)):
                        print(list(str2)[k])
                        resid_nums.append(list(str2)[k])

                                        
                    #file.write(','.join([str(x+1) for x in res_ids_n(" ".join(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')]), min_atom_id, max_atom_id, universe_name)]))
                    print('RESID NUMS')
                    print(" ".join([i for i in resid_nums]))
                    file.write(','.join([str(x+1) for x in res_ids_n(" ".join([i for i in resid_nums]), min_atom_id, max_atom_id, universe_name)]))
                    file.write('\n')
            
            resid_nums.clear()

            # RC 1 OPs to be combined
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2))+": DISTANCE ATOMS="+'C'+str(n_1)+','+'C'+str(n_2)+' NOPBC'+'\n')
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_1_r.columns)-1:
                    file.write(',')
                if i==len(OPs_1_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_1.shape[1]):
                file.write(str(RCs_1.iloc[len(time_delay)-1,i]))
                if i!=RCs_1.shape[1]-1:
                    file.write(',')
                if i==RCs_1.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # RC 2 definitions 
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_2"+" ARG=")
            
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                
                if i!=len(OPs_2_r.columns)-1:
                    file.write(',')
                if i==len(OPs_2_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs_2.shape[1]):
                file.write(str(RCs_2.iloc[len(time_delay)-1,i]))
                if i!=RCs_2.shape[1]-1:
                    file.write(',')
                if i==RCs_2.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            #Define RADIUS OF GYRATION
                        
            file.write("\nGYRATION TYPE=RADIUS ATOMS=")
            file.write(str(min_atom_id+1)+"-"+str(max_atom_id+1)+' '+'LABEL=rg '+'NOPBC'+'\n')
            
            # Define the settings of metadynamics run
            file.write("METAD ")
            file.write("LABEL=metadynamics ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG=RC_1,RC_2 ")
            file.write("SIGMA="+str(sigma_1)+","+str(sigma_2)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("GRID_MIN=-15,-15 ")
            file.write("GRID_MAX=15,15 ")
            file.write("GRID_BIN=256,256 ")
            #file.write("REWEIGHTING_NGRID=256,256 ")
            #file.write("REWEIGHTING_NHILLS=10 ")
            
            file.write("CALC_RCT ")
            file.write("RCT_USTRIDE=1 ")

            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            file.write("RC_1,RC_2,")
            for i in range(0, len(OPs_1_r.columns)):
                n_1 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('d')+len('d'):OPs_1_r.columns[i].rfind('_')])
                n_2 = int(OPs_1_r.columns[i][OPs_1_r.columns[i].find('_')+len('_'):OPs_1_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            for i in range(0, len(OPs_2_r.columns)):
                n_1 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('d')+len('d'):OPs_2_r.columns[i].rfind('_')])
                n_2 = int(OPs_2_r.columns[i][OPs_2_r.columns[i].find('_')+len('_'):OPs_2_r.columns[i].rfind('')])
                file.write(str('d'+str(n_1)+'_'+str(n_2)))
                file.write(',')
            file.write("rg,")
            file.write("metadynamics.bias,")
            file.write("metadynamics.rct ")            
            file.write("STRIDE="+str(pace)+" ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass