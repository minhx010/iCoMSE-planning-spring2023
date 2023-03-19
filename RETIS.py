# RETIS.py
# Created by PH Minh and Naomi Trampe
# Last Updated Date: 03-18-2023 
# SAMPEL Group

# Import neccessary packages
import numpy as np
import sys
import math
import copy
from math import exp
import langevin_dynamics as ld
import random
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm

def performretis(pes_type,op_type,interfacevals,nmoves,basinAloc,initcoords,basineqlen,basinlen,beta,beta_initpath):
    # RETIS Code Starts!
    # Select your PES type:
    pes_type = pes_type

    # RETIS settings
    n_moves = nmoves               # number of moves to perform for each path ensemble - path ensemble corresponding to each diff interface
    op_type = op_type              # order parameter;  1:x  2:y  3:x+y

    interfaces = interfacevals

    basinA = basinAloc          # boundary of basin A; must be less than the first interface
    basinB = interfaces[-1]     # basin B defined as last interface
    init_coords = initcoords    # initial coordinates
    init_p = [0, 0]             # initial momenta
    basineqlen = basineqlen     # basin equilibration length (steps)
    basinlen = basinlen         # basin simulation length
    sigma = 0.5                 # standard deviation of momenta perturbation

    # Langevin Dynamics Settings
    beta = beta                 # 1/kT
    gamma = 5.0                 # friction
    dt = 0.01                   # time step size

    ## Equilibration Run 
    # Declare array to store basin trajectory
    basintraj = np.zeros((basinlen + 1, 6),dtype=float)                             # matrix: each row is the next configuration and each column is x, y, px, py, fx, fy

    # Calculate initial forces
    fx,fy = ld.force(init_coords[0],init_coords[1],init_p[0],init_p[1],dt,beta,gamma,pes_type)

    # Combine positions, momenta, and forces to make an initial phase point
    init_phasepoint = init_coords + init_p + [fx,fy]
    basintrajeq = ld.vv_step(init_phasepoint,dt,beta,gamma,pes_type)                # use velocity verlet to generate the next configuration of x, y, px, py, fx, fy

    # Equilibrate in basin
    for i in range(1,basineqlen + 1):
        new_basintrajeq = ld.vv_step(basintrajeq,dt,beta,gamma,pes_type)            # array of x, y, px, py, fx, fy
        basintrajeq = new_basintrajeq
        op = ld.calc_op(op_type,basintrajeq[0],basintrajeq[1])                      #ld.calc_op(op_type, x, y of the trajectory)
        # check if trajectory reaches basin B
        if op >= basinB:
            sys.exit("Basin trajectory reached B! Exiting...")


    ## Run basin simulation 
    basintraj[0] = basintrajeq                                                      # saving the last configuration from the equilibration run as the first config for prod run
    fromBasin = False
    n_cross_flux = 0

    # Run one basin A simulation and check for first crossings
    for j in range(1,basinlen + 1):
        basintraj[j] = ld.vv_step(basintraj[j-1],dt,beta,gamma,pes_type)            # each iteration is a configuration so that basintraj at the end makes a matrix
        op = ld.calc_op(op_type,basintraj[j][0],basintraj[j][1]) 
        
        if op < basinA:
            fromBasin = True
        
        if fromBasin == True and op >= interfaces[0]: 
            n_cross_flux += 1
            fromBasin = False
        
        # Check if trajectory reaches basin B
        if op >= basinB:
            sys.exit("Basin trajectory reached B! Exiting...")

    if n_cross_flux == 0:
        sys.exit("No first crossings obtained from basin A to interface 0. Exiting...")

    ## Generate initial path at higher temp.
    # Run at higher temperature to generate an initial path
    initpath = []
    beta_initpath = beta_initpath                                          # this controls the temperature because beta = 1/kb*T
    fromBasin = False

    # Use init_phasepoint to start trajectory
    trajstep = ld.vv_step(init_phasepoint,dt,beta_initpath,gamma,pes_type)  # array of x, y, px, py, fx, fy to get the first config to generate a path
    op = ld.calc_op(op_type,trajstep[0],trajstep[1])                        # calculate the op to see where this config is at, in A or pass A and cross some lambda, or in B?

    while op < basinB:                                                      # as long as it is not in B
        nextstep = ld.vv_step(trajstep,dt,beta_initpath,gamma,pes_type)              # generate next config
        trajstep = nextstep                                                 # save this config to trajstep, so we could shoot from it again 
        op = ld.calc_op(op_type,trajstep[0],trajstep[1])

        if op < basinA:                                                     # if this then this traj is from basin A 
            fromBasin = True
            initpath = []                                                   # reset path if it returns to A - do this until it's about to leave the basin
            trajstep = np.append(trajstep,[op],axis=0)                      # array of x, y, px, py, fx, fy and op value
            initpath.append(trajstep) 
            continue

        if op >= basinA and fromBasin == True:                              # now that you have left basin A, we'll save the different configs to the list initpath
            trajstep = np.append(trajstep,[op],axis=0)
            initpath.append(trajstep)
            continue

    initpath = np.asarray(initpath)                                         # this gives all configurations of this 1 initial path
    allpaths = [[] for i in range(len(interfaces)-1)]                       # list of empty lists of size len(interfaces)-1, each empty list is a placeholder for paths in 1 interface 

    for i in range(len(interfaces)-1):                                      # loop over amt of interface starting from 0 
        allpaths[i].append(initpath)

    
    ## Perform RETIS
    # Define function to find index of a certain OP value along a path 
    def find_config_minus(path,interfaces):
        frame = np.where(path[:,6] > interfaces)[0][0]
        return frame

    def find_config_plus(path,interfaces):
        frame = np.where(path[:,6] < interfaces)[0][-1]
        return frame
    
    # DO RETIS
    for move in range(1, n_moves):
        print('Move {}'.format(move))

        # Swap or Shoot? 
        ss = np.random.ranf()
        
        # Perform swap move
        if ss < 0.5: 

            # Decides on which ensemble pair switching we're doing
            ens_switch = np.random.ranf() 
            
            if ens_switch < 0.5: 
                # print('Move = {} ; Performing Swap Move: 1-2, 3-4, 5-6, ...etc'.format(move), end='\r')
                
                # Accounting for the 1st interface, but first we need to check that it satisfies the ensemble criteria
                if np.where(allpaths[0][-1][:,6]>interfaces[0])[0].shape[0] == 0: 
                    reject_move = True
                    # print('Appending Rejected!')
                else:
                    allpaths[0].append(allpaths[0][-1])                             # append the previous 0th ensemble path to its own ensemble again because 
                                                                                    # there is no pair for this interface to swap with
                
                # For even number of interfaces:
                if (len(interfaces)-1)%2 == 0:                      
                    
                    # Accounting for the last interface 
                    if np.where(allpaths[-1][-1][:,6]>interfaces[-2])[0].shape[0] == 0: 
                        reject_move = True
                        # print('Appending Rejected!')
                    else:
                        allpaths[-1].append(allpaths[-1][-1])                       # append the previous path of the last ensemble to its own ensemble again because 
                                                                                    # there is no pair for this interface to swap with
                    
                    # Accounting for the rest of the interfaces
                    for i in np.arange(1,len(interfaces)-3,2):
                        # print('i = ', i)
                        path_1 = allpaths[i][move-1]                                # looks at the most current path of the TPE of interface i
                        path_2 = allpaths[i+1][move-1]                              # looks at the most current path of the TPE of interface i+1
                        
                        # Check if path 1 crosses path 2's interface
                        n_cross_config = np.where(path_1[:,6]>interfaces[i+1])[0].shape[0]

                        if n_cross_config != 0:                                     # if the swap is successful, do the swap
                            # print('Successful Swap!')
                            allpaths[i].append(path_2)
                            allpaths[i+1].append(path_1)
                            
                        else:                                                       # if the swap is unsuccessful, just append the same path to their own ensemble
                            # print('Unsuccessful Swap!')
                            allpaths[i].append(path_1)
                            allpaths[i+1].append(path_2)
                        
                # For odd number of interfaces:
                else:                              
                    for i in np.arange(1,len(interfaces)-1,2):
                        # print('i = ', i)
                        path_1 = allpaths[i][move-1]                                # looks at the most current path of the TPE of interface i
                        path_2 = allpaths[i+1][move-1]                              # looks at the most current path of the TPE of interface i+1
                        
                        # check if path 1 crosses path 2's interface
                        n_cross_config = np.where(path_1[:,6]>interfaces[i+1])[0].shape[0]

                        if n_cross_config != 0:                                     # if the swap is successful, do the swap
                            # print('Successful Swap!')
                            allpaths[i].append(path_2)
                            allpaths[i+1].append(path_1)
                        
                        else:                                                       # if the swap is unsuccessful, just append the same path to their own ensemble
                            # print('Unsuccessful Swap!')
                            allpaths[i].append(path_1)
                            allpaths[i+1].append(path_2)
        
            else:
                # print('Move = {} ; Performing Swap Move: 0-1, 2-3, 4-5, ...etc'.format(move), end='\r')
                
                # For even number of interfaces: 
                if (len(interfaces)-1)%2 == 0:
                    for i in range(0,len(interfaces)-1,2):
                        # print('i = ' ,i)
                        path_0 = allpaths[i][move-1]                                # looks at the most current path of the TPE of interface i
                        path_1 = allpaths[i+1][move-1]                              # looks at the most current path of the TPE of interface i+1
                        
                        # Check if path 1 crosses path 2 interface
                        n_cross_config = np.where(path_0[:,6]>interfaces[i+1])[0].shape[0]

                        if n_cross_config != 0:                                     # if the swap is successful, do the swap
                            # print('Successful Swap!')
                            allpaths[i].append(path_1)
                            allpaths[i+1].append(path_0)
                        
                        else:                                                       # if the swap is unsuccessful, just append the same path to their own ensemble
                            # print('Unsuccessful Swap!')
                            allpaths[i].append(path_0)
                            allpaths[i+1].append(path_1)
                        
                # For odd number of interfaces:
                else: 
                    
                    # Accounting for the last interface:
                    if np.where(allpaths[-1][-1][:,6]>interfaces[-2])[0].shape[0] == 0: 
                        reject_move = True
                        # print('Appending Rejected!')
                    else: 
                        allpaths[-1].append(allpaths[-1][-1])                       # append the previous path of the last ensemble to its own ensemble again 
                                                                                    # because there is no pair for this interface to swap with
                    
                    # Accounting for the rest of the interfaces:
                    for i in range(0,len(interfaces)-2,2):
                        # print('i = ' ,i)
                        path_0 = allpaths[i][move-1]                                # looks at the most current path of the TPE of interface i
                        path_1 = allpaths[i+1][move-1]                              # looks at the most current path of the TPE of interface i+1
                        
                        # Check if path 1 crosses path 2 interface
                        n_cross_config = np.where(path_0[:,6]>interfaces[i+1])[0].shape[0]

                        if n_cross_config != 0:                                     # if the swap is successful, do the swap
                            # print('Successful Swap!')
                            allpaths[i].append(path_1)
                            allpaths[i+1].append(path_0)

                        else:                                                       # if the swap is unsuccessful, just append the same path to their own ensemble
                            # print('Unsuccessful Swap!')
                            allpaths[i].append(path_0)
                            allpaths[i+1].append(path_1)
                        
        # Perform shooting move
        else: 
            # print('Move = {} ; Performing Shoot Move'.format(move), end='\r')

            for i in np.arange(len(interfaces)-1):                                  # i signifies the interface that I'm looking at
                # print('i = ', i)
                path = allpaths[i][move-1]                                          # looks at the most current path of the TPE of interface i
                reject_move = False
                lmax = round((len(path) + 1)/np.random.uniform())                   # l max = maximum path length for flexible shooting 
                
                # Pick random shooting point 
                if i == 0:
                    minus_frame = find_config_minus(path, basinA)                   # find the index corresponding to interface i-1    
                    plus_frame = find_config_plus(path, interfaces[i+1])            # find the index corresponding to interface i+1 
                    index = np.random.randint(minus_frame,plus_frame)

                # Pick random shooting point between i-1 and i+1
                else:               
                    minus_frame = find_config_minus(path, interfaces[i-1])          # find the index corresponding to interface i-1    
                    plus_frame = find_config_plus(path, interfaces[i+1])            # find the index corresponding to interface i+1  
                    index = np.random.randint(minus_frame,plus_frame)

                # Perturb the momenta
                shoot_point = copy.deepcopy(path[index])                            # gives me 1 phase space point/config along the path
                shoot_point[2] += np.random.normal(0,sigma)                         # px
                shoot_point[3] += np.random.normal(0,sigma)                         # py
                trial_path = np.asarray([shoot_point])                              # new trial path after perturbing the momenta

                if exp(-beta*(ld.kinetic_energy(shoot_point) - ld.kinetic_energy(path[index]))) < np.random.uniform():
                    # the if is testing to see if we want to use this new trial path (after momenta is perturbed) 
                    # before even shooting it
                    reject_move = True
                    # print('Shooting Move Rejected! Momenta Perturbation Failed!') 

                # Integrate backwards path if we have not rejected the move...
                if reject_move == False:
                    path_length = 1
                    trial_path[:,2] *=-1                                            # change px to the negative direction by * -1
                    trial_path[:,3] *=-1                                            # change py to the negative direction by * -1
                    trajstep = copy.deepcopy(trial_path[-1])                        # take the last config of trial path
                    op = ld.calc_op(op_type,trajstep[0],trajstep[1])                # calc op value for config

                    while op >= basinA:
                        trajstep = ld.vv_step(trajstep[:6],dt,beta,gamma,pes_type)   # ld.vv_step(x,y,px,py,fx,fy,dt,beta,gamma)
                        op = ld.calc_op(op_type,trajstep[0],trajstep[1])
                        trajstep = np.append(trajstep,[op],axis=0)
                        trial_path = np.append(trial_path,[trajstep],axis=0)
                        path_length +=1

                        # Reject if the maximum path length is exceeded
                        if path_length > lmax:
                            reject_move = True
                            # print('Shooting Move Rejected! Pathlength limit!')
                            break

                        # Reject if the backward path segment goes to B
                        if op >= basinB:
                            reject_move = True
                            # print('Shooting Move Rejected! Backward segment reaches B!')
                            break

                # Forward shooting             
                if reject_move == False: 
                    trial_path = np.flip(trial_path,axis=0)                          # flip it from the back to front of this trial path (because it was shooting backwards), 
                                                                                    # so we flip so that the last config we extract from here will be where we'll use to shoot forward
                    trial_path[:,2] *=-1
                    trial_path[:,3] *=-1
                    trajstep = copy.deepcopy(trial_path[-1])
                    op = ld.calc_op(op_type,trajstep[0],trajstep[1])
        
                    while op > basinA and op < basinB:                               # this is so that it restricts forward shooting to either end in A or B

                        nextstep = ld.vv_step(trajstep[:6],dt,beta,gamma,pes_type)
                        trajstep = nextstep
                        op = ld.calc_op(op_type,trajstep[0],trajstep[1])
                        trajstep = np.append(trajstep,[op],axis=0)
                        trial_path = np.append(trial_path,[trajstep],axis=0)
                        path_length +=1
                        
                        if path_length > lmax:
                            reject_move = True
                            # print('Shooting Move Rejected! Pathlength Limit!')
                            break

                # Final chance to reject a path (because no crossing of interface i)
                if np.where(trial_path[:,6]>interfaces[i])[0].shape[0] == 0:        # np.where(condition)[0] tells you which array or which path satisfies this condition and 
                                                                                    # the .shape[0] how many configs satisfy this, if 0, then we never cross the interface
                    reject_move = True
                    # print('Shooting Move Rejected! Did not cross its own interface!')
                
                # If we DON'T reject, then path becomes trial path
                if reject_move == False:
                    # print('Shooting Move Accepted!')
                    path = trial_path
                
                # Append current accepted path to its corresponding path ensemble
                allpaths[i].append(path)                                            # all paths contains all paths per each interface - i is to signify which interface 
                                                                                    # ensemble we're looking at. So allpaths[0] will contain all paths that cross interface lambda 0    
    # Plot sampled paths from each interface
    # Plot potential energy surface contours
    N = 100
    x_vec = np.linspace(-3.5, 3.5, N)
    y_vec = np.linspace(-3.5, 3.5, N)
    X, Y = np.meshgrid(x_vec, y_vec)
    energy = np.zeros((N, N))

    for i in range(len(x_vec)):
        for j in range(len(y_vec)):
            energy[j][i] = ld.potential(x_vec[i],y_vec[j],pes_type)
    
    # Plot the sampled paths for each interface and save the figures
    skip = 100

    for i in range(len(interfaces) - 1):
        plt.figure()

        plt.contour(x_vec,y_vec,energy,np.linspace(-3,3,20), cmap = 'jet')
            
        for j in range(1,len(allpaths[i]),skip):
            plt.plot(allpaths[i][j][:,0],allpaths[i][j][:,1])
                                                                                                # allpaths[i] tells you the interface ensemble you're in, 
                                                                                                # [j] tells you a specific path in that ensemble, [:,0/1] takes 
                                                                                                # all the x/y values of the paths in i
        # Plot basin boundaries

        if op_type == 1:
            plt.plot(np.linspace(basinA,basinA,10),np.linspace(min(y_vec),max(y_vec),10),color='r', label='Basin A')
            plt.plot(np.linspace(basinB,basinB,10),np.linspace(min(y_vec),max(y_vec),10),color='b', label='Basin B')
        elif op_type == 2:
            plt.plot(np.linspace(min(x_vec),max(x_vec),10),np.linspace(basinA,basinA,10),color='r', label='Basin A')
            plt.plot(np.linspace(min(x_vec),max(x_vec),10),np.linspace(basinB,basinB,10),color='b', label='Basin B')
        else:
            xplot = np.linspace(min(x_vec),max(x_vec),10)
            yplmax = basinB - xplot
            yplmin = basinA - xplot
            plt.plot(xplot,yplmin,color='r', label='Basin A')
            plt.plot(xplot,yplmax,color='b', label='Basin B')
            
        # Plot interfaces
        if op_type == 1:
            plt.plot(np.linspace(interfaces[i],interfaces[i],10),np.linspace(min(y_vec),max(y_vec),10), color='grey')
        elif op_type == 2:
            plt.plot(np.linspace(min(x_vec),max(x_vec),10),np.linspace(interfaces[i],interfaces[i],10), color='grey')
        else:
            xplot = np.linspace(min(x_vec),max(x_vec),10)
            yplot = interfaces[i]-xplot
            plt.plot(xplot,yplot, color='grey')    
        
        cbar = plt.colorbar(cm.ScalarMappable(cmap='jet'))
        cbar.set_ticks([])
        cbar.set_label(label = 'Energy', size=12)
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Interface Ensemble {}'.format(i))
        # plt.savefig('RETIS-paths-intf-{}'.format(i), dpi=600)
        plt.show
    return print("""\
          
            🎉 CONGRATULATIONS 🎉
                    
                """)


