from Delta_data import DatFile
import numpy as np

def load_2d_data(filename, find_troughs=True):
    """

    Args:
        filename ([type]): [description]

    Returns:
        [type]: [description]
    """
    df = DatFile(filename)

    file_data = df.get_data('gate 2 V meas', 
                            'gate 3 V meas', 
                            'sr860 x raw', 
                            'x (Right Gate Scan)',
                            'y (Left Gate Scan)')
    

    return file_data


def get_grids(file_data):
    gate1_grid = file_data.x # gets all the gate1 grid values
    gate3_grid = file_data.y # gets all the gate3 grid values
    lockin1_2D = file_data.z # gets lockin1 2D grid for all values of gate1 and Gate3
    
    return gate1_grid, gate3_grid, lockin1_2D


def get_gate1_lock1(file_data, fixedGate3Val):
    gate1 = file_data.x[:,0] 
    gate3 = file_data.y[0,:]
    # set a fixed value for gate1
    fiedGate3Val = fixedGate3Val
    idx_val = np.argmin(np.abs(gate3-fixedGate3Val)) #np.abs calculates the absolute value, np.argmin gives the indices
                                                    # of the minimums
    # print(idx_val)
    lock1 = file_data.z[:, idx_val]
    return gate1, lock1


def get_gate3_lock1(file_data, fixedGate1Val):
    gate1 = file_data.x[:,0] 
    gate3 = file_data.y[0,:]
    # set a fixed value for gate1
    fixedGate1Val = fixedGate1Val
    idx_val = np.argmin(np.abs(gate1-fixedGate1Val)) #np.abs calculates the absolute value, np.argmin gives the indices
                                                    # of the minimums
    # print(idx_val)
    lock1 = file_data.z[idx_val,:]
    return gate3, lock1

def get_lock1(file_data, fixedGate1Val, fixedGate3Val):
    '''Gives you the values of the lockin1 when both gate1 and gate3 are fixed'''
    gate1 = file_data.x[:,0] 
    gate3 = file_data.y[0,:]
    # set a fixed value for gate1
    
    idx_val_fix3 = np.argmin(np.abs(gate3-fixedGate3Val)) #np.abs calculates the absolute value, np.argmin gives the indices
                                                    # of the minimums
    
    idx_val_fix1 = np.argmin(np.abs(gate1-fixedGate1Val))
    # print(idx_val)
    lock1 = file_data.z[idx_val_fix1, idx_val_fix3]
    return  lock1

def many_lock1_values(file_data, gate1_values, gate3_values):
    lock1_values = {}
    for val1 in gate1_values:
        for val2 in gate3_values:
            lock1 = get_lock1(file_data, val1, val2)
            lock1_values[(val1, val2)] = lock1
    return lock1_values
    
        