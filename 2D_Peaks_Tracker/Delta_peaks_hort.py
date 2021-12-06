import numpy as np
from Delta_data_load import load_2d_data, get_gate1_lock1
from scipy.signal import find_peaks
from more_itertools import sort_together
import numpy as np
from numpy import transpose
import glob


def peaks_hort(gate1, lock1, find_troughs = False, sort_prominence = True, min_width=3, x_peak_max=10.0, x_peak_min=-10.0, max_peaks=2, height=0, prominence=None):
    
    peak_idx, properties = find_peaks(lock1, width=min_width, height=height, prominence=prominence)
    
    
    if find_troughs == True:
        #Multiply lock1 data by -1 so troughs become peaks, then find peaks
        negative_lock1 = [element*(-1) for element in lock1]
        peak_idx_troughs, properties_troughs = find_peaks(negative_lock1, width=min_width, height=height, prominence=prominence)
        
        #Revert trough heights to their original negative values
        properties_troughs["peak_heights"] = [peak*-1 for peak in properties_troughs["peak_heights"]]
        
        
        #combine peak and trough data
        peak_idx = np.append(peak_idx, peak_idx_troughs)
        properties = {key:np.hstack([properties[key], properties_troughs[key]]) for key in properties.keys()}
        
        
        #sort by x_value (peak_idx)
        properties_list = [properties[key] for key in properties.keys()]
        zipped_lists = zip(peak_idx, *properties_list)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        peak_idx, *properties_list = [ list(tuple) for tuple in  tuples]
        z = 0
        for key in properties:
            properties[key] = np.array(properties_list[z])
            z += 1
                                
        
    #only keep data points between x_peak_min and x_peak_max    
    boolean_array = ((gate1[peak_idx] < x_peak_max) & (gate1[peak_idx] > x_peak_min))
    x_peaks = gate1[peak_idx][boolean_array] 
    y_peaks = lock1[peak_idx][boolean_array]
    widths = properties["widths"][boolean_array] 
    prominences = properties["prominences"][boolean_array]
    
 
    
    if sort_prominence == True:
        #Sort all rows of arrays based on prominence 
        zipped_lists = zip(prominences, x_peaks, y_peaks, widths)
        sorted_pairs = sorted(zipped_lists, reverse = True)
        tuples = zip(*sorted_pairs)
        prominences, x_peaks, y_peaks, widths = [ list(tuple) for tuple in  tuples]
        
    
    
    # keep only the # of peaks = max_peaks(=2 by default) peaks
    x_peaks_low = x_peaks[0:max_peaks] 
    x_list = x_peaks_low
    
    y_peaks_low = y_peaks[0:max_peaks]
    y_list = y_peaks_low
    
    width_list = widths[0:max_peaks]
    
    prominence_list = prominences[0:max_peaks]
    
    
    #Sort remaining peaks by x_value
    zipped_lists = zip(x_list, y_list, width_list, prominence_list)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    x_list, y_list, width_list, prominence_list = [ list(tuple) for tuple in  tuples]
    
   
    return x_list, y_list, width_list, prominence_list

def find_peak_evol_hort(file_list, fixedGate3Val=0.268, min_width=1, find_troughs=False,
                   height=0, x_peak_max=10.0, x_peak_min=-10.0, max_peaks=2, prominence = None, sort_prominence = True):
   
    #goes to the file_list and gathers all peak information
 
    num_files = len(file_list) # number of files to be analyzed
    
    # np.full makes a matrix with given arguments as shapes with its elements given by the fill_value
    x_array = np.full((num_files, max_peaks), fill_value=np.nan) 
    y_array = np.full((num_files, max_peaks), fill_value=np.nan) 
    width_array = np.zeros((num_files, max_peaks))
    prominence_array = np.zeros((num_files, max_peaks))
    num_peaks_array = np.zeros(num_files, dtype=int) # no of peaks it found for each dat file
    
    valid = 0
    if find_troughs==True:
        print("Finding peaks AND troughs")
        
    if sort_prominence==True:
        print("Ordering peaks by prominence")
    
    
    for idx, filename in enumerate(file_list):
        print("\r{}".format(filename), end="")
        print(filename)
        try:
            File_1 = load_2d_data(filename) # loads the dat files
            #get the gate3 and lockin1 values from the dat file
            gate1, lock1 = get_gate1_lock1(file_data=File_1, fixedGate3Val=fixedGate3Val) 
            #gets all the information for one file
            x_list, y_list, width_list, prominence_list = peaks_hort(gate1, lock1, min_width=min_width, height=height,                                                                                       find_troughs=find_troughs, x_peak_max=x_peak_max, x_peak_min=x_peak_min, 
                                                                max_peaks=max_peaks, sort_prominence=sort_prominence, prominence=prominence)
            num_peaks = len(x_list) # gives the max numebr of peaks found
            print(num_peaks,end="")
            num_peaks_array[idx] = num_peaks
#             print(type(x_list))
            assert num_peaks <= max_peaks
            x_array[idx,:num_peaks] = x_list
            y_array[idx,:num_peaks] = y_list
            width_array[idx,:num_peaks] = width_list
            prominence_array[idx,:num_peaks] = prominence_list
                
            valid += 1
        except Exception as e:
            print(e)
    print("\nDone")
    print("Valid files #:", valid)
    
    return x_array, y_array, width_array, prominence_array, num_peaks_array

#     x_array = x_array[np.all(x_array, axis=1)]
#     y_array = y_array[np.all(y_array, axis=1)]
#     width_array = width_array[np.all(width_array, axis=1)]
#     prominence_array = prominence_array[np.all(prominence_array, axis=1)]
    
        
    return x_array, y_array, width_array, prominence_array, num_peaks_array



def find_max_2D(file_number, global_max = True, x_min=0, x_max=100, y_min = 0, y_max = 100, find_troughs=True):
    #Find and extract info of file name the user has passed
    filename = glob.glob("../*/*Sept2021_" + str(file_number) + ".dat", recursive=False)
    file_data = load_2d_data(filename[0], find_troughs=True)
    
    data = file_data.z
    
    
    if global_max == False:
        #Create a boolean array where only pixels inside desired dimensions are true   
        boolean_x = ((file_data.x < x_max) & (file_data.x > x_min))
        boolean_y = ((file_data.y < y_max) & (file_data.y > y_min))
        boolean_z = ((boolean_x) & (boolean_y))
        
        
        #Reduce data array by reducing the boolean array to only True values
        boolean_z, data = zip(*((x, y) for x, y in zip(boolean_z, data) if any(x)))        #Keep only rows with True values
        boolean_z = np.transpose(boolean_z)                                                #Transpose and repeat to keep only True columns 
        data = np.transpose(data)
        
        boolean_z, data = zip(*((x, y) for x, y in zip(boolean_z, data) if any(x)))
        boolean_z = np.transpose(boolean_z)                                                #Transpose back to original orientation
        data = np.transpose(data)
    
    
    #find max of 2D dataset
    zmax = np.amax(data)
    
    
    #If true, consider absolute value of all negatvie data as a maximum
    if find_troughs==True:
        inv_zmax = data * (-1)
        zmax_trough = -1 * np.amax(inv_zmax)
        zmax = max(zmax, zmax_trough, key = abs)
       
    #find coordinates of max
    where = np.where(file_data.z == zmax)
    coord = [index[0] for index in where]
    
    
    print("Max = ", str(zmax), " located at pixel: ", str(coord), "\nGate 2V =", str(file_data.x[coord[0], coord[1]]), "\nGate 3V =",                       str(file_data.y[coord[0], coord[1]]))
    
    return zmax, coord, file_data, data