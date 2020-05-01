import matplotlib.pyplot as plt
import numpy as np
import math 
    
def emg_threshold(emg_signal):
    noise = emg_signal[0:100]
    threshold = 3*np.std(noise)
    return threshold 

def emg_smooth(emg_signal,n):
    result = np.zeros(len(emg_signal))
    for i in range(0,len(emg_signal)):
        sum=0
        for j in range(0,n):
            if(i+j<len(emg_signal)):
                sum += emg_signal[i+j]  
        result[i] = sum*(1/n)
        i+=n
    return result

def detect_muap(original_signal,emg_signal,threshold,t):
    peaks_indices = np.array([], dtype=int)
    muap_start_index = 0
    muap_end_index = 0
    exceeds_threshold = False # to be set to true when a sample exceeds a threshold
    current_t = 0
    for i in range(0,len(emg_signal)):
        if(emg_signal[i]>threshold):
            if(exceeds_threshold==False):
                muap_start_index = i
                exceeds_threshold = True                
            current_t += 1
        else:
            if(exceeds_threshold==True): # meaning: the previous sample is an end point
                if(current_t>=t):
                    muap_end_index = i
                    current_t = 0 # reset current_t for future muap detections
                    current_muap = original_signal[muap_start_index:muap_end_index]
                    peak_index = np.argmax(current_muap)+muap_start_index
                    peaks_indices = np.append(peaks_indices,peak_index)
                exceeds_threshold = False
                current_t = 0            
    return peaks_indices

def template_matching(emg_signal,peaks_indices,diff_th,t):
    templates=[[]]
    similar_muap_timestamps=[[]]
    
    for i in range(0,len(peaks_indices)):
        current_muap = emg_signal[peaks_indices[i]-(int)(t/2):peaks_indices[i]+(int)(t/2)+1]
        template_matched = False
        if(len(templates[0])==0): # first peak, so add it to the templates & to the first array of similar_muap_timestamps
            templates[0]=current_muap.tolist()
            similar_muap_timestamps[0]=[peaks_indices[i]]
        else: 
            for j in range(0,len(templates)):
                template_asnumpy = np.asarray(templates[j])
                d = np.sum(np.square(np.subtract(current_muap,template_asnumpy)))
                if(d<diff_th): # the current_muap belongs to this template j
                    templates[j] = ((template_asnumpy + current_muap)/2).tolist()
                    similar_muap_timestamps[j] += [peaks_indices[i]]
                    template_matched = True
                    break
            if(template_matched==False):
                templates += [current_muap.tolist()]
                similar_muap_timestamps += [[peaks_indices[i]]]                        
    return templates,similar_muap_timestamps

def emg_decompose(emg_signal, t, diff_th):
    rectified_signal = np.abs(emg_signal)
    threshold = emg_threshold(rectified_signal) 
    smoothed_signal = emg_smooth(rectified_signal,t)
    peaks_indices = detect_muap(emg_signal,smoothed_signal,threshold,t)
    templates, similar_muap_timestamps = template_matching(emg_signal,peaks_indices,diff_th,t)
    
    # Plotting the templates
#     fig = plt.figure()

#     plt.subplot(2, 2, 1)
#     plt.plot(templates[0], color='green')
    
#     plt.subplot(2, 2, 2)
#     plt.plot(templates[1], color='red')

#     plt.subplot(2, 2, 3)
#     plt.plot(templates[2], color = 'magenta')
    
          
    # Plotting the emg_signal with detected peaks
    similar_muap_timestamps_sliced = [[]]*len(similar_muap_timestamps)
    for i in range(0,len(similar_muap_timestamps)):
        sim_muaps = []
        for j in range(0,len(similar_muap_timestamps[i])):
            if((similar_muap_timestamps[i][j]>=30000) and (similar_muap_timestamps[i][j]<=35000)):
                sim_muaps += [similar_muap_timestamps[i][j]]
        similar_muap_timestamps_sliced[i] = sim_muaps
    
    colors = ['g','r','m','y']
    for i in range(0,len(similar_muap_timestamps_sliced)):
        y = []
        for j in range(0,len(similar_muap_timestamps_sliced[i])):
            y += [900]
        plt.scatter(similar_muap_timestamps_sliced[i], y, marker='*', color=colors[i])

    emg_sliced = emg_signal[30000:35001]
    x_axis = []
    for i in range(30000,35001):
        x_axis += [i]
    plt.plot(x_axis,emg_sliced)
    
# Saving the figure
    figure = plt.gcf() 
    figure.set_size_inches(20, 6)
    plt.savefig("DetectedMUAP_12.5.jpg", dpi = 100)

if __name__ == "__main__":
    emg_signal = np.loadtxt(fname = "Data.txt")
    emg_decompose(emg_signal,20,math.pow(12.5,5))
    




