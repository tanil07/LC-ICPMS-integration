#Anil Timilsina & Rene Boiteau
#Change the data directory to the location of the csv files for the station you want to plot
#This script will create a pdf figure for each element in the csv files in the selected directory
#Plots will have a line for each file in the directory with the element on the y-axis and the time on the x-axis

import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
data_dir=r'C:\Users\Chemistry\Desktop\UMN Research\ICPMS_analysis\ICP_MS data\gp17-ant_env_icpms\20250206\Samples'

os.chdir(data_dir)


def baseline_subtract_linear(time, intensity, start_time, end_time, time_interval):
    """
    Performs baseline subtraction on a pandas DataFrame by linearly
    interpolating between average intensities at specified start and end times.

    Args:
        'intensity': list of intensities
        'time': list of times
        start_time (float): Start time for baseline determination.
        end_time (float): End time for baseline determination.
        time_interval (float): Time interval for averaging intensities at start and end.

    Returns:
        pandas.DataFrame: DataFrame with an added 'baseline_subtracted_intensity' column.
    """       

    
    df=pd.DataFrame({'time':time, 'intensity':intensity})

    start_avg = df[(df['time'] >= start_time - time_interval / 2) & (df['time'] <= start_time + time_interval / 2)]['intensity'].mean()
    end_avg = df[(df['time'] >= end_time - time_interval / 2) & (df['time'] <= end_time + time_interval / 2)]['intensity'].mean()

    baseline = np.interp(df['time'], [start_time, end_time], [start_avg, end_avg])
    df['baseline_subtracted_intensity'] = df['intensity'] - baseline
    return df


def integration(df, start_time, end_time):
    """
    Integrate the area under the curve of a pandas DataFrame between start_time and end_time.

    Args:
        df (pandas.DataFrame): DataFrame with 'intensity' and 'time' columns.
        start_time (float): Start time for integration.
        end_time (float): End time for integration.

    Returns:
        float: Integrated area under the curve.
    """
    df2 = df[(df['time'] > start_time) & (df['time'] <= end_time)]
    area= np.trapz(df2['baseline_subtracted_intensity'], x=df2['time'])
    return area



def plot_chromatograms(df, start_time, end_time, element):
    """
    Plots the original chromatogram, baseline-subtracted chromatogram,
    and vertical lines for start and end times.

    Args:
        df (pandas.DataFrame): DataFrame with 'time', 'intensity', and
                                'baseline_subtracted_intensity' columns.
        start_time (float): Start time for baseline determination.
        end_time (float): End time for baseline determination.
    """
    if not isinstance(df, pd.DataFrame) or not all(col in df.columns for col in ['intensity', 'time', 'baseline_subtracted_intensity']):
        raise ValueError("Input must be a pandas DataFrame with 'intensity', 'time', and 'baseline_subtracted_intensity' columns.")

    fig, ax = plt.subplots(figsize=(10, 6)) # Create figure and axes objects

    ax.plot(df['time'], df['intensity'], label="Original Chromatogram")
    ax.plot(df['time'], df['baseline_subtracted_intensity'], label="Baseline-Subtracted Chromatogram")

    ax.axvline(x=start_time, color='r', linestyle='--', label=f"Start Time ({start_time})")
    ax.axvline(x=end_time, color='g', linestyle='--', label=f"End Time ({end_time})")

    ax.set_xlabel("Time")
    ax.set_ylabel("Intensity")
    ax.legend()
    ax.set_title("Chromatogram with Baseline Subtraction:"+element)
    ax.grid(True)
    plt.show()
    return fig # return the figure object

all_files = os.listdir(data_dir)
flist = [f for f in all_files if 'subtracted.csv' in f] ###change to 'GT' if not using blank substracted files
files={}
for f in flist:
    files[f]=pd.read_csv(f)
    elements = [col for col in files[f].columns if col!= 'Number' and 'Time' not in col]


results=[]
baseline_start=300
baseline_stop=3000

#elements=['56Fe','60Ni','63Cu', '66Zn','59Co','127I','31P','79Br','55Mn','27Al']

times=[[300,3000]]
#elements=['56Fe']

plots = []

for key in files:
    print(key)
    for timerange in times:
        current={}
        current['File']=key 
        current['start_time']=timerange[0]
        current['end_time']=timerange[1]
        current['correction']=0
        for element in elements:
            time = 'Time ' + element
            df=baseline_subtract_linear(files[key][time], files[key][element], baseline_start, baseline_stop, 15)
            #plot_chromatograms(df, timerange[0], timerange[1],element)
            area=integration(df, timerange[0], timerange[1])
            current[element]=area
            
        results.append(current)


element='56Fe'

times=[[300,1000],
       [1000,2000],
       [2000,3000]]


for key in files:

    for timerange in times:
        current={}
        current['File']=key 
        current['start_time']=timerange[0]
        current['end_time']=timerange[1]
        current['correction']=0
        time = 'Time ' + element
        df=baseline_subtract_linear(files[key][time], files[key][element], baseline_start, baseline_stop, 15)
        #plot_chromatograms(df, timerange[0], timerange[1],element)
        area=integration(df, timerange[0], timerange[1])
        current[element]=area
            
        results.append(current)


element='59Co'

times=[[800,1000]]


for key in files:

    for timerange in times:
        current={}
        current['File']=key 
        current['start_time']=timerange[0]
        current['end_time']=timerange[1]
        current['correction']=0
        time = 'Time ' + element
        df=baseline_subtract_linear(files[key][time], files[key][element], baseline_start, baseline_stop, 15)
        #plot_chromatograms(df, timerange[0], timerange[1],element)
        area=integration(df, timerange[0], timerange[1])
        current[element]=area
        
        results.append(current)


results=pd.DataFrame(results)

'''cali_table = pd.read_csv('cali_table.csv')
cali_table.set_index(cali_table.columns[0], inplace=True) 

for element in elements:
    if element in results.columns:
        if element in cali_table.columns:
            if cali_table[element]['r2'] > 0.96:
                results[f'{element}_conc'] = (results[element] - cali_table[element]['intercept']) / cali_table[element]['slope']


results.to_csv('results.csv')'''
