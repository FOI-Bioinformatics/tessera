#!/usr/bin/env python3

from recan.simgen import Simgen
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from statistics import mean,median,stdev
import argparse
import shutil
import os
import sys


##### INPUTS
software_description = '''
                        Identify recombination events in a given multiple sequence alignment (MSA) file
                        Wrapper for Recan
                        '''

### Parse input arguments
## setup
parser = argparse.ArgumentParser(description=software_description)

parser.add_argument('-i','--input','--MSA','--msa',required=True,help='Path to multiple sequence alignment file')
parser.add_argument('-o','--output',required=True,help='Path to output directory')
parser.add_argument('-q','--query',required=True,help='Name of query to identiy recombination events (file name of query without extension, will try to strip the common extensions if provided)')

##/
## parse input
args = parser.parse_args()

msa_file_path = args.input
output_dir = args.output
recombination_dataset = args.query
##/
## format
if recombination_dataset.endswith('.gz'):           recombination_dataset = recombination_dataset.replace('.gz','')
for extension in ('.fna','.fasta','.fa',):
    if recombination_dataset.endswith(extension):   recombination_dataset = recombination_dataset.replace(extension,'')
##/

window_size = 1000
window_step = 100
dist_method_to_use = 'pdist'
#####/


### For INFO, define a dashed line that is terminal width. Define function to write formatted table
terminal_width = shutil.get_terminal_size().columns
separator_line = '-' * terminal_width

def print_formatted_table(data, header=None):
    num_columns = 1 # set default, if no header is passed
    # If header are not provided, generate default ones
    if header is None:
        header = ["Key"] + [f"Value {i + 1}" for i in range(num_columns)]
    else:
        num_columns = len(header)-1

    # Determine the maximum width of each column
    #key_width = max(len(str(key)) for key, _ in data) + 2  # Adding extra space for formatting
    
    # Prepare the expanded data
    expanded_data = []
    for key, value in data:
        if isinstance(value, list):
            # Split the value into sub-values and pad with None if necessary
            expanded_values = value + [None] * (num_columns - len(value))
            expanded_data.append([key] + expanded_values[:num_columns])
        else:
            expanded_data.append([key, value] + [None] * (num_columns - 1))
    
    # Calculate maximum width for each column, including header
    column_widths = [max(len(str(row[i])) for row in expanded_data if row[i] is not None) for i in range(num_columns + 1)]
    header_widths = [len(str(header[i])) for i in range(num_columns + 1)]

    # Adjust column widths to accommodate headers
    column_widths = [max(column_widths[i], header_widths[i]) + 2 for i in range(num_columns + 1)]
    
    # Define offset to use for table
    print_offset = '  '
    
    # Print the header
    print(print_offset+''.join(f"{str(header[i]):<{column_widths[i]}}" for i in range(num_columns + 1)))
    print(print_offset+'-' * sum(column_widths))  # Separator

    # Print each key-value pair
    for row in expanded_data:
        print(print_offset+''.join(f"{str(row[i]) if row[i] is not None else '':<{column_widths[i]}}" for i in range(num_columns + 1)))
###/

### Import MSA and run Recan
print('Importing MSA and executing Recan (compute dataset distance in windows over MSA)')
sim_obj = Simgen(msa_file_path)
sim_obj.get_info()

# determine index of query dataset (supply this to recan)
idx_recombination_dataset = None
for enum,record in enumerate(sim_obj.alignment_roll_window.align):
    if record.id == recombination_dataset:
        idx_recombination_dataset = enum
        break
if idx_recombination_dataset == None:
    print('FATAL: Could not determine query from MSA sequence. Please make sure to supply query name exactly as it is stated in the MSA.fasta (posted above), it should be the file-name of your query')
    print('Terminating!')
    sys.exit()
#/

sim_obj.simgen(window=window_size, shift=window_step, pot_rec=idx_recombination_dataset, dist=dist_method_to_use, )
###/

## Dump distance table
#sim_obj.save_data(out_name="sim_out") # will append .csv
##/

## Get data as dataframe (plotting) and "pythonic" format (for script)
distance_positions, dataset_distances = sim_obj.get_data(df=False)
dist_data = sim_obj.get_data()
##/

## Determine the closests datasets (to selected "pot_rec")
print(separator_line)
print('Getting distance stats...')
# for each window, get the most similar dataset
oDatasets_num_wins = {} # oDataset -> number of times declared window winner
oDatasets_num_wins_wTies = {} # --::-- with ties
for enum,pos in enumerate(distance_positions):
    # get distance for each dataset
    pos_dists = {} # dataset -> dist
    for dataset, distances in dataset_distances.items():
        dist_at_pos = distances[enum]
        pos_dists[dataset] = dist_at_pos
    #/
    # get the closest dataset
    closest_dataset_dist = sorted(pos_dists.items(),key=lambda x: x[1], reverse=True)[0][1]
    #/
    # get all datasets with distance equal to the closest dataset (i.e. allow ties)
    winners = []
    for dataset,dist in pos_dists.items():
        if dist == closest_dataset_dist:
            winners.append(dataset)
    #/
    # save
    if len(winners) == 1:
        dataset = winners[0]
        if not dataset in oDatasets_num_wins:       oDatasets_num_wins[dataset] = 0
        oDatasets_num_wins[dataset] += 1
    #/
    # save2
    for dataset in winners:
        if not dataset in oDatasets_num_wins_wTies:       oDatasets_num_wins_wTies[dataset] = 0
        oDatasets_num_wins_wTies[dataset] += 1
    #/

print('Dataset distance winners in each window (no ties):')
print_formatted_table(sorted(oDatasets_num_wins.items(),key=lambda x: x[1], reverse=True), header=['Dataset','Windows'])
print() # add empty line in-between tables
print('Dataset distance winners in each window (ties allowed):')
print_formatted_table(sorted(oDatasets_num_wins_wTies.items(),key=lambda x: x[1], reverse=True), header=['Dataset','Windows'])
print() # add empty line in-between tables
#/

# use overall distances to get closest datasets
oDatasets_dists_stats = {} # dataset -> data of ... [see arr below]
oDatasets_dists_stats_header = ["Dataset","Tot windows","Median distance","Distances >0","Distances >99%","Distances >95%","Distances >90%","Distances >80%","Distances >70%"]
for dataset, distances in dataset_distances.items():
    # compute stats
    dist_median = median(distances)
    tot_windows = len(distances)
    dists_filt = []
    dists_099 = []
    dists_095 = []
    dists_090 = []
    dists_080 = []
    dists_070 = []
    for dist in distances:
        if dist > 0:
            dists_filt.append(dist)
        if dist > 0.99:
            dists_099.append(dist)
        if dist > 0.95:
            dists_095.append(dist)
        if dist > 0.90:
            dists_090.append(dist)
        if dist > 0.80:
            dists_080.append(dist)
        if dist > 0.70:
            dists_070.append(dist)
    #/
    # save
    tmp_save = [tot_windows,round(dist_median,2),len(dists_filt),len(dists_099),len(dists_095),len(dists_090),len(dists_080),len(dists_070)]
    oDatasets_dists_stats[dataset] = tmp_save
    #/
#/
print('Dataset distance stats:')
print_formatted_table(sorted(oDatasets_dists_stats.items(),key=lambda x: (x[1][3],x[1][4],x[1][5],x[1][6]), reverse=True),header=oDatasets_dists_stats_header) # sort by 99/95/90/80/70 in order
print() # add empty line in-between tables
##/


### Output
## init outdir
print(separator_line)
print(f'Init output directory: {output_dir}')
if not os.path.exists(output_dir):          os.makedirs(output_dir)
##/

## Write distances table
tmp_output_table_path = output_dir+'/'+'distances.tsv'
print(f'Write distances table at {tmp_output_table_path}')
with open(tmp_output_table_path,'w') as nf:
    # write header
    header = oDatasets_dists_stats_header
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for dataset,colvalues in sorted(oDatasets_dists_stats.items(),key=lambda x: (x[1][3],x[1][4],x[1][5],x[1][6]), reverse=True): # sort by 99/95/90/80/70 in order
        writeArr = [dataset] + colvalues
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/

## Write number of windows where each dataset was a winner
tmp_output_table_path = output_dir+'/'+'distances_winners.tsv'
print(f'Write window winners (with ties) at {tmp_output_table_path}')
with open(tmp_output_table_path,'w') as nf:
    # write header
    header = ['Dataset','Number of windows (ties included)']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for dataset,num_windows in sorted(oDatasets_num_wins_wTies.items(),key=lambda x: x[1], reverse=True):
        writeArr = [dataset,num_windows]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
#/
##/

## Plot against top 5 queries
if 1:
    # Get top 5 candidates (from oDatasets_num_wins_wTies)
    datasets_to_include = []
    for dataset,num_wins_wTies in sorted(oDatasets_num_wins_wTies.items(),key=lambda x: x[1], reverse=True):
        datasets_to_include.append(dataset)
        if len(datasets_to_include) == 5: break
    #/
    
    # INFO-print
    print(f'Top x5 nearest datasets determined as: {", ".join(datasets_to_include)}')
    #/
    
    # Transpose the dataframe to make columns the x-axis and rows the data for plotting
    dist_data_transposed = dist_data.T
    
    # Get top x5 columns
    dist_data_transposed_extracted = dist_data_transposed[datasets_to_include]
    
    # Define line styles, markers, and colors
    line_styles = ['-', '--', '-.', ':']
    markers = ['o', 's', '^', 'X', 'D', 'P', 'v']
    
    cmap = mpl.colormaps['Set1'] # Get a colormap with enough distinct colors (it suggested tab10 but set1 works better when plotting 5 or less datasets)
    colors = cmap.colors
    
    # Combine all line styles, markers, and colors to create a unique style for each line
    unique_styles = [(ls, m, c) for ls in line_styles for m in markers for c in colors]
    np.random.shuffle(unique_styles)  # Shuffle the order to randomize
    
    # Plotting
    plt.figure(figsize=(12, 8))
    
    # Plot each row with a random marker interval
    marker_interval_tracker = set() # keep track of used marker intervals so that markers do not overlap
    for idx, row in enumerate(dist_data_transposed_extracted.columns):
        if idx >= len(unique_styles):  # In case there are more rows than unique styles
            break
        line_style, marker, color = unique_styles[idx]
        
        # Generate a random marker val and make sure the same interval is not used twice
        marker_interval = None
        iter_catch = 0
        while marker_interval == None or marker_interval in marker_interval_tracker:
            marker_interval = np.random.randint(3, 30)
            iter_catch += 1
            if iter_catch > 1000: break # try maximum this attempts to prevent endless loop for large datasets
        marker_interval_tracker.add(marker_interval)
        
        # Plot the line and markers together so both are included in the legend
        plt.plot(dist_data_transposed_extracted.index, dist_data_transposed_extracted[row], 
                 linestyle=line_style, marker=marker, color=color, markersize=4, alpha=0.5, 
                 markevery=marker_interval,  # Show markers at random intervals
                 label=row)
    
    # get mean and stdev for all plotted points. set ylim between upper 1.15 and lower "mean-stdev"
    mean_values = mean(dist_data_transposed_extracted.mean())
    std_values = mean(dist_data_transposed_extracted.std())
    plt.ylim(mean_values-std_values,1.15)
    
    # Adding labels and title
    plt.xlabel('Position (in MSA)')
    plt.ylabel('Similarity')
    plt.title('Similarity to query '+recombination_dataset)
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), title="Sequences", fontsize="small")
    
    # Adjust layout to make room for the legend
    plt.tight_layout()
    
    # Show the plot
    #plt.show()
    
    # Save plot
    plot_savepath = output_dir+'/'+'top5.pdf'
    plt.savefig(plot_savepath)
    print(f'Plot saved at: {plot_savepath}')
##/

## Plot specific
if 1:
    ## Plotting function
    def plot_pairwise(seq1_name=None,seq2_name=None,savepath=None):
        sns.set()
        
        fig_dist1 = plt.figure(figsize=(20, 8))
        
        xticks_interval = 20000
        
        plt.plot(dist_data.loc[seq1_name, : ], lw=7, alpha=0.5, label=seq1_name)
        plt.plot(dist_data.loc[seq2_name, : ], lw=7, alpha=0.5, label=seq2_name)
        
        plt.ylim(0.8, 1.15)
        plt.title(f"Similarity to {recombination_dataset}", fontsize=25)
        plt.ylabel("Similarity", fontsize=20)
        plt.xlabel("MSA position", fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xticks(np.arange(0, max(dist_data.columns), xticks_interval)) 
        
        plt.legend(prop={"size":20})
        plt.tight_layout()
        #plt.show()
        
        if savepath != None:
            plt.savefig(savepath)
    ##/
    
    ## Get top2 candidates to plot
    top2_datasets = []
    for dataset,num_wins_wTies in sorted(oDatasets_num_wins_wTies.items(),key=lambda x: x[1], reverse=True):
        top2_datasets.append(dataset)
        if len(top2_datasets) == 2: break
    
    print(f'Top x2 datasets (most likely for recombination) determined as: {" and ".join(top2_datasets)}')
    plot_savepath = output_dir+'/'+'top2.pdf'
    plot_pairwise(seq1_name=top2_datasets[0], seq2_name=top2_datasets[1], savepath=plot_savepath)
    print(f'Plot saved at: {plot_savepath}')
    ##/
###/

### Output interactive plot
if 1:
    import plotly.graph_objects as go
    
    # INFO-print
    print(f'Top x5 nearest datasets determined as: {", ".join(datasets_to_include)}')
    #/
    
    # Get top X candidates (from oDatasets_num_wins_wTies)
    datasets_to_include = []
    for dataset,num_wins_wTies in sorted(oDatasets_num_wins_wTies.items(),key=lambda x: x[1], reverse=True):
        datasets_to_include.append(dataset)
        if len(datasets_to_include) == 5: break
    #/
    # init figure
    fig = go.Figure()
    #/
    # add lines for each dataset
    for dataset in dist_data.index:
        # skip dataset if not part of topX
        if not dataset in datasets_to_include: continue
        #/
        # add to plot
        fig.add_trace(go.Scatter(x=dist_data.columns, y=dist_data.loc[dataset], mode='lines', name=dataset))
        #/
    #/
    # add layout to plot
    fig.update_layout(title=f'Similarity plot for query {recombination_dataset}', xaxis_title='Window', yaxis_title='Similarity', hovermode='closest', dragmode='zoom', showlegend=True)
    #/
    # save
    plot_savepath = output_dir+'/'+'plot.html'
    fig.write_html(plot_savepath)
    #/
    # INFO-print
    print(f'Plot saved at: {plot_savepath}')
    #/
##/
###/

### Finalize
print('All done!')
###/