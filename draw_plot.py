#!/usr/bin/python
import matplotlib.pyplot as plt

import sys
import os

global args
args = sys.argv

# constants
LABEL = 1;
VALUE = 2;

# order of colours
COLOURS = ['b','r','g','c','m','y','k']

# double check there are files
if len(args) <= 1:
  print "Please enter csv files";
  exit();

# container for y values
y_vals = []

# containers for metadata
lines=[1 for _ in range(len(args)-1)]
names=["" for _ in range(len(args)-1)]
pop_size=[0 for _ in range(len(args)-1)]
swaps=[0 for _ in range(len(args)-1)]
num_merges=[0 for _ in range(len(args)-1)]
restarts=[0 for _ in range(len(args)-1)]
cpu_times=[0 for _ in range(len(args)-1)]
wall_times=[0 for _ in range(len(args)-1)]
groups=[0 for _ in range(len(args)-1)]
max_sub_swaps=[0 for _ in range(len(args)-1)]
mip_timeout=[0 for _ in range(len(args)-1)]
mip_rel_gap=[0 for _ in range(len(args)-1)]

wall_times_idx = -1
times_idx = -1
restarts_idx = -1
skip_times = False

y_min = 999999;

for infile_name in args[1:]:
    with open(infile_name, 'r') as infile:
        y = []
        for line in infile:
            if line.startswith("#"):
                continue
            if line.startswith("!"):
                split_line = line.split();
                if split_line[LABEL] == "LINE":
                    lines[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "CPU_TIME":
                    cpu_times[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "WALL_TIME":
                    wall_times[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "MIP_REL_GAP":
                    mip_rel_gap[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "NAME":
                    names[args.index(infile_name)-1] = split_line[VALUE]
                if split_line[LABEL] == "POP_SIZE":
                    pop_size[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "TIMES":
                    times_idx = args.index(infile_name)-1
                if split_line[LABEL] == "WALL_TIMES":
                    wall_times_idx = args.index(infile_name)-1
                if split_line[LABEL] == "RESTART":
                    restarts_idx = args.index(infile_name)-1
                if split_line[LABEL] == "SWAPS":
                    swaps[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "MAX_SUB_SWAPS":
                    max_sub_swaps[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "MIP_TIMEOUT":
                    mip_timeout[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "NUM_MERGES":
                    num_merges[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "JUMP_FREQ":
                    restarts[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "GROUPS":
                    groups[args.index(infile_name)-1] = 1
                    skip_times = True
                continue

            try:
                y.append(float(line))
            except ValueError:
                continue;
	    if args.index(infile_name) == 1:
            	y_min = min(y_min, float(line))

    y_vals.append(y)

box_text = "DEAD_MOVES: " + str(swaps[0])\
    + "\nNUM_JUMPS: " + str(restarts[0])\
    
    # + "\nMAX_SUB_SWAPS: " + str(max_sub_swaps[0])
# if restarts_idx > -1:
#     box_text = box_text + "\nNUM_MERGES: " + str(int(len(y_vals[0])/pop_size[0] - y_vals[restarts_idx][-1]))\
#     + "\nRESTARTS: " + str(int(y_vals[restarts_idx][-1]))
# else:
#     box_text = box_text + "\nNUM_MERGES: " + str(int(len(y_vals[0])/pop_size[0] - restarts[0]))\
#     + "\nRESTARTS: " + str(restarts[0])

# box_text = box_text +"\nMIP_TIMEOUT: " + str(mip_timeout[0]) + " sec"\
#                     +"\nMIP_REL_GAP: " + str(mip_rel_gap[0]*100) + "%"

# if wall_times_idx > -1:
#     box_text = box_text + "\nWALL_TIME: " + str(y_vals[wall_times_idx][-1]) + " sec"
# else:
#     box_text = box_text + "\nWALL_TIME: " + str(wall_times[0]) + " sec"

# if times_idx > -1:
#     box_text = box_text + "\nCPU_TIME: " + str(y_vals[times_idx][-1]) + " sec"
# else:
#     box_text = box_text + "\nCPU_TIME: " + str(cpu_times[0]) + " sec"

fig,ax1 = plt.subplots(figsize=(12,8))

ax1.set_ylabel('Objective value')

ax2 = plt.twinx()

ax2.yaxis.set_label_position("right")
ax = ax1.twinx()

x_max = 0

ax1.set_xlabel('Iterations')
plt.title("Plot of " + names[0] + " - MIN: " + "{:.3E}".format(y_min))

for y_idx in range(len(y_vals)):
    if y_idx == wall_times_idx or y_idx == restarts_idx or (y_idx == times_idx and skip_times is True):
        continue
    line_type = ''
    line_style = 'solid'
    marker_size = 0
    if  lines[y_idx] == 1:
        line_type = 'o'
        marker_size = 4
        line_style = 'None'

    if groups[y_idx] == 1:
        ax2.plot([x * pop_size[y_idx] + pop_size[y_idx] for x in range(len(y_vals[y_idx]))],y_vals[y_idx], linewidth=lines[y_idx],color='gold', marker=line_type, mew=0, ls=line_style, ms=marker_size)

        ax2.set_ylabel('Number of groups per merge')
    elif y_idx == times_idx and skip_times is False:
         ax2.plot(range(len(y_vals[y_idx])),y_vals[y_idx], linewidth=lines[y_idx],color='gold', marker=line_type, mew=0, ls=line_style, ms=marker_size)

         ax2.set_ylabel('CPU time [sec]')
    else:
        ax1.plot(range(len(y_vals[y_idx])),y_vals[y_idx], linewidth=lines[y_idx],color=COLOURS[y_idx], marker=line_type, mew=0, ls=line_style, ms=marker_size)
        x_max = max(x_max,len(y_vals[y_idx]))

ax1.minorticks_on()
ax1.grid(True,which='major')
ax1.grid(True,which='minor', linestyle='--',alpha=0.5)

ax.yaxis.set_visible(False)
ax.text(0.96,0.96, box_text,transform=ax.transAxes, fontsize=14,
bbox={'facecolor':'white', 'alpha':0.5, 'pad':10}, ha='right', va='top')

ax1.set_xlim(right=x_max)

plt.tight_layout()
plt.rcParams["savefig.directory"] = os.getcwd()
plt.show()

filename = './tracking_data/' + names[0] + '_' + str(pop_size[0]) + '_' + str(swaps[0]) + '_' + str(num_merges[0]) + '.png'

fig.savefig(filename)
