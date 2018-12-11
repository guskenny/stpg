#!/usr/bin/python

#TODO: make separate functions for -nodisplay

import sys;
import os.path;
import Tkinter as tk;
import random as r;
import math;

# program constants
WINDOW_SIZE = 800;
SKIP_FRAMES = 10;
NODE_SIZE = 8;
END_CROSS = 80;

# physical constants
ALPHA = 1.0;
BETA = 0.0001;
K = 1.0;
DELTA_T = 0.01;
ETA = 0.99;

# test if arguments
if len(sys.argv) <= 1:
    print "No STP instance found!\nExiting...\n";
    exit();

stp_filename = str(sys.argv[1]) + ".stp";
tree_filename = str(sys.argv[1]) + ".tree";

# test if instance exists
if not os.path.isfile(stp_filename):
    print "No STP instance: " + stp_filename + " found!\nExiting...\n";
    exit();

def load_graph(fname):
    with open(fname, 'r') as infile:
        line = "";
        while not line.startswith('END'):
            line = next(infile);
            continue;
        
        for i in range(3):
            line = next(infile);
            
        split_line = line.split();
        n_nodes = int(split_line[1]);
        line = next(infile);
        split_line = line.split();
        n_edges = int(split_line[1]);

        G = [[0.0 for i in range(n_nodes)] for j in range(n_nodes)];

        line = next(infile);

        while not line.startswith('END'):
            if line.startswith('#'):
                line = next(infile);
                continue;
            split_line = line.split();
            G[int(split_line[1])-1][int(split_line[2])-1] = int(split_line[3]) * 0.01;
            G[int(split_line[2])-1][int(split_line[1])-1] = int(split_line[3]) * 0.01;
            line = next(infile);

        for i in range(3):
            line = next(infile);
         
        split_line = line.split();
        n_terminals = int(split_line[1]);

        T = [0 for i in range(n_nodes)];

        line = next(infile);

        while not line.startswith('END'):
            split_line = line.split();
            T[int(split_line[1])-1] = 1;
            line = next(infile);

        return G, T, n_nodes, n_edges, n_terminals;

def load_tree(fname,n_nodes):
    with open(fname, 'r') as infile:
        tree = [[0 for i in range(n_nodes)] for j in range(n_nodes)];
        
        for line in infile:         
            split_line = line.split();
            tree[int(split_line[0])][int(split_line[1])] = 1;
            tree[int(split_line[1])][int(split_line[0])] = 1;
            
        return tree;

G, T, n_nodes, n_edges, n_terminals = load_graph(stp_filename);

if os.path.isfile(tree_filename):
    tree = load_tree(tree_filename,n_nodes);
else:
    tree = [[0 for i in range(n_nodes)] for j in range(n_nodes)];


'''
# generate grid graph
G = [];

n_rows = 5;
n_cols = 5;

for i in xrange(n_rows * n_cols):
    ci = i % n_cols;
    ri = i / n_cols;
    dr = [];
    for j in xrange(n_rows * n_cols):
        cj = j % n_cols;
        rj = j / n_cols;
        if (ci == cj and (ri == rj-1 or ri == rj+1)) or\
                (ri == rj and (ci == cj-1 or ci == cj+1)):
            dr.append(0.01);
        else:
            dr.append(0.0);
    G.append(dr); 
n_nodes = length(G);
T = [0 for i in xrange(n_nodes));
'''

# initialise Tkinter
root = tk.Tk();

# initialise canvas
canvas = tk.Canvas(root, width=WINDOW_SIZE, height=WINDOW_SIZE, background="black");

# open canvas
canvas.pack();

# variables for position and velocity vectors
x = [];
v = [];

# variables for crossings
old_crossings = 0;
cross_count = 0;

# variables for nodes and edges
node_ids = [];
edge_ids = [];

# function to calculate the Coulomb force
def F_c(xi, xj):
    dx = xj[0] - xi[0];
    dy = xj[1] - xi[1];
    ds2 = dx*dx + dy*dy;
    ds = math.sqrt(ds2);
    const = 1.0 * BETA / (ds2*ds)
    return [-const*dx, -const*dy];

# function to calculate the Hooke's force
def F_h(xi, xj, dij):
    dx = xj[0] - xi[0];
    dy = xj[1] - xi[1];
    ds = math.sqrt(dx*dx + dy*dy);
    dl = ds - dij;
    const = K * dl/ds;
    return [const*dx, const*dy];   

# function to calculate sum of force vectors in system
def move():
    global cross_count;
    global old_crossings;

    # outer loop to speed up calculations without rendering
    for frame in xrange(SKIP_FRAMES):
        for i in xrange(n_nodes):
            F_x = 0.0;
            F_y = 0.0;
            for j in xrange(n_nodes):
                if j == i:
                    continue;
                G_ij = G[i][j];
                # if not connected, just calculate coulomb force
                if G_ij == 0.0:
                    F_ij = F_c(x[i], x[j]);
                # if connected, calculate Hooke's force
                else:
                    F_ij = F_h(x[i], x[j], G_ij);

                # sum forces
                F_x += F_ij[0];
                F_y += F_ij[1];
            
            # update velocity vectors    
            v[i][0] = (v[i][0] + ALPHA * F_x * DELTA_T) * ETA;
            v[i][1] = (v[i][1] + ALPHA * F_y * DELTA_T) * ETA;
        
        # update position vectors
        for i in xrange(n_nodes):
            x[i][0] += v[i][0] * DELTA_T;
            x[i][1] += v[i][1] * DELTA_T;

    crossings = compute_crossings();
    
    if crossings == old_crossings:
        cross_count += 1;
        print "cross count: " + str(cross_count);
    else:
        cross_count = 0;

    stop_flag = False
    if cross_count == END_CROSS:
        print "system stable"
        print "crossings: " + str(crossings);
        update_graph();
        return # return if not quitting

    old_crossings = crossings;


    if '-nodisplay' not in sys.argv:
        update_graph();
        root.after(1, move);
    else:
        move();

def update_graph():
    # render nodes and edges in canvas
    ei = 0;
    for i in xrange(n_nodes):
        #move_node(i);
        for j in xrange(i+1, n_nodes):
            if G[i][j] != 0:
                id = edge_ids[ei];
                move_edge(id, x[i], x[j]);
                ei += 1;
    
    for i in xrange(n_nodes):
        move_node(i);


# compute crossings

# TODO: FIX 3D ITERATIONS!!!!!!
#       maybe by creating an edge list with node indicies
#       that is probably best..

def compute_crossings():
    crossings = 0;
    #for i in xrange(n_edges):
    for i in xrange(n_nodes):
        for j in xrange(i+1, n_nodes):
            if G[i][j] != 0:
                #edge1 = canvas.coords(edge_ids[i]);
                edge1_i = [temp * WINDOW_SIZE for temp in x[i]];
                edge1_j = [temp * WINDOW_SIZE for temp in x[j]];
                edge1 = edge1_i + edge1_j;
        #for j in xrange(i+1,n_edges):
                for k in xrange(n_nodes):
                    for l in xrange(k+1, n_nodes):
                        if G[k][l] != 0:
                            #edge2 = canvas.coords(edge_ids[j]);
                            edge2_k = [temp * WINDOW_SIZE for temp in x[k]];
                            edge2_l = [temp * WINDOW_SIZE for temp in x[l]];
                            edge2 = edge2_k + edge2_l;
#            if check_b_boxes([edge1[0], edge1[1]],\
#                    [edge1[2], edge1[3]], [edge2[0], edge2[1]],\
#                    [edge2[2], edge2[3]]): 
                            crossings += check_crossing([edge1[0], edge1[1]],\
                                [edge1[2], edge1[3]], [edge2[0], edge2[1]],\
                                [edge2[2], edge2[3]]);

    return crossings;


# function to check if bounding boxes of line segments intersect
def check_b_boxes(a,b,c,d):
    # order coordinates
    a_x = min(a)
    a_y = max(a)
    b_x = max(b)
    b_y = min(b)
    c_x = min(c)
    c_y = max(c)
    d_x = max(d)
    d_y = min(d)

    if (a_x < d_x) and (b_x > c_x) and (a_y < d_y) and (b_y > c_y):
        return 1
    else:
        return 0

# function to check orientation of ordered points
def check_orientation(a,b,c):
    if (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0]):
        return 1;
    else:
        return 0;

# function to check if two segments cross
def check_crossing(a, b, c, d):
    if (check_orientation(a,c,d) != check_orientation(b,c,d)) and\
       (check_orientation(a,b,c) != check_orientation(a,b,d)):
        return 1;
    else:
        return 0;

# function to check if two lines cross
#def check_crossing(e_x, e_y, f_x, f_y, p_x, p_y, q_x, q_y):
#    cross1 = (f_x - e_x)*(p_y - f_y) - (f_y - e_y)*(p_x - f_x);
#    cross2 = (f_x - e_x)*(q_y - f_y) - (f_y - e_y)*(q_x - f_x);
#    if cross1 * cross2 >= 0:
#        return 0
#    else: 
#        return 1

# function to move node on screen
def move_node(i):
    newx = int(x[i][0] * WINDOW_SIZE);
    newy = int(x[i][1] * WINDOW_SIZE);
    # canvas.coords(node_ids[i], newx-NODE_SIZE, newy-NODE_SIZE,\
            # newx+NODE_SIZE, newy+NODE_SIZE);
    canvas.coords(node_ids[i], newx, newy);

# function to move edge on screen
def move_edge(i, xi, xj):
    newx1 = int(xi[0] * WINDOW_SIZE);
    newy1 = int(xi[1] * WINDOW_SIZE);
    newx2 = int(xj[0] * WINDOW_SIZE);
    newy2 = int(xj[1] * WINDOW_SIZE);
    canvas.coords(i, newx1, newy1, newx2, newy2);

# initialise nodes    
for i in xrange(n_nodes):
    xi = [r.random(), r.random()];
    x.append(xi);
    v.append([0.0, 0.0]);
    if T[i] != 1:
        # id = canvas.create_oval(0, 0, 0, 0, fill = "white");
        id = canvas.create_text(0, 0,text=str(i), fill = "white");
    else:
        # id = canvas.create_oval(0, 0, 0, 0, fill = "red");
        id = canvas.create_text(0, 0, text=str(i), fill = "red");
    node_ids.append(id);
    move_node(i);

# initialise edges    
for i in xrange(n_nodes):
    for j in xrange(i+1, n_nodes):
        if G[i][j] != 0:
            if tree[i][j] != 0:
                id = canvas.create_line(0, 0, 0, 0, fill = "blue",width=3);
            else:
                id = canvas.create_line(0, 0, 0, 0, fill = "white");
            edge_ids.append(id);
            move_edge(id, x[i], x[j]);

root.after(1000,move);   
 
root.mainloop();
