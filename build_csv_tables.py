#!/usr/bin/python
import sys, optparse, os, csv, math

data = {}
opt = {}

if len(sys.argv) < 5:
  print("Please enter csv files in order LS, MIP, GRASP, MERGE");
  exit();

for in_file in sys.argv[1:]:

	inFile = open(str(in_file))
	reader = csv.DictReader(inFile)

	for row in reader:
		if row["Instance"]:
			try:
				data[row["Instance"]][in_file].append(int(row["Best obj"]))
			except: 
				try:
					data[row["Instance"]][in_file] = [int(row["Best obj"])]
				except:
					data[row["Instance"]] = {}
					data[row["Instance"]][in_file] =[int(row["Best obj"])]

			if row["Instance"] not in opt:
				opt[row["Instance"]] = {}
				for file in sys.argv[1:]:
					opt[row["Instance"]][file] = False
				opt[row["Instance"]]["value"] = row["Optimal obj"]

			if row["Found optimal"]:
				opt[row["Instance"]][in_file] = True
				
	inFile.close()

for inst,file_data in data.items():
	avgs = {}
	square_diffs = {}
	std_devs = {}
	for file in sys.argv[1:]:
		avgs[file] = float(sum(file_data[file])/len(file_data[file]))
		square_diffs[file] = [math.pow(x - avgs[file],2) for x in file_data[file]]
		std_devs[file] = '%.2f'%(math.sqrt(float(sum(square_diffs[file]))/len(square_diffs[file])))
	# std_dev = math.sqrt(float(sum(square_diff))/len(square_diff))
	key_min = min(avgs.keys(), key=(lambda k: avgs[k]))
	min_avg = avgs[key_min]
	print("\\texttt{"+inst+"} & "+opt[inst]["value"])
	for file in sys.argv[1:]:
		if avgs[file] == min_avg:
			if opt[inst][file]:
				print("& "+"\\best{\\opt{"+str('%.2f'%(avgs[file]))+"}} & "+str(std_devs[file]),end='')
			else:
				print("& "+"\\best{"+str('%.2f'%(avgs[file]))+"} & "+str(std_devs[file]),end='')
		elif opt[inst][file]:
			print("& "+"\\opt{"+str('%.2f'%(avgs[file]))+"} & "+str(std_devs[file]),end='')
		else:
			print("& "+str('%.2f'%(avgs[file]))+" & "+str(std_devs[file]),end='')
		if file == sys.argv[-1]:
			print(" \\\\")
		else:
			print("")
	print("%")