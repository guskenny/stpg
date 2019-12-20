#!/usr/bin/python
import sys, optparse, os, csv, math

data = {}
opt = {}

in_file = str(sys.argv[1])

inFile = open(in_file)
reader = csv.DictReader(inFile)

for row in reader:
	if row["Instance"]:
		try:
			data[row["Instance"]].append(int(row["Best obj"]))
		except: 
			data[row["Instance"]] = [int(row["Best obj"])]

		if row["Found optimal"]:
			opt[row["Instance"]] = True
inFile.close()

for inst,vals in data.items():
	avg = float(sum(vals)/len(vals))
	square_diff = [math.pow(x - avg,2) for x in vals]
	std_dev = '%.2f'%(math.sqrt(float(sum(square_diff))/len(square_diff)))
	# std_dev = math.sqrt(float(sum(square_diff))/len(square_diff))
	if inst in opt:
		if std_dev == "0.00":
			print(inst + ":\t & "+"\\best{\\opt{"+str(avg)+"}} & "+str(std_dev))
		else:
			print(inst + ":\t & "+"\\opt{"+str(avg)+"} & "+str(std_dev))
	else:
		print(inst + ":\t & "+str(avg)+" & "+str(std_dev))