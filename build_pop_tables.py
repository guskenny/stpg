#!/usr/bin/python
import sys, optparse, os, csv, math

data = {}
sizes = {}

if len(sys.argv) < 2:
  print("Please enter csv file");
  exit();

for in_file in sys.argv[1:]:

	# print(in_file)

	inFile = open(str(in_file))
	reader = csv.DictReader(inFile)

	for row in reader:
		if row["Instance"]:
			if row["Population size"] == "500" or row["Population size"] == "1000":
				continue
			sizes[row["Instance"]] = int(row["Problem size"])
			try:
				data[row["Instance"]][row["Population size"]].append(int(row["Num groups"]))
			except:
				try:
					data[row["Instance"]][row["Population size"]] = [int(row["Num groups"])]
				except:
					data[row["Instance"]] = {}
					data[row["Instance"]][row["Population size"]] = [int(row["Num groups"])]
				
	inFile.close()

for inst,file_data in data.items():
	# print(file_data)
	avgs = {}
	square_diffs = {}
	std_devs = {}

	for pop_size,pop_data in file_data.items():
		avgs[pop_size] = float(sum(pop_data)/len(pop_data))
		square_diffs[pop_size] = [math.pow(x - avgs[pop_size],2) for x in pop_data]
		std_devs[pop_size] = '%.2f'%(math.sqrt(float(sum(square_diffs[pop_size]))/len(square_diffs[pop_size])))
	# std_dev = math.sqrt(float(sum(square_diff))/len(square_diff))
	key_min = min(avgs.keys(), key=(lambda k: avgs[k]))
	min_avg = avgs[key_min]
	print("\\texttt{"+inst+"} & "+str(sizes[inst]))
	avg_vals = []
	count = 1
	for pop_size in ["50","100","1500","10000"] :  
	# for pop_size,pop_data in file_data.items():
		avg_vals.append(avgs[pop_size])
		print("& "+str('%.2f'%(avgs[pop_size]))+" & "+str(std_devs[pop_size]),end='')
		if count < len(file_data):
			print("")
			count += 1
		else:
			print(" \\\\")
	print("%")