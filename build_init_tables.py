#!/usr/bin/python
import sys, optparse, os, csv, math

data = {}

if len(sys.argv) < 3:
  print("Please enter csv files in order DETERMINISTIC, RANDOM");
  exit();

for in_file in sys.argv[1:]:

	# print(in_file)

	inFile = open(str(in_file))
	reader = csv.DictReader(inFile)

	for row in reader:
		if row["Instance"]:
			try:
				data[row["Instance"]][in_file] =[]
			except:
				data[row["Instance"]] = {}
				data[row["Instance"]][in_file] =[]
			for i in range(1,21):
				data[row["Instance"]][in_file].append(int(row[str(i)]))
				
	inFile.close()

for inst,file_data in data.items():
	# print(file_data)
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
	print("\\texttt{"+inst+"}")
	avg_vals = []
	for file in sys.argv[1:]:
		avg_vals.append(avgs[file])
		print("& "+str('%.2f'%(avgs[file]))+" & "+str(std_devs[file]),end='')
		if file == sys.argv[-1]:
			print("\n& "+str('%.2f'%(100-100*(avg_vals[0]/avg_vals[1])))+" \\(\\%\\)",end='')
			print(" \\\\")
		else:
			print("")
	print("%")