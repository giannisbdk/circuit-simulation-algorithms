import matplotlib.pyplot as plt
import re
import os

plots_path = "./plots"
tr_out_path = "./"

def create_plots_dir():
	try:
		if not os.path.exists(plots_path):
			os.makedirs(plots_path)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

def plot_all(dir_path):
	filenames = []
	files = os.listdir(dir_path)

	for file in files:
		if "tr_analysis_" in file:
			filenames.append(file)
	
	# Plot everything and plot to single figure
	for file in filenames:
		plot_file(file)

	plt.title("Transient Analysis")
	plt.xlabel("Time step (seconds)")
	plt.ylabel("Value (Voltage)")
	plt.grid(True)
	plt.legend(loc='best')
	plt.savefig(plots_path + "/" + "Transient_Analysis")
	plt.show()


def plot_file(filename):
	# Get the Node name from the output file
	v_label = "V" + filename.split("_")[2]
	with open(filename ,'rb') as fd:
		lines = [x.strip("\n") for x in fd.readlines()]

	step_list = []
	val_list  = []
	pat = r"([0-9]*\.?[0-9]+)\s+([0-9]*\.?[0-9]+)"

	for line in lines[1:]:
		match = re.search(pat, line)
		if match:
			step = match.group(1)
			val  = match.group(2)
			step_list.append(step)
			val_list.append(val)

	plt.plot(step_list, val_list, label=v_label)


def tr_exist(tr_out_path):
	""" Checks if there are any transient output files in the current
	directory and returns the appropriate value"""
	files = os.listdir(tr_out_path)
	found = False
	for file in files:
		if "tr_analysis_" in file:
			found = True
	return found


def main():
	if not tr_exist(tr_out_path):
		print "Error: Cannot find any transient output files."
		return
	print "Creating plots directory...........OK"
	create_plots_dir()
	print "Plotting transient analysis........OK"
	plot_all(tr_out_path)
	print "Saving figure to plots directory...OK"


if __name__ == '__main__':
	main()