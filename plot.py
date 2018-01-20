import matplotlib.pyplot as plt
import re
import os
import sys

plots_path  = "./plots"
tr_out_path = "./"
ac_out_path = "./"


def create_plots_dirs(plot_path):
	try:
		if not os.path.exists(plot_path):
			os.makedirs(plot_path)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise


def plot_all(plot_path, out_path, prefix):
	"""
	plot_path: contains the path to the figures either plots/TRAN or plots/AC
	out_path:  contains the path to the folder which the output files reside
	prefix:    contains the prefix of the output files either transient or ac
	"""

	# Create an id from the number of figures in plots directory
	fig_id = len([fig for fig in os.listdir(plot_path) if fig.endswith(".png")]) + 1

	filenames = []
	files = os.listdir(out_path)

	for file in files:
		if prefix in file:
			filenames.append(file)

	# Plot everything and plot to single figure
	for file in filenames:
		plot_transient_file(file)

	# This is only for transient
	# TODO need to create different for AC
	# TODO create two functions/methods for transient and ac
	plt.title("Transient Analysis")
	plt.xlabel("Time step (seconds)")
	plt.ylabel("Value (Voltage)")
	plt.grid(True)
	plt.legend(loc='best')
	plt.savefig(plot_path + "/" + str(fig_id) + "_Transient_Analysis")
	plt.autoscale(True, True, True)
	plt.show()


def plot_transient_file(filename):
	# Get the Node name from the output file
	v_label = "V" + filename.split("_")[2]
	with open(filename, 'rb') as fd:
		lines = [x.strip("\n") for x in fd.readlines()]

	step_list = []
	val_list = []
	pat = r"([0-9]*\.?[0-9]+)\s+(-?[0-9]*\.?[0-9]+)"

	for line in lines[1:]:
		match = re.search(pat, line)
		if match:
			step = match.group(1)
			val = match.group(2)
			step_list.append(step)
			val_list.append(val)

	plt.plot(step_list, val_list, label=v_label)


def out_exist(out_path, prefix):
	"""
	Checks if there are any transient/ac output files in the
	out_path directory and returns the appropriate value
	"""

	files = os.listdir(out_path)
	found = False
	for file in files:
		if prefix in file:
			found = True
	return found


def get_options(analysis):
	"""
	Returns a tuple with the plot_path, out_path, prefix and type of analysis
	according to the analysis string givem from the cmd.
	"""

	if analysis == "TRAN" or analysis == "TR":
		plot_path = plots_path + "/" + "TRAN"
		out_path  = tr_out_path
		prefix    = "tr_analysis_"
		title     = "transient"
	elif analysis == "AC":
		plot_path = plots_path + "/" + "AC"
		out_path  = ac_out_path
		prefix    = "ac_analysis_"
		title     = "ac"
	else:
		print "Error: Wrong input type '{}'".format(analysis)
		print "Try again, valid options are TR,TRAN or AC, either uppercase or lowercase or mixed."
		sys.exit(0)
	return (plot_path, out_path, prefix, title)


def main():
	input_message = "\n--- Input can be either lowercase or uppercase or mixed ---\n"\
					"What would you like to plot?\n"\
					"Valid options: [TR, TRAN, AC]\n"
	try:
		analysis = str(raw_input(input_message)).upper()

		plot_path, out_path, prefix, title = get_options(analysis)

		if not out_exist(out_path, prefix):
			print "Error: Cannot find any {} output files.".format(title)
			sys.exit(0)

		print "\nCreating plots directory...........OK"
		create_plots_dirs(plot_path)

		print "Plotting {} analysis........OK".format(title)
		plot_all(plot_path, out_path, prefix)

		print "Saving figure to plots directory...OK"
	except KeyboardInterrupt:
		sys.exit(0)

if __name__ == '__main__':
	main()
