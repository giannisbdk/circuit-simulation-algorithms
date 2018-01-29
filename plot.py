import matplotlib.pyplot as plt
import re
import os
import sys

plots_path  = "./plots"
tr_out_path = "./"
ac_out_path = "./"


def create_plots_dirs(plot_path):
	"""
	Creates the required directories to save the plots.
	"""
	try:
		if not os.path.exists(plot_path):
			os.makedirs(plot_path)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise


def plot_all(plot_path, out_path, prefix, a_type):
	"""
	Plots the output files into a single figure.

	plot_path: contains the path to the figures either plots/TRAN or plots/AC
	out_path:  contains the path to the folder which the output files reside
	prefix:    contains the prefix of the output files either transient or ac
	a_type:    contains the type of the analysis either transient or ac
	"""

	tran = True if a_type == "transient" else False

	# Create an id from the number of figures in plots directory
	fig_id = len([fig for fig in os.listdir(plot_path) if fig.endswith(".png")]) + 1

	filenames = []
	files = os.listdir(out_path)

	for file in files:
		if prefix in file:
			filenames.append(file)

	# Set the appropriate fields for the figure according to the analysis
	if tran:
		fig, ax1 = plt.subplots(nrows=1, figsize=(14, 9))
		fig.suptitle("Transient Analysis")
		ax1.set_xlabel("Time step (seconds)")
		ax1.set_ylabel("Value (Voltage)")
		ax1.grid(True)
		fig_suffix = "_Transient_Analysis"
	else:
		fig, (ax2, ax1) = plt.subplots(nrows=2, figsize=(14, 9))
		fig.suptitle("AC Analysis")
		ax1.set_xlabel("Frequency (Hz)")
		ax1.set_ylabel("Phase (degrees)")
		ax2.set_xlabel("Frequency (Hz)")
		ax2.set_ylabel("Magnitude (db)")
		ax1.grid(True)
		ax2.grid(True)
		fig_suffix = "_AC_Analysis"

	# Create the path that the figure is going to be saved
	fig_path = plot_path + "/" + str(fig_id) + fig_suffix

	# Plot everything and plot to single figure
	for filename in filenames:
		if tran:
			plot_transient_file(filename, ax1)
		else:
			plot_ac_file(filename, ax1, ax2)

	# Get handles and labels from one of the subplots
	handles, labels = ax1.get_legend_handles_labels()
	fig.legend(handles, labels, loc='center right', shadow=True)
	fig.subplots_adjust(right=0.89)

	# Save figure with a lower resolution to fit in sublime when someone opens it through it
	plt.savefig(fig_path, dpi=85)
	plt.show()


def plot_transient_file(filename, ax1):
	"""
	Plot a single TRANSIENT output file into a figure containing
	a single subplot ax1.
	"""

	# Get the Node name from the output file
	v_label = filename.split("_")[2]
	with open(filename, 'rb') as fd:
		lines = [x.strip("\n") for x in fd.readlines()]

	step_list = []
	val_list = []
	pat = r"([0-9]+\.?[0-9]+)\s+(-?[0-9]+\.?[0-9]+)"

	for line in lines[1:]:
		match = re.search(pat, line)
		if match:
			step = match.group(1)
			val  = match.group(2)
			step_list.append(step)
			val_list.append(val)

	ax1.plot(step_list, val_list, label=v_label, linewidth=1.5)


def plot_ac_file(filename, ax1, ax2):
	"""
	Plot a single AC output file into a figure containing
	two subplots ax1, ax2.
	"""

	# Get the Node name and the sweep type from the output file
	tokens  = filename.split("_")
	v_label = tokens[2]
	sweep   = tokens[len(tokens)-1].split(".")[0]

	with open(filename, 'rb') as fd:
		lines = [x.strip("\n") for x in fd.readlines()]

	freq_list  = []
	magn_list  = []
	phase_list = []
	pat = r"([0-9]+\.?[0-9]+)\s+(-?[0-9]+\.?[0-9]+)\s+(-?[0-9]+\.?[0-9]+)"

	for line in lines[1:]:
		match = re.search(pat, line)
		if match:
			freq  = match.group(1)
			magn  = match.group(2)
			phase = match.group(3)
			freq_list.append(freq)
			magn_list.append(magn)
			phase_list.append(phase)

	if sweep == "LIN":
		ax1.plot(freq_list, magn_list, label=v_label, linewidth=1.5)
		ax2.plot(freq_list, phase_list, label=v_label, linewidth=1.5)
	elif sweep == "LOG":
		ax1.semilogx(freq_list, phase_list, label=v_label, linewidth=1.5)
		ax2.semilogx(freq_list, magn_list, label=v_label, linewidth=1.5)
	else:
		print "Error: Sweep type '{}' from output file is wrong.".format(sweep)
		print "Valid options (LIN, LOG)."
		sys.exit(0)


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

	analysis = analysis.upper()

	if analysis == "TRAN" or analysis == "TR":
		plot_path = plots_path + "/" + "TRAN"
		out_path  = tr_out_path
		prefix    = "tr_analysis_"
		a_type    = "transient"
	elif analysis == "AC":
		plot_path = plots_path + "/" + "AC"
		out_path  = ac_out_path
		prefix    = "ac_analysis_"
		a_type    = "ac"
	else:
		print "Error: Wrong input type '{}'".format(analysis)
		print "Try again, valid options are TR,TRAN or AC, either uppercase or lowercase or mixed."
		sys.exit(0)

	return (plot_path, out_path, prefix, a_type)


def main():
	input_message = "\n--- Input can be either lowercase or uppercase or mixed ---\n"\
					"What would you like to plot?\n"\
					"Valid options: [TR, TRAN, AC]\n"
	try:
		analysis = raw_input(input_message)

		plot_path, out_path, prefix, a_type = get_options(analysis)

		if not out_exist(out_path, prefix):
			print "Error: Cannot find any {} output files.".format(a_type)
			sys.exit(0)

		print "\nCreating plots directory...........OK"
		create_plots_dirs(plot_path)

		print "Plotting requested analysis........OK"
		plot_all(plot_path, out_path, prefix, a_type)

		print "Saving figure to plots directory...OK"
	except KeyboardInterrupt:
		sys.exit(0)

if __name__ == '__main__':
	main()
