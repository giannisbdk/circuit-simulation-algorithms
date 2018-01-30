import matplotlib.pyplot as plt
import re
import os
import sys

plots_path  = "./plots"
dc_out_path = "./"
tr_out_path = "./"
ac_out_path = "./"


def plot_dc_or_transient_file(filename, ax_tr):
	"""
	Plot a single dc or transient output file into a figure
	containing a single subplot ax.
	"""

	# Get the Node name from the output file
	if "dc" in filename:
		# v_label will contain .txt so we remove it [:-4]
		v_label = filename.split("_")[3][:-4]
	else:
		v_label = filename.split("_")[2]

	with open(filename, 'rb') as fd:
		lines = [x.strip("\n") for x in fd.readlines()]

	step_list = []
	val_list  = []
	pat = r"([0-9]+\.?[0-9]+)\s+(-?[0-9]+\.?[0-9]+)"

	for line in lines[1:]:
		match = re.search(pat, line)
		if match:
			step = match.group(1)
			val  = match.group(2)
			step_list.append(step)
			val_list.append(val)

	ax_tr.plot(step_list, val_list, label=v_label, linewidth=1.5)


def plot_ac_file(filename, ax_ac_1, ax_ac_2, sweep):
	"""
	Plot a single AC output file into a figure containing
	two subplots ax1, ax2.
	"""

	# Get the Node name and the sweep type from the output file
	v_label  = filename.split("_")[2]

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
		ax_ac_1.plot(freq_list, magn_list, label=v_label, linewidth=1.5)
		ax_ac_2.plot(freq_list, phase_list, label=v_label, linewidth=1.5)
	elif sweep == "LOG":
		ax_ac_1.semilogx(freq_list, phase_list, label=v_label, linewidth=1.5)
		ax_ac_2.semilogx(freq_list, magn_list, label=v_label, linewidth=1.5)


def plot_analyses(analyses, paths):
	"""
	Plots the output files from the analyses into a single figure and
	saves the figure to the corresponding folder, indicated by the paths
	dictionary.
	"""

	dc_files = analyses.get("DC")
	tr_files = analyses.get("TRAN")
	ac_files = analyses.get("AC")

	# AC Analysis plot/figure
	if len(ac_files):
		sweep = get_sweep(ac_files)
		fig_ac, (ax_ac_2, ax_ac_1) = plt.subplots(nrows=2, figsize=(14, 9))
		fig_ac.suptitle("AC Analysis")
		ax_ac_1.set_xlabel("Frequency (Hz)")
		ax_ac_1.set_ylabel("Phase (degrees)")
		ax_ac_2.set_xlabel("Frequency (Hz)")
		if sweep == "LIN":
			ax_ac_2.set_ylabel("Magnitude (volts)")
		else:
			ax_ac_2.set_ylabel("Magnitude (dB)")
		ax_ac_1.grid(True)
		ax_ac_2.grid(True)
		fig_ac_suffix = "_AC_Analysis"

		# Plot everything into a figure
		for file in ac_files:
			plot_ac_file(file, ax_ac_1, ax_ac_2, sweep)

		# Create and adjust the legend
		ac_handles, ac_labels = ax_ac_1.get_legend_handles_labels()
		fig_ac.legend(ac_handles, ac_labels, loc='center right', shadow=True)
		fig_ac.subplots_adjust(right=0.89)

		# Save figure to path
		ac_plot_path = paths.get("AC")
		fig_id = len([fig for fig in os.listdir(ac_plot_path) if fig.endswith(".png")]) + 1
		fig_ac.savefig(ac_plot_path + str(fig_id) + fig_ac_suffix)

	# TRAN Analysis plot/figure
	if len(tr_files):
		fig_tr, ax_tr = plt.subplots(nrows=1, figsize=(14, 9))
		fig_tr.suptitle("Transient Analysis")
		ax_tr.set_xlabel("Time step (seconds)")
		ax_tr.set_ylabel("Value (volts)")
		ax_tr.grid(True)
		fig_tr_suffix = "_Transient_Analysis"

		# Plot everything into a figure
		for file in tr_files:
			plot_dc_or_transient_file(file, ax_tr)

		# Create and adjust the legend
		tr_handles, tr_labels = ax_tr.get_legend_handles_labels()
		fig_tr.legend(tr_handles, tr_labels, loc='center right', shadow=True)
		fig_tr.subplots_adjust(right=0.89)

		# Save figure to path
		tr_plot_path = paths.get("TRAN")
		fig_id = len([fig for fig in os.listdir(tr_plot_path) if fig.endswith(".png")]) + 1
		fig_tr.savefig(tr_plot_path + str(fig_id) + fig_tr_suffix)

	# DC Analysis plot/figure
	if len(dc_files):
		fig_dc, ax_dc = plt.subplots(nrows=1, figsize=(14, 9))
		fig_dc.suptitle("DC Analysis")
		ax_dc.set_xlabel("Voltage step (volts)")
		ax_dc.set_ylabel("Value (volts)")
		ax_dc.grid(True)
		fig_dc_suffix = "_DC_Analysis"

		# Plot everything into a figure
		for file in dc_files:
			plot_dc_or_transient_file(file, ax_dc)

		# Create and adjust the legend
		dc_handles, dc_labels = ax_dc.get_legend_handles_labels()
		fig_dc.legend(dc_handles, dc_labels, loc='center right', shadow=True)
		fig_dc.subplots_adjust(right=0.89)

		# Save figure to path
		dc_plot_path = paths.get("DC")
		fig_id = len([fig for fig in os.listdir(dc_plot_path) if fig.endswith(".png")]) + 1
		fig_dc.savefig(dc_plot_path + str(fig_id) + fig_dc_suffix)

	# Show all the created figures
	plt.show()


def get_analyses():
	"""
	Checks the output path for every analysis DC, TRAN or AC for output files
	and returns a dictionary that includes these files for every analysis.
	"""

	analyses = {}
	dc_files = []
	tr_files = []
	ac_files = []

	dc_prefix = "dc_analysis_"
	tr_prefix = "tr_analysis_"
	ac_prefix = "ac_analysis_"

	dc_path_files = os.listdir(dc_out_path)
	tr_path_files = os.listdir(tr_out_path)
	ac_path_files = os.listdir(ac_out_path)

	for file in dc_path_files:
		if dc_prefix in file:
			dc_files.append(file)

	for file in tr_path_files:
		if tr_prefix in file:
			tr_files.append(file)

	for file in ac_path_files:
		if ac_prefix in file:
			ac_files.append(file)

	analyses["DC"]   = dc_files
	analyses["TRAN"] = tr_files
	analyses["AC"]   = ac_files

	return analyses

def check_analyses(analyses):
	"""
	Given the analyses dictionary checks its values which are lists
	for output files. In case all values/lists are empty i.e. False
	it terminates the program.
	"""

	if not any(analyses.values()):
		print "Error: Can't find any output files either for DC, TRAN or AC analysis."
		sys.exit(0)


def get_sweep(ac_files):
	"""
	Given the ac_files checks and returns the type of the sweep, either LIN or LOG.
	"""

	sweep = None
	for file in ac_files:
		if "LIN" in file:
			sweep = "LIN"
			break
		elif "LOG" in file:
			sweep = "LOG"
			break

	if sweep is None:
		print "Error: Sweep type from output file is wrong."
		print "Valid options (LIN, LOG)."
		sys.exit(0)

	return sweep


def get_paths_create_folders(analyses):
	"""
	Creates the required folders to store the figures after the plot is done.
	And returns the path dictionary that contains the paths to that folders.
	"""

	paths   = {}
	lists   = []
	dc_list = analyses.get("DC")
	tr_list = analyses.get("TRAN")
	ac_list = analyses.get("AC")

	if len(dc_list):
		paths["DC"] = plots_path + "/DC/"

	if len(tr_list):
		paths["TRAN"] = plots_path + "/TRAN/"

	if len(ac_list):
		paths["AC"] = plots_path + "/AC/"

	flag = False
	for path in paths.values():
		if not os.path.exists(path):
			flag = True
			os.makedirs(path)

	if flag:
		print "\nCreating plot directories..........OK",

	return paths


def main():

	# Get the analyses dictionary
	analyses = get_analyses()

	# Check for errors
	check_analyses(analyses)

	paths = get_paths_create_folders(analyses)

	print "\nPlotting requested analyses........OK",

	plot_analyses(analyses, paths)

	print "\nSaving figures to directories......OK"

if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(0)
