#!/user/bin/env python

import os.path
import os
import sys
from pans_labyrinth import loggingFunctions

"""
Files and directory functions go here.
"""

def get_files(file_or_directory):
	"""
	Creates a list of files from either the given file, or all files within the
	directory specified (where each file name is its absolute path).
	:param file_or_directory (str): file or directory name given on commandline
	:return: a list of all the files found.
	"""

	files_list = []
	if file_or_directory:
		if os.path.isdir(file_or_directory):
			# Create a list containing the file names
			for root, dirs, files in os.walk(file_or_directory):
				for filename in files:
					files_list.append(os.path.join(root, filename))
		# check if input is concatenated file locations
		elif ',' in file_or_directory:
			for filename in file_or_directory.split(','):
				files_list.append(os.path.abspath(filename))
		else:
			files_list.append(os.path.abspath(file_or_directory))

	if not files_list:
		exit("No files were found")

	sorted_files = sorted(files_list)
	return sorted_files


def walkdir(folder):
	'''
	Walk through each files in a directory and yeild all the paths in the dir
	:param folder: the folder containing all the files to walk through
	'''
	LOG = loggingFunctions.create_logger()

	if len(os.listdir(folder)) == 0:
		LOG.critical("Empty directory")
		sys.exit()
	for dirpath, dirs, files in os.walk(folder):
		for filename in files:
			yield os.path.abspath(os.path.join(dirpath, filename))
