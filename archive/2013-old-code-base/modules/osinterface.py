"""Module to interface with OS
Just collection of functions since an object would not hold any data
- main functions: directory traversal, checking/creating paths
- determining system load, estimate free resources
"""

import os as os
import shutil as shu
import time as ti
import collections as col
import fnmatch as fn
import multiprocessing as mp
import math as ma

def check_or_create_dir(dirpath):
	"""Check existence of dirpath or
		try to create if necessary
	"""
	if not os.path.exists(dirpath):
		try:
			os.makedirs(dirpath)
		except OSError:
			raise Exception("OSInterface: Could not create directory %s" %dirpath)
	return None
def collect_files(path, pattern):
	"""Create a list of the files compliant to 'pattern'
		in the subtree starting in 'path'
		Returned list might be empty
		Returned list contains filenames with full path
	"""
	assert os.path.exists(path), "OSInterface.collect_files(): path {} is not accessible".format(path)
	all_files = []
	for root,subdir,filenames in os.walk(path,topdown=True,onerror=None,followlinks=False):
		for f in filenames:
			all_files.append(os.path.join(root,f))
	red_files = fn.filter(all_files,pattern)
	return red_files
def copy_file_to_dir(absfilename,directory):
	"""Copy a file (absolute location) to a new location
		Returns the absolute path of the file in that new location
	"""
	assert os.path.exists(absfilename), "OSInterace.copy_file_to_dir(): could not find file {}".format(absfilename)
	if not os.path.exists(directory):
		try:
			ret = check_or_create_dir(directory)
		except:
			raise Exception("OSInterface: Cannot create directory %s\n" %directory)
	try:
		ret = shu.copy(absfilename,directory)
		newlocation = os.path.join(directory,os.path.basename(absfilename))
		assert os.path.exists(newlocation), "OSInterface.copy_file_to_dir(): copying file to {} failed".format(newlocation)
		return newlocation
	except Exception as e:
		raise Exception("OSInterface.copy_file_to_dir: Caught error {} from shutil.copy".format(e))
def move_file_to_newfile(oldabsname, newabsname):
	"""Move a file to another location (it is renamed at least)
	"""
	assert os.path.exists(oldabsname), "OSIinterface.move_file_to_newfile(): could not find file {}".format(oldabsname)
	if not os.path.exists(os.path.dirname(newabsname)):
		try:
			ret = check_or_create_dir(os.path.dirname(newabsname))
		except:
			raise Exception("OSInterface: cannot create directory %s" %os.path.dirname(newabsname))
	try:
		ret = shu.move(oldabsname, newabsname)
		assert os.path.exists(newabsname), "OSInterface.move_file_to_newfile(): shutil.move to {} failed".format(newabsname)
		return newabsname
	except Exception as e:
		raise Exception("OSInterface.move_file_to_newfile(): Caught error {} from shutil.move".format(e))
def collect_all_files(path):
	"""Same functionality as above, but ignore file extension
	"""
	assert os.path.exists(path), "OSInterface.collect_all_files(): path {} is not accessible".format(path)
	all_files = []
	for root,subdir,filenames in os.walk(path,topdown=True,onerror=None,followlinks=False):
		for f in filenames:
			all_files.append(os.path.join(root,f))
	return all_files
def estimate_free_cpus():
	"""Estimate number of free CPUs for parallelizing tasks
	"""
	try:
		cpus = mp.cpu_count()
		# tuple[1] is avg over last 5 minutes
		avgload = int(ma.ceil((os.getloadavg()[1])))
		if cpus > 4 > avgload - cpus:
			# assuming server is/was under high load
			# so wait a little bit and see what happens
			ti.sleep(60)
			avgload = int(ma.ceil((os.getloadavg()[1])))
		freecpus = 1 if 1 >= cpus - avgload else cpus - avgload - 1
		return freecpus
	except OSError:
		# this means load average could not be obtained
		# return 1 as number of free CPUs
		return 1
	except Exception as e:
		raise Exception("Estimating free CPUs failed with exception {}".format(e))
