#!./vpy-3.2.5
# -*- coding: utf-8 -*-

import sys as sys
import os as os
import random as rand
import configparser as cp
import time as ti
import logging as log
import traceback as trb
import multiprocessing as mp

import modules.controller as ctrl

def main():
	'''
	main() only used to check config file, set up logging
	principal control over pipeline is delegated to controller
	'''
	try:
		# everything lands on stderr until log file is set up
		log.basicConfig(stream=sys.stderr, level=log.DEBUG)
		cfgfile = None
		usecase = None
		try:
			cfgfile = sys.argv[1]
			usecase = sys.argv[2]
			assert os.path.exists(cfgfile), "No config file found at location {}\n".format(cfgfile)
			controller = ctrl.Controller(cfgfile, usecase)
		except IndexError:
			log.error("Got {} as config and {} as usecase - could not proceed\nInvoke: ./run_pipeline CFGFILE USECASE".format(cfgfile, usecase))
			return 1
		except AssertionError as ae:
			log.error(str(ae))
			return 1
		except Exception as e:
			log.error("Caught exception before log file was ready or while setup: " + str(e) + "\n")
			return 1
		else:
			logger = log.getLogger(__name__)
			rand.seed(os.urandom(1024))
			# controller gets control over pipeline execution
			logger.info("Controller init successfull, starting computation")
			status = controller.run()
			if not status:
				logger.warning("Controller returned with bad status - check previous log messages")
			log.shutdown()
			return 0
	except Exception as e:
		# make sure logging system shuts down properly
		logger = log.getLogger(__name__)
		logger.error(str(e))
		logger.error(trb.format_exc(limit=3))
		log.shutdown()
		return 1
	return 0

if __name__ == "__main__":
	sys.exit(main())
