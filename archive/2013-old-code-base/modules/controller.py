"""Controller module
Exerts principal control over pipeline execution
Parses configuration file and finalizes init of logging
system, reports ready to main() and then receives kick-off command
"""

import os as os
import sys as sys
import configparser as cp
import logging as log

import modules.osinterface as osif
import modules.dbmanager as dbman
import modules.toolbox as tlbx

class Controller(object):
	"""Controller class
	Runs appropriate part of pipeline
	Needs path to config file and a usecase
	"""
	def __init__(self, cfgfile, usecase):
		self.cfgfile = cfgfile
		self.usecase = usecase
		self.config, self.dbconfig, self.tlbxconfig, logdict = self._parse_config(cfgfile, usecase)
		ham = self._reset_root_logger(logdict)
		self.logger = log.getLogger(__name__)
		self.logger.debug("Initial configuration set")    
	def _parse_config(self, cfgfile, usecase):
		"""Use Python's configparser to read section "usecase"
		of given INI file, returns configuration as dict
		data structure (w/o coercing types)
		"""
		cfgparse = cp.ConfigParser()
		cfgparse.read(cfgfile)
        # use of raw is necessary because of format strings for logging
		cfgdict = dict(cfgparse.items(usecase, raw=True))
		logdict = dict(cfgparse.items("Logging", raw=True))
		tbxdict = dict(cfgparse.items("Toolbox", raw=True))
		dbdict = dict(cfgparse.items("Database", raw=True))
		assert cfgdict and logdict and tbxdict and dbdict, "_parse_config(): got empty configurations"
		return cfgdict, dbdict, tbxdict, logdict
	def _reset_root_logger(self, logdict):
		"""Set up logging system with log file and
		delete StreamHandler to stderr (all log messages go
		to a single file for simplicity)
		"""
		osif.check_or_create_dir(logdict['logdir'])
		# we may now assume path is accessible
		fhdlr = log.FileHandler(os.path.join(logdict['logdir'], logdict['logfile']))
		form = log.Formatter(logdict['format'])
		fhdlr.setFormatter(form)
		rootlogger = log.getLogger('')
		rootlogger.addHandler(fhdlr)
		# if at some point the logging system is changed, this probably needs to be changed, too
		eggs = {'DEBUG':log.DEBUG, 'ERROR':log.ERROR, 'INFO':log.INFO}
		rootlogger.setLevel(eggs[logdict['level']])
		allhandlers = rootlogger.handlers
		# there should be 2 handlers now, 1 from basicConfig and 1 from here
		assert len(allhandlers) == 2, "_reset_root_logger(): too few or too many handlers"
		if isinstance(allhandlers[0], log.StreamHandler):
			rootlogger.removeHandler(allhandlers[0])
		else:
			assert isinstance(allhandlers[1], log.StreamHandler), "_reset_root_logger(): no StreamHandler found"
			rootlogger.removeHandler(allhandlers[1])
		return None
	def run(self):
		"""Kick-off method, import appropriate module,
		start execution for usecase
		"""
		self.logger.debug("Kick-off for pipeline initiated")
		self.logger.info("Running on machine {}".format(os.uname()[1]))
		# before we can do anything, check if DB dir exists
		assert os.path.exists(self.dbconfig['dbdir']), "run(): DB directory {} appears not to be accessible".format(self.dbconfig['dbdir'])
		dbmanager = dbman.DBManager(self.dbconfig)
		toolbox = tlbx.Toolbox(self.tlbxconfig)
		# guess there is a more elegant way for the following
		if self.usecase == "Test":
			self.logger.debug("Usecase TEST: do nothing")
		elif self.usecase == "CreateBlank":
			mod = __import__("mapping", globals(), locals(), [], -1)
			obj = mod.Mapping(self.config)
			status = obj.prepare_mapping(dbmanager, toolbox)
		elif self.usecase == "ProcessBlank":
			mod = __import__("mapping", globals(), locals(), [], -1)
			obj = mod.Mapping(self.config)
			status = obj.process_mapping(dbmanager, toolbox)
		elif self.usecase == "MapFiles":
			mod = __import__("mapping", globals(), locals(), [], -1)
			obj = mod.Mapping(self.config)
			status = obj.map_signal_tracks(dbmanager, toolbox)
		elif self.usecase == "PrepareDatasets":
			mod = __import__("prediction", globals(), locals(), [], -1)
			obj = mod.Prediction(self.config)
			status = obj.prepare_datasets(dbmanager, toolbox)
		elif self.usecase == "PrepareExpression":
			mod = __import__("prediction", globals(), locals(), [], -1)
			obj = mod.Prediction(self.config)
			status = obj.prepare_expression_data(dbmanager, toolbox)
		elif self.usecase == "ComputeFeatures":
			mod = __import__("prediction", globals(), locals(), [], -1)
			obj = mod.Prediction(self.config)
			status = obj.compute_features(dbmanager, toolbox)
		elif self.usecase == "Toolbox":
			status = toolbox.stand_alone()
		else:
			self.logger.error("Usecase {} unknown".format(self.usecase))
			return False
		return status
