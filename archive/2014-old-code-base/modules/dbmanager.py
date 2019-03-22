"""DBmanager module
Handles requests to the database files
Reads a file on demand (no bottleneck, since
files are very small), converts it into a dict
of dicts and returns that dictionary
"""

import collections as col
import logging as log
import os as os

import modules.osinterface as osif

class DBManager(object):
    def __init__(self, configuration):
        self.config = configuration
        self.logger = log.getLogger(__name__)
        self.dbfiles = osif.collect_files(self.config['dbdir'], "*.tsv")
        assert self.dbfiles, "DBManager: No database files found in {}".format(dbdir)
        self.available_DBs = self._prepare_DB_access()
        assert self.available_DBs, "DBManager: Found database files but could not initiate available DBs"
    def _prepare_DB_access(self):
        resdict = {}
        for k,v in self.config.items():
            for f in self.dbfiles:
                if v in f:
                    resdict[k] = f
        return resdict          
    def create_DBdict(self, whichDB):
        """Generate a list of dicts for DB file
        Raises an Exception if the file does not exist
        """
        try:
            dbfile = self.available_DBs[whichDB]
            dbdict = {}
            with open(dbfile, "r") as dbf:
                header = dbf.readline().strip().split("\t")
                for line in dbf:
                    if not line: continue # in case of newlines at end of file
                    dbentry = line.strip().split("\t")
                    d = dict(zip(header, dbentry))
                    dbdict[d['key']] = d
            if not dbdict:
                self.logger.warning("DBdict created from {} is empty".format(dbfile))
            return dbdict
        except KeyError:
            raise Exception("Requested DB {} does not exist".format(whichDB))
    def update_DBfile_make_id(self, whichDB, dbdict, newentries):
        """Open DB file in 'a' mode and write new entries to DB file
        Assume that these new entries do not have a key yet - determine
        largest key already in DB file
        Getting the original dbdict as input simplifies determining the currently
        largest key in the database

        dbdict parameter is deprecated

        """
        try:
            dbfile = self.available_DBs[whichDB]
            dbdict = self.create_DBdict(whichDB)
            key = max(map(int, dbdict.keys()))
            order = []
            with open(dbfile, "r") as dbf:
                header = dbf.readline().strip()
                order = header.split("\t")
            with open(dbfile, "a") as dbf:
                for entry in sorted(newentries, key=lambda d: int(d['key'])):
                    key += 1
                    entry['key'] = str(key)
                    assert set(entry.keys()) == set(order), "New entry is missing key - please provide reasonable default\n{}\n{}".format(entry.keys(), order)
                    newline = "\t".join( [ entry[i] for i in order ] )
                    dbf.write(newline)
                    dbf.write("\n")
            return
        except KeyError:
            raise Exception("Add entry: requested DB {} no longer exists".format(whichDB))
    def update_DBfile_no_id(self, whichDB, newentries):
        """Same as above but entries are assumed to have already the correct id
        """
        try:
            dbfile = self.available_DBs[whichDB]
            order = []
            with open(dbfile, "r") as dbf:
                header = dbf.readline().strip()
                order = header.strip().split("\t")
            with open(dbfile, "a") as dbf:
                for entry in sorted(newentries, key=lambda d: int(d["key"])): # keep DBfile sorted... currently only sugar
                    assert set(entry.keys()) == set(order), "New entry is missing key - please provide reasonable default\n{}\n{}".format(entry.keys(), order)
                    newline = "\t".join( [ entry[i] for i in order ] )
                    dbf.write(newline)
                    dbf.write("\n")
            return
        except KeyError:
            raise Exception("Add entry: requested DB {} no longer exists".format(whichDB))
    def update_DBentry(self, entry, whichDB):
        """Read DB file, replace changed entry and (over-)write DB file
        """
        try:
            db = self.create_DBdict(whichDB)
            db[entry['key']] = entry
            dbfile = self.available_DBs[whichDB]
            header = None
            backup = None
            with open(dbfile, "r") as dbf:
                backup = dbf.read()
                dbf.seek(0)
                header = dbf.readline()
            with open(dbfile, "w") as dbf:
                dbf.write(header)
                order = header.strip().split("\t")
                for i in sorted(db.keys(), key=int):
                    dbentry = db[i]
                    outline = "\t".join( [dbentry[k] for k in order] )
                    dbf.write(outline + "\n")
            return          
        except KeyError:
            raise Exception("Update entry: requested DB {} no longer exists".format(whichDB))
        except Exception: # this is important to make sure that the DB is dumped to disk if anything happens
            backout = open("/tmp/Py_DB_backup.txt", "w")
            backout.write(backup)
            backout.close()
            raise Exception("Updating DB entry failed - DB dump written to /tmp/Py_DB_backup.txt")
    def get_next_id(self, whichDB):
        """Return the next id for a new entry for a DB file
        """
        dbdict = self.create_DBdict(whichDB)
        if not dbdict:
            return 1
        keys = map(int, dbdict.keys())
        return (max(keys) + 1)
