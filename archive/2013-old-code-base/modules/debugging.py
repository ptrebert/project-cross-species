"""
Debugging module to monitor memory usage
Each module level function return a simple
string to be passed to self.logger.debug
"""

import resource as rsrc
import gc as gc

def print_mem_io_load(unit='GB'):
    """Return a combined string for memory
    and I/O load
    """
    gc.collect(2)
    units = {'GB':1024*1024, 'MB':1024, 'KB':1, 'TB':1024*1024*1024}
    normalize = units.get(unit,1)
    normunit = unit if unit in units else 'KB'
    resobj = rsrc.getrusage(rsrc.RUSAGE_SELF)
    mem = resobj.ru_maxrss / normalize
    smem = resobj.ru_ixrss / normalize
    unsmem = resobj.ru_idrss / normalize
    ioin = resobj.ru_inblock
    ioout = resobj.ru_oublock
    logstr = 'Process RSS {normunit}: {mem}\n'\
            'Process shared memory {normunit}: {smem}\n'\
            'Process unshared memory {normunit}: {unsmem}\n'\
            'I/O block input ops: {ioin}\n'\
            'I/O block output ops: {ioout}'.format(**locals())
    return logstr
