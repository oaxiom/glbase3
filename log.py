"""

log.py

sets up a logging environment for the rest of glbase

"""

import logging

class log:
    LOG_FILENAME = '/tmp/glbase_log.txt'
    handle = logging.getLogger('glbase_log')
    handle.setLevel(logging.DEBUG)


