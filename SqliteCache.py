## taken from https://gist.github.com/leonjza/b5043b5602e9222e6760
## original author leonjza


import os
import errno
import sqlite3
import sys
from time import time
try:
    from cPickle import loads, dumps  # Python 2
    try:
        buffer_func = buffer
    except NameError:
        buffer_func = lambda x: x
except ImportError:
    from pickle import loads, dumps   # Python 3
    buffer_func = lambda x: x

import logging

logger = logging.getLogger(__name__)


class SqliteCache:
    """
        SqliteCache
        Ripped heavily from: http://flask.pocoo.org/snippets/87/
        This implementation is a simple Sqlite based cache that
        supports cache timers too. Not specifying a timeout will
        mean that the value will exist forever.
    """

    # prepared queries for cache operations
    _create_sql = (
        'CREATE TABLE IF NOT EXISTS entries '
        '( key TEXT PRIMARY KEY, val BLOB, exp FLOAT )'
    )
    _create_index = 'CREATE INDEX IF NOT EXISTS keyname_index ON entries (key)'
    _get_sql = 'SELECT val, exp FROM entries WHERE key = ?'
    _del_sql = 'DELETE FROM entries WHERE key = ?'
    _set_sql = 'REPLACE INTO entries (key, val, exp) VALUES (?, ?, ?)'
    _add_sql = 'INSERT INTO entries (key, val, exp) VALUES (?, ?, ?)'
    _clear_sql = 'DELETE FROM entries'

    # other properties
    connection = None

    def __init__(self, path):

        """ Inits a new SqliteCache instance """

        self.path = os.path.abspath(path)
        logger.debug('Instantiated with cache_db path as {path}'.format(path=self.path))

        # prepare the directory for the cache sqlite db
        try:

            os.mkdir(self.path)
            logger.debug('Successfully created the storage path for {path}'.format(
                path=self.path))

        except OSError as e:

            if e.errno != errno.EEXIST or not os.path.isdir(self.path):
                raise

    def _get_conn(self):

        """ Returns a Sqlite connection """

        if self.connection:
            return self.connection

        # specify where we want the cache db to live
        cache_db_path = os.path.join(self.path, "cache.sqlite")

        # setup the connection
        conn = sqlite3.Connection(cache_db_path, timeout=60)
        logger.debug('Connected to {path}'.format(path=cache_db_path))

        # ensure that the table schema is available. The
        # 'IF NOT EXISTS' in the create_sql should be
        # pretty self explanitory
        with conn:
            conn.execute(self._create_sql)
            conn.execute(self._create_index)
            logger.debug('Ran the create table && index SQL.')

        # set the connection property
        self.connection = conn

        # return the connection
        return self.connection

    def get(self, key):

        """ Retreive a value from the Cache """

        return_value = None

        # get a connection to run the lookup query with
        with self._get_conn() as conn:

            # loop the response rows looking for a result
            # that is not expired
            for row in conn.execute(self._get_sql, (key,)):

                expire = row[1]
                if expire == 0 or expire > time():
                    return_value = loads(str(row[0]))
                    # TODO: Delete the value that is expired?

                break

        return return_value

    def delete(self, key):

        """ Delete a cache entry """

        with self._get_conn() as conn:
            conn.execute(self._del_sql, (key,))

    def update(self, key, value, timeout=None):

        """ Sets a k,v pair with a optional timeout """

        # if no timeout is specified, then we will
        # leave it as a non-expiring value. Other-
        # wise, we add the timeout in seconds to
        # the current time
        expire = 0 if not timeout else time() + timeout

        # serialize the value with protocol 2
        # ref: https://docs.python.org/2/library/pickle.html#data-stream-format
        data = buffer_func(dumps(value, 2))

        # write the updated value to the db
        with self._get_conn() as conn:
            conn.execute(self._set_sql, (key, data, expire))

    def set(self, key, value, timeout=None):

        """ Adds a k,v pair with a optional timeout """

        # if no timeout is specified, then we will
        # leave it as a non-expiring value. Other-
        # wise, we add the timeout in seconds to
        # the current time
        expire = 0 if not timeout else time() + timeout

        # serialize the value with protocol 2
        # ref: https://docs.python.org/2/library/pickle.html#data-stream-format
        data = buffer_func(dumps(value, 2))

        # adding a new entry that may cause a duplicate key
        # error if they key already exists. In this case
        # we will fall back to the update method.
        with self._get_conn() as conn:

            try:

                conn.execute(self._add_sql, (key, data, expire))

            except sqlite3.IntegrityError:

                # call the update method as fallback
                logger.warning(
                    'Attempting to set an existing key {k}. Falling back to update method.'.format(
                        k=key))
                self.update(key, value, timeout)
                pass

    def clear(self):

        """ Clear a cache """

        with self._get_conn() as conn:
            conn.execute(self._clear_sql, )

    def __del__(self):

        """ Cleans up the object by destroying the sqlite connection """

        if self.connection:
            self.connection.close()


# allow this module to be used to clear the cache
if __name__ == '__main__':

    # check args
    if len(sys.argv) != 3 or sys.argv[1] != 'clear':
        print('[!] Error: You have to specify the clear with `python %s clear <path to cache.sqlite>`' % sys.argv[0])
        sys.exit(1)

    if not os.path.isdir(sys.argv[2]):
        print('[!] Error: %s does not seem to be a path!' % sys.argv[2])
        sys.exit(1)

    # setup the cache instance and clear it.
    c = SqliteCache(sys.argv[2])
    c.clear()
    print(' * Cache cleared')