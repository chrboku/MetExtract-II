"""
ParquetCache - A Parquet-based caching implementation
Converted from SqliteCache to use Polars DataFrames and Parquet files
"""

import os
import errno
import polars as pl
import sys
from time import time

try:
    from cPickle import loads, dumps  # Python 2
except ImportError:
    from pickle import loads, dumps  # Python 3

import logging
import base64

logger = logging.getLogger(__name__)


class ParquetCache:
    """
    ParquetCache
    A simple Parquet-based cache that supports cache timers.
    Not specifying a timeout will mean that the value will exist forever.
    Converted from SqliteCache to use Polars DataFrames.
    """

    # cache DataFrame
    cache_df = None
    cache_path = None

    def __init__(self, path):
        """Inits a new ParquetCache instance"""

        self.path = os.path.abspath(path)
        logger.debug("Instantiated with cache_db path as {path}".format(path=self.path))

        # prepare the directory for the cache parquet file
        try:
            os.mkdir(self.path)
            logger.debug("Successfully created the storage path for {path}".format(path=self.path))

        except OSError as e:
            if e.errno != errno.EEXIST or not os.path.isdir(self.path):
                raise

        # specify where we want the cache file to live
        self.cache_path = os.path.join(self.path, "cache.parquet")

    def _load_cache(self):
        """Load the cache DataFrame from disk"""

        if self.cache_df is not None:
            return self.cache_df

        if os.path.exists(self.cache_path):
            try:
                self.cache_df = pl.read_parquet(self.cache_path)
                logger.debug("Loaded cache from {path}".format(path=self.cache_path))
            except Exception as e:
                logger.warning(f"Could not load cache file: {e}. Creating new cache.")
                self.cache_df = None

        if self.cache_df is None:
            # Create empty DataFrame with correct schema
            self.cache_df = pl.DataFrame({"key": [], "val": [], "exp": []}, schema={"key": pl.Utf8, "val": pl.Utf8, "exp": pl.Float64})
            logger.debug("Created new cache DataFrame")

        return self.cache_df

    def _save_cache(self):
        """Save the cache DataFrame to disk"""

        if self.cache_df is not None:
            self.cache_df.write_parquet(self.cache_path)
            logger.debug("Saved cache to {path}".format(path=self.cache_path))

    def get(self, key):
        """Retrieve a value from the Cache"""

        return_value = None

        # Load cache
        cache_df = self._load_cache()

        # Filter for the key
        filtered = cache_df.filter(pl.col("key") == key)

        # Check if key exists and is not expired
        for row_dict in filtered.to_dicts():
            expire = row_dict["exp"]
            if expire == 0 or expire > time():
                # Decode the base64-encoded pickled value
                return_value = loads(base64.b64decode(row_dict["val"].encode("utf-8")))
            else:
                # Delete expired entry
                self.delete(key)
            break

        return return_value

    def delete(self, key):
        """Delete a cache entry"""

        # Load cache
        cache_df = self._load_cache()

        # Remove the key
        self.cache_df = cache_df.filter(pl.col("key") != key)

        # Save back
        self._save_cache()

    def update(self, key, value, timeout=None):
        """Sets a k,v pair with an optional timeout"""

        # if no timeout is specified, then we will
        # leave it as a non-expiring value. Other-
        # wise, we add the timeout in seconds to
        # the current time
        expire = 0 if not timeout else time() + timeout

        # serialize the value and encode as base64 string
        data = base64.b64encode(dumps(value, 2)).decode("utf-8")

        # Load cache
        cache_df = self._load_cache()

        # Remove old entry if it exists
        cache_df = cache_df.filter(pl.col("key") != key)

        # Add new entry
        new_row = pl.DataFrame({"key": [key], "val": [data], "exp": [expire]})

        self.cache_df = pl.concat([cache_df, new_row])

        # Save back
        self._save_cache()

    def set(self, key, value, timeout=None):
        """Adds a k,v pair with an optional timeout"""

        # if no timeout is specified, then we will
        # leave it as a non-expiring value. Other-
        # wise, we add the timeout in seconds to
        # the current time
        expire = 0 if not timeout else time() + timeout

        # serialize the value and encode as base64 string
        data = base64.b64encode(dumps(value, 2)).decode("utf-8")

        # Load cache
        cache_df = self._load_cache()

        # Check if key already exists
        existing = cache_df.filter(pl.col("key") == key)

        if len(existing) > 0:
            # Key exists, use update method
            logger.warning("Attempting to set an existing key {k}. Falling back to update method.".format(k=key))
            self.update(key, value, timeout)
        else:
            # Add new entry
            new_row = pl.DataFrame({"key": [key], "val": [data], "exp": [expire]})

            self.cache_df = pl.concat([cache_df, new_row])

            # Save back
            self._save_cache()

    def clear(self):
        """Clear a cache"""

        # Create empty DataFrame
        self.cache_df = pl.DataFrame({"key": [], "val": [], "exp": []}, schema={"key": pl.Utf8, "val": pl.Utf8, "exp": pl.Float64})

        # Save empty cache
        self._save_cache()

    def __del__(self):
        """Cleans up the object"""
        # No connection to close with Parquet files
        pass


# allow this module to be used to clear the cache
if __name__ == "__main__":
    # check args
    if len(sys.argv) != 3 or sys.argv[1] != "clear":
        print("[!] Error: You have to specify the clear with `python %s clear <path to cache.parquet>`" % sys.argv[0])
        sys.exit(1)

    if not os.path.isdir(sys.argv[2]):
        print("[!] Error: %s does not seem to be a path!" % sys.argv[2])
        sys.exit(1)

    # setup the cache instance and clear it.
    c = ParquetCache(sys.argv[2])
    c.clear()
    print(" * Cache cleared")
