import logging
import os

# Handle both relative imports (when used as module) and absolute imports (when run directly)
try:
    from .DesignPatterns.SingletonDecorator import Singleton
except ImportError:
    # When run directly, use absolute imports
    import sys
    import os as os_module

    # Handle case where __file__ might not be defined (e.g., when using exec())
    if "__file__" in globals():
        sys.path.insert(0, os_module.path.dirname(os_module.path.abspath(__file__)))
    else:
        sys.path.insert(0, os_module.path.dirname(os_module.path.abspath(os_module.getcwd() + "/src")))
    from DesignPatterns.SingletonDecorator import Singleton


@Singleton
class LoggingSetup:
    def __init__(self):
        self.location = None
        pass

    def initLogging(self, location=None):
        from .utils import get_main_dir

        while len(logging.getLogger().handlers) > 0:
            handler = logging.getLogger().handlers[0]
            logging.getLogger().removeHandler(handler)

        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        ## default log to PyMetExtract's directory
        fileHandler = logging.FileHandler(os.environ["LOCALAPPDATA"] + "/MetExtractII/Log.txt")
        fileHandler.setFormatter(formatter)
        logging.getLogger().addHandler(fileHandler)

        if location is not None:
            fileHandler = logging.FileHandler(location + "/Log.txt")
            fileHandler.setFormatter(formatter)
            logging.getLogger().addHandler(fileHandler)

        consoleHandler = logging.StreamHandler()
        logging.getLogger().addHandler(consoleHandler)

        logging.getLogger().setLevel(logging.INFO)
