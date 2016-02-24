import logging
from DesignPatterns.SingletonDecorator import Singleton

@Singleton
class LoggingSetup:
    def __init__(self):
        self.initialized=False

    def initLogging(self):

        if not self.initialized:
            from utils import get_main_dir
            for handler in logging.getLogger().handlers:
                logging.getLogger().removeHandler(handler)

            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fileHandler=logging.FileHandler(get_main_dir()+"/Log.txt")
            fileHandler.setFormatter(formatter)
            logging.getLogger().addHandler(fileHandler)

            if len(logging.getLogger().handlers)<2:     ## ugly, but it works
                consoleHandler = logging.StreamHandler()
                logging.getLogger().addHandler(consoleHandler)

            logging.getLogger().setLevel(logging.INFO)

            self.initialized=True