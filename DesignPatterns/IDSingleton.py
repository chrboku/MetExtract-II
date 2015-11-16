from SingletonDecorator import Singleton

@Singleton
class IDProvider:
    def __init__(self):
        self.curID=1

    def getNextID(self):
        ret=self.curID
        self.curID += 1
        return ret