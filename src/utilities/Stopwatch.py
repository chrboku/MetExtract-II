from time import time


class Stopwatch:
    def __init__(self):
        self._start = 0
        self._stop = 0

    def start(self):
        self._start = time()

    def stop(self):
        self._stop = time()

    def getDurationInSeconds(self):
        return self._stop - self._start


if __name__ == "__main__":
    sw = Stopwatch()

    sw.start()

    sum = 0
    for i in range(100000000):
        sum = sum + i

    sw.stop()

    print("it took %.1f seconds to sum up the values from 0 to 10.000.000" % (sw.getDurationInSeconds()))
