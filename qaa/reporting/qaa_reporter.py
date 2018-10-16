import sys
from abc import ABC, abstractmethod

class QAAReporter(ABC):
    def __init__(self, indir, header=list(), out=sys.stdout):
        self.indir = indir
        self.outstream = out
        self.header = header
        super().__init__()
    @abstractmethod
    def generateReport(self):
        pass
