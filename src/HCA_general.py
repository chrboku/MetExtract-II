from math import pow, sqrt

from abc import ABCMeta, abstractmethod

from .DesignPatterns.IDSingleton import IDProvider

import logging


class HCANode:
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def getLeaves(self):
        pass

    @abstractmethod
    def getDistance(self):
        pass


class HCALeaf(HCANode):
    def __init__(self, obj, _id, _ind):
        self.obj = obj
        self._id = _id

        self._ind = _ind

    def getObj(self):
        return self.obj

    def getID(self):
        return self._id

    def getLeaves(self):
        return [self]

    def getDistance(self):
        return 0.0


class HCAComposite(HCANode):
    def __init__(self, kidL, kidR, distance, _id=None):
        self.kids = [kidL, kidR]
        self.distance = distance
        if _id is None:
            _id = "_%d" % IDProvider.Instance().getNextID()
        self._id = _id

    def getLeftKid(self):
        return self.kids[0]

    def getRightKid(self):
        return self.kids[1]

    def getLeaves(self):
        leaves = []
        for kid in self.kids:
            leaves.extend(kid.getLeaves())
        return leaves

    def getDistance(self):
        return self.distance

    def getID(self):
        return self._id


def euclidean(a, b):
    assert len(a) == len(b)

    su = sum([pow(a[i] - b[i], 2) for i in range(len(a))])
    return sqrt(su)


def average(clust):
    kids = clust.getLeaves()

    if len(kids) == 0:
        raise Exception("Strange Error in average of HCA. No kids in node.")

    sums = [0 for i in kids[0].getObj()]
    for kid in kids:
        for i in range(len(kid.getObj())):
            sums[i] = sums[i] + kid.getObj()[i]
    for i in range(len(sums)):
        sums[i] = sums[i] * 1.0 / len(kids)

    return sums


import numpy as np


class HCA_generic:
    def __init__(self, dist=euclidean, link=average):
        self.dist = dist
        self.link = link

        self.links = {}
        self.dists = {}

    def generateTree(self, objs, ids=None):
        # logging.info("  .. started generating tree for %d objects"%(len(objs)))
        if ids is None:
            ids = range(len(objs))
        clusts = [HCALeaf(objs[i], _id=ids[i], _ind=i) for i in range(len(objs))]
        nClusts = len(clusts)

        self.dists = np.ones((nClusts, nClusts), dtype=float) * 1e6
        self._updateDist(clusts)

        # logging.info("  .. calculated initial distances")

        while nClusts > 1:
            pos = np.argmin(self.dists)
            i = pos // self.dists.shape[1]  # Use integer division
            j = pos % self.dists.shape[1]
            distance = self.dists[i, j]

            c1 = min(i, j)
            c2 = max(i, j)

            comp = HCAComposite(kidL=clusts[c1], kidR=clusts[c2], distance=distance)
            self.getLinkFor(comp)
            clusts[i] = comp
            clusts[j] = None

            self._updateDist(clusts, updateInd=i)
            self.dists[j, :] = 1e6
            self.dists[:, j] = 1e6

            nClusts = nClusts - 1

        return clusts[0]

    def _updateDist(self, clusts, updateInd=None):
        nClusts = len(clusts)

        for i in range(nClusts - 1) if updateInd is None else [updateInd]:
            for j in range(i + 1, nClusts) if updateInd is None else range(nClusts - 1):
                cl1 = clusts[i]
                cl2 = clusts[j]

                if i != j and cl1 is not None and cl2 is not None:
                    dist = self.dist(self.getLinkFor(cl1), self.getLinkFor(cl2))
                else:
                    dist = 1e6

                self.dists[i, j] = dist
                self.dists[j, i] = dist

    def getLinkFor(self, cl):
        clID = cl.getID()
        if clID in self.links.keys():
            return self.links[clID]
        else:
            link = self.link(cl)
            self.links[clID] = link
            return link

    def getIndsFor(self, c1):
        return [leaf._ind for leaf in c1.getLeaves()]

    def plotTree(self, tree, intend="", newIntend="   "):
        if isinstance(tree, HCALeaf):
            print(intend, "Leaf ID:", tree.getID(), "Data:", tree.getObj())

        elif isinstance(tree, HCAComposite):
            self.plotTree(tree.getLeftKid(), intend=intend + newIntend, newIntend=newIntend)
            print(
                intend,
                "Composite + ID:",
                tree.getID(),
                "Distance:",
                tree.getDistance(),
                "Data:",
                self.getLinkFor(tree),
            )  # self.link(tree)
            self.plotTree(tree.getRightKid(), intend=intend + newIntend, newIntend=newIntend)

        else:
            raise Exception("Strange error in plotTree of HCA_generic.")

    def splitTreeWithCallback(self, tree, callBackFunction, recursive=True):
        if callBackFunction(tree, hca=self):
            if isinstance(tree, HCAComposite):
                subClusts = []
                if recursive:
                    subClusts.extend(self.splitTreeWithCallback(tree=tree.getLeftKid(), callBackFunction=callBackFunction))
                    subClusts.extend(self.splitTreeWithCallback(tree=tree.getRightKid(), callBackFunction=callBackFunction))
                else:
                    subClusts.append(tree.getLeftKid())
                    subClusts.append(tree.getRightKid())
                return subClusts
            else:
                return [tree]
        else:
            return [tree]

    def getObjsOrderInTree(self, tree):
        if isinstance(tree, HCALeaf):
            return [tree.getID()]

        elif isinstance(tree, HCAComposite):
            return self.getObjsOrderInTree(tree.getLeftKid()) + self.getObjsOrderInTree(tree.getRightKid())

        else:
            raise Exception("Strange error in getObjsOrderInTree of HCA_generic.")


if __name__ == "__main__":
    obs = [
        [1, 1, 1, 1],
        [1, 1.1, 2.2, 3.1],
        [2, 2, 2, 2],
        [3, 3, 3, 3],
        [1, 2.1, 3.1, 4],
        [1.1, 1.1, 1.1, 1.1],
    ]

    hc = HCA_generic()
    tree = hc.generateTree(obs)
    hc.plotTree(tree)
    print(hc.getObjsOrderInTree(tree))
