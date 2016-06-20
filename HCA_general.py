from math import pow, sqrt

from abc import ABCMeta, abstractmethod

from DesignPatterns.IDSingleton import IDProvider


class HCANode:
    __metaclass__=ABCMeta

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
    def __init__(self, obj, _id):
        self.obj=obj
        self._id=_id

    def getObj(self):
        return self.obj

    def getID(self):
        return self._id

    def getLeaves(self):
        return [self]

    def getDistance(self):
        return 0.


class HCAComposite(HCANode):
    def __init__(self, kidL, kidR, distance, _id=None):
        self.kids=[kidL, kidR]
        self.distance=distance
        if _id is None:
            _id="_%d"%IDProvider.Instance().getNextID()
        self._id=_id

    def getLeftKid(self):
        return self.kids[0]
    def getRightKid(self):
        return self.kids[1]

    def getLeaves(self):
        leaves=[]
        for kid in self.kids:
            leaves.extend(kid.getLeaves())
        return leaves

    def getDistance(self):
        return self.distance

    def getID(self):
        return self._id


def euclidean(a, b):
    assert len(a)==len(b)

    su=sum([pow(a[i]-b[i],2) for i in range(len(a))])
    return sqrt(su)

def average(clust):
    kids=clust.getLeaves()

    if len(kids)==0:
        raise Exception("Strange Error in average of HCA. No kids in node.")

    sums=[0 for i in kids[0].getObj()]
    for kid in kids:
        for i in range(len(kid.getObj())):
            sums[i]=sums[i]+kid.getObj()[i]
    for i in range(len(sums)):
        sums[i]=sums[i]*1./len(kids)

    return sums

class HCA_generic:
    def __init__(self, dist=euclidean, link=average):
        self.dist=dist
        self.link=link

        self.links={}
        self.dists={}

    def generateTree(self, objs, ids=None):
        if ids is None:
            ids=range(len(objs))
        clusts=[HCALeaf(objs[i], _id=ids[i]) for i in range(len(objs))]

        while len(clusts)>1:
            (i, j), distance=self.HCA_findNearestClusts(clusts)
            c1=min(i,j)
            c2=max(i,j)
            c2=clusts.pop(c2)
            c1=clusts.pop(c1)

            comp=HCAComposite(kidL=c1, kidR=c2, distance=distance)
            clusts.append(comp)

        return clusts[0]

    def HCA_findNearestClusts(self, clusts):
        nearestDist=1e6
        nearestPair=(0,0)
        for i in range(len(clusts)-1):
            for j in range(i+1, len(clusts)):
                d=self.getDistFor(clusts[i], clusts[j])     #d=self.dist(self.link(clusts[i]), self.link(clusts[j]))
                if d<nearestDist:
                    nearestDist=d
                    nearestPair=(i,j)
        return nearestPair, nearestDist

    def getDistFor(self, cl1, cl2):
        cl1ID=cl1.getID()
        cl2ID=cl2.getID()
        if cl1ID in self.dists.keys() and cl2ID in self.dists[cl1ID].keys():
            return self.dists[cl1ID][cl2ID]
        else:
            dist=self.dist(self.getLinkFor(cl1), self.getLinkFor(cl2))
            if cl1ID not in self.dists.keys():
                self.dists[cl1ID]={}
            if cl2ID not in self.dists.keys():
                self.dists[cl2ID]={}
            self.dists[cl1ID][cl2ID]=dist
            self.dists[cl2ID][cl1ID]=dist
            return dist

    def getLinkFor(self, cl):
        clID=cl.getID()
        if clID in self.links.keys():
            return self.links[clID]
        else:
            link=self.link(cl)
            self.links[clID]=link
            return link

    def plotTree(self, tree, intend="", newIntend="   "):
        if isinstance(tree, HCALeaf):
            print intend, "Leaf ID:" ,tree.getID(), "Data:", tree.getObj()

        elif isinstance(tree, HCAComposite):
            self.plotTree(tree.getLeftKid(), intend=intend+newIntend, newIntend=newIntend)
            print intend, "Composite + ID:" ,tree.getID(), "Distance:", tree.getDistance(), "Data:", self.getLinkFor(tree)#self.link(tree)
            self.plotTree(tree.getRightKid(), intend=intend+newIntend, newIntend=newIntend)

        else:
            raise Exception("Strange error in plotTree of HCA_generic.")

    def splitTreeWithCallback(self, tree, callBackFunction, recursive=True):
        if callBackFunction(tree, hca=self):
            if isinstance(tree, HCAComposite):
                subClusts=[]
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


if __name__=="__main__":
    obs=[[1,1,1,1],
         [1,1.1,2.2,3.1],
         [2,2,2,2],
         [3,3,3,3],
         [1,2.1,3.1,4],
         [1.1,1.1,1.1,1.1]]

    hc=HCA_generic()
    tree=hc.generateTree(obs)
    hc.plotTree(tree)
    print hc.getObjsOrderInTree(tree)
