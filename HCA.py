# generic hierarchical clustering algorithm (HCA)

import functools


# Generic HCA object used to store a component the tree that is either a leaf or a composite
class HCNode(object):
    # initialise the componet with a certain value
    def __init__(self, val=lambda x: x):
        self.val = val

    # returns value of the component
    def getValue(self):
        raise Exception("Abstract: Not implemented")

    # returns the object
    def getObject(self):
        raise Exception("Abstract: Not implemented")

    # returns the number of kid components
    def getKidsCount(self):
        raise Exception("Abstract: Not implemented")

    # returns all values of the kid components
    def getKidsValues(self):
        raise Exception("Abstract: Not implemented")

    # returns the kids themselves
    def getKids(self):
        raise Exception("Abstract: Not implemented")


# HCA leaf object
class HCLeaf(HCNode):
    def __init__(self, data, val=lambda x: x):
        super(HCLeaf, self).__init__(val=val)
        self.data = data

    def getValue(self):
        return self.val(self.data)

    def getObject(self):
        return self.data

    def getKidsCount(self):
        return 1

    def getKidsValues(self):
        return self.val(self.data)

    def getKids(self):
        return [self]


# HCA composite object
# implemented as a binary composite for a binary tree structure. It has one left and one right neighbour only
class HCComposite(HCNode):
    def __init__(self, left, right, val=lambda x: x, mean=None, add=None):
        assert mean is not None and add is not None

        super(HCComposite, self).__init__(val=val)
        self.mean = mean
        self.add = add
        self.left = left
        self.right = right

    def getValue(self):
        return self.mean(self.getKidsValues(), self.getKidsCount())

    def getObject(self):
        raise Exception("Composite: Use forbidden")

    def getKidsCount(self):
        return self.left.getKidsCount() + self.right.getKidsCount()

    def getKidsValues(self):
        return self.add(self.left.getKidsValues(), self.right.getKidsValues())

    # flatten all kids to an array
    def getKids(self):
        ret = []
        for kid in self.left.getKids():
            ret.append(kid)
        for kid in self.right.getKids():
            ret.append(kid)
        return ret


# generic method to compare different distances between HCA components
def distCmp(x, y, dist):
    d = dist(x, y)
    if d > 0:
        return 1
    if d < 0:
        return -1
    return 0


# performs the hierarchical clustering analysis and stores the results using HCLeaf and HCComposite
class HierarchicalClustering:
    def __init__(self, dataP, dist, val=lambda x: x, mean=None, add=None):
        assert mean is not None and add is not None

        if len(dataP) == 0:
            raise Exception("HC needs more than zero elements")

        # create a HCLeaf object for each provided data
        data = [HCLeaf(d, val=val) for d in dataP]

        # iteratively merge the two closest two HCNode objects, merge them into a new HCComposite object
        # and place it back to the clustering. Stop, if only one HCNode object remains
        while len(data) > 1:
            data.sort(key=functools.cmp_to_key(lambda x, y: distCmp(x, y, dist)))
            diff = [dist(data[x + 1], data[x]) for x in range(0, len(data) - 1)]
            minindex, minvalue = min(enumerate(diff), key=lambda x: x[1])
            d = HCComposite(
                data[minindex], data[minindex + 1], val=val, mean=mean, add=add
            )
            data.pop(minindex)
            data.pop(minindex)
            data.append(d)

        self.tree = data[0]

    def getTree(self):
        return self.tree


# HELPER METHOD that returns the HCLeaf object with the lowest value
def _getMinKid(node, minK=100000):
    if isinstance(node, HCLeaf):
        return min(minK, node.getValue())
    elif isinstance(node, HCComposite):
        return min(minK, _getMinKid(node.left, minK), _getMinKid(node.right, minK))
    else:
        raise Exception("HCA Clustering error")


# HELPER METHOD that returns the HCLeaf object with the highest value
def _getMaxKid(node, maxK=-1):
    if isinstance(node, HCLeaf):
        return max(maxK, node.getValue())
    elif isinstance(node, HCComposite):
        return max(maxK, _getMaxKid(node.left, maxK), _getMaxKid(node.right, maxK))
    else:
        raise Exception("HCA Clustering error")


# cut HCA tree according to a given maximal sub-cluster size (i.e. maximal ppm deviations
# between lowest and highest mz value)
def cutTreeSized(node, ppm):
    if isinstance(node, HCLeaf):
        return [node]
    if isinstance(node, HCComposite):
        minV = _getMinKid(node)
        maxV = _getMaxKid(node)

        if (1000000.0 * (maxV - minV) / node.getValue()) < ppm:
            return [node]
        else:
            ret = []
            for x in cutTreeSized(node.left, ppm):
                ret.append(x)
            for x in cutTreeSized(node.right, ppm):
                ret.append(x)
            return ret
    else:
        raise Exception("HCA Clustering error")


# print(HCA tree (recursive with intend))
def printTree(node, inlet=""):
    if isinstance(node, HCLeaf):
        print(inlet + str(node.getValue()))
    elif isinstance(node, HCComposite):
        printTree(node.left, inlet + "  ")
        print(inlet + ":" + str(node.getValue()))
        printTree(node.right, inlet + "  ")
    else:
        print("Error")


if __name__ == "__main__":
    data = [
        1,
        3,
        5,
        6,
        7,
        22,
        23,
        25,
        26,
        33,
        39,
        43,
        45,
        46,
        65,
        66,
        68,
        69,
        71,
        72,
        82,
        83,
        86,
        95,
        96,
        99,
    ]
    hc = HierarchicalClustering(
        data,
        dist=lambda x, y: abs(x.getValue() - y.getValue()),
        val=lambda x: x,
        mean=lambda x, y: x / y,
        add=lambda x, y: x + y,
    )

    printTree(hc.getTree())
