## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
==============
Tree Utilities
==============

Build and log trees, get clades and node heights, convert to NEWICK format and
so on. 

Unless explicitly specified, any tree is assumed to be Ultrametric
(tips are contemporaneous).

"""

import operator, sys, os.path
from genericutils import fileFromName

# Bio.Nexus.Tree stuff

from Bio.Nexus import Trees, Nodes

__all__ = ["TreeBuilder", "TreeLogger", "getClade", "getTreeClades",
           "getCommonAncesstor", "countNexusTrees", "toNewick", "nodeHeights",
           "nodeHeight", "treeHeight", "setLabels", "convertDemographics",
           "coalLogLike"]


class TreeBuilder(object) :
  """ A basic helper for building BioPython trees.

  Use:
   - tb = TreeBuilder()
   - Create some terminals via leaf = tb.createLeaf(name)
   - use mergeNodes to successively merge nodes.
   - Call tb.finalize(root-node) to get the tree
  """
  
  def __init__(self, weight=1.0, rooted = True, name='') :
    self.t = Trees.Tree(weight=weight, rooted = rooted, name=name)

  def createLeaf(self, name) :
    """ Create a new leaf.

    :param name: taxon name
    :returns: node"""
    
    nd = Trees.NodeData()
    nd.taxon = name
    leaf = Nodes.Node(nd)
    self.t.add(leaf, None)
    return leaf

  def newNode(self) :
    nd = Trees.NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)
    return node
  
  def mergeNodes(self, subtrees) :
    """ Join subtrees in subtrees into one sub-tree.

    Each element in subtrees is a (node,branch) pair, where branch is the length
    of the branch connecting the node to its parent. 

    :param subtrees: sequence of [node,branch]
    :returns: node"""
    
    nd = Trees.NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)

    for n1,h1 in subtrees:
      n1.set_prev(node.id)
      n1.data.branchlength = h1

    node.add_succ([x.id for x,h in subtrees])
    return node

  def finalize(self, rootNode) :
    """ Convert a node to a proper tree with root rootNode.
    
    :param rootNode: node
    :returns: biopython tree
    """
    t = self.t
    rr = t.node(t.root)
    rr.set_succ(rootNode.succ)
    for p in rootNode.succ :
      t.node(p).set_prev(t.root)
    if hasattr(rootNode.data,"attributes") :
      rr.data.attributes = rootNode.data.attributes
    t.kill(rootNode.id)
    return t

class TreeLogger(object) :
  def __init__(self, outName = None, argv = None,
               version = None, overwrite = False) :
    self.outName = outName
    if outName is not None:
      if not overwrite and os.path.isfile(outName) \
             and os.path.getsize(outName) > 0 :
        raise RuntimeError("not overwriting existing file " +  outName)
      
      self.outFile = outFile = file(outName, "w")
      self.count = 0
      if argv is not None :
        print >> outFile, "[Generated by %s%s]" % \
              (os.path.basename(argv[0]),
               (", version " + version) if version is not None else "")
        print >> outFile, "[%s]" % " ".join(argv)
        
      print >> outFile, "#NEXUS"
      print >> outFile, "begin trees;"
      
  def outTree(self, tree, treeAttributes = {'R' : None}, name = None) :
    c = ""
    if treeAttributes is not None :
      c = " ".join(["[&%s %s]" % (k,v) if v else "[&%s]" % k
                    for k,v in treeAttributes.items()]) + " "
    if self.outName :
      if not name :
        name = "tree_%d" % self.count
      print >> self.outFile, "tree %s = %s%s ;" % (name,c,tree)
      self.count += 1
    else :
      print "%s%s" % (c,tree)
      
  def close(self) :
    if self.outName :
      print >> self.outFile, "end;"
      self.outFile.close()

def getClade(tree, nodeId) :
  n = tree.node(nodeId)
  if n.data.taxon :
    return [nodeId,]
  return reduce(operator.add, [getClade(tree, x) for x in n.succ])

def getCommonAncesstor(tree, taxaIds) :
  return reduce(tree.common_ancestor, taxaIds)

def _getTreeClades_i(tree, nodeID) :
  node = tree.node(nodeID)
  if node.data.taxon :
    return ([node.data.taxon], [])

  cl = [_getTreeClades_i(tree, ch) for ch in node.succ]
  allt = apply(operator.add, [x[0] for x in cl])
  clades = [(allt,node)]
  for cl1 in [x[1] for x in cl if x[1]] :
    clades.extend(cl1)
  return (allt, clades)

def getTreeClades(tree, trivialClades = False):
  """ Clades of subtree as a list of (taxa-list, tree-node).

  taxa-list is a list of strings.
  """

  c = _getTreeClades_i(tree, tree.root)[1]

  if trivialClades :
    for n in tree.get_terminals():
      nd = tree.node(n)
      c.append(([nd.data.taxon], nd))
    
  return c


def countNexusTrees(nexFileName) :
  """ Number of trees in a nexus file."""
  nexFile = fileFromName(nexFileName)
  c = 0
  for l in nexFile:
    if l.startswith("tree ") :
      c += 1
  return c

def toNewick(tree, nodeId = None, topologyOnly = False, attributes = None) :
  """ BioPython tree or sub-tree to unique NEWICK format.

  Children nodes are sorted (via text), so representation is unique and does not
  depend on arbitrary children ordering.
  """
  
  if not nodeId :
    nodeId = tree.root
  node = tree.node(nodeId)
  data = node.data
  if data.taxon :
    rep = data.taxon
  else :
    reps = [toNewick(tree, n, topologyOnly, attributes) for n in node.succ]
    reps.sort()
    rep = "(" + ",".join(reps) + ")"

  if attributes is not None :
    attrs = getattr(data, attributes, None)
    if attrs is not None and len(attrs) :
      s = '[&'
      for a in attrs:
        s += (a+'='+str(attrs[a])+",")
      s = s[:-1] + ']'
      rep += s
    
  if not topologyOnly and nodeId != tree.root and data.branchlength is not None:
    rep = rep + (":%r" % data.branchlength)
  return rep


def _getNodeHeight(tree, n, heights, w = 0):
  if n.id in heights:
    return heights[n.id]
  
  if len(n.succ) == 0 :
    heights[n.id] = 0.0
    return 0.0

  i = n.succ[w]
  if i not in heights:
    if n.succ[1-w] in heights:
      w = 1-w
      i = n.succ[w]
      
  c = tree.node(i)
  h = _getNodeHeight(tree, c, heights, 1-w) + c.data.branchlength
  
  heights[n.id] = h
  return h

def _collectNodeHeights(tree, nid, heights) :
  node = tree.node(nid)
  
  if len(node.succ) == 0 :
    heights[node.id] = 0.0
    return (node.data.branchlength, [node.id])

  (hs0,tips0), (hs1,tips1) = [_collectNodeHeights(tree, c, heights)
                               for c in node.succ]
  if hs0 != hs1 :
    if hs0 < hs1 :
      h = hs1
      dx = hs1-hs0
      # add hs1-hs0 to all left (0) side
      for n in tips0 :
        heights[n] += dx
    else :
      # add hs0-hs1 to all right (1) side 
      h = hs0
      dx = hs0-hs1
      for n in tips1 :
        heights[n] += dx
  else :
    h = hs0
  heights[node.id] = h
  return (h + node.data.branchlength, [nid] + tips0 + tips1)

def _collectNodeHeightsSimple(tree, nid, heights) :
  node = tree.node(nid)
  
  if len(node.succ) == 0 :
    heights[node.id] = 0.0
    return node.data.branchlength

  h = max([_collectNodeHeightsSimple(tree, c, heights) for c in node.succ])
  heights[node.id] = h
  return h + node.data.branchlength

def nodeHeights(tree, nids = None, allTipsZero = True) :
  """ Return a mapping from node ids to node heights.
  Without nids - for all nodes.

  The mapping may contain heights for other nodes as well.

  With !allTipsZero, handle non-ultrametric trees as well.
  """

  heights = dict()

  if allTipsZero :
    if nids == None :
      _collectNodeHeightsSimple(tree, tree.root, heights)
    else :
      w = 0
      for nid in nids :
        _getNodeHeight(tree, tree.node(nid), heights, w)
        w = 1-w
  else :
      # have to scan all tree to account for tip heights, ignore nids if
      # present
      _collectNodeHeights(tree, tree.root, heights)
     
  return heights

def nodeHeight(tree, nid) :
  """ Height of node. """
  
  node = tree.node(nid)
  if not node.succ :
    h = 0
  else :
    h = max([tree.node(c).data.branchlength + nodeHeight(tree, c)
             for c in node.succ])
    
  return h

def treeHeight(tree) :
  """ Height of tree. """
  return _getNodeHeight(tree, tree.node(tree.root), dict())


def setLabels(trees) :
  """ Convert labels attribute to tree node member (python list).

  Return boolean array per tree indicating if tree has complete meta-data.
"""
  hasAll = []
  for tree in trees :
    has = True 
    for i in tree.get_terminals() :
      data = tree.node(i).data
      ok = False
      if hasattr(data, "attributes") :
        if "labels" in data.attributes:
          l = data.attributes["labels"].split(' ')
          data.labels = l
          ok = len(l) > 0
      if not ok :
        has = False
    hasAll.append(has)

  return hasAll

import demographic
LPP = demographic.LinearPiecewisePopulation

def _tod(xt, yt, b = None) :
  xt = [float(x) for x in xt.split(',') if len(x)]
  yt = [float(x) for x in yt.split(',')]
  if b is not None and len(xt) + 2 == len(yt) :
    xt.append(b)
    
  return LPP(yt, xt)

def _toDemog(dtxt) :
  return _tod(*dtxt.split('|'))

def _toDemog1(dmt, dmv, branch) :
  return _tod(dmt if dmt is not None else "", dmv, b = branch)

def convertDemographics(tree, formatUnited = "dmf",
                        formatSeparated = ("dmv", "dmt"),
                        dattr = "demographic") :
  """ Convert demographic function stored in BEAST trees attributes to biopy
  demographic.

  Support old-style dmf attribute and new-style dmv,dmt
  """
  
  dmf = formatUnited
  dmv,dmt = formatSeparated

  missing = 0
  for i in tree.all_ids() :
    data = tree.node(i).data
    d = None
    if hasattr(data, "attributes") :
      if dmf in data.attributes:
        dtxt = data.attributes[dmf]
        d = _toDemog(dtxt)
      elif dmv in data.attributes:
        d = _toDemog1(data.attributes.get(dmt), data.attributes.get(dmv),
                      data.branchlength)
    if d is not None :
      setattr(data, dattr, d)
    else :
      missing += 1
  return missing


def _toAttrText(vals) :
  if len(vals) > 1 :
    return "{" + ",".join(["%f" % x for x in vals]) + "}"
  return "%f" % vals[0]
   
def revertDemographics(tree, dattr = "demographic",
                       formatSeparated = ("dmv", "dmt")) :
  dmv, dmt = formatSeparated
  
  for i in tree.all_ids() :
    data = tree.node(i).data

    d = getattr(data, dattr, None)
    if d :
      if not hasattr(data, "attributes") :
        data.attributes = dict()
      else :
        for l in formatSeparated :
          if l in data.attributes:
            del data.attributes[l]
            
      if isinstance(d, demographic.ConstantPopulation) :
        data.attributes[dmv] = d.population(0)
      elif isinstance(d, demographic.LinearPiecewisePopulation) :
        data.attributes[dmv] = _toAttrText(d.vals)
        if len(d.xvals) :
          data.attributes[dmt] = _toAttrText(d.xvals)
      elif isinstance(d, demographic.StepFunctionPopulation) :
        data.attributes[dmv] = _toAttrText(d.vals)
        data.attributes[dmt] = _toAttrText([0] + d.xvals)
      else :
        raise RuntimeError("unsupported")

import coalescent

def coalLogLike(tree, demog, condOnTree = False) :
  nh = nodeHeights(tree, allTipsZero = False)  
  terms = tree.get_terminals()
  
  times = sorted([(t,nid not in terms) for nid,t in nh.items()])

  return coalescent.coalLogLike(demog, times, condOnTree)
