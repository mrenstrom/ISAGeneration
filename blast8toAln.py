#!/usr/bin/python3.4
#
# read blast8 output, reduce to translation table of seq to align (or multi)
#
#
# fasta:
#
#>MAPRIS:1:10:100
#CACCATGACCAAGTCGGATTTAGCTGCCGA.......................
#
# Blast Tabular Format:
# 1.	 qseqid	  - query (e.g., gene) sequence id
# 2.	 sseqid	  - subject (e.g., reference genome) sequence id
# 3.	 pident	  - percentage of identical matches
# 4.	 length	  - alignment length
# 5.	 mismatch - number of mismatches
# 6.	 gapopen  - number of gap openings
# 7.	 qstart	  - start of alignment in query
# 8.	 qend	  - end of alignment in query
# 9.	 sstart	  - start of alignment in subject
# 10.	 send	  - end of alignment in subject
# 11.	 evalue	  - expect value
# 12.	 bitscore - bit score
#
#             ID         chr    match%  len     mismatch  gap open      qStart  qEnd  chr_st    chr_end        ev       bit score
#--------------------------------------------------------------------------------------------------------------------------------
#MAPRIS:1:10:100	chr7	100.00	39	0  	  0		1	39    135229908	135229870	4.9e-14	76.0
# .
# .
# .
#
#
from sys import argv
import sys
import csv
import re
from decimal import Decimal
from collections import defaultdict
from collections import namedtuple
import operator
   
#---------------------------------------------------------------------------------------------------------------
#
# Alignment Class - used to store alignment values 
#
#
#
#---------------------------------------------------------------------------------------------------------------
class Alignment():
  def __init__(self,cChr,cStart,cEnd):
    self._id       = ''
    self._fileID   = ''
    self._cloneID  = ''
    self._class    = 'S'
    self._aChr     = cChr
    self._aMatch   = 0.0
    self._aLen     = 0
    self._aStart   = 0
    self._aEnd     = 0
    self._chrStart = cStart
    self._chrEnd   = cEnd
    self._chrDir   = '+'
    self._q        = 0.0
    self._count    = 0
    self._span     = 0
    self.calcDir()
  def calcDir(self):
    if (self._chrStart > self._chrEnd): 
      self._chrDir = '-'
    else:
      self._chrDir = '+'
  #
  # any alignment match
  #
  def alignMatch(self,a):
    if (self._aChr == a._aChr) and (self._chrDir == a._chrDir) and (abs(self._chrStart - a._chrStart) < 100):
      return(True)
    return(False)
  #
  # merge
  #
  def merge(self,a):
    self._count += a._count
    a._class = 'X'
  #
  # print representation
  #
  def __repr__(self):
    return(self._id + " " + self._class + " " + self._aChr + " " + self._chrDir + " " + str(self._chrStart))
  
  def __str__(self):
    return(self._id + " " + self._class + " " + self._aChr + " " + self._chrDir + " " + str(self._chrStart))
#
#
#
#--------------------------------------------------------------------------------------------------------------
#
# Main program entry 
#
#--------------------------------------------------------------------------------------------------------------
print("---Blast8toTran.py---")
print(argv)
script ,blast8, outfile = argv
#
# store alignments in dict based on ID
#
#   each entry will be a dict with key = seq id and value = list of alignments
#
invalid = []
alDict = {}
failDict = {}
mergeLog = []
with open(blast8) as fh:
    k = 1
    if (k % 10000) == 0:
       print("Loop = " + str(k))
    k += 1

    for line in csv.reader(fh, delimiter='\t'):
      #
      # parse ID,COUNT,SPAN from fasta ID
      #       MAPRIS:file:clone:count
      #  ie.  MARRIS:22:100:44
      #i            0  1   2  
      idPos     =  [pos for pos, char in enumerate(line[0]) if char == ':']
      aFileID   = int(line[0][idPos[0]+1:idPos[1]])
      aCloneID  = int(line[0][idPos[1]+1:idPos[2]])
      aCount    = int(line[0][idPos[2]+1:])
      # 
      # std blast8 format
      #
      aLen     = int(line[3])
      aStart   = int(line[6])
      aChr     = line[1]
      chrStart = int(line[8])
      chrEnd   = int(line[9])
      #
      # validate
      #
      if aLen < 30:
        invalid.append((id,"Length"))
        continue
      #
      # start can be very displaced for some seq?
      #
      if aStart > 10:
        invalid.append((id,"start"))
        continue
      #
      # create new object
      #
      a = Alignment(aChr,chrStart,chrEnd)
      #
      # build object and add to list
      #
      a._id       = line[0]
      a._fileID   = aFileID
      a._cloneID  = aCloneID
      a._class    = 'S'   # single alignment - default
      a._count    = aCount
      a._span     = 0
      a._aMatch   = float(line[2])
      a._aLen     = aLen
      a._aStart   = aStart
      a._aEnd     = int(line[7])
      a._q        = float(line[11])
      a._mergeID  = "X"
      #
      # id already in dict?
      #
      try:
          # ID already in dict
          alDict[a._id].append(a)
      except KeyError:
          #ID is new
          alDict[a._id] = [a]
print("Read ", len(alDict), "good alignments")

#---------------------------------------------------------------------
#
# search top alignments
#
#---------------------------------------------------------------------
topDict = {}
midList = []
idSort = sorted(alDict.items(),key=lambda al: al[0])  # smallest first
for id,aList in idSort:
  #print(id + " #alignments = " + str(len(aList)))
  #
  # sort on score
  #
  alSort = sorted(aList,key=lambda al: al._q,reverse=True)
  #
  # r is the multi-alignment ratio, 0 = unique, 1 = total duplication
  #
  r = 0
  if (len(alSort) >= 2):
    r = alSort[1]._q / alSort[0]._q
  #
  # add multi-alignments >= 0.98 for multi-alignment search match
  #
  k = 0
  top = alSort[0]
  topList = [top]
  for al in alSort:
    if k != 0:
       q = al._q / top._q 
       if q > 0.95:
         topList.append(al)
    k += 1
  #
  # add to next storage struct
  #     
  midList.append((id,r,topList))
#
# debug print - sort on size
#
sort1 =  sorted(midList,key=lambda l: len(l[2]),reverse=True)
#for (k,r,l) in sort1:
#   print("{0:s}     {1:2.2f}     {2:n}     {3:n}".format(k,r,len(l),l[0]._count))
#
#---------------------------------------------------------------------
#
# write output
#
#---------------------------------------------------------------------
with open(outfile, 'w') as out:
  out.write('id\tcount\talignStr\tclass\tAlign\tchromo\tstart\tdir\tspan\tAlignments\tr\trefseq\n')
  for (id,r,alList) in sort1:
    for al in alList:
      cloneID  = al._id  #"MAPRIS:{0:2n}:{1:05n}".format(al._fileID,al._cloneID)
      alCount  = al._count 
      alSpan   = al._span
      startPos = al._chrStart
      endPos   = al._chrEnd
      dirStr   = al._chrDir
      alStr    = al._aChr + ":" + str(startPos) + ":" + dirStr

      out.write("{0:16s}\t{1:5n}\t{2:16s}\t".format(cloneID,alCount,al._aChr))
      out.write("{0:9n}\t".format(startPos))
      out.write("{0:s}\t".format(dirStr))
      out.write("{0:3.2f}\t".format(r))
      out.write("{0:3.2f}\n".format(al._q))
