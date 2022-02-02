#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 08:36:23 2022

@author: bennett
"""

from numba import float64, int32, deferred_type, jit
from numba.experimental import jitclass
from datetime import datetime
from toto import toto

toto_type = deferred_type()
toto_type.define(toto.class_type.instance_type)

spec = [('a',float64),('b',float64),('c',toto_type)]

@jitclass(spec)
class fifi(object):
  def __init__(self, combis):
    self.a = combis
    self.b = 2
    self.c = toto(combis)
    self.d = "onetwo3"

  def mySqrt(self,x):
    s = x
    for i in range(self.a):
      s = (s + x/s) / 2.0
    return s





@jit(nopython = True)
def run(n,results):
  for i in range(n):
    q = fifi(200)
    results[i+1] = q.mySqrt(i + 1)

if __name__ == '__main__':
  n = int(1e6)
  results = [0.0] * (n+1)
  starttime = datetime.now()
  run(n,results)
  endtime = datetime.now()

  print("Script running time: %s"%str(endtime-starttime))
  print("Sqrt of 144 is %f"%results[145])