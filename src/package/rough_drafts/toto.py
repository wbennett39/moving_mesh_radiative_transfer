#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 08:39:58 2022

@author: bennett
"""
from numba import int32
from numba.experimental import jitclass

spec2 = [('n',int32)]

@jitclass(spec2)
class toto(object):
  def __init__(self,n):
    self.n = 42 + n

  def work(self,y):
    return y + self.n
