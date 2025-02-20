#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 13:13:06 2025

@author: cox
"""
from sympy import *
init_printing(pretty_print=True)

T, s2a, s2e = symbols("T, sigma^2_\\alpha, sigma^2_\\epsilon", positive=True)

omega = 1-sqrt(s2e)/sqrt(Ti*s2a + s2e)


vary  = ((1-omega)**2 * s2a + s2e + omega**2 * (s2e/Ti) - 2*omega*s2e/Ti)/s2e
simplify(vary)







