import sympy as sp
import os
import sys


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)
from symbols.symbol import *


if __name__ == "__main__":
    
    a, b, c, d = symbols('a, b, c, d')
    
    