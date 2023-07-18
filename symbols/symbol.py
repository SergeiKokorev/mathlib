from __future__ import annotations
from typing import List


class Expression:
    
    def __init__(self, exp: object):
        self.__builder = exp


class Symbol:
    
    symbols = {}

    def __init__(self, symbol: str, **kwd):
        self.__symbol = symbol
        self.__commutative = kwd.get('commutative', True)

    @property
    def s(self):
        return self.__symbol
    
    @property
    def commutative(self):
        return self.__commutative
    
    @commutative.setter
    def commutative(self, commutative: bool):
        if isinstance(commutative, bool):
            self.__commutative = commutative

    def __add__(self, s: Symbol) -> Expression:
        pass


def symbols(sym: str, *args, **kwds) -> List[Symbol]:

    res = []

    if not isinstance(sym, str):
        raise TypeError(f'Symbols must be string type not {type(sym)}')

    return tuple(Symbol(s.strip()) for s in sym.split(','))
        

