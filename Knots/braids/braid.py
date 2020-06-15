from MakeBraid import *
import os

globalvars = {}       # We will store the calculator's variables here



def mult(x,y):
    elems = y.split("_")
    MyBraid[elems[0]].append([MyBraid['idx'],int(elems[1])])
    MyBraid['idx']=MyBraid['idx']+1
    return x+"*"+y
    
def div(x,y):
    elems = y.split("_")
    MyBraid[elems[0]+'i'].append([MyBraid['idx'],int(elems[1])])
    MyBraid['idx']=MyBraid['idx']+1
    return x+"/"+y

def lookup(map, name):
    for x,v in map:  
        if x==name: return v
    if name not in globalvars.keys(): print 'Undefined:', name
    return globalvars.get(name, 0)


from string import *
import re
from yappsrt import *

class sbraidScanner(Scanner):
    patterns = [
        ('"in"', re.compile('in')),
        ('"="', re.compile('=')),
        ('"let"', re.compile('let')),
        ('r"\\)"', re.compile('\\)')),
        ('"\\("', re.compile('\(')),
        ('"_"', re.compile('_')),
        ('"/"', re.compile('/')),
        ('"[*]"', re.compile('[*]')),
        ('"set"', re.compile('set')),
        ('[ \r\t\n]+', re.compile('[ \r\t\n]+')),
        ('END', re.compile('$')),
        ('GEN', re.compile('[abc]')),
        ('NUM', re.compile('[0-9]+')),
        ('VAR', re.compile('[d-zA-Z_]+')),
    ]
    def __init__(self, str):
        Scanner.__init__(self,None,['[ \r\t\n]+'],str)

class sbraid(Parser):
    def goal(self):
        _token_ = self._peek('"set"', 'GEN', 'NUM', 'VAR', '"\\("', '"let"')
        if _token_ != '"set"':
            expr = self.expr([])
            END = self._scan('END')
            print '=', expr
            return expr
        else:# == '"set"'
            self._scan('"set"')
            VAR = self._scan('VAR')
            expr = self.expr([])
            END = self._scan('END')
            globalvars[VAR] = expr
            print VAR, '=', expr
            return expr

    def expr(self, V):
        elem = self.elem(V)
        v = elem
        while self._peek('"[*]"', '"/"', 'END', 'r"\\)"', '"in"', '"_"') in ['"[*]"', '"/"']:
            _token_ = self._peek('"[*]"', '"/"')
            if _token_ == '"[*]"':
                self._scan('"[*]"')
                elem = self.elem(V)
                v = mult(v,elem)
            else:# == '"/"'
                self._scan('"/"')
                elem = self.elem(V)
                v = div(v,elem)
        return v

    def elem(self, V):
        term = self.term(V)
        n = term
        while self._peek('"_"', '"[*]"', '"/"', 'END', 'r"\\)"', '"in"') == '"_"':
            self._scan('"_"')
            term = self.term(V)
            if MyBraid['idx']==2:
                MyBraid[n].append([2,int(term)])
                MyBraid['idx']=3
            n = n+"_"+term                
        return n

    def term(self, V):
        _token_ = self._peek('GEN', 'NUM', 'VAR', '"\\("', '"let"')
        if _token_ == 'GEN':
            GEN = self._scan('GEN')
            return GEN
        elif _token_ == 'NUM':
            NUM = self._scan('NUM')
            return NUM
        elif _token_ == 'VAR':
            VAR = self._scan('VAR')
            return lookup(V, VAR)
        elif _token_ == '"\\("':
            self._scan('"\\("')
            expr = self.expr(V)
            self._scan('r"\\)"')
            return expr
        else:# == '"let"'
            self._scan('"let"')
            VAR = self._scan('VAR')
            self._scan('"="')
            expr = self.expr(V)
            V = [(VAR, expr)] + V
            self._scan('"in"')
            expr = self.expr(V)
            return expr


def parse(rule, text):
    P = sbraid(sbraidScanner(text))
    return wrap_error_reporter(P, rule)



if __name__=='__main__':
    print '?'*31
    print 'Welcome to the sBraid Calculator'
    print 'enter group elements in the form a_i'
    print '?'*31

    # We could have put this loop into the parser, by making the
    # `goal' rule use (expr | set var expr)*, but by putting the
    # loop into Python code, we can make it interactive (i.e., enter
    # one expression, get the result, enter another expression, etc.)
    while 1:
        MyBraid = {'a':[],
                   'ai':[],
                   'b':[],
                   'bi':[],
                   'c':[],
                   'ci':[],
                   'idx':2}
        try: s = raw_input('>>> ')
	except EOFError: break
        if not strip(s): break
        parse('goal', s)
        print MyBraid
        gens = {'magbraid':MyBraid['a'],
                'magibraid':MyBraid['ai'],
                'magbbraid':MyBraid['b'],
                'magibbraid':MyBraid['bi'],
                'magnot':MyBraid['c'],
                'maginot':MyBraid['ci']
                }

        j=0
        for i in MyBraid.keys():
            if i not in ['idx']:
                try:
                    for k in MyBraid[i]:
                        if k[1]>j:
                            j=k[1]
                except ValueError:
                    j = j
        print 'j', j
        j = j+1
        ids = GetIdsM(MyBraid['idx'],j,gens)

        gens.update({'magident':ids})
        eBraid = 'temp'
        MakeBraidM("./from_BraidCalc/"+eBraid+".ps",gens)

    print 'Bye.'
    print MyBraid


