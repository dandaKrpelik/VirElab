# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="darkgrid")

eps = 1e-10

class UC:
    def __init__(self, name, x, Dx = None, dx = None):
        self.name = name
        self._x = x
        self._Dx = Dx
        self._dx = dx
        
        if Dx is None and dx is None:
            self._Dx = 0
            self._dx = 0
    
    def sample(self):
        u = 2*np.random.rand()-1
        return self.x + u*self.Dx
    
    def __str__(self):
        out = '({:.4g}, {:.4g}, {:.4g} %) is {}'.format( self.x,self.Dx,100*self.dx, self.name)
        return out
    
    @property
    def x(self):
        return self._x
    
    def set_x(self,x):
        self.x = x
    
    @property
    def Dx(self):
        if self._Dx is None:
            self._recalc_Dx()
        return self._Dx
    @property
    def dx(self):
        if self._dx is None:
            self._recalc_dx()
        return self._dx
    @property
    def dxp(self):
        return self.dx * 100
    
    def _recalc_Dx(self):
        assert not self._dx is None
        self._Dx = abs(self.x) * self._dx
    def _recalc_dx(self):
        assert not self._Dx is None
        self._dx = self._Dx / abs(self.x)
    

    
    def __add__(self, o):
        if type(o) is UC:
            return UC( '({})+({})'.format(self.name, o.name), self.x + o.x , self.Dx + o.Dx )
        else:
            return UC( '({})+({:.3g})'.format(self.name, o), self.x + o , self.Dx )
    def __mul__(self,o):
        if type(o) is UC:
            return UC( '({})*({})'.format(self.name, o.name), self.x * o.x , dx = self.dx + o.dx )
        else:
            return UC( '({})*({:.3g})'.format(self.name, o), self.x * o , dx = self.dx )
    def __pow__(self,o):    
        if type(o) is UC:
            x = self.x ** o.x
            Dx = self.Dx * o * self.x**(o.x-1) + o.Dx * x * np.log( self.x )
            return UC( '({})^({})'.format(self.name, o.name), x, Dx = abs(Dx))
        else:
            x = self.x ** o
            Dx = self.Dx * o * self.x**(o-1)
            return UC( '({})^({:.3g})'.format(self.name, o), x , Dx = abs( Dx ))
        
    def __sub__(self,o):
        return self + ((-1)*o)
    def __truediv__(self,o):
        return self * (o ** (-1))
    
    
    
    def __radd__(self, o):
        return self + o
    def __rmul__(self, o):
        return self * o
    def __rsub__(self,o):
        return o + (-1) * self
    def __rtruediv__(self,o):
        return self **(-1) * o
    
    def sqrt(self):
        return UC( 'sqrt({})'.format(self.name), np.sqrt(self.x), dx = self.dx/2 )

def print_var_table(rvs, labels = None, sigdig = 4, rows = None, figformat = 'f'):
    if rows is None:
        rows = ['X','X_M (j)', 'DXM (j)', 'dXM (%)']
        
    no_format = '[:.{}{}]'.format(sigdig,figformat)
    no_format = no_format.replace('[','{')
    no_format = no_format.replace(']','}')
    
    n = len(rvs)
    out = '------------------------------\n'
    out += '---\t\t UQ RV table\n'
    
    if labels is None:
        labels = [rv.name for rv in rvs]
        
    out += rows[0] + '\t|'
    for i in range(n):
        out += ' \t|'+labels[i]
    out += '||\n'
    
    out += rows[1] + '\t|'
    for i in range(n):
        out += ' \t|'+no_format.format( rvs[i].x )
    out += '||\n'    
    out += rows[2] + '\t|'
    for i in range(n):
        out += ' \t|'+no_format.format( rvs[i].Dx )
    out += '||\n'    
    out += rows[3] + '\t|'
    for i in range(n):
        out += ' \t|'+no_format.format( rvs[i].dxp )
    out += '||\n'
    out += '---------------------------------------\n'
    out += '---------------------------------------\n'
    print(out)
    
def print_var_ltable(rvs, labels = None, sigdig = 4, rows = None, figformat = 'g'):
    if rows is None:
        rows = ['$X$','$X_M$ (j)', '$\Delta_{XM}$  (j)', '$\delta_{XM}$ (\%)']
    
    no_format = '[:.{}{}]'.format(sigdig,figformat)
    no_format = no_format.replace('[','{')
    no_format = no_format.replace(']','}')
    
    n = len(rvs)
    out = '\\begin{tabular}{|r|'
    for i in range(n):
        out += ' c '
    out += '|}\n'
    
    if labels is None:
        labels = [rv.name for rv in rvs]
        
    out += '\\hline\n'+rows[0]
    for i in range(n):
        out += ' & ' + labels[i]
    out += '\\\\\n'
    
    out += '\\hline\n'+rows[1]
    for i in range(n):
        out += ' & '+no_format.format( rvs[i].x )
    out += '\\\\\n'    
    out += rows[2]
    for i in range(n):
        out += ' & '+no_format.format( rvs[i].Dx )
    out += '\\\\\n'    
    out += rows[3]
    for i in range(n):
        out += ' & '+no_format.format( rvs[i].dxp )
    out += '\\\\\n\\hline\n'
        
    out += '\\end{tabular}\n'
    print(out)
    
    
############# PRISTROJE    

class Merak:
    def __init__(self, name, ranges = [0.5,1,2,5,10,50,150], errs = (1/100, 1/100, 1) , unit = None):
        ## errs (a,b,c) : (of rdg (%/100), of rng (%/100), lsd (count))
        ## chyby v rel. jednotkach (ne procentech!)
        self.name = name
        
        if not type(ranges) in [list, np.ndarray, tuple]:
            ranges = [ ranges ]
        self.ranges = np.array( ranges )
        
        self.unit = unit if unit else '-'
        self.set_unit4latex()
        
        if type(errs) is tuple:     ## stejna chyba pro vsechny rozsahy
            self.errs = errs
        elif type(errs) is list:     
            if type(errs[0]) is tuple:      ## list vsech chyb pro kazdy rozsah
                self.errs = errs
            else:   
                self.errs = [(0,tp,1) for tp in errs]   ## list s tridami presnosti, dig=1 kvuli chybe cteni
        else:
            self.errs = (0,errs,1)          ## jedina TP - stejna pro rozsahy

    def read(self, Xs):
        out = UC( 'read('+self.name+')', Xs, self.err(Xs) )
        return out
    def rread(self, Xs):
        ## adds epistemic noise to readings
        out = UC( self.name + '_sample', Xs, self.err(Xs) )
        out = self.read(out.sample())
        return out

    def set_unit4latex(self):
        dick = {'OHM':'$\\Omega$'}        ## replace units for latex
        
        out = self.unit
        for key in dick:
            out = out.replace(key, dick[key])
            
        self.unit4latex = out


    def get_range_i(self, X):
        ## X is float
        rngs = self.ranges
        assert abs(X) < rngs[-1] , 'measured value is above the highest range'   ## X < max_range
        
        n = len(rngs)
        return n - len(rngs[ rngs > abs(X) ])
        
    def get_range(self, X):
        ## X is float
        rng_i = self.get_range_i(X)
        return self.ranges[ rng_i ]
    
    def get_lsd(self, X):
        ## X is float
        return 0.
    
    def get_errs(self, i):
        return self.errs[i] if type(self.errs) is list else self.errs
    
    def err(self, Xin):
        ## X is float
        X = abs(Xin)
        rng_i = self.get_range_i(X)
        rng = self.ranges[rng_i]
        lsd = self.get_lsd(X)
        
        errs = self.get_errs(rng_i)
        return X * errs[0] + rng * errs[1] + lsd * errs[2]

        
    def __str__(self):
        return '( {} ) '.format(self.unit)+self.name
    
    def strlatex(self):
        return '( {} ) '.format(self.unit4latex)+self.name
    
    @property
    def summary(self):
        skip = '\t|\t'
        
        out = '---------------------------------\n'
        out += '----------------------- ! MP !--\n'
        out += '== '+str(self)+'\n'
        out += '-----------------------\n'
        out += skip + ' range \t\t| lsd (Km) \t| err (/rdg, /rng, *lsd) \n'
        for ri in range(len(self.ranges)):
            rng = self.ranges[ri]
            rng_ = rng-eps
            
            errs = self.get_errs(ri)
            err_txt = '{:3.4g} (%), {:3.4g} (%), {:2d} (dig)'.format(100*errs[0],100*errs[1], errs[2])
            out += skip + '{:.3e} \t| {:.3e} \t| {}\n'.format(rng, self.get_lsd(rng_), err_txt )
        out += '---------------------------------\n'
        return out
    @property
    def lsummary(self):
        
        out = '---------------------------------\n\n'
        out += '\\begin{tabular}{|r| l | c c c |}\n'
        out += '\\hline\n\\multicolumn{5}{c}{'+self.strlatex()+'} \\\\\n\\hline\n'
            
        out += '\\hline\n\t range ({0}) & lsd ({0}) & $/rdg$ (\%) & $/rng$ (\%) & dig \\\\\n\t\\hline\n'.format(self.unit4latex)
        for ri in range(len(self.ranges)):
            rng = self.ranges[ri]
            rng_ = rng-eps
            
            out += '\t '+'{:.3e} & {:.3e} '.format(rng, self.get_lsd(rng_) )        
            
            errs = self.get_errs(ri)
            out += '& {:3.4g} & {:3.4g} & {:2d} \\\\\n'.format(100*errs[0],100*errs[1], errs[2])
        
        out += '\t\\hline\n'
        out += '\\end{tabular}\n\n'
        out += '---------------------------------\n'
        return out


class anMerak(Merak):
    def __init__(self, name, ndilku, rngs, TPs, cast_dilku = 1, unit = None):
        super().__init__(name, ranges = rngs, errs = TPs , unit = unit)
        self.ndilku = ndilku
        self.cast_dilku = cast_dilku    ## rozliseni uzivatele: cd = 3 ->> zaokrouhlujeme na tretinu dilku
    
    def read(self, Xs):
        rng = self.get_range(Xs)
        K = rng / self.ndilku
        alpha = Xs / K
        alpha = round(alpha * self.cast_dilku) / self.cast_dilku
        return super().read(K * alpha)       

    def get_lsd(self, X):
        rng = self.get_range(X)
        return rng / ( self.ndilku * self.cast_dilku) ## rozliseni uzivatele: cd = 3 ->> zaokrouhlujeme na tretinu dilku

    def __str__(self):
        return super().__str__()+' Stupn.:'+str(self.ndilku)

class digMerak(Merak):
    
    def __init__(self, name, dd = (5,3), errs = (1/100, 2/100, 3), rngs = None, unit = None):
        super().__init__( name, rngs ,errs, unit = unit)
        self.dig_disp = dd ## (5,3) ### 4999
        
        if rngs is None:
            ranges = []
            for i in range(dd[1]+1):
                ranges += [ 10 ** i ]
            
            if type(errs) is list: assert len(errs) == len(ranges), 'jinak dlouhe tabulky'
            self.ranges = dd[0] * np.array( ranges )
        
    
    def get_lsd(self, X):
        dd = self.dig_disp
        t = int(np.log10(abs(X)/dd[0])-dd[1])   ### lsd = 10^t
        return 10**t

    def read(self, X):
        exp = int(np.log10( self.get_lsd(X) ))
        Xr = round(X, -exp)                 ## zaokrouhluje na rozliseni displeje v danem rozsahu
        return super().read( Xr )

    def __str__(self):
        return super().__str__()+' Disp:'+str(self.dig_disp)
    
