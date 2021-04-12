import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as mani
from enum import Enum

import os	## for relative imports
script_dir = os.path.dirname(__file__)
to_path = lambda p : os.path.relpath(p)

import time
T0 = time.time()
PI = 3.1415926535

class Trig(Enum):
    Ascending = 0
    Descending = 1

class Ring:
    def __init__(self, n = 1):
        self.x = []
        self.i = -1
        self.full = False
        self.n = n
        
    def set_n(self, n):
        if n > self.n: self.full = False
        if n < self.n:
            if self.full or self.i > n:
                self.x = x[:n]
                
            self.full = self.full or self.i <= n
            self.i = self.i % n
        self.n = n
    
    def put(self, x):
        self.i += 1
        if self.i == self.n:
            self.full = True
            self.i %= self.n
        
        if self.full:
            self.x[self.i] = x
        else:
            self.x += [x]
        
    def __iter__(self):
        self.iter_counter = 0
        return self
    def __next__(self):
        i = self.iter_counter
        if i >= len(self.x):
            raise StopIteration
        
        self.iter_counter += 1
        return self.x[i]
        
            

# figure preparation
#fig, ax = plt.subplots(1, 1, figsize = (8, 6))

class Signal:
    def __init__(self, offset = 0, spectrum = [], noise = 0, square = [], trap = []):
        self.offset = offset
        self.spectrum = spectrum    ## [ [f, A, phi] , ...    ]
        self.square = square        ## [ [f, A, phi, duty], ... ]
        self.trap = trap					## [ [f, A, phi, a,b,c,d,e] ]				## a first max, b end first max, c zero, d first min, e end min
        self.noise = noise
        
    
    def __call__(self, t):
        if type(t) in [list,tuple,np.ndarray]:
            out = np.zeros_like(t, dtype=float)        
            out += self.noise * np.random.randn( len(out) )
        else:
            out = 0
            out += self.noise * np.random.randn()

        for s in self.spectrum:
            p = 0 if len(s) < 3 else s[2]
            out += s[1] * np.sin( 2*PI*s[0]*(t+p) )       
        for s in self.square:
            T = 1./s[0]
            p = 0 if len(s) < 3 else s[2]
            tau = ( t % T - p) % T
            d = 0.5*T if len(s) < 4 else s[3]*T
            out += s[1] * (tau < d)
        
        for s in self.trap:
            T = 1./s[0]
            A = s[1]
            phi = s[2]
            
            a,b,c,d,e = s[3:]
            
            tau = ( t % T - phi) % T
            tau /= T
            if type(t) in [list,tuple,np.ndarray]: 
                out += A* tau * 1/a
                out[tau > a] += A*(tau[tau > a] -a )* (-1/a)
                out[tau > b] += A*(tau[tau > b] -b )* (-1/(c-b))
                out[tau > c] += A*(tau[tau > c] -c )* (1/(c-b) -1/(d-c))
                out[tau > d] += A*(tau[tau > d] -d )* (1/(d-c))
                out[tau > e] += A*(tau[tau > e] -e )* (1/(1-e))
            else:
                out += A* tau * 1/a
                if tau > a: out += A*(tau-a )* (-1/a)
                if tau > b: out += A*(tau-b )* (-1/(c-b))
                if tau > c: out += A*(tau-c )* (1/(c-b) -1/(d-c))
                if tau > d: out += A*(tau-d )* (1/(d-c))
                if tau > e: out += A*(tau-e )* (1/(1-e))
				
        return out + self.offset
    


class Channel:
    ID_CNT = 0
    colors = ['yellow',  'deepskyblue',  'orangered', 'lawngreen',  'magenta', 'dodgerblue']
    def __init__(self, signal = Signal()):
        self.id = Channel.ID_CNT
        Channel.ID_CNT += 1
        self.name = 'ch'+str(self.id)
        self.signal = signal
        self.active = True
        
        self.voltdiv = 2e-7
        
        self.trig = 0
        self.trig_edge = Trig.Ascending
        
        self.dv = 0
        self.dh = 0
        
        self.invert = False
        self.AC = False
        
    def set_id(self, i):
        self.id = i
        self.name = 'ch'+str(i)
        
    def color(self):
        return Channel.colors[self.id % len(Channel.colors)]
        
    def call(self, t):
        out = self.signal(t)
        if self.AC:
            out -= self.signal.offset
            for sq in self.signal.square:
                d = 0.5 if len(sq) < 4 else sq[3]
                out -= sq[1]*d
        return out
    
    def __call__(self, t):
        out = self.call(t)
                
        if self.invert: out *= -1
        return out

class DiffChannel(Channel):
    def __init__(self, a, b):
        super().__init__()
        self.a = a
        self.b = b
        self.voltdiv = max(a.voltdiv, b.voltdiv)
        
    def call(self, t):
        x = self.a(t)
        y = self.b(t)
        return x - y

class Oscilloscope:
    ## basic noise due to ele grid
    def __init__(self, noise = Signal(spectrum = [[50,1e-6]], noise = 1e-7)):
        self.noise = Channel(noise)
        self.noise.set_id(0)
        self.noise.voltdiv = 3* max( [x[1] for x in noise.spectrum] + [x[1] for x in noise.square] )
        self.channels = {}#ch1.name:ch1}
        self.trig_channel = self.noise
    
        self.samples = Ring(5)
    
        self.hdiv = 12
        self.vdiv = 8
    
        self.secdiv = 5e-3
        #self.voltdiv = 2e-7
        
        self.divsamples = 150
    
        try:
            self.bkg = plt.imread(to_path("lab_meta/OSCI_unit_bkg.png"))
        except:
            self.bkg = None
        
        self.mode = 'TY'   ## for x-y set it a channel list , eg.[ch1, ch2]
        

    def add_channel(self, ch, trig = True):
        n = len(self.channels) + 1
        ch.set_id( n )
        
        self.channels[ch.name] = ch
        if trig: self.trig_channel = ch
        
        echo = 'added channel '+str(n)
        if trig: echo += ' and set it as trigger signal'
        print(echo)
        
    def find_Trig(self, t0 = None, dt = None):
        #print('looking for trig, side is ', edge)
        source = self.trig_channel
        val = source.trig
        edge = source.trig_edge

        if t0 is None:
            t0 = time.time()-T0
        if dt is None:
            dt = 10* self.secdiv / self.divsamples
            
        t = t0
        x = source(t) + self.noise(t)
        
        
        LIM = t+10    ## 10 sec max search
        while(t < LIM):
            t += dt
            y = source(t) + self.noise(t)
            if edge == Trig.Ascending and y > val and x < val:
                return t
            elif edge == Trig.Descending and y < val and x > val:
                return t
            
            x = y
        return -1

    def sample(self, t0 = None):
        trig_time = self.find_Trig(t0)
        if trig_time < 0:
            print('trig time not found')
            return
        t = np.linspace( trig_time - self.secdiv * self.hdiv/2 , 
                        trig_time + self.secdiv * self.hdiv/2,
                        self.hdiv * self.divsamples )
        
        noise = self.noise(t)
        
        scans = {}
        for chname in self.channels:
            ch = self.channels[chname]
            if ch.active:
                x = ch(t)
                scans[ch.name] = (x + noise) 
        self.samples.put( [t - trig_time,scans] )
        # plot and set axes limits
        
    def init_fig(self):
        fig, ax = self.new_fig()
        self.fig = fig; self.ax = ax
        self.set_fig(self.fig, self.ax)        
    
        
    def set_fig(self, fig , ax ):
        #self.ax.margins(0,0)
        fig.subplots_adjust(left=0, bottom=0,right=1,top=1,wspace=0,hspace=0)
        ax.tick_params(left=False,
                bottom=False,
                labelleft=False,
                labelbottom=False)
        ax.set_facecolor('black')
        
        if self.bkg is None:
            self.plotax = ax
        else:
            ax.imshow(self.bkg, extent = [0 , 1200 , 0, 542 ], zorder = 1)
        
            menuax = fig.add_axes( [ 180/1200 , 72/542 , (670-180)/1200 , (121-72)/542  ] , label = 'menu')
            menuax.set_facecolor('firebrick')
            menuax.tick_params(left=False,
						bottom=False,
						labelleft=False,
						labelbottom=False)       
        
            self.menuax = menuax
            self.plotax = fig.add_axes( [ 180/1200 , 121/542 , (670-180)/1200 , (445-121)/542  ] , label = 'plot')
        
            ax = self.menuax
            ax.cla()
            ax.set_facecolor('firebrick')
        
            draw = []
            if type(self.mode) is list:
                draw = self.mode                    
            else:
                for chname in self.channels:
                    ch = self.channels[chname]
                    if ch.active:
                        draw += [ch]
                        
            cnt = 0
            for ch in draw:
                x = cnt*0.25
                y = 0.6
                ax.text( x,y, ch.name, zorder = 2, fontsize = 11, color = 'black', fontweight = 'bold')
                text = '{:3g} V/d'.format(ch.voltdiv)
                ax.text( x+0.08 ,y, text,
                        zorder = 2, fontsize = 11,
                        backgroundcolor = 'black',
                        color = ch.color())
                if ch.invert:
                    ax.text( x+0.06, y, 'v', fontsize=12, zorder = 4, color = self.noise.color(), fontweight='bold' )
                    
                ax.text(x+0.05, 0.1, '({:2g},{:2g})'.format(ch.dh,ch.dv),
                        zorder = 2, fontsize = 9,
                        backgroundcolor = 'black',
                        color = 'white')
                cnt += 1 

            text2 = '|dT {:.0e} s/d'.format(self.secdiv)
            ax.text(0.75,0.4, text2, zorder = 3, fontsize = 10, color = 'black', fontweight = 'bold')

        
    def step(self, n = -1):
        if n < 0: n = self.samples.n
        for i in range(n):
            self.sample()
        
    def clear(self):
        self.samples = Ring(5)

    def animation(self, frame = 0):
        self.step(1)
    
        ax = self.plotax
        ax.cla()
    
        ax.set_facecolor('black')

        ax.tick_params(left=False,
                    bottom=False,
                    labelleft=False,
                    labelbottom=False)

        if type(self.mode) is list:
            ch = self.mode
            
            for a in self.samples:
                t, xs = a
                for chname in xs:
                    if chname == ch[0].name:
                        x = xs[chname]
                    if chname == ch[1].name:
                        y = xs[chname]
                
                ax.plot( x/ch[0].voltdiv, y/ch[1].voltdiv , color = ch[0].color(),zorder = 1)
                ax.scatter( x/ch[0].voltdiv, y/ch[1].voltdiv , color = ch[1].color(), s = 15,zorder = 2)

                
            ax.set_xlim( [ -self.hdiv / 2 , self.hdiv / 2] )
            ax.set_ylim( [ -self.vdiv / 2 , self.vdiv / 2] )
            
            ax.set_xticks( [ - self.hdiv / 2 + i for i in range(self.hdiv) ] )
            ax.set_yticks( [ - self.vdiv / 2 + i for i in range(self.vdiv) ] )           
        else:
            t_base = None      
            for a in self.samples:
                t, xs = a
                if t_base is None:
                    t_base = t
                for chname in xs:
                    x = xs[chname]
                    ch = self.channels[chname]
                    ax.plot(t + ch.dh, (x + ch.dv) / ch.voltdiv , color = ch.color() , lw = 1)
                    ax.scatter([t_base[0]],[ch.dv/ch.voltdiv], marker = '>', s = 550, color = ch.color(),zorder=2)
                
            trig_y = (self.trig_channel.trig + self.trig_channel.dv) / self.trig_channel.voltdiv
            #ax.axhline(trig_y, 0,1 , linestyle='-.', color = self.trig_channel.color())
            ax.axhline(trig_y, 0,1 , linestyle='-.', color = self.noise.color(), lw = 0.5)
            ax.scatter([t_base[0]],[trig_y], marker = '>', s = 250, color = self.noise.color(), zorder=3)
            
            mid = (t_base[-1]+t_base[0])/2
            ax.scatter([ self.trig_channel.dh + mid ],[4], marker = 'v', s = 250, color = self.noise.color())
            
                
            ax.set_xlim( [ t_base[0] , t_base[-1]] )
            ax.set_ylim( [ -self.vdiv / 2 , self.vdiv / 2] )
            
            ax.set_xticks( [ t_base[0] + i*self.secdiv for i in range(self.hdiv) ] )
            ax.set_yticks( [ - self.vdiv / 2 + i for i in range(self.vdiv) ] )
         

        
        ax.grid(True, color = self.noise.color(), alpha=0.8, linestyle=':')
        ax.axhline(0, 0,1, lw = 1, color = self.noise.color(), zorder = 2)
        ax.axvline(0, 0,1, lw=1, color = self.noise.color(), zorder = 2)
        
    def new_fig(self):    
        if self.bkg is None:
            fig, ax = plt.subplots(1, 1, figsize = (12, 9))   ## 1:1 aspect
        else:
            fig, ax = plt.subplots(1, 1, figsize = (12, 4.1))
        
        return fig,ax
            
    def show(self):
        fig, ax = self.new_fig()
        self.set_fig(fig,ax)
        self.animation()
        plt.show(block=False)
    
    def screen(self):
        fig, ax = plt.subplots(1, 1, figsize = (12, 9))
        self.plotax = ax
        self.animation()
        plt.show(block=False)
            
    
    def prompt(self):
        line = input()
        while(not line in ['exit','quit']):
            print(line)
            line = input()        
    
    def run(self):
        self.init_fig()
        anim = mani.FuncAnimation(self.fig, self.animation, interval = 250)
        plt.show(block = False)
        self.prompt()

if __name__=='__main__':
    osci = Oscilloscope()
    #osci.trig_edge=Trig.Descending
    #osci.trig = 3e-7
    
    osci.init_fig()
    osci.animation(0)
    plt.show(block=False)
        
    
    ch1 = Channel(Signal(spectrum = [[15,2],[30,3.5],[45,1.5]],noise = 1e-4))
    ch1.voltdiv = 2
    ch1.trig = 3
    
    ch2 = Channel( Signal(spectrum = [[50,2.5]], noise = 1e-4) )
    ch2.voltdiv = 1
    
    osci.add_channel( ch1 )
    osci.add_channel( ch2 , trig = False)
    osci.show()
    
    ch3 = Channel(Signal(offset = 0, square=[[50, 2, 3e-3, 0.5]], noise = 0.05))
    ch3.voltdiv = 0.5
    ch1.active = False
    ch2.active = False
    ch4 = Channel(Signal(offset = -1.5, square=[[50, 2.5, -2e-3, 0.5]], noise = 0.05))
    ch4.voltdiv = 0.5
    
    osci.add_channel(ch3, trig = False)
    osci.add_channel(ch4)
    
    osci.step(25)
    osci.show()
    
    ch3.invert = True
    ch3.AC = True
    ch4.AC = True
    
    osci.step(25)
    osci.show()
    #osci.run()
    
    
    ch1.active = True
    osci.trig_channel = ch1
    
    ch3.AC = False
    osci.step(25)
    osci.show()
    
    osci.mode = [ch1, ch3]
    osci.show()
    
    osci.mode = [ch4, ch1]
    osci.show()
    
    osci.mode = [ch4, ch3]
    osci.show()

    pass
