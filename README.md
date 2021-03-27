# VirElab
Virtual labs for labs in electrical measurements

```python
import numpy as np
import matplotlib.pyplot as plt
import osci
```


```python
## channel color scheme; yellow for no channel
osci.Channel.colors = ['yellow',  'deepskyblue',  'orangered', 'lawngreen',  'magenta', 'dodgerblue']
```


```python
o = osci.Oscilloscope()

## no channels connected - no signal shown
o.show()


## null signal to show base (cofigurable) background noise
null_sig = osci.Signal()
c1 = osci.Channel(null_sig)  ## noise is additional channel signal noise; 
                                                       ##ind of osci noise
o.add_channel(c1)
o.show()
```


    
![png](Demo/output_2_0.png)
    


    added channel 1 and set it as trigger signal



    
![png](Demo/output_2_2.png)
    



```python
o = osci.Oscilloscope()

#### harmonic signal
## osci.Signal is the basic class for representing signal
# spectrum is a list of [freq., Amp] pairs
# noise is gaussian std
s1 = osci.Signal(spectrum = [[50, 2.5]], noise = 0.05)

## osci.Channel represents input channel. it has extra attribs for displays and triggers
c1 = osci.Channel(s1)
c1.voltdiv = 1

o.add_channel(c1)

o.step()
o.show()
```

    added channel 1 and set it as trigger signal



    
![png](Demo/output_3_1.png)
    



```python
## trigger
c1.trig = 1.5
o.step()
o.show()
```


    
![png](Demo/output_4_0.png)
    



```python
## shift
c1.dh = 0.005
c1.dv = 1
o.step()
o.show()
```


    
![png](Demo/output_5_0.png)
    



```python
o = osci.Oscilloscope()

## two channels, TY mode
s1 = osci.Signal(spectrum = [[50, 2.5]], noise = 0.05)
s2 = osci.Signal(spectrum = [[100, 1.5], [20, 2]], noise = 0.1)

c1 = osci.Channel(s1)
c1.voltdiv = 1
c2 = osci.Channel(s2)
c2.voltdiv = 1

o.add_channel(c1)
o.add_channel(c2, main = False)  ## c1 remains triggering

o.step()
o.show()
```

    added channel 1 and set it as trigger signal
    added channel 2



    
![png](Demo/output_6_1.png)
    



```python
## take just the screen; 1:1 aspect
o.screen()
```


    
![png](Demo/output_7_0.png)
    



```python
o = osci.Oscilloscope()

## square signal in channel 2
s1 = osci.Signal(offset = 1, spectrum = [[50, 3]], noise = 0.05)
s2 = osci.Signal(square = [[125, 2, 0.7], [20, 2, 0.2]])

c1 = osci.Channel(s1)
c1.voltdiv = 1
c2 = osci.Channel(s2)
c2.voltdiv = 1

o.add_channel(c1)
o.add_channel(c2, main = False)

o.step()
o.show()


## only one sample instead of default 5
o.clear()  
c1.dh = 0.01
c2.dv = -1
o.show()
```

    added channel 1 and set it as trigger signal
    added channel 2



    
![png](Demo/output_8_1.png)
    



    
![png](Demo/output_8_2.png)
    



```python
o.trig_channel = c2   ## c2 is now triggering
c2.trig = 3
c2.dh = 0.005
o.step()
o.show()
```


    
![png](Demo/output_9_0.png)
    



```python
o = osci.Oscilloscope()

## x-y mode
s1 = osci.Signal(spectrum = [[50, 2.5]], noise = 0.05)
s2 = osci.Signal(spectrum = [[200, 1.5]], noise = 0.05)

c1 = osci.Channel(s1)
c1.voltdiv = 0.5
c2 = osci.Channel(s2)
c2.voltdiv = 0.5

o.add_channel(c1)
o.add_channel(c2, main = False)

o.mode = [c1,c2] ### set this to list of the x-y channels. any other shows all TYs

o.step()
o.show()
```

    added channel 1 and set it as trigger signal
    added channel 2



    
![png](Demo/output_10_1.png)
    



```python
s3 = osci.Signal(spectrum = [[150, 1.5]], noise = 0.01)
c3 = osci.Channel(s3)
c3.voltdiv = 0.5

o.add_channel(c3, main = False)
o.mode = [c1, c3]

o.step()
o.show()
```

    added channel 3



    
![png](Demo/output_11_1.png)
    



```python
## 3 channels
o.mode = 'for TY, this may be whatever-just not a list of [chX,chY]'
o.secdiv = 2e-3
c1.trig = 1
c1.dh = - 6 * o.secdiv   ## put trigger to origin
o.step()
o.show()
```


    
![png](Demo/output_12_0.png)
    



```python
o = osci.Oscilloscope()

## x-y with squares
s1 = osci.Signal(offset = -1, square = [[50, 2.5, 0.3]], noise = 0.05)
s2 = osci.Signal(square = [[200, 1.5,0.7], [25, -3, 0.2]], noise = 0.05)

c1 = osci.Channel(s1)
c1.voltdiv = 0.5
c2 = osci.Channel(s2)
c2.voltdiv = 0.5

o.add_channel(c1)
o.add_channel(c2, main = False)

o.mode = [c1,c2] ### set this to list of the x-y channels. any other shows all TYs

o.step()
o.show()
```

    added channel 1 and set it as trigger signal
    added channel 2



    
![png](Demo/output_13_1.png)
    



```python

```
