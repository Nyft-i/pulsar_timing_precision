
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cbook import get_sample_data
from matplotlib.offsetbox import (AnnotationBbox, DrawingArea, OffsetImage,
                                  TextArea)
from matplotlib.patches import Circle

f0 = 4
f1 = -1.8e-12
x0 = 60000 * 24 * 3600
xg = 60000 * 24 * 3600
df0 = 1e-7
df1 = -4.2e-15
dfd = 3.23e-8
taud = 50*24*3600

def annotate_dim(ax,xyfrom,xyto,text=None):
  if text is None:
    text = str(np.sqrt( (xyfrom[0]-xyto[0])**2 + (xyfrom[1]-xyto[1])**2 ))

  ax.annotate("",xyfrom,xyto,arrowprops=dict(arrowstyle='<->'))
  ax.text((xyto[0]+xyfrom[0])/2,(xyto[1]+xyfrom[1])/2,text,fontsize=9)


# define a line decreasing with gradient m
def nu(x):
  x = x * 24 * 3600
 
  if x < x0:
    return 0
  else:
    return (#f0 + 
            #f1*(x-x0) + 
            df0 + 
            ((df1)*((x-xg))) +
            dfd * np.exp(-(x-xg)/taud)
            )
    
def rec_line(x):
  x = x * 24 * 3600
  if x < x0:
    return 0
  else:
    return df0 + df1*(x-xg)
    
  
fig = plt.figure()
ax = fig.add_subplot(111)
# change shape of the figure
fig.set_size_inches(7, 3)

#time
t_before = np.linspace(59950, 60000, 100)
t = np.linspace(60000, 60200, 10000)


# pre-glitch
ax.plot(t_before, np.zeros(len(t_before)), color='navy', label=f"$\\nu(t)$, pre-glitch model subtracted.")

# post-glitch
nut = np.zeros(len(t))
for enu, x in enumerate(t):
  nut[enu] = nu(x)
  
recovery_line = np.zeros(len(t))
for enu, x in enumerate(t):
  recovery_line[enu] = rec_line(x)
ax.plot(t, nut, color='navy')

#glitc epoch
ax.axvline(x=60000, color='pink', linestyle='--')

# post-recovery line
ax.plot(t, recovery_line, color='gray', linestyle=':')

# x/y labels
ax.set_xlabel('Time (MJD)')
# latex formatting
ax.set_ylabel(f'$\\nu(t)$ $\minus$ $\\nu_{{\\text{{model}}}}(t)$')

plt.legend(loc='lower right',
          ncol=3, fancybox=True)
plt.savefig('glitch_structure.png', dpi=400, bbox_inches='tight')
plt.show()