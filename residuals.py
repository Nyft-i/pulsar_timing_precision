
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
dph = 0
df0 = -1e-7
df1 = +4.2e-15
dfd = -3.23e-8
taud = 50*24*3600

def annotate_dim(ax,xyfrom,xyto,text=None):
  if text is None:
    text = str(np.sqrt( (xyfrom[0]-xyto[0])**2 + (xyfrom[1]-xyto[1])**2 ))

  ax.annotate("",xyfrom,xyto,arrowprops=dict(arrowstyle='<->'))
  ax.text((xyto[0]+xyfrom[0])/2,(xyto[1]+xyfrom[1])/2,text,fontsize=9)


# define a line decreasing with gradient m
def phi(x):
  x = x * 24 * 3600
 
  if x < x0:
    return 0
  else:
    return df0*(x-xg)+ 0.5 * df1 * (x-xg)**2+dfd * taud * (1-np.exp(-(x-xg)/taud))
    
def rec_line(x):
  x = x * 24 * 3600
  if x < x0:
    return 0
  else:
    return (df0 + df1*(x-xg)) * f0
    
  
fig = plt.figure()
ax = fig.add_subplot(111)
# change shape of the figure
fig.set_size_inches(5, 3)

#time
t_before = np.linspace(59900, 60000, 100)
t = np.linspace(60000, 60500, 10000)


nobs= 20
toas = np.linspace(59910, 60490, nobs)
# reascale to 59900-65000
errs = 0.05 *2
noise = np.random.normal(0,0.05,nobs)



# pre-glitch
ax.plot(t_before, np.zeros(len(t_before)), color='k', label=f"$\\phi(t)$")

# post-glitch
nut = np.zeros(len(t))
for enu, x in enumerate(t):
  nut[enu] = phi(x)
  
recovery_line = np.zeros(len(t))
for enu, x in enumerate(t):
  recovery_line[enu] = rec_line(x)
ax.plot(t, nut, color='k')

# datapoins
data_points = np.zeros(len(toas))
for enu, x in enumerate(toas):
  data_points[enu] = phi(x)
plt.errorbar(toas, data_points+noise, yerr=errs, linestyle='',marker=",",c="magenta")
plt.scatter(toas, data_points+noise, marker="s",s=7,c="magenta", label="TOA residuals")

#glitc epoch
ax.axvline(x=60000, color='pink', linestyle='--', label=r"$t_g=60000$")

# post-recovery line
#ax.plot(t, recovery_line, color='gray', linestyle=':')

# x/y labels
ax.set_xlabel('Time (MJD)')
# latex formatting
ax.set_ylabel(f'Residuals')
ax.set_xlim(59900, 60500)

plt.legend(loc='upper right',
          fancybox=True)
plt.savefig('residuals.png', dpi=400, bbox_inches='tight')
plt.show()