# toolkit =====================================================================

# get path of this script
path = os.path.dirname(os.path.abspath(__file__))

# calculation =================================================================

# check whether is an array/list or scalar
hasattr(N, "__len__")

# check whether is string
isinstance(test_string, str)

# drawing =====================================================================

# change default matplotlib fonts (window layouts)
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 13

# change default matplotlib fonts (Jupyterlab layouts)
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 32
plt.rcParams['legend.fontsize'] = 26

# draw a circle
ax.add_artist(plt.Circle((x,y),r,facecolor=RGB,edgecolor='none'))

# add a panel with no margins
ax = fig.add_axes([0,0,1,1])

# panel with no axis lines, ticks, or labels
ax.axis('off')

# remove tick labels
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

# set background colors
fig = plt.figure()
fig.patch.set_facecolor('black') # background that surrounds the panel 
ax.set_facecolor('k') # background in the panel

# change axis colors
ax.spines['bottom'].set_color('#dddddd') # axis lines
ax.spines['top'].set_color('#dddddd') 
ax.spines['right'].set_color('red')
ax.spines['left'].set_color('red')
ax.tick_params(axis='x',colors='red',which="both") # ticks 
ax.yaxis.label.set_color('red') # label
ax.xaxis.label.set_color('red')

# title color
ax.title.set_color('red') 

# change tick positions & labels
ticks = ax.get_xticks()
ax.set_yticks(ticks)
ax.set_yticklabels(labels)

# value to color
plt.get_cmap('GnBu_r')(value)

# make a gif with moviepy -----------------------------------------------------

from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

fig = plt.figure(figsize=(4,6))
ax = plt.gca()
dot = plt.scatter(.2,y0/1e3,fc='w',ec='k',s=60,zorder=49)

def make_frame(t1):
    dot.set_offsets(new_coord)
    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame,duration=duration)
animation.write_gif('fall.gif',fps=fps)

# interactive plot with sliders -----------------------------------------------
from matplotlib.widgets import Slider, Button

fig = plt.figure(figsize=(fig_x,fig_y))

# plot pannel
ax1 = plt.axes([fx_p1,fy_p,fw_p,fh_p])
l1,=ax1.plot([0,x1[0],x2[0]],[0,y1[0],y2[0]],color='C0')

# lamda slider
ax_s = plt.axes([fx_s,fy_s,fw_s,fh_s])
sld = Slider(ax=ax_s,label=r'$\lambda=\omega^2l/g$',
             valmin=0,valmax=20,valinit=lamda,orientation='horizontal')
sld.vline._linewidth = 0  # no initial line

# slider ticks
tick_s = [lamda_cr1,lamda_cr2,lamda_cr3,lamda_cr4]
ticklabel_s = [r'$\lambda_1$',r'$\lambda_2$',r'$\lambda_3$',r'$\lambda_4$']
sld.ax.xaxis.set_ticks(tick_s)
sld.ax.xaxis.set_ticklabels(ticklabel_s)

# The function to be called anytime a slider's value changes
def update(val,ax1=ax1,ax2=ax2,ax3=ax3,ax4=ax4):
    # new slider values
    lamda = sld.val 

    # upload plots
    l1.set_xdata([0,x1[0],x2[0]])
    l1.set_ydata([0,y1[0],y2[0]])
    fig.canvas.draw_idle()

# register the update function with each slider
sld.on_changed(update)

