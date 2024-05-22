# change default matplotlib fonts
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 13

# change default matplotlib fonts (Jupyter notebook)
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 32
plt.rcParams['legend.fontsize'] = 26

# get path of this script
path = os.path.dirname(os.path.abspath(__file__))

# draw a circle
ax = plt.gca()
ax.add_artist(plt.Circle((x,y),r,facecolor=RGB,edgecolor='none'))


# add a panel with no border
ax = fig.add_axes([0,0,1,1])


# panel with no axis lines, ticks, or labels
ax.axis('off')


# set background colors
fig = plt.figure()
fig.patch.set_facecolor('black') # background that surrounds the panel 
ax.set_facecolor('k') # background in the panel


# change axis colors
ax.spines['bottom'].set_color('#dddddd') # axis lines
ax.spines['top'].set_color('#dddddd') 
ax.spines['right'].set_color('red')
ax.spines['left'].set_color('red')

ax.tick_params(axis='x',colors='red') # ticks, which="both" changes major&minor
ax.tick_params(axis='y',colors='red')

ax.yaxis.label.set_color('red') # label
ax.xaxis.label.set_color('red')

ax.title.set_color('red') # title

# remove tick labels
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])


# change tick positions & labels
ticks = ax.get_xticks()
ax.set_yticks(ticks)
ax.set_yticklabels(labels)


# value to color
plasma = plt.get_cmap('GnBu_r')
plasma(value)


# check whether is an array/list or scalar
hasattr(N, "__len__")


# check whether is string
isinstance(test_string, str)

# make a movie with moviepy ...................................................

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


