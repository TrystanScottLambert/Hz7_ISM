#### Testing interactive plots ####

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])

coords = []
def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    plt.scatter(event.xdata, event.ydata,color='r')
    fig.canvas.draw()
    coords.append((event.xdata,event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()