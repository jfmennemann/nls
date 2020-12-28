"""      
self.fig = plt.figure(self.fig_name, facecolor="white", constrained_layout=False)



window = self.fig.canvas.window()

window.findChild(QtWidgets.QToolBar).setVisible(False)
window.statusBar().setVisible(False)
"""

"""
if screen_size == 'small':
    
    resolution = '1920x1080'
    
elif screen_size == 'large':

    resolution = '2560x1440'


if resolution == '2560x1440':
    
    n_pixels_x = 1000
    n_pixels_y = 1000
    
    pos_x = 2560 - n_pixels_x
    pos_y = 0
    
    plt.rcParams.update({'font.size': 10})
    

if resolution == '1920x1080':
    
    n_pixels_x = 800
    n_pixels_y = 600
    
    pos_x = 1920 - n_pixels_x
    pos_y = 0
    
    plt.rcParams.update({'font.size': 6})


window.setGeometry(pos_x, pos_y, n_pixels_x, n_pixels_y)
"""




