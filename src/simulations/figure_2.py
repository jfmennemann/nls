import matplotlib.pyplot as plt

from PyQt5 import QtWidgets 


import numpy as np


color_gridlines_major = '#666666'
color_gridlines_minor = '#999999'


class Figure2(object):
    
    
    def __init__(self, x, density_min, density_max, V_min, V_max, V, u_ref):
        
        Jx = x.size
        
        index_center_x = Jx // 2
        
        assert(np.abs(x[index_center_x])<1e-15)
        
        
        
        
        self.fig_name = "figure_1"
        
        width  = 6
        height = 4

        self.fig = plt.figure(self.fig_name, figsize=(width, height), facecolor="white", constrained_layout=False)
        
        window = self.fig.canvas.window()
        
        window.findChild(QtWidgets.QToolBar).setVisible(False)
        window.statusBar().setVisible(False)
        
        
        
        self.gridspec = self.fig.add_gridspec(ncols=1, nrows=3, left=0.125, bottom=0.125, right=0.9, top=0.98, wspace=0.40, hspace=0.65, width_ratios=[1], height_ratios=[1, 1, 1])
        
        
        
        #==========================================================================================
        ax_00 = self.fig.add_subplot(self.gridspec[0, 0])
        ax_10 = self.fig.add_subplot(self.gridspec[1, 0])
        ax_20 = self.fig.add_subplot(self.gridspec[2, 0])
        #==========================================================================================
        
        #==========================================================================================
        self.line_u_abs_squared, = ax_00.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='-', color='k',         label='u')
        
        if u_ref is not None:
        
            self.line_u_ref_abs_squared, = ax_00.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='--', color='tab:green', label='u_ref')
        
        if V is not None:
            
            ax_V = ax_00.twinx()
            
            self.line_V, = ax_V.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='-', color='red', label='V')
            
            ax_V.set_ylim(V_min, V_max)
            
            ax_V.set_ylabel('V')
            
            
        ax_00.set_ylim(density_min, density_max)
        
        ax_00.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_00.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_00.set_xlabel('x')
        ax_00.set_ylabel('|u|^2')
        
        # ax_00.legend(loc='upper right', fancybox=False, ncol=3)
        #==========================================================================================
        
        #==========================================================================================
        self.line_u_real, = ax_10.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='-', color='k', label='u')
        
        if u_ref is not None:
            self.line_u_ref_real, = ax_10.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='--', color='tab:green', label='u_ref')
        
        tmp = max(abs(density_min), abs(density_max))
        
        ax_10.set_ylim(-tmp, tmp)
        
        ax_10.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_10.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_10.set_xlabel('x')
        ax_10.set_ylabel('Re u')
        
        # ax_10.legend(loc='upper right', fancybox=False, ncol=2)
        #==========================================================================================
        
        #==========================================================================================
        self.line_u_imag, = ax_20.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='-', color='k', label='u')
        
        if u_ref is not None:
            self.line_u_ref_imag, = ax_20.plot(x, np.zeros_like(x), linewidth=1.0, linestyle='--', color='tab:green', label='u_ref')
        
        tmp = max(abs(density_min), abs(density_max))
        
        ax_20.set_ylim(-tmp, tmp)
        
        ax_20.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_20.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_20.set_xlabel('x')
        ax_20.set_ylabel('Im u')
        
        # ax_20.legend(loc='upper right', fancybox=False, ncol=2)
        #==========================================================================================
           
        plt.ion()
        
        plt.draw()
        plt.pause(0.001)
        
    
    def update_u(self, u_complete, u_ref=None):
        
        self.line_u_abs_squared.set_ydata(np.abs(u_complete)**2)
        
        self.line_u_real.set_ydata(np.imag(u_complete))
        self.line_u_imag.set_ydata(np.real(u_complete))
        
        if u_ref is not None:
            
            self.line_u_ref_abs_squared.set_ydata(np.abs(u_ref)**2)
            
            self.line_u_ref_real.set_ydata(np.imag(u_ref))
            self.line_u_ref_imag.set_ydata(np.real(u_ref))
    
    
    def update_V(self, V):
        
        self.line_V.set_ydata(V)
        
    
    def redraw(self):
        
        plt.figure(self.fig_name)
        
        plt.draw()
        
        self.fig.canvas.start_event_loop(0.001)
        
    