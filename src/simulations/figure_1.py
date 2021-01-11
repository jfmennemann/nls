import matplotlib.pyplot as plt

from PyQt5 import QtWidgets 


import numpy as np


color_gridlines_major = '#666666'
color_gridlines_minor = '#999999'


class Figure1(object):
    
    
    def __init__(self, x, density_max, V_max, V, u_ref):
        
        Jx = x.size
        
        # integer division
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
            
            ax_V.set_ylim(0, 1.2*abs(V_max))
            
            ax_V.set_ylabel('V')
            
            
        ax_00.set_ylim(-0.1*density_max, 1.1*density_max)
        
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
        
        ax_10.set_ylim(-1.2*abs(density_max), 1.2*abs(density_max))
        
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
        
        ax_20.set_ylim(-1.2*abs(density_max), 1.2*abs(density_max))
        
        ax_20.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_20.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_20.set_xlabel('x')
        ax_20.set_ylabel('Im u')
        
        # ax_20.legend(loc='upper right', fancybox=False, ncol=2)
        #==========================================================================================
        
        #==========================================================================================
        """
        self.line_defect_of_mass, = ax_30.semilogy([], [], linewidth=1.0, linestyle='-', color='k')
        
        ax_30.set_xlim(0, T)
        # ax_30.set_ylim(-0.1, 2.1)
        ax_30.set_ylim(1e-15, 1)
        
        ax_30.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_30.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_30.set_xlabel('t')
        ax_30.set_ylabel('defect of mass')
        """
        #==========================================================================================
        
        #==========================================================================================
        """
        self.line_rel_error, = ax_40.semilogy([], [], linewidth=1.0, linestyle='-', color='k')
        
        ax_40.set_xlim(0, T)
        # ax_30.set_ylim(-0.1, 2.1)
        ax_40.set_ylim(1e-6, 1e1)
        
        ax_40.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5)
        ax_40.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5)
        
        ax_40.set_xlabel('t')
        ax_40.set_ylabel('rel_error')
        """
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
        
    
    """
    def update_rel_error(self, rel_error_of_times_analysis, times_analysis, nr_times_analysis):
        
        self.line_rel_error.set_xdata(times_analysis[0:nr_times_analysis])
        self.line_rel_error.set_ydata(rel_error_of_times_analysis[0:nr_times_analysis])
    """
    
    """ 
    def update_defect_of_mass(self, defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis):
        
        self.line_defect_of_mass.set_xdata(times_analysis[0:nr_times_analysis])
        self.line_defect_of_mass.set_ydata(defect_of_mass_of_times_analysis[0:nr_times_analysis])
    """
    
    
    def redraw(self):
        
        plt.figure(self.fig_name)
        
        plt.draw()
        
        self.fig.canvas.start_event_loop(0.001)
        
        
        
        
    
    
    
        
    
    