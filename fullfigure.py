#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 22:15:36 2021

@author: mkravche

Wrapper figure class to call the animations in ways "just my way". Everyone seems to have their own collection of matplotlib scripts that
are ungeneralizable. If you want it, sure come and have it. I left my treasure in on place.
Ironically this is like the fifth wrapper function that I have made for just my way. It seems to change every time. This is disposable code.
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['toolbar']='toolbar2'

class full_figure:
    """
    This is my custom figure class that i import for making animations easier.
    """
    def __init__(self,nrows=2,ncols=3,super_title=None,  fullscreen=True,
                                                         subplot_kw=dict(frameon=False,
                                                                         aspect='equal',
                                                                         animated=True,
                                                                         xticks=[],
                                                                         yticks=[]),
                                                         
                                                         **fig_kw,):
        self.nrows=nrows
        self.ncols=ncols
        self.fig,self.axes = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw=subplot_kw,**fig_kw)
        self.manager = plt.get_current_fig_manager()
        if fullscreen:
            mpl.rcParams['toolbar'] = 'None'
            self.fig.canvas.window().statusBar().setVisible(False) 
            self.manager.full_screen_toggle()
        else:
            mpl.rcParams['toolbar']='toolbar2'
    
        if super_title:
            self.fig.suptitle(super_title)
        plt.rcParams.update({'font.size': 15})
        self.fig.set_figheight(9)
        self.fig.set_figwidth(16)
        
        # self._untightened=True
        # self._unset_colorbar = True
    
        
        
        
        
    def _set_aspects_equal(self):
        for ax in self.axes.reshape(-1):
            ax.set_aspect('equal')
            
    def draw_all_artists(self):
        # Refresh all the artists. The animated = True leaves the work of what to update to the programmer (me)
        for ax in self.axes.reshape(-1):
            self.fig.draw_artist(ax)
            
    #https://matplotlib.org/stable/api/_as_gen/matplotlib.artist.Artist.set_animated.html#matplotlib.artist.Artist.set_animated
    #https://matplotlib.org/stable/tutorials/advanced/blitting.html
    #https://www.geeksforgeeks.org/matplotlib-axes-axes-add_artist-in-python/#:~:text=The%20Axes.add_artist%20%28%29%20function%20in%20axes%20module%20of,.%20Return%20value%3A%20This%20method%20returns%20the%20artist.
