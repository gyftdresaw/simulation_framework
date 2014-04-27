
# consolidating various plotting routines to clarify
# simulation vs analysis script components
#
# plotting routines are tailored to displaying development on lattice

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle,Circle
import saving

class PlotHelper:
    # simulation data consists of:
    #  - timepoints
    #  - data per species per cell
    #  - cell positions
    #  - simulation class
    def __init__(self,t,cdata,centers,sim):
        self.t = t
        self.cdata = cdata 
        self.centers = centers
        self.sim = sim
        self.NCells = centers.shape[0]

    # we'll have it accept these arguments in a dictionary to make save/load easier
    @classmethod
    def from_pickle(cls,pickle_data):
        t = pickle_data['t']
        cdata = pickle_data['cdata']
        centers = pickle_data['centers']
        sim = pickle_data['sim']
        return cls(t,cdata,centers,sim)

    # straight from pickle file
    @classmethod
    def from_file(cls,filename):
        return cls.from_pickle(saving.load(filename))

    ## TYPICAL PLOTS ##
    # plot 2D picture
    def plot_2D(self,species_name,time):
        tindex = np.abs(self.t-time).argmin()
        astatus = np.zeros(self.NCells)
        for i in xrange(self.NCells):
            im_id = self.sim.get_IM(i)
            astatus[i] = self.cdata[i][tindex,self.sim.IMs[im_id].get_node_id(species_name)]

        # plot using drawn rectangles
        fig = plt.figure()
        
        gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        # set aspect ratios
        ax1.set(aspect='equal')
        ax2.set(aspect=15)

        # set axis limits
        x,y = self.centers[:,0],self.centers[:,1]
        ax1.axis([min(x)-1., max(x)+1., min(y)-1., max(y)+1.])

        # data normalizer for color
        data_norm = mpl.colors.Normalize(vmin=min(astatus),vmax=max(astatus))

        dx = [0.8]*self.NCells # cells are roughly 1 unit apart
        for x,y,c,h in zip(self.centers[:,0],self.centers[:,1],plt.cm.jet(data_norm(astatus)),dx):
            # ax1.add_artist(Rectangle(xy=(x-h/2.,y-h/2.),facecolor=c,
            #                          width=h,height=h,edgecolor='black'))
            # circles look nicer, radius assumed to be 0.5
            ax1.add_artist(Circle(xy=(x,y),radius = 0.5,facecolor=c))

        # now plot colorbar
        mpl.colorbar.ColorbarBase(ax2,cmap=plt.cm.jet,norm=data_norm)
        
        # plt.scatter(self.centers[:,0],self.centers[:,1],
        #             c=astatus,s=100,marker='s',edgecolors='none')
        plt.show()

    # we need some better functions for plotting 2D

    def plot_horizontal(self,species_name,time):
        tindex = np.abs(self.t-time).argmin()
        astatus = np.zeros(self.NCells)
        for i in xrange(self.NCells):
            im_id = self.sim.get_IM(i)
            astatus[i] = self.cdata[i][tindex,self.sim.IMs[im_id].get_node_id(species_name)]

        plt.scatter(self.centers[:,1],astatus,s=50,edgecolors='none')
        plt.title(species_name,fontsize=24)
        plt.show()

    # plot species at various times
    def plot_species(self,species_name,times,cids=None):
        if cids == None:
            cids = range(self.NCells)
        NCells = len(cids)
        tindices = [np.abs(self.t-v).argmin() for v in times]
        astatus = np.zeros((NCells,len(tindices)))
        for i in xrange(len(cids)):
            for j in xrange(len(tindices)):
                im_id = self.sim.get_IM(i)
                astatus[i,j] = self.cdata[cids[i]][tindices[j],self.sim.IMs[im_id].get_node_id(species_name)]
        '''
        plt.figure()
        plt.scatter(x_coord,np.zeros(NCells),
                    c=astatus,s=50,marker='s',edgecolors='none')
        plt.show()
        '''

        # plot along x axis
        plt.figure()
        colors = cm.Dark2(np.linspace(0, 1, len(tindices)))
        x_coord = np.linspace(1,NCells,NCells)
        for j in xrange(len(tindices)):
            plt.scatter(x_coord,astatus[:,j],color=colors[j],s=50)
        plt.legend(['%.1f'% time for time in times],loc='best')
        plt.title(species_name,fontsize=24)
        plt.show()

    def plot_all_species(self,times):
        # only plot species that everybody has
        IM_nodes = [set(IM.node_names.keys()) for IM in self.sim.IMs]
        to_plot = reduce(lambda x,y: x & y,IM_nodes)
        for n in to_plot:
            self.plot_species(n,times)

    def plot_cell(self,cid,times):
        tstart,tend = [np.abs(self.t-v).argmin() for v in times]
        plt.figure()
        plt.plot(self.t[tstart:tend],self.cdata[cid][tstart:tend,:])
        plt.legend(self.sim.IMs[self.sim.get_IM(cid)].get_node_names())
        plt.show()

    # differentiate cell types
    def get_r8s(self,time):
        tindex = np.abs(self.t-time).argmin()
        return [i for i in xrange(self.NCells) if self.cdata[i][tindex,self.sim.IMs[self.sim.get_IM(i)].get_node_id('u')] > 0.010]

    def get_r25s(self,time):
        r25_indices = [ [r8-1,r8+1] for r8 in self.get_r8s(time)]
        # not the most readable just flattens the list
        return [c for pair in r25_indices for c in pair if c >= 0 and c < self.NCells]

    def get_progenitors(self,time):
        r8s = self.get_r8s(time)
        r25s = self.get_r25s(time)
        return [c for c in xrange(self.NCells) if not c in r8s and not c in r25s]

    ## plot based on these categorizations

    def plot_cell_species(self,species_name,time):
        NCells = self.NCells
        tindex = np.abs(self.t-time).argmin()

        astatus = np.zeros((NCells,1))
        for i in xrange(NCells):
            im_id = self.sim.get_IM(i)
            astatus[i] = self.cdata[i][tindex,self.sim.IMs[im_id].get_node_id(species_name)]

        r8s = self.get_r8s(time)
        r25s = self.get_r25s(time)
        progenitors = self.get_progenitors(time)
        cell_types = [r8s,r25s,progenitors]

        # plot along x axis
        x_coord = np.linspace(1,NCells,NCells)
        plt.figure()
        gs = gridspec.GridSpec(2,1,height_ratios=[19,1])
        # first plot values 
        ax1 = plt.subplot(gs[0])
        ax1.xaxis.grid(True,linewidth=0.1,linestyle='-',color='0.4');
        colors = cm.Set1(np.linspace(0, 1, 10))
        for j in xrange(len(cell_types)):
            plt.scatter(x_coord[cell_types[j]],astatus[cell_types[j]],color=colors[j],s=50)
        plt.legend(['R8','R2/5','MP'],loc='best')
        for j in xrange(len(cell_types)):
            plt.plot(x_coord[cell_types[j]],astatus[cell_types[j]],'--',color=colors[j])
        plt.title(species_name,fontsize=24)
        plt.xlabel('Cell',labelpad=0)

        # now plot stripes
        ax2 = plt.subplot(gs[1],sharex=ax1)
        xcorners = np.tile(np.linspace(0,NCells,NCells+1)+0.5,(2,1))
        ycorners = np.tile(np.array([[0],[1]]),(1,NCells+1))
        plt.pcolor(xcorners,ycorners,astatus.T,edgecolors='k')
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        plt.show()

    def plot_all_cell_species(self,time):
        # only plot species that everybody has
        IM_nodes = [set(IM.node_names.keys()) for IM in self.sim.IMs]
        to_plot = reduce(lambda x,y: x & y,IM_nodes)
        for n in to_plot:
            self.plot_cell_species(n,time)
