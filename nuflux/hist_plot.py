import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.colors as colors
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm

from . import fluxMC

def histogram1D(plotname, DATA, TMIN, TMAX,  XLABEL, TITLE, nbins):
	
	fsize = 11
	
	x1 = DATA[0]
	w1 = DATA[1]
	I1 = DATA[2]
	
	rc('text', usetex=True)
	params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
					'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(params)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	hist1 = np.histogram(x1, weights=w1, bins=nbins, density = False, range = (TMIN,TMAX) )

	ans0 = hist1[1][:nbins]
	ans1 = hist1[0]#/(ans0[1]-ans0[0])

	# comb1 = (ans1 + 0*8*case2[2]/case1[2] * ans2)
	# comb2 = (ans3*0 + 8*case4[2]/case3[2] * ans4)

	comb1 = ans1

	# comb1 = comb1/np.sum(comb1) #* I1  #/(ans0[1]-ans0[0])

	ax.bar(ans0,comb1, ans0[1]-ans0[0], label=r"PDF",\
			ec=None, fc='indigo', alpha=0.4, align='edge', lw = 0.0)	

	ax.step(np.append(ans0,10e10), np.append(comb1, 0.0), where='post',
				c='indigo', lw = 2.0)

	ax.set_title(TITLE, fontsize=fsize)

	plt.legend(loc="upper left", frameon=False, fontsize=fsize)
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(r"PDF",fontsize=fsize)

	ax.set_xlim(TMIN,TMAX)
	ax.set_ylim(0,ax.get_ylim()[1]*1.5)

	# plt.show()
	plt.savefig(plotname)
	plt.close()

def histogram2D(plotname, DATACOHX, DATACOHY, XMIN, XMAX, YMIN, YMAX,  XLABEL,  YLABEL, TITLE, NBINS):
	
	fsize = 9
	
	x1 = DATACOHX[0]
	y1 = DATACOHY[0]
	w1 = DATACOHY[1]
	I1 = DATACOHX[2]

	rc('text', usetex=True)
	params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
					'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(params)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	bar = ax.hist2d(x1,y1, bins=NBINS, weights=w1, range=[[XMIN,XMAX],[YMIN,YMAX]],cmap="Blues",normed=True)

	ax.set_title(TITLE, fontsize=fsize)
	cbar_R = fig.colorbar(bar[3],ax=ax)
	cbar_R.ax.set_ylabel(r'a.u.', rotation=90)

	plt.legend(loc="upper left", frameon=False, fontsize=fsize)
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	# ax.set_xlim(XMIN,XMAX)
	# ax.set_ylim(YMIN,YMAX)

	# plt.show()
	plt.savefig(plotname)
	plt.close()

def scatter2D(plotname, DATACOHX, DATACOHY, XMIN, XMAX, YMIN, YMAX,  XLABEL,  YLABEL, TITLE, NBINS, log_plot_x=False,log_plot_y=False):
	
	fsize = 9
	
	x1 = DATACOHX[0]
	y1 = DATACOHY[0]
	w1 = DATACOHY[1]
	I1 = DATACOHX[2]

	rc('text', usetex=True)
	params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
					'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(params)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	bar = ax.scatter(x1,y1, color='purple', s=0.2, edgecolor=None,alpha=0.6)

	ax.set_title(TITLE, fontsize=fsize)

	plt.legend(loc="upper left", frameon=False, fontsize=fsize)
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	if log_plot_x:
		ax.set_xscale("log")
	if log_plot_y:
		ax.set_yscale("log")
	ax.set_xlim(np.min([x1]),np.max([x1]))
	ax.set_ylim(np.min([y1]),np.max([y1]))


	plt.savefig(plotname)
	plt.close()

def histogram2DLOG(plotname, DATACOHX, DATACOHY, XMIN, XMAX, YMIN, YMAX,  XLABEL,  YLABEL, TITLE, NBINS):
	
	fsize = 11
	
	x1 = DATACOHX[0]
	y1 = DATACOHY[0]
	w1 = DATACOHY[1]
	I1 = DATACOHX[2]

	rc('text', usetex=True)
	params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
					'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(params)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)
			
	x1 = np.log10(x1)
	y1 = np.log10(y1)


	bar = ax.hist2d(x1,y1, bins=NBINS, weights=w1, range=[[np.log10(XMIN),np.log10(XMAX)],[np.log10(YMIN),np.log10(YMAX)]],cmap="Blues",normed=True)
	# hist[1][:nbins] = 10**hist[1][:nbins]

	ax.set_title(TITLE, fontsize=fsize)
	cbar_R = fig.colorbar(bar[3],ax=ax)
	cbar_R.ax.set_ylabel(r'a.u.', rotation=90)

	plt.legend(loc="upper left", frameon=False, fontsize=fsize)
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	# ax.set_xlim(XMIN,XMAX)
	# ax.set_ylim(YMIN,YMAX)

	# ax.set_xscale("log")
	# ax.set_yscale("log")
	# plt.show()
	plt.savefig(plotname)
	plt.close()


def data_plot(plotname, X, BINW, MODEL, DATA, ERRORLOW, ERRORUP, XMIN, XMAX, XLABEL,  YLABEL, TITLE):
	
	fsize = 11
	
	rc('text', usetex=True)
	params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
					'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(params)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	ax.bar(X, MODEL, BINW, align='center',\
				facecolor='lightblue', lw = 0.0)
	ax.step(np.append(X-BINW/2.0, X[-1]+BINW[-1]/2.0), np.append(MODEL,0.0), where='post',\
				color='dodgerblue', lw = 1.0)

	ax.errorbar(X, DATA, yerr= np.array([ERRORLOW,ERRORUP]), xerr = BINW/2.0, \
													marker="o", markeredgewidth=0.5, capsize=2.0,markerfacecolor="white",\
													markeredgecolor="black",ms=3, color='black', lw = 0.0, elinewidth=1.0, zorder=100)



	ax.set_title(TITLE, fontsize=fsize)

	plt.legend(loc="upper left", frameon=False, fontsize=fsize)
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	# ax.set_xlim(XMIN,XMAX)
	# ax.set_ylim(YMIN,YMAX)

	# plt.show()
	plt.savefig(plotname)
	plt.close()


def cut(true, noZvtx=False):
	##########################
	def cut_in_ptoverpz_min(Zvt):
	    return 0.1/(600-400)*(Zvt-400.0)+0.1
	def cut_in_ptoverpz_max(Zvt):
	    return (0.41-0.128)/(600.-200.)*(Zvt-200.)+0.128
	def cut_in_pi0Energy_min(Zvt):
	    return (0.2-1.0)/(600.-200.)*(Zvt-200)+1.0
	cuts={}

	cuts['total_photon_energy_min'] = 0.650
	cuts['photon1_energy_min'] = 0.1
	cuts['photon2_energy_min'] = 0.1
	cuts['photon1_energy_max'] = 2
	cuts['photon2_energy_max'] = 2
	cuts['ratio_egamma_min'] = 0.2
	cuts['photon_radius_max'] = 85 # cm
	cuts['photon_xy_min'] = 15 # cm
	cuts['interphoton_distance_min'] = 30.0 # cm
	cuts['angle_XY_2gammas_max'] = 150.0*np.pi/180.0
	cuts['etheta_min'] = 2.5 # GeV degrees
	cuts['ptoverpz_min'] = 0.1 # absolute minimum
	cuts['ptoverpz_functionZvtx_min'] =  cut_in_ptoverpz_min #
	cuts['ptoverpz_functionZvtx_max'] = cut_in_ptoverpz_max # 
	cuts['pi0Energy_functionZvtx_min'] = cut_in_pi0Energy_min #
	cuts['Zvtx_max'] = 490.0 # cm
	cuts['Zvtx_min'] = 316.0 # cm
	if noZvtx==True:
		cuts['Zvtx_max'] = 1e100 # cm
		cuts['Zvtx_min'] = 0.0 # cm
	cut_events, eff = analysis.analysis_cuts(true,cuts)

	flow = {'cut_events': cut_events,
		   'effs' : eff}

	return flow

def cut_flow(true,labels=None):
	##########################
	# CUT FLOW 0
	def cut_in_ptoverpz_min(Zvt):
	    return 0#0.1/(600-400)*(Zvt-400.0)+0.1
	def cut_in_ptoverpz_max(Zvt):
	    return 1e100#(0.41-0.128)/(600.-200.)*(Zvt-200.)+0.128
	def cut_in_pi0Energy_min(Zvt):
	    return 0#(0.2-1.0)/(600.-200.)*(Zvt-200)+1.0
	cuts={}
	cuts['total_photon_energy_min'] = 0.650*0
	cuts['photon1_energy_min'] = 0.1*0
	cuts['photon2_energy_min'] = 0.1*0
	cuts['photon1_energy_max'] = 2*1e100
	cuts['photon2_energy_max'] = 2*1e100
	cuts['photon_radius_max'] = 1e2 # cm
	cuts['photon_xy_min'] = 10 # cm
	cuts['interphoton_distance_min'] = 0*30.0 # cm
	cuts['ratio_egamma_min'] = 0*0.2
	cuts['angle_XY_2gammas_max'] = 1e100*150.0*np.pi/180.0
	cuts['etheta_min'] = 0*2.5 # GeV degrees
	cuts['ptoverpz_min'] = 0*0.1 # absolute minimum
	cuts['ptoverpz_functionZvtx_min'] =  cut_in_ptoverpz_min #
	cuts['ptoverpz_functionZvtx_max'] = cut_in_ptoverpz_max # 
	cuts['pi0Energy_functionZvtx_min'] = cut_in_pi0Energy_min # 
	cuts['Zvtx_max'] = 800 # cm
	cuts['Zvtx_min'] = 0.0 # cm
	cut_events_0, eff_0 = analysis.analysis_cuts(true,cuts)

	##########################
	# CUT FLOW 1
	cuts['total_photon_energy_min'] = 0.650
	cuts['photon1_energy_min'] = 0.1
	cuts['photon2_energy_min'] = 0.1
	cuts['photon1_energy_max'] = 2
	cuts['photon2_energy_max'] = 2
	cuts['ratio_egamma_min'] = 0.2
	cut_events_1, eff_1 = analysis.analysis_cuts(true,cuts)

	##########################
	# CUT FLOW 2
	cuts['photon_radius_max'] = 85 # cm
	cuts['photon_xy_min'] = 15 # cm
	cuts['interphoton_distance_min'] = 30.0 # cm
	cuts['angle_XY_2gammas_max'] = 150.0*np.pi/180.0
	cut_events_2, eff_2 = analysis.analysis_cuts(true,cuts)

	##########################
	# CUT FLOW 3
	def cut_in_ptoverpz_min(Zvt):
	    return 0.1/(600-400)*(Zvt-400.0)+0.1
	def cut_in_ptoverpz_max(Zvt):
	    return (0.41-0.128)/(600.-200.)*(Zvt-200.)+0.128
	def cut_in_pi0Energy_min(Zvt):
	    return (0.2-1.0)/(600.-200.)*(Zvt-200)+1.0
	cuts['etheta_min'] = 2.5 # GeV degrees
	cuts['ptoverpz_min'] = 0.1 # absolute minimum
	cuts['ptoverpz_functionZvtx_min'] =  cut_in_ptoverpz_min #
	cuts['ptoverpz_functionZvtx_max'] = cut_in_ptoverpz_max # 
	cuts['pi0Energy_functionZvtx_min'] = cut_in_pi0Energy_min # 
	cut_events_3, eff_3 = analysis.analysis_cuts(true,cuts)

	##########################
	# CUT FLOW 4
	cuts['Zvtx_max'] = 490.0 # cm
	cuts['Zvtx_min'] = 316.0 # cm
	cut_events_4, eff_4 = analysis.analysis_cuts(true,cuts)

	cut_events = [	cut_events_0,
					cut_events_1,
					cut_events_2,
					cut_events_3,
					cut_events_4]
	eff = [	eff_0,
			eff_1,
			eff_2,
			eff_3,
			eff_4]


	flow = {'cut_events': cut_events,
		   'effs' : eff,
		   'labels' : labels}

	return flow


def flow_plot(flow,
	true_var=None,reco_var='Pt_pion_rec',filename=None,
	units=1e3,xmin=0.0,xmax=300,nbins=30, 
	title='title', XLABEL=r'x', YLABEL=r"pdf",
	fsize=10, ymin=1e-6,ymax=100,
	plot_signal=False):


	##############################
	# FIGURE
	##############################
	rc('text', usetex=True)
	paramsrc={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
	                'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(paramsrc)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	if true_var==None:
		#############################
		# RECO_0
		hist1 = np.histogram(flow['cut_events'][0][reco_var]*units, 
								weights=flow['cut_events'][0]['w'], 
								bins=nbins, 
								density = False, 
								range = (xmin,xmax) )
		HISTNORM=np.sum(hist1[0])*(hist1[1][1]-[hist1[1][0]])
		#########
	else: 
		#########
		# TRUE 
		hist1 = np.histogram(flow['cut_events'][0][true_var]*units,
								 weights=flow['cut_events'][0]['w'], 
								 bins=nbins, 
								 density = False, 
								 range = (xmin,xmax) )
		x0 = hist1[1][:nbins]
		HISTNORM=np.sum(hist1[0])*(hist1[1][1]-[hist1[1][0]])
		y0 = hist1[0]/HISTNORM
		# ax.bar(x0,y0, x0[1]-x0[0],\
		#         ec="black", fc='black', alpha=0.4, align='edge', lw = 0.0)	
		ax.step(np.append(x0,10e10), np.append(y0,0.0), where='post',
		        c='black', dashes=(5,0), lw = 1.2,\
		label=r'T '+true_var.replace("_", "\_"))
		# ax.bar([0,0],[0,0],0,ec='black', fc=[0,0,0,0.4], align='edge', ls='--', lw = .5)

	###############
	# Reco+cuts
	def my_hist(cut_events,eff,color,label,alpha=0.4):

	    hist1 = np.histogram(cut_events[reco_var]*units, 
		    					weights=cut_events['w'], 
		    					bins=nbins, 
		    					density = False, 
		    					range = (xmin,xmax))
	    
	    x0 = hist1[1][:nbins]
	    y0 = hist1[0]/HISTNORM
	    ax.bar(x0,y0, x0[1]-x0[0],\
	            ec=color, fc=color, alpha=alpha, align='edge', lw = 0.0)
	    ax.step(np.append(x0,10e10), np.append(y0,0.0), 
	    	where='post', c='black', lw = 0.4)
	    # for correct legend
	    ax.bar([0,0],[0,0],0,ec='black', fc=np.array(color)*[1,1,1,alpha], align='edge', lw = 0.4,label=label)

	NTOT=np.shape(flow['effs'])[0]
	cmap = cm.BuPu(np.linspace(0.1,1,int(NTOT)))
	cmap = colors.ListedColormap(cmap)	

	###########################
	# PLOT THE FLOW HISTOGRAMS
	for i in range(NTOT):
		my_hist(flow['cut_events'][i],flow['effs'][i],cmap(float(i)/NTOT),alpha=0.9,label=flow['labels'][i])

	# set limits now
	ax.set_xlim(xmin,xmax)

	if ax.get_ylim()[1]>ymax:
		ax.set_ylim(ymin,ax.get_ylim()[1]*1000)
	else:
		ax.set_ylim(ymin,ymax)

	
	#############################
	# PLOT THE KOTO EVENTS 
	if plot_signal:
		pi0Z, pi0pt = np.genfromtxt("digitized/KOTO_data.dat",unpack=True)
		pi0pt*=1e3
		ax.vlines(pi0pt*units,ax.get_ylim()[0], ax.get_ylim()[1], linestyle=':', linewidth=0.8)
		for i in range(np.size(pi0pt)):
		    if i!=0:
		        if i==2:
		            ax.annotate(r'', xy=(pi0pt[i],ax.get_ylim()[0]), xytext=(pi0pt[i],0.3e-5), 
		                        arrowprops=dict(width=1,headwidth=4,headlength=5,facecolor='red',edgecolor='black',lw=0.4))
		        else:
		            ax.annotate(r'', xy=(pi0pt[i],ax.get_ylim()[0]), xytext=(pi0pt[i],0.3e-5), 
		                        arrowprops=dict(width=1,headwidth=4,headlength=5,facecolor='gold',edgecolor='black',lw=0.4))

		ax.vlines(250,ax.get_ylim()[0], ax.get_ylim()[1], linestyle='--', linewidth=0.7)
		ax.vlines(130,ax.get_ylim()[0], ax.get_ylim()[1], linestyle='--', linewidth=0.7)
		ax.annotate(r'', xy=(130+20,20e-2), xytext=(128,20e-2), 
		                        arrowprops=dict(arrowstyle="-|>",color='black',lw=0.6,mutation_scale=6))
		ax.annotate(r'', xy=(250-20,20e-2), xytext=(252,20e-2), 
		                        arrowprops=dict(arrowstyle="-|>",color='black',lw=0.6,mutation_scale=6))
		ax.annotate(r'\noindent\quad KOTO\\signal region', xy=(158,10e-2),color='black',fontsize=0.8*fsize)



	################
	# PLOTS STYLE
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	plt.legend(loc='upper right',frameon=False,fontsize=0.6*fsize)


	ax.set_title(title, fontsize=0.8*fsize)
	ax.set_yscale('log')

	if filename!=None:
		fig.savefig(filename)
	
	plt.close()

	return ax,fig



def hist_2D(true,
	var_x='Zvtx',var_y='Pt_pion_rec',filename=None,
	unitsx=1e3,	unitsy=10,
	xmin=1.5e3,xmax=6.5e3,ymin=0,ymax=300,nbins=20, 
	title='title', XLABEL=r'x', YLABEL=r"y",
	fsize=10,
	plot_signal=False):
	def cut_in_ptoverpz_min(Zvt):
	    return 0.1/(600-400)*(Zvt-400.0)+0.1
	def cut_in_ptoverpz_max(Zvt):
	    return (0.41-0.128)/(600.-200.)*(Zvt-200.)+0.128
	def cut_in_pi0Energy_min(Zvt):
	    return (0.2-1.0)/(600.-200.)*(Zvt-200)+1.0
	cuts={}
	cuts['total_photon_energy_min'] = 0.650
	cuts['photon1_energy_min'] = 0.1
	cuts['photon2_energy_min'] = 0.1
	cuts['photon1_energy_max'] = 2
	cuts['photon2_energy_max'] = 2
	cuts['photon_radius_max'] = 85 # cm
	cuts['photon_xy_min'] = 15 # cm
	cuts['interphoton_distance_min'] = 30.0 # cm
	cuts['ratio_egamma_min'] = 0.2
	cuts['angle_XY_2gammas_max'] = 150.0*np.pi/180.0
	cuts['etheta_min'] = 2.5 # GeV degrees
	cuts['ptoverpz_min'] = 0.1 # absolute minimum
	cuts['ptoverpz_functionZvtx_min'] =  cut_in_ptoverpz_min #
	cuts['ptoverpz_functionZvtx_max'] = cut_in_ptoverpz_max # 
	cuts['pi0Energy_functionZvtx_min'] = cut_in_pi0Energy_min # 
	cuts['Zvtx_max'] = 1e100*490.0 # cm
	cuts['Zvtx_min'] = 0*316.0 # cm
	cut_events, eff = analysis.analysis_cuts(true,cuts)

	y1 = cut_events[var_y]*unitsy
	x1 = cut_events[var_x]*unitsx
	w1 = cut_events['w']

	rc('text', usetex=True)
	rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
	                'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(rcparams)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	bar = ax.hist2d(x1,y1, bins=nbins, weights=w1/np.max(w1), 
	                range=[[xmin,xmax],[ymin,ymax]],cmap="Blues",normed=False,norm=LogNorm())

	ax.set_title(title, fontsize=fsize)
	cbar_R = fig.colorbar(bar[3],ax=ax)
	cbar_R.ax.set_ylabel(r'a.u.', rotation=90)


	if plot_signal:
		pts = np.array([[3160,130], [3950,130], [4900,160], [4900,250],[3160,250]])
		p = Polygon(pts, closed=True,facecolor=None,edgecolor='red', linestyle='--',fill=False)
		ax = plt.gca()
		ax.add_patch(p)

		pi0Z, pi0pt = np.genfromtxt("digitized/KOTO_data.dat",unpack=True)
		pi0pt*=1e3;pi0Z*=1e3
		for i in range(np.size(pi0pt)):
		    if i!=0:
			    if i == 2 :
			        ax.scatter(pi0Z[i],pi0pt[i],color='black',marker='+',lw=2.)
			        ax.scatter(pi0Z[i],pi0pt[i],color='orangered',marker='+',lw=1.,s=20)
			    else:
			        ax.scatter(pi0Z[i],pi0pt[i],color='black',marker='+',lw=2.)
			        ax.scatter(pi0Z[i],pi0pt[i],color='lightgreen',marker='+',lw=1,s=20)

	################
	# PLOTS STYLE
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)
	
	ax.set_title(title, fontsize=0.8*fsize)
	# ax.set_yscale('log')

	if filename!=None:
		fig.savefig(filename)
	
	plt.close()

	return ax,fig



import matplotlib.tri as tri
def scatter2D(bag,
	var_x='Zvtx',var_y='Pt_pion_rec',filename=None,
	unitsx=1e3,	unitsy=10,
	xmin=1.5e3,xmax=6.5e3,ymin=0,ymax=300,nbins=20, 
	title='title', XLABEL=r'x', YLABEL=r"y",
	fsize=10,
	plot_signal=False):

	rc('text', usetex=True)
	rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
	                'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(rcparams)
	axes_form  = [0.15,0.15,0.78,0.74/2]
	axes_form2  = [0.15,0.15+0.74/2,0.78,0.74/2]

	fig = plt.figure()
	ax1 = fig.add_axes(axes_form)
	ax2 = fig.add_axes(axes_form2)

	colors=['black','lightblue']
	axes=[ax1,ax2]
	
	j=0
	for cut_events in bag:
		ax = axes[j]
		yy1 = cut_events['cut_events'][var_y]*unitsy
		xx1 = cut_events['cut_events'][var_x]*unitsx
		ww1 = cut_events['cut_events']['w']

		TOT_EVENTS=int(1e4)
		AllEntries = np.array(range(np.shape(ww1)[0]))
		AccEntries = np.random.choice(AllEntries, size=TOT_EVENTS, replace=True, p=ww1/np.sum(ww1))

		x1,y1,w1 = xx1[AccEntries], yy1[AccEntries], ww1[AccEntries]
		w1=np.ones(np.size(AccEntries))

		# bar = ax.hist2d(x1,y1, bins=nbins, weights=w1/np.max(w1), 
		#                 range=[[xmin,xmax],[ymin,ymax]],cmap="Blues",normed=False)#,norm=LogNorm(vmin=np.min(w1/np.max(w1)),vmax=np.max(w1/np.max(w1))))
		# cbar_R = fig.colorbar(bar[3],ax=ax)
		# cbar_R.ax.set_ylabel(r'a.u.', rotation=90)
		
		# bar = ax.scatter(x1,y1, c=w1 , facecolor=colors[i], s=3, lw=0.0, marker = 'o', edgecolor='none',alpha=1,rasterized=True)
		# bar = ax.scatter(xx1,yy1, c=ww1/np.max(ww1), s=3, lw=0.0, marker = 'o', edgecolor='none',alpha=1,rasterized=True)
		bar = ax.hist2d(xx1,yy1, (40,20),norm=LogNorm(),cmap='RdBu_r')
		
		if plot_signal:
			pts = np.array([[3160,130], [3950,130], [4900,160], [4900,250],[3160,250]])
			p = Polygon(pts, closed=True,facecolor=None,edgecolor='red', linestyle='--',fill=False)
			ax.add_patch(p)

			pi0Z, pi0pt = np.genfromtxt("digitized/KOTO_data.dat",unpack=True)
			pi0pt*=1e3;pi0Z*=1e3
			for i in range(np.size(pi0pt)):
				if i!=0:
					if i == 2 :
						ax.scatter(pi0Z[i],pi0pt[i],color='black',marker='+',lw=2.)
						ax.scatter(pi0Z[i],pi0pt[i],color='orangered',marker='+',lw=1.,s=20)
					else:
						ax.scatter(pi0Z[i],pi0pt[i],color='black',marker='+',lw=2.)
						ax.scatter(pi0Z[i],pi0pt[i],color='gold',marker='+',lw=1,s=20)


		j+=1
	
	ax1.set_xlim(xmin,xmax)
	ax1.set_ylim(ymin,ymax)
	ax2.set_xlim(xmin,xmax)
	ax2.set_ylim(ymin,ymax)
		
	ax2.set_xticks([])

	# ax1.set_title(title, fontsize=fsize)

	################
	# PLOTS STYLE
	ax1.set_xlabel(XLABEL,fontsize=fsize)
	ax1.set_ylabel(YLABEL,fontsize=fsize)
	ax2.set_ylabel(YLABEL,fontsize=fsize)

	ax2.set_title(title, fontsize=0.8*fsize)
	# ax.set_yscale('log')

	if filename!=None:
		fig.savefig(filename,dpi=400)
	
	plt.close()

	return axes,fig


########################################################
# 1D FLOW PLOTS
########################################################
def batch_plot_1D(dic,PATH,title,params=None,NBINS=30,YMIN=1e-6,YMAX=1,MODEL=''):
########################################################

	# Pion pT
	filename=PATH+'flow_pT_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var="Pt_pion",reco_var='Pt_pion_rec',
	filename=filename,
	units=1e3,xmin=0.0,xmax=350,nbins=NBINS, 
	title=title, XLABEL=r'$|\vec{p}_\pi^{\,T}|$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=True)

	## Pion pz
	filename=PATH+'flow_pz_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='Pz_pion_rec',
	filename=filename,
	units=1e3,xmin=0.0,xmax=5e3,nbins=NBINS, 
	title=title, XLABEL=r'$p_\pi^z$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion pt/pz
	filename=PATH+'flow_pToverpZ_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='PtoverPz_pion_rec',
	filename=filename,
	units=1,xmin=0.0,xmax=2,nbins=NBINS, 
	title=title, XLABEL=r'$|\vec{p}_\pi^{\,T}|/p_\pi^z$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosT
	filename=PATH+'flow_cosT_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='cosT',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{T}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosgamma1
	filename=PATH+'flow_cosgamma1_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var='CthetaGamma1',reco_var='CthetaGamma1_rec',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosgamma2
	filename=PATH+'flow_cosgamma2_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var='CthetaGamma2',reco_var='CthetaGamma2_rec',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_2}}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion Zvtx
	filename=PATH+'flow_Zvtx_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var='Zvtx_true',reco_var='Zvtx',
	filename=filename,
	units=10,xmin=0,xmax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion Xvtx
	filename=PATH+'flow_Xvtx_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='Xvtx_true',
	filename=filename,
	units=10,xmin=-200,xmax=200,nbins=NBINS, 
	title=title, XLABEL=r'$X_{\rm vtx}^T$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion Yvtx
	filename=PATH+'flow_Yvtx_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='Yvtx_true',
	filename=filename,
	units=10,xmin=-200,xmax=200,nbins=NBINS, 
	title=title, XLABEL=r'$Y_{\rm vtx}^T$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	## Pion theta_XY_2gamma_rec
	filename=PATH+'flow_thetaXY2gamma_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='theta_XY_2gamma_rec',
	filename=filename,
	units=180/np.pi,xmin=0,xmax=180,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma\gamma}^{\rm XY}$ ($^\circ$)', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion r1
	filename=PATH+'flow_r1_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='r1',
	filename=filename,
	units=10,xmin=0,xmax=1e3,nbins=NBINS, 
	title=title, XLABEL=r'$r_1$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion r2
	filename=PATH+'flow_r2_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='r2',
	filename=filename,
	units=10,xmin=0,xmax=1e3,nbins=NBINS, 
	title=title, XLABEL=r'$r_2$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion Epion
	filename=PATH+'flow_Epion_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var=None,reco_var='Epion_rec',
	filename=filename,
	units=1e3,xmin=0,xmax=3e3,nbins=NBINS, 
	title=title, XLABEL=r'$E_\pi$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# theta12
	filename=PATH+'flow_theta12_'+MODEL+'.pdf'
	_ = flow_plot(dic,
	true_var='thetaGamma12',reco_var='thetaGamma12_rec',
	filename=filename,
	units=180/np.pi,xmin=0,xmax=180,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma \gamma}/^\circ$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)


########################################################
# 2D PLOTS
########################################################
def batch_plot_2D(true,
					PATH,
					title,
					params=None,
					NBINS=20,
					MODEL=''):
########################################################

	# Pion signal
	filename=PATH+'2Dhist_pT_Zvtx_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Zvtx',var_y='Pt_pion_rec',
	filename=filename,
	unitsy=1e3,unitsx=10,
	xmin=1.5e3,xmax=6.5e3,ymin=0,ymax=400,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}$/mm', 
	YLABEL=r"$|\vec{p}_\pi^{\,T}|$/MeV",
	fsize=10, plot_signal=True)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_ZvtxR_ZvtxT_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Zvtx',var_y='Zvtx_true',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=0,xmax=7e3,ymin=0,ymax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}^{\rm RECO}$/mm', 
	YLABEL=r"$Z_{\rm vtx}^{\rm TRUE}$/mm",
	fsize=10, plot_signal=False)
	
	# Pion xk vs yk
	filename=PATH+'2Dhist_Xvtx_Yvtx_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Xvtx_true',var_y='Yvtx_true',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=-100.0,xmax=100.0,ymin=-100.0,ymax=100.0,nbins=NBINS, 
	title=title, XLABEL=r'$X_{\rm vtx}^T$/mm', 
	YLABEL=r"$Y_{\rm vtx}^T$/mm",
	fsize=10, plot_signal=False)

	# Pion xk vs Zvtx_true
	filename=PATH+'2Dhist_Xvtx_ZvtxT_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Xvtx_true',var_y='Zvtx_true',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=-100.0,xmax=100.0,ymin=0,ymax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$X_{\rm vtx}^T$/mm', 
	YLABEL=r"$Z_{\rm vtx}^T$/mm",
	fsize=10, plot_signal=False)	

	# Pion xk vs Zvtx
	filename=PATH+'2Dhist_Xvtx_ZvtxR_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Xvtx_true',var_y='Zvtx',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=-100.0,xmax=100.0,ymin=0,ymax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$X_{\rm vtx}^T$/mm', 
	YLABEL=r"$Z_{\rm vtx}^R$/mm",
	fsize=10, plot_signal=False)


	# Pion CTHETA1_RECO vs CTHETA2_RECO
	filename=PATH+'2Dhist_costheta12_rec_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1_rec',var_y='CthetaGamma2_rec',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.9,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}^{\rm REC}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_2}^{\rm REC}}$",
	fsize=10, plot_signal=False)

	# Pion CTHETA1_TRUE vs CTHETA2_TRUE
	filename=PATH+'2Dhist_costheta12_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1',var_y='CthetaGamma2',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.9,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_2}}$",
	fsize=10, plot_signal=False)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_costheta1_smear_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1_rec',var_y='CthetaGamma1',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.95,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}^{\rm REC}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_1}}$",
	fsize=10, plot_signal=False)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_theta12_smear_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='thetaGamma12_rec',var_y='thetaGamma12',
	filename=filename,
	unitsy=180.0/np.pi,unitsx=180.0/np.pi,
	xmin=0.,xmax=60,ymin=0.,ymax=60,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma \gamma}^{\rm REC}/^\circ$', 
	YLABEL=r"${\theta_{\gamma \gamma}}/^\circ$",
	fsize=10, plot_signal=False)


########################################################
# 1D FLOW PLOTS
########################################################
def batch_plot_1D_dipole(dic,PATH,title,params=None,NBINS=30,YMIN=1e-6,YMAX=1,MODEL=''):
########################################################

	# Pion pT
	filename=PATH+'flow_pT_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var="pt_gamma1",reco_var='Pt_pion_rec',
	filename=filename,
	units=1e3,xmin=0.0,xmax=400,nbins=NBINS, 
	title=title, XLABEL=r'$|\overrightarrow{p_\pi}^T|$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=True)

	## Pion pz
	filename=PATH+'flow_pz_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='Pz_pion_rec',
	filename=filename,
	units=1e3,xmin=0.0,xmax=5e3,nbins=NBINS, 
	title=title, XLABEL=r'$p_\pi^z$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion pt/pz
	filename=PATH+'flow_pToverpZ_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='PtoverPz_pion_rec',
	filename=filename,
	units=1,xmin=0.0,xmax=2,nbins=NBINS, 
	title=title, XLABEL=r'$|\overrightarrow{p_\pi}^T|/p_\pi^z$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosT
	filename=PATH+'flow_cosT_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='cosT',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{T}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosgamma1
	filename=PATH+'flow_cosgamma1_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='CthetaGamma1',reco_var='CthetaGamma1_rec',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion cosgamma2
	filename=PATH+'flow_cosgamma2_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='CthetaGamma2',reco_var='CthetaGamma2_rec',
	filename=filename,
	units=1,xmin=0,xmax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_2}}$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)


	# Pion Zvtx
	filename=PATH+'flow_Zvtx_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='Zvtx2',reco_var='Zvtx',
	filename=filename,
	units=10,xmin=0,xmax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)


	## Pion theta_XY_2gamma_rec
	filename=PATH+'flow_thetaXY2gamma_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='theta_XY_2gamma_rec',
	filename=filename,
	units=180/np.pi,xmin=0,xmax=180,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma\gamma}^{\rm XY}$ ($^\circ$)', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion r1
	filename=PATH+'flow_r1_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='r1',
	filename=filename,
	units=10,xmin=0,xmax=1e3,nbins=NBINS, 
	title=title, XLABEL=r'$r_1$/mm', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# Pion Epion
	filename=PATH+'flow_Epion_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var=None,reco_var='Epion_rec',
	filename=filename,
	units=1e3,xmin=0,xmax=3e3,nbins=NBINS, 
	title=title, XLABEL=r'$E_\pi$/MeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# theta12
	filename=PATH+'flow_theta12_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='thetaGamma12',reco_var='thetaGamma12_rec',
	filename=filename,
	units=180/np.pi,xmin=0,xmax=180,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma \gamma}/^\circ$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# PX2
	filename=PATH+'flow_pX2_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='pX2',reco_var='pX2',
	filename=filename,
	units=1,xmin=0,xmax=5,nbins=NBINS, 
	title=title, XLABEL=r'$\vec{p}_{X2}/$GeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# # PX2prime
	filename=PATH+'flow_pX2prime_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='pX2prime',reco_var='pX2prime',
	filename=filename,
	units=1,xmin=0,xmax=5,nbins=NBINS, 
	title=title, XLABEL=r'$\vec{p}_{X2}^\prime/$GeV', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)

	# PX2prime
	filename=PATH+'flow_thetaX2X2_dipole_%.0f_%.0f.pdf'%(params.m2*1e3,params.m1*1e3)
	_ = flow_plot(dic,
	true_var='thetaX2X2',reco_var='thetaX2X2',
	filename=filename,
	units=180/np.pi,xmin=0,xmax=90,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{X2 X2}/^\circ$', YLABEL=r"pdf",
	fsize=10, ymin=YMIN,ymax=YMAX,
	plot_signal=False)


########################################################
# 2D PLOTS
########################################################
def batch_plot_2D_dipole(true,
					PATH,
					title,
					params=None,
					NBINS=20,
					MODEL=''):
########################################################

	# Pion signal
	filename=PATH+'2Dhist_pT_Zvtx_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Zvtx',var_y='Pt_pion_rec',
	filename=filename,
	unitsy=1e3,unitsx=10,
	xmin=1.5e3,xmax=6.5e3,ymin=0,ymax=400,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}$/mm', 
	YLABEL=r"$|\vec{p}_\pi^{\,T}|$/MeV",
	fsize=10, plot_signal=True)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_Zvtx2_Zvtx_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Zvtx',var_y='Zvtx2',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=0,xmax=7e3,ymin=0,ymax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}^{\rm RECO}$/mm', 
	YLABEL=r"$Z_{\rm vtx}^{X_2}$/mm",
	fsize=10, plot_signal=False)
	
	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_Zvtx2prime_Zvtx_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='Zvtx',var_y='Zvtx2prime',
	filename=filename,
	unitsy=10,unitsx=10,
	xmin=0,xmax=7e3,ymin=0,ymax=7e3,nbins=NBINS, 
	title=title, XLABEL=r'$Z_{\rm vtx}$/mm', 
	YLABEL=r"$Z_{\rm vtx}^{X_2^\prime}$/mm",
	fsize=10, plot_signal=False)

	# Pion CTHETA1_RECO vs CTHETA2_RECO
	filename=PATH+'2Dhist_costheta12_rec_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1_rec',var_y='CthetaGamma2_rec',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.9,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}^{\rm REC}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_2}^{\rm REC}}$",
	fsize=10, plot_signal=False)

	# Pion CTHETA1_TRUE vs CTHETA2_TRUE
	filename=PATH+'2Dhist_costheta12_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1',var_y='CthetaGamma2',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.9,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_2}}$",
	fsize=10, plot_signal=False)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_costheta1_smear_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='CthetaGamma1_rec',var_y='CthetaGamma1',
	filename=filename,
	unitsy=1,unitsx=1,
	xmin=0.9,xmax=1,ymin=0.95,ymax=1,nbins=NBINS, 
	title=title, XLABEL=r'$\cos{\theta_{\gamma_1}^{\rm REC}}$', 
	YLABEL=r"$\cos{\theta_{\gamma_1}}$",
	fsize=10, plot_signal=False)

	# Pion zk vs Zvtx
	filename=PATH+'2Dhist_theta12_smear_'+MODEL+'.pdf'
	_ = hist_2D(true,
	var_x='thetaGamma12_rec',var_y='thetaGamma12',
	filename=filename,
	unitsy=180.0/np.pi,unitsx=180.0/np.pi,
	xmin=0.,xmax=60,ymin=0.,ymax=60,nbins=NBINS, 
	title=title, XLABEL=r'$\theta_{\gamma \gamma}^{\rm REC}/^\circ$', 
	YLABEL=r"${\theta_{\gamma \gamma}}/^\circ$",
	fsize=10, plot_signal=False)


def pt_plot(flow,
	true_var=None,reco_var='Pt_pion_rec',filename=None,
	units=1e3,xmin=0.0,xmax=300,nbins=30, 
	title='title', XLABEL=r'x', YLABEL=r"pdf",
	fsize=10, ymin=1e-6,ymax=100,
	plot_signal=False):


	##############################
	# FIGURE
	##############################
	rc('text', usetex=True)
	paramsrc={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
	                'figure.figsize':(1.2*3.7,1.4*2.3617)	}
	rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
	rcParams.update(paramsrc)
	axes_form  = [0.15,0.15,0.78,0.74]

	fig = plt.figure()
	ax = fig.add_axes(axes_form)

	if true_var==None:
		#############################
		# RECO_0
		hist1 = np.histogram(flow['cut_events'][0][reco_var]*units, 
								weights=flow['cut_events'][0]['w'], 
								bins=nbins, 
								density = False, 
								range = (xmin,xmax) )
		HISTNORM=np.sum(hist1[0])*(hist1[1][1]-[hist1[1][0]])
		#########
	else: 
		#########
		# TRUE 
		hist1 = np.histogram(flow['cut_events'][0][true_var]*units,
								 weights=flow['cut_events'][0]['w'], 
								 bins=nbins, 
								 density = False, 
								 range = (xmin,xmax) )
		x0 = hist1[1][:nbins]
		HISTNORM=np.sum(hist1[0])*(hist1[1][1]-[hist1[1][0]])
		y0 = hist1[0]/HISTNORM
		# ax.bar(x0,y0, x0[1]-x0[0],\
		#         ec="black", fc='black', alpha=0.4, align='edge', lw = 0.0)	
		ax.step(np.append(x0,10e10), np.append(y0,0.0), where='post',
		        c='black', dashes=(5,0), lw = 1.2,\
		label=r'T '+true_var.replace("_", "\_"))
		# ax.bar([0,0],[0,0],0,ec='black', fc=[0,0,0,0.4], align='edge', ls='--', lw = .5)

	###############
	# Reco+cuts
	def my_hist(cut_events,eff,color,label,alpha=0.4):

	    hist1 = np.histogram(cut_events[reco_var]*units, 
		    					weights=cut_events['w'], 
		    					bins=nbins, 
		    					density = False, 
		    					range = (xmin,xmax))
	    
	    x0 = hist1[1][:nbins]
	    y0 = hist1[0]/HISTNORM
	    ax.bar(x0,y0, x0[1]-x0[0],\
	            ec=color, fc=color, alpha=alpha, align='edge', lw = 0.0)
	    ax.step(np.append(x0,10e10), np.append(y0,0.0), 
	    	where='post', c='black', lw = 0.4)
	    # for correct legend
	    ax.bar([0,0],[0,0],0,ec='black', fc=np.array(color)*[1,1,1,alpha], align='edge', lw = 0.4,label=label)

	NTOT=np.shape(flow['effs'])[0]
	cmap = cm.BuPu(np.linspace(0.1,1,int(NTOT)))
	cmap = colors.ListedColormap(cmap)	

	###########################
	# PLOT THE FLOW HISTOGRAMS
	for i in range(NTOT):
		my_hist(flow['cut_events'][i],flow['effs'][i],cmap(float(i)/NTOT),alpha=0.9,label=flow['labels'][i])

	# set limits now
	ax.set_xlim(xmin,xmax)

	if ax.get_ylim()[1]>ymax:
		ax.set_ylim(ymin,ax.get_ylim()[1]*1000)
	else:
		ax.set_ylim(ymin,ymax)

	
	#############################
	# PLOT THE KOTO EVENTS 
	if plot_signal:
		pi0Z, pi0pt = np.genfromtxt("digitized/KOTO_data.dat",unpack=True)
		pi0pt*=1e3
		ax.vlines(pi0pt*units,ax.get_ylim()[0], ax.get_ylim()[1], linestyle=':', linewidth=0.8)
		for i in range(np.size(pi0pt)):
		    if i!=0:
		        if i==2:
		            ax.annotate(r'', xy=(pi0pt[i],ax.get_ylim()[0]), xytext=(pi0pt[i],0.3e-5), 
		                        arrowprops=dict(width=1,headwidth=4,headlength=5,facecolor='red',edgecolor='black',lw=0.4))
		        else:
		            ax.annotate(r'', xy=(pi0pt[i],ax.get_ylim()[0]), xytext=(pi0pt[i],0.3e-5), 
		                        arrowprops=dict(width=1,headwidth=4,headlength=5,facecolor='gold',edgecolor='black',lw=0.4))

		ax.vlines(250,ax.get_ylim()[0], ax.get_ylim()[1], linestyle='--', linewidth=0.7)
		ax.vlines(130,ax.get_ylim()[0], ax.get_ylim()[1], linestyle='--', linewidth=0.7)
		ax.annotate(r'', xy=(130+20,20e-2), xytext=(128,20e-2), 
		                        arrowprops=dict(arrowstyle="-|>",color='black',lw=0.6,mutation_scale=6))
		ax.annotate(r'', xy=(250-20,20e-2), xytext=(252,20e-2), 
		                        arrowprops=dict(arrowstyle="-|>",color='black',lw=0.6,mutation_scale=6))
		ax.annotate(r'\noindent\quad KOTO\\signal region', xy=(158,10e-2),color='black',fontsize=0.8*fsize)



	################
	# PLOTS STYLE
	ax.set_xlabel(XLABEL,fontsize=fsize)
	ax.set_ylabel(YLABEL,fontsize=fsize)

	plt.legend(loc='upper right',frameon=False,fontsize=0.6*fsize)


	ax.set_title(title, fontsize=0.8*fsize)
	# ax.set_yscale('log')

	if filename!=None:
		fig.savefig(filename)
	
	plt.close()

	return ax,fig
