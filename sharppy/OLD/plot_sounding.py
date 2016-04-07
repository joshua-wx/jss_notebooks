import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.params as params
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils

import numpy as np
import matplotlib.pyplot as plt

from StringIO import StringIO

from matplotlib.axes import Axes
import matplotlib.transforms as transforms
import matplotlib.axis as maxis
import matplotlib.spines as mspines
import matplotlib.path as mpath
from matplotlib.projections import register_projection

import sndfilter as sndfilter

def parseSPC(spc_file):
    ## read in the file
    data = np.array([l.strip() for l in spc_file.split('\n')])

    ## necessary index points
    title_idx = np.where( data == '%TITLE%')[0][0]
    start_idx = np.where( data == '%RAW%' )[0] + 1
    finish_idx = np.where( data == '%END%')[0]

    ## create the plot title
    data_header = data[title_idx + 1].split()
    location = data_header[0]
    time = data_header[1][:11]

    ## put it all together for StringIO
    full_data = '\n'.join(data[start_idx : finish_idx][:])
    sound_data = StringIO( full_data )

    ## read the data into arrays
    pres, h, temp, dewp, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments="%", unpack=True )
    #p, h, T, Td, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments="%", unpack=True )

    ## filter sounding
    #dict_out = sndfilter.soundingFilter(temp,dewp,pres,wspd,wdir)
    #print dict_out
    raw_input()
    return p, h, T, Td, wdir, wspd

# The sole purpose of this class is to look at the upper, lower, or total
# interval as appropriate and see what parts of the tick to draw, if any.
class SkewXTick(maxis.XTick):
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group(self.__name__)

        lower_interval = self.axes.xaxis.lower_interval
        upper_interval = self.axes.xaxis.upper_interval

        if self.gridOn and transforms.interval_contains(
                self.axes.xaxis.get_view_interval(), self.get_loc()):
            self.gridline.draw(renderer)

        if transforms.interval_contains(lower_interval, self.get_loc()):
            if self.tick1On:
                self.tick1line.draw(renderer)
            if self.label1On:
                self.label1.draw(renderer)

        if transforms.interval_contains(upper_interval, self.get_loc()):
            if self.tick2On:
                self.tick2line.draw(renderer)
            if self.label2On:
                self.label2.draw(renderer)

        renderer.close_group(self.__name__)


# This class exists to provide two separate sets of intervals to the tick,
# as well as create instances of the custom tick
class SkewXAxis(maxis.XAxis):
    def __init__(self, *args, **kwargs):
        maxis.XAxis.__init__(self, *args, **kwargs)
        self.upper_interval = 0.0, 1.0

    def _get_tick(self, major):
        return SkewXTick(self.axes, 0, '', major=major)

    @property
    def lower_interval(self):
        return self.axes.viewLim.intervalx

    def get_view_interval(self):
        return self.upper_interval[0], self.axes.viewLim.intervalx[1]


# This class exists to calculate the separate data range of the
# upper X-axis and draw the spine there. It also provides this range
# to the X-axis artist for ticking and gridlines
class SkewSpine(mspines.Spine):
    def _adjust_location(self):
        trans = self.axes.transDataToAxes.inverted()
        if self.spine_type == 'top':
            yloc = 1.0
        else:
            yloc = 0.0
        left = trans.transform_point((0.0, yloc))[0]
        right = trans.transform_point((1.0, yloc))[0]

        pts  = self._path.vertices
        pts[0, 0] = left
        pts[1, 0] = right
        self.axis.upper_interval = (left, right)


# This class handles registration of the skew-xaxes as a projection as well
# as setting up the appropriate transformations. It also overrides standard
# spines and axes instances as appropriate.
class SkewXAxes(Axes):
    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e. ``subplot(111,
    # projection='skewx')``.
    name = 'skewx'

    def _init_axis(self):
        #Taken from Axes and modified to use our modified X-axis
        self.xaxis = SkewXAxis(self)
        self.spines['top'].register_axis(self.xaxis)
        self.spines['bottom'].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines['left'].register_axis(self.yaxis)
        self.spines['right'].register_axis(self.yaxis)

    def _gen_axes_spines(self):
        spines = {'top':SkewSpine.linear_spine(self, 'top'),
                  'bottom':mspines.Spine.linear_spine(self, 'bottom'),
                  'left':mspines.Spine.linear_spine(self, 'left'),
                  'right':mspines.Spine.linear_spine(self, 'right')}
        return spines

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        rot = 30

        #Get the standard transform setup from the Axes base class
        Axes._set_lim_and_transforms(self)

        # Need to put the skew in the middle, after the scale and limits,
        # but before the transAxes. This way, the skew is done in Axes
        # coordinates thus performing the transform around the proper origin
        # We keep the pre-transAxes transform around for other users, like the
        # spines for finding bounds
        self.transDataToAxes = self.transScale + (self.transLimits +
                transforms.Affine2D().skew_deg(rot, 0))

        # Create the full transform from Data to Pixels
        self.transData = self.transDataToAxes + self.transAxes

        # Blended transforms like this need to have the skewing applied using
        # both axes, in axes coords like before.
        self._xaxis_transform = (transforms.blended_transform_factory(
                    self.transScale + self.transLimits,
                    transforms.IdentityTransform()) +
                transforms.Affine2D().skew_deg(rot, 0)) + self.transAxes

#read in data
spc_file = open('201411270400_27-11-14_04Z_BAP.oax', 'r').read()

#prase and create profile object
pres, hght, tmpc, dwpc, wdir, wspd = parseSPC(spc_file)
prof = profile.create_profile(profile='default', pres=pres, hght=hght, tmpc=tmpc, \
                                    dwpc=dwpc, wspd=wspd, wdir=wdir, missing=-9999, strictQC=True)

#create parcel pbjects
sfcpcl = params.parcelx( prof, flag=1 ) # Surface Parcel
fcstpcl = params.parcelx( prof, flag=2 ) # Forecast Parcel
mupcl = params.parcelx( prof, flag=3 ) # Most-Unstable Parcel
mlpcl = params.parcelx( prof, flag=4 ) # 100 mb Mean Layer Parcel
meipcl = params.parcelx( prof, flag=6 ) # Mean Effective Layer Parcel

# Now register the projection with matplotlib so the user can select
# it.
register_projection(SkewXAxes)

#calculate dynamic stats
sfc = prof.pres[prof.sfc]
p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
mean_3km = winds.mean_wind(prof, pbot=sfc, ptop=p3km)
sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
sfc_3km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
#print "0-3 km Pressure-Weighted Mean Wind (kt):", utils.comp2vec(mean_3km[0], mean_3km[1])[1]
#print "0-6 km Shear (kt):", utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1]
srwind = params.bunkers_storm_motion(prof)
#print "Bunker's Storm Motion (right-mover) [deg,kts]:", utils.comp2vec(srwind[0], srwind[1])
#print "Bunker's Storm Motion (left-mover) [deg,kts]:", utils.comp2vec(srwind[2], srwind[3])
srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])
srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])
#print "0-3 km Storm Relative Helicity [m2/s2]:",srh3km[0]

#effective inflow depth stats
stp_fixed = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km[0], utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])
ship = params.ship(prof)
eff_inflow = params.effective_inflow_layer(prof)
ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
#print "Effective Inflow Layer Bottom Height (m AGL):", ebot_hght
#print "Effective Inflow Layer Top Height (m AGL):", etop_hght
effective_srh = winds.helicity(prof, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])
#print "Effective Inflow Layer SRH (m2/s2):", effective_srh[0]
ebwd = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
ebwspd = utils.mag( ebwd[0], ebwd[1] )
#print "Effective Bulk Wind Difference:", ebwspd
scp = params.scp(mupcl.bplus, effective_srh[0], ebwspd)
stp_cin = params.stp_cin(mlpcl.bplus, effective_srh[0], ebwspd, mlpcl.lclhght, mlpcl.bminus)
#print "Supercell Composite Parameter:", scp
#print "Significant Tornado Parameter (w/CIN):", stp_cin
#print "Significant Tornado Parameter (fixed):", stp_fixed

#PLOTTING
indices = {'SB CAPE': [int(sfcpcl.bplus), 'J/kg'],\
           'SB CIN': [int(sfcpcl.bminus), 'J/kg'],\
           'SB LCL': [int(sfcpcl.lclhght), 'm AGL'],\
           'SB LFC': [int(sfcpcl.lfchght), 'm AGL'],\
           'SB EL': [int(sfcpcl.elhght), 'm AGL'],\
           'SB LI': [int(sfcpcl.li5), 'C'],\
           'MEI CAPE': [int(meipcl.bplus), 'J/kg'],\
           'MEI CIN': [int(meipcl.bminus), 'J/kg'],\
           'MEI LCL': [int(meipcl.lclhght), 'm AGL'],\
           'MEI LFC': [meipcl.lfchght, 'm AGL'],\
           'MEI EL': [int(meipcl.elhght), 'm AGL'],\
           'MEI LI': [int(meipcl.li5), 'C'],\
           'EI SRH': [int(effective_srh[0]), 'm2/s2'],\
           'EI Base':  [int(ebot_hght), 'm AGL'],\
           'EI Top':   [int(etop_hght), 'm AGL'],\
           'EI Dpth': [int(etop_hght-ebot_hght), 'm AGL']}

# Set the parcel trace to be plotted as the Mean Effective Layer Parcel
pcl = meipcl

# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='skewx')
ax.grid(True)

pmax = 1000
pmin = 10
dp = -10
presvals = np.arange(int(pmax), int(pmin)+dp, dp)

# plot the moist-adiabats
for t in np.arange(-10,45,5):
    tw = []
    for p in presvals:
        tw.append(thermo.wetlift(1000., t, p))
    ax.semilogy(tw, presvals, 'k-', alpha=.2)

def thetas(theta, presvals):
    return ((theta + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK

# plot the dry adiabats
for t in np.arange(-50,110,10):
    ax.semilogy(thetas(t, presvals), presvals, 'r-', alpha=.2)

plt.title('Test Sounding', fontsize=12, loc='left')
# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dicatated by the typical meteorological plot
ax.semilogy(prof.tmpc, prof.pres, 'r', lw=2) # Plot the temperature profile
ax.semilogy(prof.vtmp, prof.pres, 'r', lw=1) # Plot the wetbulb profile
ax.semilogy(prof.dwpc, prof.pres, 'g', lw=2) # plot the dewpoint profile
ax.semilogy(pcl.ttrace, pcl.ptrace, 'k', lw=2) # plot the parcel trace 

# Plot the effective inflow layer using blue horizontal lines
ax.axhline(eff_inflow[0], color='b')
ax.axhline(eff_inflow[1], color='b')

#plt.barbs(10*np.ones(len(prof.pres)), prof.pres, prof.u, prof.v)
# Disables the log-formatting that comes with semilogy
ax.yaxis.set_major_formatter(plt.ScalarFormatter())
ax.set_yticks(np.linspace(100,1000,10))
ax.set_ylim(1050,100)
ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.set_xlim(-50,50)

# List the indices within the indices dictionary on the side of the plot.
string = ''
for key in np.sort(indices.keys()):
    string = string + key + ': ' + str(indices[key][0]) + ' ' + indices[key][1] + '\n'
plt.text(.75, .95, string, verticalalignment='top', transform=plt.gca().transAxes)

# Draw the hodograph on the Skew-T.
# TAS 2015-4-16: hodograph doesn't plot for some reason ...
#ax2 = plt.axes([.625,.625,.25,.25])
#below_12km = np.where(interp.to_agl(prof, prof.hght) < 12000)[0]
#u_prof = prof.u[below_12km]
#v_prof = prof.v[below_12km]
#ax2.plot(u_prof[~u_prof.mask], v_prof[~u_prof.mask], 'k-', lw=2)
#ax2.get_xaxis().set_visible(False)
#ax2.get_yaxis().set_visible(False)
#for i in range(10,90,10):
#    # Draw the range rings around the hodograph.
#    circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
#    ax2.add_artist(circle)
#ax2.plot(srwind[0], srwind[1], 'ro') # Plot Bunker's Storm motion right mover as a red dot
#ax2.plot(srwind[2], srwind[3], 'bo') # Plot Bunker's Storm motion left mover as a blue dot
#ax2.set_xlim(-60,60)
#ax2.set_ylim(-60,60)
#ax2.axhline(y=0, color='k')
#ax2.axvline(x=0, color='k')

plt.show()











