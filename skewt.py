#coding = utf-8
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rc('font', family='Arial')
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units
import metpy.interpolate as mpi
import datetime

def get_pressure_level_index(p, plevel, reverse=False):
    if reverse:
        idx = 0
    else:
        idx = -1
    return np.where(p.magnitude>=plevel)[0][idx]

def cal_tt(t850, td850, t500):
    return t850 + td850 - 2 * t500

def execute_sweat(tt, td8, u8, v8, u5, v5):
    s8 = np.sqrt(u8 * u8 + v8 * v8)
    s5 = np.sqrt(u5 * u5 + v5 * v5)
    s = s8 * s5
    z_mask = s == 0
    nz_mask = s != 0
    sdir = np.ndarray(tt.shape, np.float32)
    sdir[z_mask] = 0
    sdir[nz_mask] = (u5[nz_mask] * v8[nz_mask] - v5[nz_mask] * u8[nz_mask]) / s[nz_mask]
    tt49 = np.ndarray(tt.shape, np.float32)
    tt49_mask = tt >= 49
    tt49[tt < 49] = 0
    tt49[tt49_mask] = tt[tt49_mask] - 49.0
    result = 12 * td8
    result += 20 * tt49
    result += 2 * 1.944 * s8
    result += s5 * 1.944
    result += 125 * (sdir + 0.2)
    return result

def c_sweat(t850, td850, t500, u850, v850, u500, v500):
    tt = cal_tt(t850, td850, t500)
    return execute_sweat(tt, td850, u850, v850, u500, v500)

def showalter_index(t850, td850, t500):
    plcl, tlcl = mpcalc.lcl(850 * units('hPa'), t850, td850)
    p = np.array([plcl.magnitude, 500]) * units('hPa')
    out = mpcalc.moist_lapse(p, tlcl)
    return t500.magnitude - out[1].magnitude

def lifted_index(tsfc, tdsfc, psfc, t500):
    plcl, tlcl = mpcalc.lcl(psfc, tsfc, tdsfc)
    p = np.array([plcl.magnitude, 500]) * units('hPa')
    out = mpcalc.moist_lapse(p, tlcl)
    return t500.magnitude - out[1].magnitude

def delta_height(p1, p2):
    h1 = mpcalc.pressure_to_height_std(p1)
    h2 = mpcalc.pressure_to_height_std(p2)
    return h2 - h1

class Skew_T(SkewT):
    def __init__(self, pres, temp, rh=None, td=None, u=None, v=None, wspd=None,
                 wdir=None, alt=None, station=None, time=None, fig=None, **kwargs):
        if not fig:
            fig = plt.figure(figsize=(9, 9), dpi=200)
        super().__init__(fig=fig, rotation=30)
        # Parameter conversion
        nonetype = type(None)
        if isinstance(td, nonetype):
            td = mpcalc.dewpoint_rh(temp, rh * units.percent).magnitude
        if isinstance(u, nonetype) or isinstance(v, nonetype):
            u, v = mpcalc.wind_components(wspd * units('m/s'), wdir * units.degree)
            u = u.magnitude
            v = v.magnitude
        self.kw = kwargs
        # Interpolate Nans
        xi = np.arange(0, len(pres), 1)
        self.p_i = mpi.interpolate_nans_1d(xi, pres) * units('hPa')
        self.t_i = mpi.interpolate_nans_1d(self.p_i, temp) * units.degC
        self.td_i = mpi.interpolate_nans_1d(self.p_i, td) * units.degC
        self.u_i = mpi.interpolate_nans_1d(self.p_i, u) * units('m/s')
        self.v_i = mpi.interpolate_nans_1d(self.p_i, v) * units('m/s')
        self.alt = mpi.interpolate_nans_1d(self.p_i, alt) * units('m')
        self.st = station
        self.time = time
        self.dp_idx = np.where(~np.isnan(td))[0][-1]

        self.process_skewt()
        
    def process_skewt(self):
        # Calculation
        index_p100 = get_pressure_level_index(self.p_i, 100)
        lcl_p, lcl_t = mpcalc.lcl(self.p_i[0], self.t_i[0], self.td_i[0])
        lfc_p, lfc_t = mpcalc.lfc(self.p_i, self.t_i, self.td_i)
        el_p, el_t = mpcalc.el(self.p_i, self.t_i, self.td_i)
        prof = mpcalc.parcel_profile(self.p_i, self.t_i[0], self.td_i[0]).to('degC')
        cape, cin = mpcalc.cape_cin(self.p_i, self.t_i, self.td_i, prof)
        mucape, mucin = mpcalc.most_unstable_cape_cin(self.p_i, self.t_i, self.td_i)
        pwat = mpcalc.precipitable_water(self.td_i, self.p_i)
        i8 = get_pressure_level_index(self.p_i, 850)
        i7 = get_pressure_level_index(self.p_i, 700)
        i5 = get_pressure_level_index(self.p_i, 500)
        theta850 = mpcalc.equivalent_potential_temperature(850 * units('hPa'), self.t_i[i8], self.td_i[i5])
        theta500 = mpcalc.equivalent_potential_temperature(500 * units('hPa'), self.t_i[i5], self.td_i[i5])
        thetadiff = theta850 - theta500
        k = self.t_i[i8] - self.t_i[i5] + self.td_i[i8] - (self.t_i[i7] - self.td_i[i7])
        a = ((self.t_i[i8] - self.t_i[i5]) - (self.t_i[i8] - self.td_i[i5]) -
            (self.t_i[i7] - self.td_i[i7]) - (self.t_i[i5] - self.td_i[i5]))
        sw = c_sweat(np.array(self.t_i[i8].magnitude), np.array(self.td_i[i8].magnitude),
                     np.array(self.t_i[i5].magnitude), np.array(self.u_i[i8].magnitude),
                     np.array(self.v_i[i8].magnitude), np.array(self.u_i[i5].magnitude),
                     np.array(self.v_i[i5].magnitude))
        si = showalter_index(self.t_i[i8], self.td_i[i8], self.t_i[i5])
        li = lifted_index(self.t_i[0], self.td_i[0], self.p_i[0], self.t_i[i5])
        srh_pos, srh_neg, srh_tot = mpcalc.storm_relative_helicity(self.u_i, self.v_i, self.alt, 1000 * units('m'))
        sbcape, sbcin = mpcalc.surface_based_cape_cin(self.p_i, self.t_i, self.td_i)
        shr6km = mpcalc.bulk_shear(self.p_i, self.u_i, self.v_i, heights=self.alt, depth=6000 * units('m'))
        wshr6km = mpcalc.wind_speed(*shr6km)
        sigtor = mpcalc.significant_tornado(sbcape, delta_height(self.p_i[0], lcl_p), srh_tot, wshr6km)[0]
        # Plotting
        self.ax.set_ylim(1050, 100)
        self.ax.set_xlim(-40, 50)
        self.plot(self.p_i, self.t_i, 'r', linewidth=1)
        self.plot(self.p_i[:self.dp_idx], self.td_i[:self.dp_idx], 'g', linewidth=1)
        self.plot_barbs(self.p_i[:index_p100], self.u_i[:index_p100] * 1.94, self.v_i[:index_p100] * 1.94)
        self.plot(lcl_p, lcl_t, 'ko', markerfacecolor='black')
        self.plot(self.p_i, prof, 'k', linewidth=2)
        if cin.magnitude < 0:
            chi = -1 * cin.magnitude
            skew.shade_cin(self.p_i, self.t_i, prof)
        elif cin.magnitude > 0:
            chi = cin.magnitude
            skew.shade_cin(self.p_i, self.t_i, prof)
        else:
            chi = 0.
        self.shade_cape(self.p_i, self.t_i, prof)
        self.plot_dry_adiabats(linewidth=0.5)
        self.plot_moist_adiabats(linewidth=0.5)
        self.plot_mixing_lines(linewidth=0.5)
        plt.title('Skew-T Plot \nStation: {} Time: {}'.format(self.st, self.time.strftime('%Y.%m.%d %H:%M')), fontsize=14, loc='left')
        # Add hodograph
        ax = self._fig.add_axes([0.95, 0.71, 0.17, 0.17])
        h = Hodograph(ax, component_range=50)
        h.add_grid(increment=20)
        h.plot_colormapped(self.u_i[:index_p100], self.v_i[:index_p100], self.alt[:index_p100], linewidth=1.2)
        # Annotate parameters
        # Annotate names
        namelist = ['CAPE', 'CIN', 'MUCAPE', 'PWAT', 'K', 'A', 'SWEAT', 'LCL', 'LFC', 'EL', 'SI', 'LI', 'T850-500',
                    'θse850-500', 'SRH', 'STP']
        xcor = -50
        ycor = -90
        spacing = -9
        for nm in namelist:
            ax.text(xcor, ycor, '{}: '.format(nm), fontsize=10)
            ycor += spacing
        # Annotate values
        varlist = [cape, chi, mucape, pwat, k, a, sw, lcl_p, lfc_p, el_p, si, li, self.t_i[i8] - self.t_i[i5], thetadiff,
                   srh_tot, sigtor]
        xcor = 10
        ycor = -90
        for v in varlist:
            if hasattr(v, 'magnitude'):
                v = v.magnitude
            ax.text(xcor, ycor, str(np.round_(v, 2)), fontsize=10)
            ycor += spacing
        # Annotate units
        unitlist = ['J/kg', 'J/kg', 'J/kg', 'mm', '°C', '°C', '', 'hPa', 'hPa', 'hPa', '°C', '°C', '°C', '°C']
        xcor = 45
        ycor = -90
        for u in unitlist:
            ax.text(xcor, ycor, ' {}'.format(u), fontsize=10)
            ycor += spacing
    
    def save(self, fpath=None):
        if not fpath:
            fpath = 'D:\\SKEWT_{}_{}.png'.format(self.time.strftime('%Y%m%d%H%M%S'), self.st)
        plt.savefig(fpath, bbox_inches='tight')
