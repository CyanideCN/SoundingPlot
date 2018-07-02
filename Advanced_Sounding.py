# -*- coding: UTF-8 -*-
#Version: 1.4.3
#Updated at 2018.05.19
#增加KY指数计算

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units
from matplotlib.font_manager import FontProperties
import os
from tkinter import filedialog

deg_rad=np.pi/180
font = FontProperties(fname=r"C:\\WINDOWS\\Fonts\\Dengl.ttf")
font2 = FontProperties(fname=r"C:\\WINDOWS\\Fonts\\Consola.ttf")

def get_layer(level, array, counts=False, refpres=True, refarray=[]):
    if str(level)=='nan':
        return None
    if refpres==True:
        ref=p1
    elif refpres==False:
        ref=refarray
    count=0
    if np.max(ref)>level:
        while count<len(ref):
            if ref[count]==level:
                tar=array[count]
            elif ref[count]<level:
                det1=ref[count-1]-level
                det2=level-ref[count]
                if det1<det2:
                    tar=array[count-1]
                elif det1>det2:
                    tar=array[count]
                elif det1==det2:
                    tar=array[count]
                break
            count=count+1
    elif np.max(ref)<level:
        tar=None
    if counts==False:
        return tar
    if counts==True:
        return count

file = filedialog.askopenfilename(filetypes=[('text Files','*.txt')],title='Open L-band Radar Files')
full_name = os.path.split(file)[-1]
(file_name, file_type) = os.path.splitext(full_name)
station=file_name.split('_', 7)[3]
time=file_name.split('_', 7)[4]
year=time[:4]
month=time[4:6]
date=time[6:8]
hour=time[8:10]
min=time[10:12]
sec=time[12:14]

f=open(file)
f.readlines(1)
info=f.readlines(1)[0][:-1].split(' ', 4)
stationlon=float(info[1])
stationlat=float(info[2])

mnames=['temp','pres', 'hum', 'E', 'F', 'G', 'H', 'I', 'dir', 'wind', 'alt', 'x']
table=pd.read_table(file ,sep=' ',names=mnames, header=None, skiprows=6, skipfooter=150, engine="python", encoding='utf-8')

temp=table['temp']
temp1=np.array(temp)
temp1[temp1 =='/////']=None
t=np.array(temp1, dtype=float)
pres=table['pres']
pres1=np.array(pres)
pres1[pres1 =='//////']=None
p=np.array(pres1, dtype=float)
hum=table['hum']
hum1=np.array(hum)
hum1[hum1 =='///']=None
h=np.array(hum1, dtype=float)
dir=table['dir']
dir1=np.array(dir)
dir1[dir1 =='///']=None
d=np.array(dir1, dtype=float)
wind=table['wind']
wnd=np.array(wind)
wnd[wnd =='///']=None
wnd1=3.6*np.array(wnd, dtype=float)
alt=table['alt']
alt=np.array(alt)

a=np.arange(0, len(pres), 1)
p1=mpcalc.interpolate_nans(a, p)
p2=p1* units.hPa
t1=mpcalc.interpolate_nans(p1, t)
t2=t1* units.degC
h1=mpcalc.interpolate_nans(p1, h)
h2=h1* units.percent
d1=mpcalc.interpolate_nans(p1, d)
d2=d1* units.degree
wd1=mpcalc.interpolate_nans(p1, wnd1)/3.6
wd2=wd1* units.kph
del temp, pres, hum, dir, wind
#Dew point
td=t1-((14.55+0.114*t1)*(1-0.01*h1)+np.power((2.5+0.007*t1)*(1-0.01*h1),3)+(15.9+0.37*t1)*np.power((1-0.01*h1),14))
td2=td* units.degC

fig = plt.figure(figsize=(9, 9), dpi=200)
skew = SkewT(fig, rotation=30)

skew.plot(p2, t2, 'r', linewidth=1)
skew.plot(p2, td2, 'g', linewidth=1)
u, v = mpcalc.get_wind_components(wd2, d2)
s=slice(0, len(u), 100)
u1=u[s]
v1=v[s]
p3=p2[s]
count=0
while count<len(p3):
    if p3.magnitude[count]<100:
        break
    count=count+1
u1=u1[0:count]
v1=v1[0:count]
p3=p3[0:count]
skew.plot_barbs(p3, u1, v1)
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 50)
# Calculate LCL height and plot as black dot
lcl_pressure, lcl_temperature = mpcalc.lcl(p2[0], t2[0], td2[0])
lfc_pressure, lfc_temperature = mpcalc.lfc(p2, t2, td2)
el_pressure, el_temperature = mpcalc.el(p2, t2, td2)
lclh=get_layer(lcl_pressure.magnitude, alt)
lfch=get_layer(lfc_pressure.magnitude, alt)
elh=get_layer(el_pressure.magnitude, alt)
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
prof = mpcalc.parcel_profile(p2, t2[0], td2[0]).to('degC')
skew.plot(p2, prof, 'k', linewidth=2)
cape, cin=mpcalc.cape_cin(p2, t2, td2, prof)
# Shade areas of CAPE and CIN
if cin.magnitude<0:
    chi=-1*cin.magnitude
    skew.shade_cin(p2, t2, prof)
elif cin.magnitude>0:
    chi=cin.magnitude
    skew.shade_cin(p2, t2, prof)
elif cin.magnitude==0:
    chi=cin.magnitude
if cape.magnitude<0:
    cpe=-1*cape.magnitude
else:
    cpe=cape.magnitude
skew.shade_cape(p2, t2, prof)
# Add the relevant special lines
skew.plot_dry_adiabats(linewidth=0.5)
skew.plot_moist_adiabats(linewidth=0.5)
skew.plot_mixing_lines(linewidth=0.5)

#Annotate sounding index
mucape, mucin=mpcalc.most_unstable_cape_cin(p2, t2, td2)
s2=slice(0, 5001, 25)
td3=td2[s2]
p4=p2[s2]
pwat=mpcalc.precipitable_water(td3, p4)

t850=get_layer(850, t1)
td850=get_layer(850, td)
t700=get_layer(700, t1)
td700=get_layer(700, td)
t500=get_layer(500, t1)
td500=get_layer(500, td)
wf850=get_layer(850, wd1)
wf500=get_layer(500, wd1)
wd850=get_layer(850, d1)
wd500=get_layer(500, d1)
wf200=get_layer(200, wd1)
wd200=get_layer(200, d1)
h0=get_layer(0, alt, refpres=False, refarray=t1)
hn20=get_layer(-20, alt, refpres=False, refarray=t1)
u850, v850=mpcalc.get_wind_components(wd850*deg_rad, wf850)
u200, v200=mpcalc.get_wind_components(wd200*deg_rad, wf200)
shear200=np.sqrt((u200-u850)**2+(v200-v850)**2)
theta850=mpcalc.equivalent_potential_temperature(850*units.hectopascal, t850*units.degC, td850*units.degC)
theta500=mpcalc.equivalent_potential_temperature(500*units.hectopascal, t500*units.degC, td500*units.degC)
thetadiff=theta850-theta500
#重写
def get_sweat():
    a=12*td850
    b=20*(TT-49)
    c=4*wf850
    d=2*wf500
    e=125*(np.sin((wd500-wd850)/deg_rad)+0.2)
    if a<0:
        a=0
    if b<0:
        b=0
    if c<0:
        c=0
    if d<0:
        d=0
    if e<0:
        e=0
    if wd850<130 or wd850>250 or wd500<210 or wd500>310 or wd500<wd850 or wf850<7.5 or wf500<7.5:
        e=0
    sweat=a+b+c+d+e
    return sweat
ki=None
ai=None
TT=None
sweat=None
try:
    ki=t850-t500+td850-(t700-td700)
    ai=(t850-t500)-(t850-td850)-(t700-td700)-(t500-td500)
    TT=t850+td850-2*t500
    sweat=get_sweat()
except TypeError:
    pass

lvl850=get_layer(850, p1, counts=True)
lvl500=get_layer(500, p1, counts=True)
lvllfc=get_layer(lfc_pressure.magnitude, p1, counts=True)
prof2 = mpcalc.parcel_profile(p1[lvl850:]*units.hectopascal, t850*units.degC, td850*units.degC).to('degC')
pf500=prof2[get_layer(500, p1, counts=True)-lvl850]
si=t500-pf500.magnitude

li=t500-mpcalc.moist_lapse(p1[lvllfc:]*units.hectopascal, lfc_temperature).magnitude[lvl500-lvllfc]
#KY
omega=7.29e-5
rdx=2.87e-3
kts=(-2*omega*np.sin(stationlat*deg_rad))/(np.log(500/850)*rdx)
ta=kts*wf500*wf850*np.sin((wd500-wd850)*deg_rad)
alpha=1
beta=1
gamma=0
ky=((beta*ta)-si+gamma)/(alpha+(t850-td850))
if ta<si:
    ky=0
plt.text(-93.3, 97, 'L-band Radar Skew-T Plot @HCl' +'\nStation: '+station+' Time (UTC): '+year+'.'+month+'.'
            +date+' '+hour+':'+min+':'+sec, fontsize=14)

ax = fig.add_axes([0.95, 0.71, 0.17, 0.17])
h = Hodograph(ax, component_range=50)
h.add_grid(increment=20)
height=get_layer(100, p1, counts=True)
h.plot_colormapped(u[0:height], v[0:height], alt[0:height], linewidth=1.2)
spacing=-9

ax.text(-50, -90, 'CAPE: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing, 'CIN: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*2, 'MUCAPE: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*3, 'PWAT: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*4, 'KI: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*5, 'AI: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*6, 'SWEAT: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*7, 'LCL: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*8, 'LFC: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*9, 'EL: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*10, '0°C lev: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*11, '-20°C lev: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*12, 'T850-500: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*13, 'θse850-500: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*14, 'SI: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*15, 'KY: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*16, 'LI: ', fontsize=14, FontProperties=font2)
ax.text(-50, -90+spacing*17, 'SHR850-200: ', fontsize=14, FontProperties=font2)

ax.text(10, -90, str(np.round_(cpe, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing, str(np.round_(chi, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*2, str(np.round_(mucape.magnitude, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*3, str(np.round_(pwat.magnitude, 2)), fontsize=14, FontProperties=font2)
if ki==None:
    k=''
else:
    k=str(np.round_(ki, 2))
if ai==None:
    a=''
else:
    a=str(np.round_(ai, 2))
if sweat==None:
    sw=''
else:
    sw=str(np.round_(sweat, 2))
ax.text(10, -90+spacing*4, k, fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*5, a, fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*6, sw, fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*7, str(int(lclh)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*8, str(int(lfch)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*9, str(int(elh)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*10, str(h0), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*11, str(hn20), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*12, str(np.round_(t850-t500, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*13, str(np.round_(thetadiff.magnitude, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*14, str(np.round_(si, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*15, str(np.round_(ky, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*16, str(np.round_(li, 2)), fontsize=14, FontProperties=font2)
ax.text(10, -90+spacing*17, str(np.round_(shear200, 2)), fontsize=14, FontProperties=font2)

ax.text(45, -90, ' J/kg', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing, ' J/kg', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*2, ' J/kg', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*3, ' mm', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*7, ' m', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*8, ' m', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*9, ' m', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*10, ' m', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*11, ' m', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*12, ' °C', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*13, ' °C', fontsize=14, FontProperties=font2)
ax.text(45, -90+spacing*17, ' m/s', fontsize=14, FontProperties=font2)

plt.savefig('D:\\Meteorology\\L波段探空\\Skew-T_'+station+'_'+time+'.png', bbox_inches='tight')
#