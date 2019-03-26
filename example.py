from diamond5 import Diamond5
from skewt import Skew_T
from tkinter.filedialog import askopenfilename
import datetime
from pathlib import Path

fp = askopenfilename()
date = str(Path(fp)).split('\\')[-1].split('.')[0]
dtime = datetime.datetime.strptime(date, '%Y%m%d%H%M%S')
d = Diamond5(fp)
data = d.get_data(56187)

skewt = Skew_T(data.pres.values, data.temp.values, td=data.dwpt.values,
               wdir=data.wdir.values, wspd=data.wspd.values, alt=data.hght.values * 10, station=56187, time=dtime)
skewt.save('D:\\1.png')