import numpy as np
import pandas as pd

class Diamond5(object):
    def __init__(self, filepath):
        self.f = open(filepath, 'r')

    def _get_data_section(self, station_id):
        if isinstance(station_id, int):
            station_id = str(station_id)
        self.f.seek(0)
        content = self.f.read()
        try:
            pos = content.index(station_id)
        except Exception as e:
            return None, None
        self.f.seek(pos)
        out = list()
        info = self.f.readline()
        while 1:
            line = self.f.readline()
            if line.startswith(' '):
                out.append(line.strip().split('  '))
            else:
                self.f.close()
                return info.strip().split(' '), out
    
    def get_data(self, station_id):
        info, data = self._get_data_section(station_id)
        data = np.array(data, dtype=float)
        data[data == 9999] = np.nan
        df = pd.DataFrame(data, columns=['pres', 'hght', 'temp', 'dwpt', 'wdir', 'wspd'])
        df_new = df.drop_duplicates('pres', 'last')
        return df_new