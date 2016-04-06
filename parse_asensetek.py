import json
import spectrum

def parse(fname):
    lamp = spectrum.Spectrum()
    lamp.zero()
    handle=open(fname,'r')
    data = json.load(handle)
    handle.close()
    for sample in data['data_points'][0]['spectrumPoints']:
        nm = int(str(sample.keys()[0]))
        if nm >= lamp.start and nm <= lamp.end:
            lamp.data[nm - lamp.start] = sample[sample.keys()[0]]
    return lamp        

