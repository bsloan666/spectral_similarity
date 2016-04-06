import spectrum

def export(spec, fname):
    handle = open(fname, 'w')
    handle.write("wavelength (nanometers),power (peak normalized)\n")
    for i, val in enumerate(spec.data):
        handle.write("%d, %f\n"%(i+spec.start, val))
    handle.close()    
