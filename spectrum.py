import os
import re
import sys
import math
import numpy
import json

class Spectrum(object):
    """
    Python implementation of ACP core spectrum class
    """
    def __init__(self):
        self.start = 370
        self.end = 730 
        self.data = []
        for i in range(self.start,self.end+1):
            self.data.append(0.5)

    def sample(self, wavelen):
        """
        return the power for a given wavelength
        """
        if wavelen >= self.start and wavelen <= self.end:
            return self.data[wavelen-self.start]
        elif wavelen > self.end:
            return self.data[-1]
        else: 
            return self.data[0]

    def to_file(self, whandle):
        """
        serialization
        """
        whandle.write(json.dumps(self.__dict__))

    def from_file(self, rhandle):
        """
        construction from file 
        """
        temp = rhandle.read()
        obj = json.loads(temp)
        self.start =    obj['start']
        self.end =      obj['end']
        self.data =     obj['data']

    def flood_fill(self):
        dist = 0
        nm = self.start
        while  self.data[nm - self.start] < 0.0:
            nm=nm+1
        interval_start = nm
        nm=nm+1;
        while  self.data[nm - self.start] < 0.0:
            nm=nm+1
        interval_end = nm

        last_val = 0.0001
        next_val = 0.0001
        interval_len = interval_end - interval_start

        for nm in range(self.start,self.end+1):
            if self.data[nm - self.start] < 0.0:
                inv_weight = float(dist)/float(interval_len)
                weight = 1.0 - inv_weight
                self.data[nm - self.start] = last_val * weight + next_val * inv_weight
                dist=dist+1
            else:
                last_val = self.data[nm - self.start]
                index = (nm - self.start) + interval_len
                if index < (self.end - self.start):
                    next_val = self.data[index]
                dist = 0


    def import_csv(self, fname):
        """
        parser for certain csv archives
        """
        handle = open(fname, 'r')
        lines = handle.readlines()
        handle.close()
        start = False
        self.fill(-1.0)
        for line in lines:
            if 'wavelength' not in line: 
                if not start:
                    continue
                else: 
                    nm, power = re.split(',', line )
                    if int(nm) >= 370 and int(nm) <= 730:
                        self.data[int(nm) - self.start] = float(power)
            else:
                start=True

        self.flood_fill()

        

    def normalize(self):
        """
        make the power sum to unity
        """
        pow = self.power();
        if pow > 0.0:
            for nm in range(self.start, self.end+1):
                self.data[nm - self.start] = self.data[nm - self.start]/pow

    def scale(self, scalar):
        """
        scale every value by scalar
        """
        for nm in range(self.start, self.end+1):
            self.data[nm - self.start] = self.data[nm - self.start] * scalar

    def max(self):
        """
        return the power of the highest band
        """
        m_max = 0.0
        for i in self.data:
            if i > m_max:
                m_max = i
        return m_max        
         
    def power(self):
        """
        return the sum of the powers at all wavelengths
        """
        accum = 0.0
        for nm in range(self.start, self.end+1):
            accum = accum + self.data[nm - self.start]
        return accum

    def zero(self):
        """
        set to zero
        """
        for nm in range(self.start, self.end+1):
            self.data[nm - self.start] = 0.0

    def fill(self, val):
        """
        set to val 
        """
        for nm in range(self.start, self.end+1):
            self.data[nm - self.start] = val 

    def add_gaussian(self, center, amplitude, sigma):
        """
        create a lobe centered at a given wavelength,
        with a given amplitude (heigiht) and sigma (width)
        """
        e = 2.71828
        for nm in range(self.start, self.end+1):
            self.data[nm - self.start] += ( 
                amplitude * math.pow( e, -(math.pow(
                float(nm-center), 2.0)/(2.0 * math.pow(sigma,2.0)))))

    def trapezoid_bin_10nm(self):
        """
        An SSL SSI-specific subsampling of spectral data
        """
        result = []
        for n in range(self.start, self.end+1, 10):
            result.append(
            self.sample(n-5)/2 + self.sample(n-4) + self.sample(n-3) + 
            self.sample(n-2) + self.sample(n-1) + self.sample(n) + 
            self.sample(n+1) + self.sample(n+2) + self.sample(n+3) +
            self.sample(n+4) + self.sample(n+5)/2)
        return result 

def mult4(  a1, a2,  a3, a4,  result):
    """
    convenience for obtaining the product of four spectra
    """
    for nm in range(a1.start, a1.end +1):
        result.data[nm - a1.start] =  ( 
            a1.data[nm - a1.start] * 
            a2.data[nm - a1.start] *
            a3.data[nm - a1.start] *
            a4.data[nm - a1.start])

def mult2(  a1, a2, result):
    """
    multiply two spectra
    """
    for nm in range(a1.start, a1.end +1):
        result.data[nm - a1.start] =  ( 
            a1.data[nm - a1.start] * 
            a2.data[nm - a1.start] )

def normdiff2(a1, a2, result):
    """
    compute (A - B)/A
    """
    for nm in range(a1.start, a1.end +1):
        result.data[nm - a1.start] =  (( 
            a1.data[nm - a1.start] - 
            a2.data[nm - a1.start] )/a1.data[nm - a1.start])

def diff2(a1, a2, result):
    """
    compute A - B
    """
    for nm in range(a1.start, a1.end +1):
        result.data[nm - a1.start] =  ( 
            a1.data[nm - a1.start] - 
            a2.data[nm - a1.start] )

def div2(  a1, a2, result):
    """
    compute A/B
    """
    for nm in range(a1.start, a1.end +1):
        result.data[nm - a1.start] =  ( 
            a1.data[nm - a1.start] / 
            a2.data[nm - a1.start] )
