#!/usr/bin/env python

"""
    test_illumimant.py

    Test program to determine to apply the Solid-State Lighting Committee's 
    Spectral Similarity Index. 

    Jack Holm's written description of the computation [since slightly revised]:
    The following is a description of the steps in the weighted SF metric:

        1. The spectral power measurement data are binned into 10 nm samples with a
        half-width of 10 nm from 380 to 670 nm. At the meeting we discussed trying
        both tetrahedral and triangular binning functions. For example 1 nm
        measurements would be binned as follows:
        a. tetrahedral: p[n-5]/2 + p[n-4] + p[n-3] + p[n-2] + p[n-1] + p[n] + p[n+1]
        + p[n+2] + p[n+3] + p[n+4] + p[n+5]/2
        b. triangular: p[n-9]/10 + p[n-8]/5 + 3 p[n-7]/10 + 2 p[n-6]/5 + p[n-5]/2 +
        3 p[n-4]/5 + 7 p[n-3]/10 + 4 p[n-2]/5 + 9 p[n-1]/10 + p[n] + 9 p[n+1]/10 + 4
        p[n+2]/5 + 7 p[n+3]/10 + 3 p[n+4]/5 + p[n+5]/2 + 2 p[n+6]/5 + 3 p[n+7]/10 +
        p[n+8]/5 + p[n+9]/10
        where p[n] is the measured power at wavelength n.

        2. Each 10 nm sample is divided by the sum of all the 10 nm samples for the
        illuminant, to normalize the total power to unity. This is done for both the
        reference and test illuminants.

        3. The reference illuminant vector is subtracted from the test illuminant
        vector to create a difference vector.

        4. The difference vector is divided by the normalized reference illuminant
        vector to create a relative difference vector.

        5. The relative difference vector is multiplied by the spectral weighting
        vector which rolls off the short and long wavelength limits similarly to ISO
        7589. The spectral weighting vector is as follows: {12/45, 22/45, 32/45,
        40/45, 44/45, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 11/15, 3/15}

        6. A discrete Fourier transform is applied to the weighted relative
        difference vector to get a vector of Fourier coefficients.

        7. A vector consisting of the first 15 Fourier coefficients is multiplied by
        a frequency weighting vector to get a vector of weighted coefficients. The
        frequency weighting vector is {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3,
        2, 1}

        8. The vector of weighted coefficients is summed to get an error value.

        9. Some function TBD is applied to the error value to get the index value.
        The starting function will be to multiply the error value by a constant and
        subtract it from 100. The constant will be determined as required to get
        reasonable index values based on experience with different illumination
        sources. If necessary a more complex function could be used to convert the
        error values to index values.
"""
# for sys.argv
import sys

# for array operations and FFT
import numpy

# internal Spectrum class
import spectrum

# useful arrays
import constant_arrays

# read the test spectrum from a file
test = spectrum.Spectrum()
testfile = open(sys.argv[1], 'r')
test.from_file(testfile)
testfile.close()

# set the ref to 7589 tungsten
ref_array = constant_arrays.tungsten_7589 

# normalize  ref
p1 = numpy.sum(ref_array)
ref_array = numpy.divide(ref_array, float(p1))

# trapezoid subsample test spectrum
test_array = test.trapezoid_bin_10nm()[0:36]

# normalize  test
p1 = numpy.sum(test_array)
test_array = numpy.divide(test_array, p1)

# compute the normalized difference
diff_array = numpy.subtract(test_array, ref_array)

diff_array = numpy.divide(diff_array, ref_array)

# falloff data
falloff_array = constant_arrays.falloff_array[0:36] 

# multiply normdiff by the falloff spectrum
diff_array = numpy.multiply(diff_array, falloff_array) 

# compute real valued fft
freq_array = numpy.fft.rfft(diff_array)

# create an array of the first 16 fft components
first_sixteen = numpy.absolute(freq_array[0:16])

# multiply coefficients by weights
weighted_freqs = numpy.multiply(first_sixteen, constant_arrays.freq_weights)

# sum 
error = numpy.sum(weighted_freqs)

# voila
print error


