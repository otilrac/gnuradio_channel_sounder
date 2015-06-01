#!/usr/bin/env python
#
# Copyright 2011-2013 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

from gnuradio import gr, digital, filter
from gnuradio import blocks
from gnuradio import channels
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser
import sys

try:
    import scipy
except ImportError:
    print "Error: could not import scipy (http://www.scipy.org/)"
    sys.exit(1)

try:
    import pylab
except ImportError:
    print "Error: could not import pylab (http://matplotlib.sourceforge.net/)"
    sys.exit(1)

from scipy import fftpack

class example_timing(gr.top_block):
    def __init__(self, N, sps, rolloff, ntaps, bw, noise,
                 foffset, toffset, poffset):
        gr.top_block.__init__(self)

        rrc_taps = filter.firdes.root_raised_cosine(
            sps, sps, 1.0, rolloff, ntaps)

        gain = bw
        nfilts = 32
        rrc_taps_rx = filter.firdes.root_raised_cosine(
            nfilts, sps*nfilts, 1.0, rolloff, ntaps*nfilts)

        data0 = 2.0*scipy.random.randint(0, 2, N) - 1.0
        data0 = scipy.exp(1j*poffset) * data0
        data1 = 2.0*scipy.random.randint(0, 2, N) - 1.0
        data1 = scipy.exp(1j*poffset) * data1

        # MIMO channels
        data01 = data0 + 1j*data1
        data10 = data1 + 1j*data0

        self.src0 = blocks.vector_source_c(data01.tolist(), False)
        self.src1 = blocks.vector_source_c(data10.tolist(), False)
        self.rrc0 = filter.interp_fir_filter_ccf(sps, rrc_taps)
        self.rrc1 = filter.interp_fir_filter_ccf(sps, rrc_taps)
        self.chn0 = channels.channel_model(noise, foffset, toffset)
        self.chn1 = channels.channel_model(noise, foffset, toffset)
        self.off0 = filter.fractional_resampler_cc(0.20, 1.0)
        self.off1 = filter.fractional_resampler_cc(0.20, 1.0)

        self.clk = digital.pfb_clock_sync_ccf(sps, gain, rrc_taps_rx,
                                                  nfilts, nfilts//2, 1)
        self.taps = self.clk.taps()
        self.dtaps = self.clk.diff_taps()

        self.delay = int(scipy.ceil(((len(rrc_taps)-1)/2 +
                                     (len(self.taps[0])-1)/2)/float(sps))) + 1


        self.vsnk_err = blocks.vector_sink_f()
        self.vsnk_rat = blocks.vector_sink_f()
        self.vsnk_phs = blocks.vector_sink_f()

        self.connect((self.clk,1), blocks.null_sink(gr.sizeof_gr_complex))
        self.connect((self.clk,2), self.vsnk_err)
        self.connect((self.clk,3), self.vsnk_rat)
        self.connect((self.clk,4), self.vsnk_phs)


        self.vsnk_src0 = blocks.vector_sink_c()
        self.vsnk_src1 = blocks.vector_sink_c()
        self.vsnk_clk0 = blocks.vector_sink_c()
        self.vsnk_clk1 = blocks.vector_sink_c()

        self.connect(self.src0, self.rrc0, self.chn0, self.off0, (self.clk,0), self.vsnk_clk0)
        self.connect(self.src1, self.rrc1, self.chn1, self.off1, (self.clk,1), self.vsnk_clk1)
        self.connect(self.src0, self.vsnk_src0)
        self.connect(self.src1, self.vsnk_src1)


def main():
    parser = OptionParser(option_class=eng_option, conflict_handler="resolve")
    parser.add_option("-N", "--nsamples", type="int", default=2000,
                      help="Set the number of samples to process [default=%default]")
    parser.add_option("-S", "--sps", type="int", default=4,
                      help="Set the samples per symbol [default=%default]")
    parser.add_option("-r", "--rolloff", type="eng_float", default=0.35,
                      help="Set the rolloff factor [default=%default]")
    parser.add_option("-W", "--bandwidth", type="eng_float", default=2*scipy.pi/100.0,
                      help="Set the loop bandwidth (PFB) or gain (M&M) [default=%default]")
    parser.add_option("-n", "--ntaps", type="int", default=45,
                      help="Set the number of taps in the filters [default=%default]")
    parser.add_option("", "--noise", type="eng_float", default=0.0,
                      help="Set the simulation noise voltage [default=%default]")
    parser.add_option("-f", "--foffset", type="eng_float", default=0.0,
                      help="Set the simulation's normalized frequency offset (in Hz) [default=%default]")
    parser.add_option("-t", "--toffset", type="eng_float", default=1.0,
                      help="Set the simulation's timing offset [default=%default]")
    parser.add_option("-p", "--poffset", type="eng_float", default=0.0,
                      help="Set the simulation's phase offset [default=%default]")
    (options, args) = parser.parse_args ()

    # Adjust N for the interpolation by sps
    options.nsamples = options.nsamples // options.sps

    # Set up the program-under-test
    put = example_timing(options.nsamples, options.sps, options.rolloff,
                         options.ntaps, options.bandwidth, options.noise,
                         options.foffset, options.toffset, options.poffset)
    put.run()

    data_src0 = scipy.array(put.vsnk_src0.data()[20:])
    data_src1 = scipy.array(put.vsnk_src1.data()[20:])
    data_clk0 = scipy.array(put.vsnk_clk0.data()[20:])
    data_clk1 = scipy.array(put.vsnk_clk1.data()[20:])

    data_err = scipy.array(put.vsnk_err.data()[20:])
    data_rat = scipy.array(put.vsnk_rat.data()[20:])
    data_phs = scipy.array(put.vsnk_phs.data()[20:])

    f1 = pylab.figure(1, figsize=(12,10), facecolor='w')

    # Plot the IQ symbols
    s1 = f1.add_subplot(3,2,1)
    s1.plot(data_src0.real, data_src0.imag, "bo")
    s1.plot(data_clk0.real, data_clk0.imag, "ro")
    s1.set_title("IQ 0")
    s1.set_xlabel("Real part")
    s1.set_ylabel("Imag part")
    s1.set_xlim([-2, 2])
    s1.set_ylim([-2, 2])

    # Plot the symbols in time
    delay = put.delay
    m = len(data_clk0.real)
    s2 = f1.add_subplot(3,2,2)
    s2.plot(data_src0.real, "bs", markersize=10, label="Input")
    s2.plot(data_clk0.real[delay:], "ro", label="Recovered")
    s2.set_title("Symbols 0")
    s2.set_xlabel("Samples")
    s2.set_ylabel("Real Part of Signals")
    s2.legend()

    # Plot the IQ symbols
    s1 = f1.add_subplot(3,2,3)
    s1.plot(data_src1.real, data_src1.imag, "bo")
    s1.plot(data_clk1.real, data_clk1.imag, "ro")
    s1.set_title("IQ 1")
    s1.set_xlabel("Real part")
    s1.set_ylabel("Imag part")
    s1.set_xlim([-2, 2])
    s1.set_ylim([-2, 2])

    # Plot the symbols in time
    delay = put.delay
    m = len(data_clk1.real)
    s2 = f1.add_subplot(3,2,4)
    s2.plot(data_src1.real, "bs", markersize=10, label="Input")
    s2.plot(data_clk1.real[delay:], "ro", label="Recovered")
    s2.set_title("Symbols 1")
    s2.set_xlabel("Samples")
    s2.set_ylabel("Real Part of Signals")
    s2.legend()

    # Plot the clock recovery loop's error
    s3 = f1.add_subplot(3,2,5)
    s3.plot(data_err, label="Error")
    s3.plot(data_rat, 'r', label="Update rate")
    s3.set_title("Clock Recovery Loop Error")
    s3.set_xlabel("Samples")
    s3.set_ylabel("Error")
    s3.set_ylim([-0.5, 0.5])
    s3.legend()

    # Plot the clock recovery loop's error
    s4 = f1.add_subplot(3,2,6)
    s4.plot(data_phs)
    s4.set_title("Clock Recovery Loop Filter Phase")
    s4.set_xlabel("Samples")
    s4.set_ylabel("Filter Phase")


    diff_taps = put.dtaps
    ntaps = len(diff_taps[0])
    nfilts = len(diff_taps)
    t = scipy.arange(0, ntaps*nfilts)

    f3 = pylab.figure(3, figsize=(12,10), facecolor='w')
    s31 = f3.add_subplot(2,1,1)
    s32 = f3.add_subplot(2,1,2)
    s31.set_title("Differential Filters")
    s32.set_title("FFT of Differential Filters")

    for i,d in enumerate(diff_taps):
        D = 20.0*scipy.log10(1e-20+abs(fftpack.fftshift(fftpack.fft(d, 10000))))
        s31.plot(t[i::nfilts].real, d, "-o")
        s32.plot(D)
    s32.set_ylim([-120, 10])


    pylab.show()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass

