import time
import numpy as np
import os
from copy import deepcopy
import pymzml
from pymzml.utils.utils import index_gzip
from scipy.interpolate import interp1d
import pymzml.obo
import math
import re
from collections import defaultdict
from scipy.interpolate import griddata
import numpy as np 
from scipy import signal
from bisect import bisect_left

def compute_bin_indices(mzs, intensities, min_mz, bin_width):
    indices = ((mzs - min_mz) / bin_width).astype(np.int32)
    return indices, mzs, intensities

class IndexedScan:
    """
    A class to store an in individual scan with the peaks indexed into bins by m/z.
    """
    def __init__(self, scan, rt, scan_num, min_mz, bin_width):
        self.scan = scan
        self.retention_time = rt
        self.scan_num = scan_num
        self.min_mz = min_mz
        self.bin_width = bin_width
        self.indexed_peaks = self.index_peaks(self.min_mz, self.bin_width)


    def index_peaks(self, min_mz, bin_width):
        mzs = np.array(self.scan[:,0])
        intensities = np.array(self.scan[:,1])
        indices, mzs, intensities = compute_bin_indices(mzs, intensities, min_mz, bin_width)

        indexed_peaks = defaultdict(list)
        for mz, intensity, idx in zip(mzs, intensities, indices):
            indexed_peaks[idx].append([mz, intensity])

        return indexed_peaks

    def get_intensity(self, mz, mz_tolerance):
        """
        Calculates the total intensity within :param mz_tolerance
        around :param mz within the scan.
        """
 
        min_index = math.trunc((mz - mz_tolerance - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tolerance - self.min_mz) / self.bin_width)
        if min_index == max_index:
            sum_intensity = 0
            #sum the intensities of all peaks in the bin within the tolerance of the mz
            if self.indexed_peaks.get(min_index, None) is None:
                return 0
            else:
                bin_peaks = self.indexed_peaks.get(min_index)
                for i in range(len(bin_peaks)):
                    if abs(bin_peaks[i][0] - mz) <= mz_tolerance:
                        sum_intensity += bin_peaks[i][1]
            return sum_intensity
        else:
            #check peaks from both bins..
            sum_intensity = 0
            indices = range(min_index, max_index + 1, 1)
            for i in range(len(indices)):
                if self.indexed_peaks.get(indices[i], None) is None:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(indices[i])
                    for j in range(len(bin_peaks)):
                        if abs(bin_peaks[j][0] - mz) <= mz_tolerance:
                            sum_intensity += bin_peaks[j][1]
            return sum_intensity

    def get_exp_mz(self, mz, mz_tol):
        """
        Takes an m/z and an absolute width and finds the nearest experimental m/z in the scan
        within the given tolerance to the query m/z.

        :param mz: the target m/z
        :param mz_tol: the width (in m/z space) within which
        """
        min_index = math.trunc((mz - mz_tol - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tol - self.min_mz) / self.bin_width)
        if min_index == max_index:
            peaks_sorted = sorted(self.indexed_peaks.get(min_index, None), key=lambda x:x[1], reverse=True)
        else:
            indices = range(min_index, max_index+1, 1)
            peaks_within_tol = []
            for s in indices:
                if self.indexed_peaks.get(s, None) is None:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(s)
                    for p in bin_peaks:
                        if abs(p[0] - mz) <= mz_tol:
                            peaks_within_tol.append(p)
            peaks_sorted = sorted(peaks_within_tol, key=lambda x:x[1], reverse=True)

        if peaks_sorted is None:
            return None
        else:
            return peaks_sorted[0][0]

class IndexedFile:
    """
    A class to store a set of indexed scans originating from a single mzML file.
    """
    def __init__(self):
        self.min_mz = 0
        self.bin_width = 1
        self.indexed_scans = []

    def extract_xic(self, mz, mz_width, rt_range=None):
        """
        Takes an mz and an mz width (and optionally an rt range)
        and produces an eic (scans, retention times, and intensities)
        """
        rts = []
        scans = []
        intensities = []

        t_min = -1
        t_max = -1
        if rt_range is not None:
            t_min = rt_range[0]
            t_max = rt_range[1]

        for i in range(len(self.indexed_scans)):
            if rt_range is not None and (self.indexed_scans[i].retention_time < t_min or self.indexed_scans[i].retention_time > t_max):
                continue
            elif rt_range is None or (t_min < self.indexed_scans[i].retention_time < t_max):
                intensity = self.indexed_scans[i].get_intensity(mz, mz_width)
                scans.append(self.indexed_scans[i].scan_num)
                rts.append(self.indexed_scans[i].retention_time)
                intensities.append(intensity)
        return np.transpose([rts, intensities])

    def get_indexed_spectrum_atRt(self, rt):
        """
        Grabs the spectrum with an RT nearest to the provided rt (minutes)
        """
        index_guess = int(round(len(self.indexed_scans)/2))
        finding_spec = True
        min_error = self.indexed_scans[index_guess].retention_time - rt
        while finding_spec:
            current_rt = self.indexed_scans[index_guess].retention_time
            current_error = current_rt - rt

            if abs(current_error) <= abs(min_error):
                min_error = current_error

                if min_error < 0:
                    index_guess += 1
                if min_error > 0:
                    index_guess -= 1
                if min_error == 0:
                    finding_spec = False

            else:
                finding_spec = False


        return self.indexed_scans[index_guess]


def nonlinear_axis(start, end, res):
    """
    Creates a nonlinear axis with the m/z values spaced with a defined and constant resolution.
    :param start: Minimum m/z value
    :param end: Maximum m/z value
    :param res: Resolution of the axis ( m / delta m)
    :return: One dimensional array of the nonlinear axis.
    """
    axis = []
    i = start
    axis.append(i)
    i += i / float(res)
    while i < end:
        axis.append(i)
        i += i / float(res)
    return np.array(axis)

def mergedata(data1, data2):
    """
    Uses a 1D intepolation to merge the intensity of data2 into the axis of data1
    :param data1: Reference data (N x 2 array)
    :param data2: Data to merge (M x 2 array)
    :return: data2 interpolated into the same axis as data1 (N x 2)
    """
    f = interp1d(data2[:, 0], data2[:, 1], bounds_error=False, fill_value=0)
    intx = data1[:, 0]
    inty = f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat

def gzip_files(mzml_path, out_path):
    """
    Create and indexed gzip mzML file from a plain mzML.
    """
    with open(mzml_path) as fin:
        fin.seek(0, 2)
        max_offset_len = fin.tell()
        max_spec_no = pymzml.run.Reader(mzml_path).get_spectrum_count() + 32

    index_gzip(
        mzml_path, out_path, max_idx=max_spec_no, idx_len=len(str(max_offset_len))
    )


def auto_gzip(mzml_path):
    out_path = mzml_path + ".gz"
    gzip_files(mzml_path, out_path)
    return out_path


def get_data_from_spectrum(spectrum, threshold=-1):
    if spectrum == None:
        return
    impdat = np.transpose([spectrum.mz, spectrum.i])
    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 1] > threshold]
    return impdat


def get_im_data_from_spectrum(spectrum, threshold=-1):
    array_params = spectrum._get_encoding_parameters("raw ion mobility array")
    dtarray = spectrum._decode(*array_params)

    impdat = np.transpose([spectrum.mz, dtarray, spectrum.i])

    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 2] > threshold]
    return impdat


def get_inj_time(spectrum):
    element = spectrum.element
    it = 1
    for child in element.iter():
        if 'name' in child.attrib:
            if child.attrib['name'] == 'ion injection time':
                it = child.attrib['value']
                try:
                    it = float(it)
                except:
                    it = 1
    return it


def search_by_id(obo, id):
    key = "MS:{0}".format(id)
    return_value = ""
    for lookup in obo.lookups:
        if key in lookup:
            if obo.MS_tag_regex.match(key):
                for fn in FIELDNAMES:
                    if fn in lookup[key].keys():
                        return_value += "{0}\n".format(lookup[key][fn])
    return return_value


FIELDNAMES = ["id", "name", "def", "is_a"]

def mergedata2d(x1, y1, x2, y2, z2, method="linear"):
    """
    Uses a 2D interpolation to merge a 2D grid of data from one set of axes to another.

    Requires rectangular grids of all inputs such as might be created by numpy.meshgrid.
    :param x1: Reference x-value grid (N x M)
    :param y1: Reference y-value grid (N x M)
    :param x2: Data x-value grid (N2 x M2)
    :param y2: Data y-value grid (N2 x M2)
    :param z2: Data intensity grdi (N2 x M2)
    :return: z2 interpolated onto x1 y1 2D grid (N x M)
    """
    oldx = np.ravel(x2)
    oldy = np.ravel(y2)
    newx = np.ravel(x1)
    newy = np.ravel(y1)
    zout = griddata(np.transpose([oldx, oldy]), np.ravel(z2), (newx, newy), method=method, fill_value=0)
    return zout

def get_longest_index(datalist):
    lengths = [len(x) for x in datalist]
    return np.argmax(lengths)

def get_resolution_im(data):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    testdata = deepcopy(data)
    testdata = testdata[testdata[:, 2] > 0]
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    b1 = diffs[:,1] > 0
    diffs = diffs[b1]
    resolutions = safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)

def merge_im_spectra(datalist, mzbins=None, type="Integrate"):
    """
    Merge together a list of ion mobility data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    if len(datalist) > 1:
        concat = np.concatenate(datalist)
    else:
        concat = np.array(datalist[0])
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.
    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution_im(datalist[maxlenpos])
        mzaxis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        mzaxis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))

    # For drift time, use just unique drift time values. May need to make this fancier.
    dtaxis = np.sort(np.unique(concat[:, 1]))

    # Create the mesh grid from the new axes
    X, Y = np.meshgrid(mzaxis, dtaxis, indexing="ij")

    template = np.transpose([np.ravel(X), np.ravel(Y), np.ravel(np.zeros_like(X))])
    print("Shape merge axis:", X.shape)
    xbins = deepcopy(mzaxis)
    xbins[1:] -= np.diff(xbins)
    xbins = np.append(xbins, xbins[-1] + np.diff(xbins)[-1])
    ybins = deepcopy(dtaxis)
    ybins[1:] -= np.diff(ybins) / 2.
    ybins = np.append(ybins, ybins[-1] + np.diff(ybins)[-1])
    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = mergedata2d(template[:, 0], template[:, 1], d[:, 0], d[:, 1], d[:, 2])
            elif type == "Integrate":
                newdat, xedges, yedges = np.histogram2d(d[:, 0], d[:, 1], bins=[xbins, ybins], weights=d[:, 2])
            else:
                print("ERROR: unrecognized merge spectra type:", type)
                continue
            template[:, 2] += np.ravel(newdat)
    return template

def get_autocorr_ratio(data):
    # evaluate a 1 dimensional merged array ratio to check if a dataset is centroided
    # If ratio >= 0.5, it is more than likely NOT centroided
    # otherwise it probably is
    corry = signal.fftconvolve(data[:, 1], data[:, 1][::-1], mode='same')

    if np.amax(corry) != 0:
        corry /= np.amax(corry)

    maxindex = np.argmax(corry)
    ratio = corry[maxindex+1]
    return ratio

def isempty(thing):
    """
    Checks wether the thing, a list or array, is empty
    :param thing: Object to check
    :return: Boolean, True if not empty
    """
    try:
        if np.asarray(thing, dtype=object).size == 0 or thing is None:
            out = True
        else:
            out = False
    except (TypeError, ValueError, AttributeError):
        print("Error testing emptiness")
        out = False
    return out

def linear_interpolation(x1, x2, x):
    """
    :param x1:
    :param x2:
    :param x:
    :return: float(x - x1) / float(x2 - x1)
    """
    return float(x - x1) / float(x2 - x1)

def nearest(array, target):
    """
    In a sorted array, quickly find the position of the element closest to the target.
    :param array: Array
    :param target: Value
    :return: np.argmin(np.abs(array - target))
    """
    i = bisect_left(array, target)
    if i <= 0:
        return 0
    elif i >= len(array) - 1:
        return len(array) - 1
    if np.abs(array[i] - target) > np.abs(array[i - 1] - target):
        i -= 1
    return i

def lintegrate(datatop, intx, fastmode=False):
    """
    Linearize x-axis by integration.

    Each intensity value in the old data gets proportionally added into the new x-axis.

    The total sum of the intensity values should be constant.
    :param datatop: Data array
    :param intx: New x-axis for data
    :param fastmode: If True, uses a faster but less accurate method of integration that just picks the nearest point.
    :return: Integration of intensity from original data onto the new x-axis.
    Same shape as the old data but new length.
    """
    length = len(datatop)
    l2 = len(intx)
    inty = np.zeros_like(intx)
    for i in range(0, length):
        x = datatop[i, 0]
        y = datatop[i, 1]
        if intx[0] < x < intx[len(intx) - 1]:
            index = nearest(intx, x)
            if fastmode:
                inty[index] += y
            else:
                if intx[index] == x:
                    inty[index] += y
                elif intx[index] < x and index < l2 - 1:
                    index2 = index + 1
                    interpos = linear_interpolation(intx[index], intx[index2], x)
                    inty[index] += (1 - interpos) * y
                    inty[index2] += interpos * y
                elif intx[index] > x and index > 0:
                    index2 = index - 1
                    interpos = linear_interpolation(intx[index], intx[index2], x)
                    inty[index] += (1 - interpos) * y
                    inty[index2] += interpos * y
    newdat = np.column_stack((intx, inty))
    return newdat

def merge_spectra(datalist, mzbins=None, type="Interpolate"):
    """
    Merge together a list of data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Filter out junk spectra that are empty
    datalist = [x for x in datalist if len(x) > 0]
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    # print "Median Resolution:", resolution
    # axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    # for d in datalist:
    #    print(d)
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.

    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution(datalist[maxlenpos])
        if resolution < 0:
            print("ERROR with auto resolution:", resolution, maxlenpos, datalist[maxlenpos])
            print("Using ABS")
            resolution = np.abs(resolution)
        elif resolution == 0:
            print("ERROR, resolution is 0, using 20000.", maxlenpos, datalist[maxlenpos])
            resolution = 20000

        axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        axis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))
    template = np.transpose([axis, np.zeros_like(axis)])

    print("Length merge axis:", len(template))

    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = mergedata(template, d)
            elif type == "Integrate":
                newdat = lintegrate(d, axis)
            else:
                print("ERROR: unrecognized trdtrmerge spectra type:", type)
                continue

            template[:, 1] += newdat[:, 1]

    # Trying to catch if the data is backwards
    try:
        if template[1, 0] < template[0, 0]:
            template = template[::-1]
    except:
        pass
    return template

def reduce_spectrum_density(spectrum_data, max_points=5000, method="decimate"):
    """
    Reduce the density of a spectrum while preserving important features.
    
    Args:
        spectrum_data (numpy.ndarray): 2D array with columns [m/z, intensity]
        max_points (int): Maximum points to keep
        method (str): "decimate" for every nth point, "binning" for intensity binning
    
    Returns:
        numpy.ndarray: Reduced density spectrum
    """
    if spectrum_data is None or len(spectrum_data) == 0:
        return spectrum_data
    
    if len(spectrum_data) <= max_points:
        return spectrum_data
    
    if method == "decimate":
        # Simple decimation - take every nth point
        reduction_factor = max(1, len(spectrum_data) // max_points)
        reduced_data = spectrum_data[::reduction_factor]
        
    elif method == "binning":
        # Bin the data to reduce points while preserving signal
        mz = spectrum_data[:, 0]
        intensity = spectrum_data[:, 1]
        
        # Create bins
        n_bins = min(max_points, len(spectrum_data) // 2)
        bin_edges = np.linspace(np.min(mz), np.max(mz), n_bins + 1)
        
        # Bin the data (sum intensities in each bin)
        binned_intensity, _ = np.histogram(mz, bins=bin_edges, weights=intensity)
        binned_mz = (bin_edges[:-1] + bin_edges[1:]) / 2  # Bin centers
        
        # Remove zero intensity bins
        non_zero_mask = binned_intensity > 0
        binned_mz = binned_mz[non_zero_mask]
        binned_intensity = binned_intensity[non_zero_mask]
        
        reduced_data = np.column_stack((binned_mz, binned_intensity))
    
    else:
        # Default to decimation
        reduction_factor = max(1, len(spectrum_data) // max_points)
        reduced_data = spectrum_data[::reduction_factor]
    
    print(f"Spectrum density reduction: {len(spectrum_data)} → {len(reduced_data)} points")
    return reduced_data


class Importer:
    def __init__(self, file_path, **kwargs):
        self._file_path = file_path
        self._params = kwargs
        self.ext = os.path.splitext(file_path)[1]
        self.filesize = os.stat(file_path).st_size
        self.scans = None
        self.times = None
        self.levels = None
        self.centroided = False
        self.polarity = "Positive"
        self.scan_range = None # Note, we will use inclusive scan range such that the last scan is self.scan_range[1]
        self.data = None
        self.immsdata = None
        self.cdms_support = False
        self.imms_support = False
        self.chrom_support = False
        self.centroid_threshold = 0.8
        self.thermo_support = False

        self.indexed_file = None


    def get_polarity(self, scan=None):
        return self.polarity

    def get_ms_order(self, scan=1):
        if self.levels is not None:
            index = self.get_scan_index(scan)
            return self.levels[index]
        else:
            return 1

    # list of all scans
    def get_all_scans(self):
        pass

    # averaged scans
    def get_avg_scan(self, scan_range=None, time_range=None):
        pass

    # Single scan data
    def get_single_scan(self, scan):
        pass

    def get_max_time(self):
        return self.times[-1]

    def get_max_scan(self):
        return self.scans[-1]

    def get_scans_from_times(self, time_range):
        mins = self.get_time_scan(time_range[0])
        if time_range[1] - time_range[0] > 0:
            maxs = self.get_time_scan(time_range[1])
            return [mins, maxs]
        else:
            return [mins, mins]

    def get_times_from_scans(self, scan_range):
        mint = self.get_scan_time(scan_range[0])
        if scan_range[1] - scan_range[0] > 1:
            maxt = self.get_scan_time(scan_range[1])
            avgt = (mint + maxt) / 2.
            return [mint, avgt, maxt]
        elif scan_range[1] - scan_range[0] == 0 or scan_range[1] - scan_range[0] == 1:
            maxt = self.get_scan_time(scan_range[1]+1)
            avgt = (mint + maxt) / 2.
            return [mint, avgt, maxt]
        else:
            # Assume it's a single scan
            return [mint, mint, mint]

    def get_tic(self):
        if not self.chrom_support:
            print("TIC data not supported for this file type:", self._file_path)
            raise Exception
        pass

    def index_scans(self, min_mz, bin_width):
        if not self.chrom_support:
            print("Indexing not supported for this file type:", self._file_path)
            raise Exception
        else:
            self.indexed_file = IndexedFile()
            scans = self.get_all_scans()
            for scan_number in self.scans:
                if self.get_ms_order(scan_number) == 1:
                    scan = scans[self.get_scan_index(scan_number)]
                    rt = self.get_scan_time(scan_number)
                    ind_scan = IndexedScan(scan, rt, scan_number, min_mz, bin_width)
                    self.indexed_file.indexed_scans.append(ind_scan)
        return

    def get_eic(self, mass, mz_tol, rt_range=None):
        if not self.chrom_support:
            print("EIC data not supported for this file type:", self._file_path)
            raise Exception
        else:
            if self.indexed_file is None:
                self.index_scans(0, 1)
            return self.indexed_file.extract_xic(mass, mz_tol, rt_range)

    def get_scan_time(self, scan):
        # Find index in self.scans
        index = self.get_scan_index(scan)
        t = self.times[index]
        return t

    def get_scan_index(self, scan):
        if scan in self.scans:
            index = np.where(np.array(self.scans) == scan)[0][0]
        elif scan < self.scans[0]:
            index = 0
        elif scan > self.scans[-1]:
            index = len(self.scans) - 1
        else:
            index = np.searchsorted(self.scans, scan)
        return index

    def get_time_scan(self, time):
        index = np.argmin(np.abs(self.times - time))
        scan = self.scans[index]
        return scan

    def check_centroided(self):
        scan_data = self.get_single_scan(self.scan_range[0])
        ratio = get_autocorr_ratio(scan_data)
        print(f"Autocorr ratio {round(ratio, 3)}")
        if ratio < self.centroid_threshold:
            self.centroided = True
        else:
            self.centroided = False
        return self.centroided

    def scan_range_from_inputs(self, scan_range=None, time_range=None):
        if scan_range is None and time_range is None:
            scan_range = self.scan_range
        elif time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        scan_range = np.array(scan_range)
        if scan_range[0] < np.amin(self.scans):
            scan_range[0] = np.amin(self.scans)
        if scan_range[0] > np.amax(self.scans):
            scan_range[0] = np.amax(self.scans)
        if scan_range[1] > np.amax(self.scans) or scan_range[1] < scan_range[0]:
            scan_range[1] = np.amax(self.scans)

        print("Scan Range:", scan_range)
        return scan_range

    def avg_fast(self, scan_range=None, time_range=None):
        if self.data is None:
            self.get_all_scans()

        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        data = deepcopy(self.data)
        if scan_range is not None:
            startindex = self.get_scan_index(scan_range[0])
            endindex = self.get_scan_index(scan_range[1])
            data = data[startindex:endindex + 1]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans))

        if data is None or isempty(data):
            print("Error: Empty Data Object")
            return None
        elif len(data) > 1:
            try:
                data = merge_spectra(data)
            except Exception as e:
                print("Merge Spectra Error 2:", e)
                print(data)
        else:
            data = data[0]

        return data


    def get_cdms_data(self, scan_range=None):
        if not self.cdms_support:
            print("CDMS data not supported for this file type:", self._file_path)
            raise Exception

        raw_dat = self.get_all_scans()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])
        scans = np.ones(len(mz))  # Update this
        it = np.ones(len(mz))

        data_array = np.transpose([mz, intensity, scans, it])

        return data_array

    def get_imms_avg_scan(self, scan_range=None, time_range=None, mzbins=1):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception
        start_time = time.perf_counter()
        if self.immsdata is None:
            self.get_all_imms_scans()

        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        data = deepcopy(self.immsdata)
        if scan_range is not None:
            startindex = self.get_scan_index(scan_range[0])
            endindex = self.get_scan_index(scan_range[1])
            data = data[startindex:endindex + 1]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans))

        if data is None or isempty(data):
            print("Error: Empty Data Object")
            return None

        # Need to merge to get 2D from sparse array
        data = merge_im_spectra(data, mzbins=mzbins)
        # try:
        #     data = merge_im_spectra(data, mzbins=mzbins)
        # except Exception as e:
        #     concat = np.concatenate(data)
        #     sort = concat[concat[:, 0].argsort()]
        #     data = removeduplicates(sort)
        #     print("Mark2", e)

        print("Import Time:", time.perf_counter() - start_time, "s")
        return data

    def get_imms_scan(self, s):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception
        start_time = time.perf_counter()
        if self.immsdata is None:
            self.get_all_imms_scans()

        index = self.get_scan_index(s)
        data = self.immsdata[index]
        print("Import Time:", time.perf_counter() - start_time, "s")
        return data

    def get_all_imms_scans(self):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception

    def close(self):
        pass

def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)

def safedivide(a, b):
    """
    Divide numpy arrays b from a. Define a/0=0 to avoid zero division error.
    """
    c = deepcopy(b)
    # Ensure 'a' and 'c' are float arrays for division
    a = np.asarray(a, dtype=float)
    c = np.asarray(c, dtype=float)
    
    non_zero_mask = (b != 0)
    c[non_zero_mask] = a[non_zero_mask] / b[non_zero_mask]
    return c

def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    """
    if len(testdata) < 2:
        return 20000 # Default resolution if not enough data points
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)

def nonlinear_axis(start, end, res):
    """
    Creates a nonlinear axis with the m/z values spaced with a defined and constant resolution.
    """
    start = float(start)
    end = float(end)
    res = float(res)
    if res == 0:
        res = 20000  # Avoid division by zero

    # The original implementation built this with a Python while-loop that
    # appended one float at a time (i += i/res). That recurrence is a simple
    # geometric series: i_n = start * (1 + 1/res)**n. For high-resolution data
    # over a wide m/z range this loop could run into the hundreds of thousands
    # or millions of iterations in pure Python, which is a major source of the
    # GUI hanging during averaging. Compute it directly with numpy instead.
    if start <= 0 or end <= start:
        return np.array([start])

    ratio = 1.0 + 1.0 / res
    n_points = int(np.log(end / start) / np.log(ratio)) + 1
    n_points = max(n_points, 1)
    axis = start * np.power(ratio, np.arange(n_points))
    return axis

class MZMLImporter(Importer):
    """
    Imports mzML data files.

    Note: Some mzML files are 1 indexed, so the first scan is 1, not 0. Others are not, which is frustrating...
    """
    def __init__(self, path, *args, **kwargs):
        super().__init__(path, **kwargs)

        self.msrun = pymzml.run.Reader(path) #, build_index_from_scratch=True)
        self.levels = None
        self.init_scans()

        self.cdms_support = True
        self.imms_support = True
        self.chrom_support = True
        print("Reading mzML file:", path)

    def init_scans(self):
        get_time = lambda s: s.scan_time_in_minutes()

        times, scans, levels = [], [], []
        for spectrum in self.msrun:
            level = getattr(spectrum, 'ms_level', None)
            if level is None:
                continue
            try:
                times.append(float(get_time(spectrum)))
                scans.append(spectrum.ID)
                levels.append(level)
            except Exception:
                continue

        self.times = np.array(times, dtype=np.float32)
        self.scans = np.array(scans)
        self.levels = np.array(levels)
        self.scan_range = [int(self.scans.min()), int(self.scans.max())]
        self.reset_reader()

    def get_single_scan(self, scan):
        try:
            data = get_data_from_spectrum(self.msrun[scan])
        except Exception as e:
            print("Error in get_single_scan:", e)
            data = None
        return data

    def avg_safe(self, scan_range=None, time_range=None, max_points=5000):
        """
        Optimized version of avg_safe that reduces data density to prevent GUI freezing.
        
        Args:
            scan_range: Range of scans to average
            time_range: Time range to average
            max_points (int): Maximum points to keep in final spectrum
        
        Returns:
            numpy.ndarray: Averaged spectrum with controlled data density
        """
        scan_range = self.scan_range_from_inputs(scan_range, time_range)
        
        print(f"avg_safe: Processing scan range {scan_range}")
        
        # Collect all spectra first (like merge_spectra does)
        spectra_to_merge = []
        self.reset_reader()
        
        # Process each spectrum in the range
        for spectrum in self.msrun:
            if scan_range[0] <= spectrum.ID <= scan_range[1]:
                try:
                    # Get raw spectrum data
                    data = get_data_from_spectrum(spectrum)
                    
                    if data is not None and len(data) > 0:
                        # Apply data reduction if spectrum is too dense
                        if len(data) > max_points:
                            # Use every nth point to reduce density
                            reduction_factor = max(1, len(data) // max_points)
                            data = data[::reduction_factor]
                            print(f"Reduced spectrum {spectrum.ID}: {len(data) * reduction_factor} → {len(data)} points")
                        
                        spectra_to_merge.append(data)
                        
                except Exception as e:
                    print(f"Error processing spectrum {spectrum.ID}: {e}")
                    continue
        
        self.reset_reader()
        
        if not spectra_to_merge:
            print("No spectra found in range")
            return np.array([]).reshape(0, 2)
        
        print(f"Collected {len(spectra_to_merge)} spectra for averaging")
        
        # Use the existing merge_spectra function which handles density properly
        try:
            merged_result = merge_spectra(spectra_to_merge)
            print(f"Final merged spectrum has {len(merged_result)} points")
            return merged_result
            
        except Exception as e:
            print(f"Error in merge_spectra: {e}")
            # Fallback: simple averaging if merge_spectra fails
            return self._simple_average_fallback(spectra_to_merge, max_points)

    def _simple_average_fallback(self, spectra_list, max_points=5000):
        """
        Fallback method for simple averaging when merge_spectra fails.
        
        Args:
            spectra_list: List of spectrum arrays
            max_points (int): Maximum points in final spectrum
        
        Returns:
            numpy.ndarray: Simply averaged spectrum
        """
        if not spectra_list:
            return np.array([]).reshape(0, 2)
        
        # Find the spectrum with most reasonable length for template
        lengths = [len(spec) for spec in spectra_list]
        median_length = np.median(lengths)
        
        # Choose a spectrum close to median length as template
        best_idx = np.argmin([abs(len(spec) - median_length) for spec in spectra_list])
        template_spectrum = spectra_list[best_idx].copy()
        
        # Reduce template if too dense
        if len(template_spectrum) > max_points:
            reduction_factor = max(1, len(template_spectrum) // max_points)
            template_spectrum = template_spectrum[::reduction_factor]
        
        # Create template for averaging
        mz_axis = template_spectrum[:, 0]
        intensity_sum = np.zeros_like(template_spectrum[:, 1])
        count = 0
        
        # Average all spectra onto the template
        for spectrum in spectra_list:
            if len(spectrum) > 0:
                try:
                    # Interpolate spectrum onto template axis
                    f = interp1d(spectrum[:, 0], spectrum[:, 1], 
                            bounds_error=False, fill_value=0, kind='linear')
                    interpolated_intensity = f(mz_axis)
                    intensity_sum += interpolated_intensity
                    count += 1
                except Exception as e:
                    print(f"Error interpolating spectrum: {e}")
                    continue
        
        if count > 0:
            averaged_intensity = intensity_sum / count
            result = np.column_stack((mz_axis, averaged_intensity))
            print(f"Fallback averaging: {count} spectra → {len(result)} points")
            return result
        else:
            return np.array([]).reshape(0, 2)



    def reset_reader(self):
        self.msrun.close()
        self.msrun = pymzml.run.Reader(self._file_path)

    def get_all_scans(self, threshold=-1):
        self.reset_reader()
        newtimes = []
        newids = []
        self.data = []
        for n, spec in enumerate(self.msrun):
            if spec.ID in self.scans:
                try:
                    impdat = get_data_from_spectrum(spec, threshold=threshold)
                    self.data.append(impdat)
                    newtimes.append(self.times[n])
                    newids.append(spec.ID)
                except Exception as e:
                    print("Failed to get all scans:", e)
        self.times = np.array(newtimes)
        self.scans = np.array(newids)
        self.reset_reader()
        return self.data

    def get_avg_scan(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if self.filesize > 1e9 and self.data is None:
            data = self.avg_safe(scan_range, time_range)
        else:
            data = self.avg_fast(scan_range, time_range)
        return data

    def get_tic(self):
        try:
            tic = self.msrun["TIC"]
            ticdat = np.transpose([tic.time, tic.i])
            if len(ticdat) != len(self.scans):
                print("TIC wrong size. Likely extra non-MS scans", len(ticdat), len(self.scans))
                raise Exception
            if not np.isclose(ticdat[-1,0],  self.times[-1]):
                print("TIC end time not equal to scan end time")
                raise Exception
        except:
            print("Error getting TIC in mzML; trying to make it...")
            tic = []
            self.get_all_scans()
            print("Imported Data. Constructing TIC")
            for d in self.data:
                try:
                    tot = np.sum(d[:, 1])
                except:
                    tot = 0
                tic.append(tot)
            t = self.times
            ticdat = np.transpose([t, tic])
            print("Done")
        return ticdat

    def get_property(self, s, name):
        element = self.msrun[s].element
        it = 1
        for child in element.iter():
            if 'name' in child.attrib:
                if child.attrib['name'] == name:
                    it = child.attrib['value']
                    try:
                        it = float(it)
                    except:
                        it = 1
        return it

    def get_inj_time_array(self):
        self.reset_reader()
        its = []
        for i, s in enumerate(self.scans):
            it = get_inj_time(self.msrun[s])
            try:
                it = float(it)
            except:
                print("Error in scan header:", i, s, it)
                it = 1
            its.append(it)
        self.reset_reader()
        return np.array(its)

    def get_all_imms_scans(self):
        self.reset_reader()
        newtimes = []
        newids = []
        self.immsdata = []
        for n, spec in enumerate(self.msrun):
            if spec.ID in self.scans:
                try:
                    impdat = get_im_data_from_spectrum(spec)
                    self.immsdata.append(impdat)
                    newtimes.append(self.times[n])
                    newids.append(spec.ID)
                except Exception as e:
                    print(f"Error processing spectrum ID {spec.ID}: {e}")

        self.times = np.array(newtimes)
        self.scans = np.array(newids)
        self.reset_reader()
        return self.immsdata

    # grab pol info from single scan. Similar to old impl but use direct idxing
    def get_polarity(self, scan=1):
        if scan == -1:
            scan = self.scans[0]
        # Directly access the scan at the provided index
        try:
            spec = self.msrun[scan]
        except Exception as e:
            try:
                print(e)
                print("Error getting scan", scan, "trying to get scan", self.scans[1])
                scan = int(self.scans[1])
                spec = self.msrun[scan]
            except Exception as e:
                print("Error getting scan", scan, e)
                print("Defaulting to positive polarity")
                return "Positive"  # Default to positive if error occurs
        # to string gets the raw byte xml format of the scan
        comp = ""
        xml = spec.to_string()
        # convert the byte xml to a string
        for i in xml:
            comp += chr(i)
        # if we want to look at the raw xml data
        # look for the string telling us what polarity is
        if "negative scan" in comp:
            print("Polarity: Negative")
            self.polarity = "Negative"
            return "Negative"
        if "positive scan" in comp:
            print("Polarity: Positive")
            self.polarity = "Positive"
            return "Positive"
        return None


    def get_isolation_mz_width(self, s):
        scan = self.msrun[s]
        precursors = scan.selected_precursors
        if len(precursors) == 0:
            return None, None
        else:
            precursor = precursors[0]
            mz = precursor['mz']
            width = 3
            return mz, width

    def get_cdms_data(self, scan_range=None):
        raw_dat = self.get_all_scans(threshold=0)

        it = 1. / self.get_inj_time_array()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        scans = np.concatenate([s * np.ones(len(raw_dat[i])) for i, s in enumerate(self.scans)])
        try:
            intensity = np.concatenate([d[:, 1] * it[i] / 1000. for i, d in enumerate(raw_dat)])
        except Exception as e:
            print("Mark1:", e, it)
            intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])
        try:
            it = np.concatenate([it * np.ones(len(raw_dat[i])) for i, it in enumerate(it)])
        except Exception as e:
            print("Error with injection time correction:", e)

        data_array = np.transpose([mz, intensity, scans, it])
        return data_array


    def close(self):
        self.msrun.close()

if __name__ == "__main__":
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    mzml_file_path = r"C:\Users\suneel\Desktop\folders\IgF1 std.mzML"
    print(f"Loading file: {mzml_file_path}...")

    try:
        importer = MZMLImporter(mzml_file_path)
        print("File loaded successfully.")
    except Exception as e:
        print(f"FATAL ERROR: Could not load file. {e}")
        exit()

    try:
        while True:
            print("\n--- Main Menu ---")
            print("1: Plot Total Ion Chromatogram (TIC)")
            print("2: Plot Mass Spectrum at a specific time")
            print("3. Plot Averaged Mass Spectrum in a time range")
            print("4: Exit")
            choice = input("Enter your choice (1-4): ")

            if choice == '1':
                print("Calculating and plotting TIC...")
                tic_data = importer.get_tic()
                time = tic_data[:, 0]
                intensity = tic_data[:, 1]
                
                plt.figure(figsize=(12, 6))
                plt.plot(time, intensity, color='black')
                plt.title(f"Total Ion Chromatogram (TIC)\n{os.path.basename(mzml_file_path)}")
                plt.xlabel("Retention Time (minutes)")
                plt.ylabel("Intensity")
                plt.grid(True, linestyle='--', alpha=0.6)
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                plt.show()

            elif choice == '2':
                try:
                    user_time = float(input("Enter time point in minutes (e.g., 25.4): "))
                    time_diff = np.abs(importer.times - user_time)
                    closest_idx = time_diff.argmin()
                    actual_time = importer.times[closest_idx]
                    scan_id = importer.scans[closest_idx]
                    
                    print(f"Finding spectrum closest to {user_time} min... Found scan at {actual_time:.3f} min (Scan ID: {scan_id}).")
                    
                    spec_data = importer.get_single_scan(scan_id)
                    if spec_data is not None and len(spec_data) > 0:
                        mz = spec_data[:, 0]
                        intensity = spec_data[:, 1]
                        
                        plt.figure(figsize=(12, 6))
                        plt.plot(mz, intensity, color='navy')
                        plt.title(f"Mass Spectrum at {actual_time:.3f} min")
                        plt.xlabel("m/z")
                        plt.ylabel("Intensity")
                        plt.grid(True, linestyle='--', alpha=0.6)
                        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                        plt.show()
                    else:
                        print("Could not retrieve spectrum data.")

                except ValueError:
                    print("Invalid input. Please enter a number.")
                except Exception as e:
                    print(f"An error occurred: {e}")


            elif choice == '3':
                try:
                    start_time = float(input("Enter start time in minutes: "))
                    end_time = float(input("Enter end time in minutes: "))
                    if start_time >= end_time:
                        print("Error: Start time must be less than end time.")
                        continue
                        
                    print(f"Averaging spectra from {start_time:.2f} to {end_time:.2f} min...")
                    avg_spec_data = importer.get_avg_scan(time_range=[start_time, end_time])
                    
                    if avg_spec_data is not None and len(avg_spec_data) > 0:
                        mz = avg_spec_data[:, 0]
                        intensity = avg_spec_data[:, 1]

                        plt.figure(figsize=(12, 6))
                        plt.plot(mz, intensity, color='darkred')
                        plt.title(f"Averaged Mass Spectrum ({start_time:.2f}-{end_time:.2f} min)")
                        plt.xlabel("m/z")
                        plt.ylabel("Averaged Intensity")
                        plt.grid(True, linestyle='--', alpha=0.6)
                        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                        plt.show()
                    else:
                        print("Could not retrieve averaged spectrum. The time range might be empty.")

                except ValueError:
                    print("Invalid input. Please enter a number.")
                except Exception as e:
                    print(f"An error occurred: {e}")
                    
            # --- Option 4: Exit ---
            elif choice == '4':
                print("Exiting.")
                break
            
            else:
                print("Invalid choice. Please enter a number between 1 and 4.")

    finally:
        importer.close()
        print("File handle closed.")