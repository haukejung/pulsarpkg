import multiprocessing as mp
from .multiprocessing_helper_functions import *
from astropy.io import fits

import functions


def sort_dict_by_key(dictionary):
    """
    sorts a dictionary by key.

    :param dictionary:
    :return:
    """
    keys = []
    values = []
    for k, v in sorted(dictionary.items()):
        keys.append(k)
        values.append(v)
    return keys, values


def arr_normalize_axis(arr, axis=None, mask=None):
    """
    normalizes array with respect to the given axis, making each row/column
    have the same mean within a specified region
    :param arr: the array to be normalized
    :param axis: 'y' or 'x', will specify which axis to take as the chunks
           to be normalized
    :param mask: indicates which components of the given axis are priority.
           if mask is left None or set to an array of all the same value, then all
           points will be considered equal priority. If mask is set to
           ones for the first quarter, then zeroes for the rest, then
           the array will be normalized such that the means of the rows/columns
           will be equal within the region specified by the array. The mask can
           also be other numbers; you could specify a list
           [1/float(i)**2 for i in range(arr_len)], for instance.
    :return:
    """
    if axis is None:
        return arr
    elif axis is 'y':
        if mask is None:
            mask = [1.]*len(arr[0])
        # not used:
        # else:
        #     for e in mask:
        #         e = np.float(e)
        mean = 0.
        for row in arr:
            mean += np.mean(mask * row)
        mean /= float(len(arr))
        for row in arr:
            this_mean = np.mean(mask * row)
            row *= mean / this_mean
        return arr
    elif axis is 'x':
        return arr_normalize_axis(arr.T, axis='y', mask=mask).T
    else:
        raise Exception("invalid axis specification.")


def repl_nonvals_wmed(array):
    """
    Replaces Non-values like NaN (and 0, because it makes trouble) with the median
    Median is computed only with values other than 0 and NaN
    :param array: Numpy array
    :return tuple of the array and the median
    """
    index = np.where(np.logical_or(array < 0., array > 0.))
    median = np.median(array[index])
    for row in range(len(array)):
        for col in range(len(array[row])):
            if array[row][col] == 0 or np.isnan(array[row][col]):
                array[row][col] = median
    return array, median


class Dynamic:
    def __init__(self, data, db_header, filename=None, rotate=False):
        """
        Initialize the dynamic spectrum class with a header and the corresponding dynamic spectrum
        :param data: the raw data
        :param rotate: rotate when handling local fits files
        """
        functions.check_object_type(data, fits.HDUList)

        self.hdu_header = data[0].header
        self.db_header = db_header
        self.filename = filename if filename else self.db_header['filename']
        self.rotate = rotate
        self.dyn = self.get_dynamic_spectrum(data[0].data)

        self.dyni = Indexed2D(data=self.dyn, axes=self.get_dyn_axes())

    def get_dynamic_spectrum(self, data, normalize_frequency=False, normalize_time=True, outliers_sigma=9):
        """
        returns a numpy array containing the dynamic spectrum.
        :param data: plain numpy array
        :param normalize_frequency: if set to True, then will normalize the dynamic spectrum
               to compensate for uneven power in different frequency bins (horizontally stripey)
        :param normalize_time: if set to True, then will normalize the dynamic spectrum
               to compensate for uneven power in different time bins (vertically stripey)
        :param outliers_sigma:
        :return: the dynamic spectrum as a numpy array
        """
        # dynamic = data
        dynamic = np.rot90(data) if self.rotate else data
        # dyn_mean = np.mean(dynamic)
        # dyn_std = np.std(dynamic)
        dynamic, dyn_median = repl_nonvals_wmed(dynamic)

        dyn_med_std = np.std(dynamic - dyn_median)
        # sets values 9 SDs above the mean and values less than 0 to 0
        # index = np.where(np.logical_or(dynamic >= dyn_median + (float(outliers_sigma) * dyn_med_std), dynamic < 0.))
        index = np.where(dynamic >= dyn_median + (float(outliers_sigma) * dyn_med_std))
        dynamic[index] = dyn_median

        if normalize_frequency:
            dynamic = arr_normalize_axis(dynamic, 'y')
        if normalize_time:
            dynamic = arr_normalize_axis(dynamic, 'x')
        return dynamic

    def get_dyn_axes(self):
        """
        finds the axes that correspond to the dynamic spectrum plot.
        :return: a list of tuples ([conjugate frequency axis values],[conjugate time axis values])
        """
        t_int = float(self.hdu_header["T_INT"])  # total integration time
        naxis1 = int(self.hdu_header["NAXIS1"])
        naxis2 = int(self.hdu_header["NAXIS2"])
        nsubs = naxis2 if self.rotate else naxis1  # number of time subintegrations
        freq = float(self.hdu_header["FREQ"])  # centre frequency
        BW = abs(float(self.hdu_header["BW"]))  # bandwidth
        nchans = naxis1 if self.rotate else naxis2  # using NAXIS1 instead of NCHANS, because NCHANS isn't always set

        frequency = list(np.linspace(freq-BW/2, freq+BW/2, nchans))
        time = list(np.linspace(0, t_int, nsubs))
        return frequency, time

    def get_dyn_y_axis(self):
        return self.dyni.y_axis

    def get_dyn_x_axis(self):
        return self.dyni.x_axis


"""
EVERYTHING IS INDEXED FIRST WITH RESPECT TO Y, THEN WITH RESPECT TO X
EVERYTHING IS INDEXED FIRST WITH RESPECT TO Y, THEN WITH RESPECT TO X
EVERYTHING IS INDEXED FIRST WITH RESPECT TO Y, THEN WITH RESPECT TO X
"""


class Indexed2D:
    """
    here's a pretty schweet class that keeps track of both a numpy array and its
    associated axes. It looks long, but it's actually really simple - .data contains
    the 2D numpy array with the data, and .x_axis and .y_axis hold the x and y axes
    respectively (.axes holds both as a tuple). It's useful because it does type checking
    and stuff for you, so instead of lugging around variables like dyn_data dyn_xaxis
    dyn_yaxis all the way through your program, you can just have an Indexed2D object
    that holds all the information for you, guaranteed to work in imshow or something.

    oh, another cool thing about this - you can access the values in the Indexed2D by
    the value of its axes. For instance, if you have an array that has 200 elements from
    -1 to 1 on the x axis, and 100 elements from 0 to 1 on the y axis, you can access the
    value of a point at y=-0.32 and x=0.57 by saying my_Indexed2D[-0.32,0.57]
    """
    def __init__(self, data=None, axes=None, dtype=float):
        self._is_data_set = False
        self._are_axes_set = False
        if data is None:
            self.data = np.array([[]])
        else:
            self.set_data(data, dtype)
        if axes is None:
            self.axes = ([], [])
            self.y_axis = []
            self.x_axis = []
        else:
            self.set_axes(axes)
        return

    def __getitem__(self, tup):
        y = tup[0]
        x = tup[1]
        y_index = self.__get_y_index(y)
        x_index = self.__get_x_index(x)
        return Indexed2D(data=self.data[y_index, x_index], axes=(self.y_axis[y_index], self.x_axis[x_index]))

    def __get_y_index(self, value):
        if type(value) == slice:
            if value.start is not None:
                if value.start < min(self.y_axis) or value.start > max(self.y_axis):
                    raise IndexError('y axis index out of bounds: ' + str(value.start))
                start_index = list(np.absolute([p - value.start for p in self.y_axis]))
                start_index = start_index.index(min(start_index))
            else:
                start_index = 0
            if value.stop is not None:
                if value.stop < min(self.y_axis) or value.stop > max(self.y_axis):
                    raise IndexError('y axis index out of bounds: ' + str(value.stop))
                stop_index = list(np.absolute([p - value.stop for p in self.y_axis]))
                stop_index = stop_index.index(min(stop_index))
            else:
                stop_index = len(self.y_axis) - 1
            return slice(start_index, stop_index + 1)
        else:
            index = [p - value for p in self.y_axis]
            index = index.index(min(index))
            return index

    def __get_x_index(self, value):
        if type(value) == slice:
            if value.start is not None:
                if value.start < min(self.x_axis) or value.start > max(self.x_axis):
                    raise IndexError('x axis index out of bounds: ' + str(value.start))
                start_index = list(np.absolute([p - value.start for p in self.x_axis]))
                start_index = start_index.index(min(start_index))
            else:
                start_index = 0
            if value.stop is not None:
                if value.stop < min(self.x_axis) or value.stop > max(self.x_axis):
                    raise IndexError('x axis index out of bounds: ' + str(value.stop))
                stop_index = list(np.absolute([p - value.stop for p in self.x_axis]))
                stop_index = stop_index.index(min(stop_index))
            else:
                stop_index = len(self.x_axis) - 1
            return slice(start_index, stop_index + 1)
        else:
            index = [p - value for p in self.x_axis]
            index = index.index(min(index))
            return index

    def set_data(self, data, dtype=float):
        if type(data) is not list and type(data) is not np.ndarray:
            raise TypeError('Data does not have the right type.')
        for d in data:
            if len(d) != len(data[0]):
                raise IndexError('Data must be rectangular in shape.')
        for d in data:
            if type(d) is not list and type(d) is not np.ndarray:
                raise TypeError('Data does not have the right type.')
            for i in range(len(d)):
                if d[i] is not dtype:
                    try:
                        d[i] = dtype(d[i])
                    except Exception as e:
                        print('your data could not be casted to ' + str(dtype) + ': ')
                        raise e
        if self._are_axes_set:
            y_axis_matching = len(data) == len(self.y_axis)
            x_axis_matching = len(data[0]) == len(self.x_axis)
            if y_axis_matching and x_axis_matching:
                self._is_data_set = True
                self.data = np.array(data)
            else:
                raise IndexError('Data must have dimensions as axes')
        else:
            self._is_data_set = True
            self.data = np.array(data)
        return

    def set_axes(self, axes):
        if type(axes) is not tuple:
            raise TypeError('Axes argument should be a tuple (y_axis,x_axis)')
        y_axis = axes[0]
        x_axis = axes[1]
        if type(y_axis) is not list and type(y_axis) is not np.ndarray:
            raise TypeError('The axes should be specified with a list or numpy array')
        if type(x_axis) is not list and type(x_axis) is not np.ndarray:
            raise TypeError('The axes should be specified with a list or numpy array')
        for i in range(len(y_axis)):
            if type(y_axis[i]) is not float:
                y_axis[i] = float(y_axis[i])
        for i in range(len(x_axis)):
            if type(x_axis[i]) is not float:
                x_axis[i] = float(x_axis[i])
        if self._is_data_set:
            y_axis_matching = len(y_axis) == len(self.data)
            x_axis_matching = len(x_axis) == len(self.data[0])
            if y_axis_matching and x_axis_matching:
                self._are_axes_set = True
                self.axes = (y_axis, x_axis)
                self.x_axis = x_axis
                self.y_axis = y_axis
            else:
                raise IndexError('Axes must have dimensions as data (y: {0}!={1}, x: {2}!={3})'.format(
                        len(y_axis), len(self.data), len(x_axis), len(self.data[0])))
        else:
            self._are_axes_set = True
            self.x_axis = x_axis
            self.y_axis = y_axis
            self.axes = (y_axis, x_axis)
        return

    def get_data(self):
        return np.array(self.data)

    def get_axes(self):
        return self.axes

    def get_x_axis(self):
        return self.x_axis

    def get_y_axis(self):
        return self.y_axis


# this class constructs, contains, and displays secondary spectra.
class Secondary(Dynamic):  # Secondary inherits the Dynamic class
    def __init__(self, data, db_header, filename=None, rotate=False, hand=None):
        """
        initialize me with an dynamic object
        :param data:
        :param hand:
        :return:
        """
        functions.check_object_type(data, fits.HDUList)

        Dynamic.__init__(self, data, db_header, filename, rotate)
        data = self.get_secondary_spectrum()
        axes = self.get_sec_axes()

        self.sec = Indexed2D(data=data, axes=axes)
        self.hand = hand
        self.made_1D = False
        self.parabola_power = {}
        self.observation_name = filename
        self.band = str(self.observation_name)

    def __getitem__(self, value):
        """
        the secondary object can be accessed like a list - sec[5,5] will return
        :param value:
        :return:
        """
        return self.sec[value]

    def get_secondary_spectrum(self, subtract_secondary_background=True, normalize_frequency=True, normalize_time=True,
                               cut_off_bottom=True, xscale=1., yscale=1.):
        """
        returns a secondary spectrum given a dynamic spectrum.
        :param subtract_secondary_background: whether or not to subtract the background
               to improve contrast. Points without signal will have a mean of 0,
               the points with meaningful data will have a mean higher than 0.
        :param normalize_frequency: whether to try to get rid of stripeyness in the frequency axis
        :param normalize_time: whether to try to get rid of stripeyness in the time axis
        :param cut_off_bottom: whether to cut off the "mirror image" bottom half of the secondary spectrum
        :param xscale: the multiplicative scale by which to cut down x
        :param yscale: the multiplicative scale by which to cut down y
        :return: a numpy array containing the secondary spectrum
        """
        dynamic = self.dyn
        dynamic = dynamic - np.mean(dynamic)
        secondary = np.fft.fftn(dynamic)

        secondary, sec_median = repl_nonvals_wmed(secondary)

        # secondary /= secondary.max()
        secondary = np.abs(np.fft.fftshift(secondary))**2
        secondary = 10. * np.log10(secondary/np.max(secondary))  # in decibels

        if normalize_frequency:
            mask = [1. if i < len(secondary[0]) / 4. or i > 3 * len(secondary[0]) / 4. else 0. for i in
                    range(len(secondary[0]))]
            secondary = arr_normalize_axis(secondary, 'y', mask)
        if normalize_time:
            mask = [1. if i < len(secondary) / 4. else 0. for i in range(len(secondary))]
            secondary = arr_normalize_axis(secondary, 'x', mask)
        if subtract_secondary_background:
            # secondary_background = np.mean(secondary[:len(secondary) / 4][:len(secondary[0]) / 4])
            # secondary = secondary - secondary_background
            nbins = 25
            histSec = np.histogram(secondary, bins=nbins)
            binsize = (np.max(secondary)-np.min(secondary))//nbins
            maxindex = np.where(histSec[0] == np.max(histSec[0]))  # where frequency of occurences is highest
            xVal = int(maxindex[0])  # position of peak in noise
            xValDb = np.min(secondary)+binsize*xVal+3.  # value of peak in noise offset up by 3Db
            index = np.where(secondary < xValDb)
            if index != -1:  # if there are values less than threshold value
                secondary[index] = xValDb

        ysize = secondary.shape[0]
        xsize = secondary.shape[1]

        xmin = int(xsize / 2. - xsize / (2. * xscale))
        xmax = int(xsize / 2. + xsize / (2. * xscale))

        ymin = int(ysize / 2. - ysize / (2. * yscale))

        if cut_off_bottom:
            ymax = int(ysize / 2.)
        else:
            ymax = int(ysize / 2. + ysize / (2. * yscale))

        return secondary[ymin:ymax, xmin:xmax]

    def get_sec_axes(self):
        """
        finds the axes that correspond to the secondary spectrum plot.
        :return: a list of tuples ([conjugate frequency axis values],[conjugate time axis values])
        """
        t_int = float(self.hdu_header["T_INT"])  # total integration time
        naxis1_ = int(self.hdu_header["NAXIS1"])
        naxis2_ = int(self.hdu_header["NAXIS2"])
        naxis2 = naxis2_ if self.rotate else naxis1_
        BW = abs(float(self.hdu_header["BW"]))  # bandwidth
        nchans = naxis1_ if self.rotate else naxis2_

        nyq_t = 1000. / (2. * t_int)  # nyquist frequency for the delay axis of the secondary spectrum
        nyq_f = nchans / (2. * BW)  # nyquist frequency for the fringe frequency axis of the secondary spectrum
        fringe = list(np.linspace(-nyq_t, nyq_t, naxis2))
        delay = list(reversed(np.linspace(0, nyq_f, nchans / 2.)))
        return delay, fringe

    def get(self, value):
        return self.sec.get_data().item(tuple(value))

    def get_y_axis(self):
        """
        gives the y axis of the secondary spectrum
        :return:
        """
        return self.sec.y_axis

    def get_x_axis(self):
        """
        gives the x axis of the secondary spectrum
        :return:
        """
        return self.sec.x_axis

    def crop_percent(self, y_scale, x_scale):
        """
        crops the secondary spectrum by y_scale and x_scale percent
        :param y_scale:
        :param x_scale:
        :return:
        """
        y_scale = float(y_scale)
        x_scale = float(x_scale)
        if y_scale < 0 or x_scale < 0 or y_scale > 1 or x_scale > 1:
            raise ValueError('x_scale and y_scale must be between 0 and 1.')
        y_max = max(self.sec.get_y_axis())
        x_max = max(self.sec.get_x_axis())
        self.sec = self.sec[y_max * y_scale:, -x_max * x_scale:x_max * x_scale]
        return

    def crop(self, y_lim, x_lim):
        """
        crops the secondary spectrum to the tuples specified by x_lim
        and y_lim - i.e. if the secondary spectrum goes from -10 to 10
        on the x axis and 0 to 5 on the y axis, you couls specify
        crop( (-2.5,2.5) , (0,3) ) to crop the secondary spectrum to those ranges.
        :param y_lim:
        :param x_lim:
        :return:
        """
        self.sec = self.sec[float(y_lim[0]):float(y_lim[1]), float(x_lim[0]):float(x_lim[1])]

    def get_sec(self):
        """
        gives sec as a numpy 2D array
        :return:
        """
        return self.sec.get_data()


def distparab(x, y, a):  # Calculate distance away from parabola
    x = -abs(x)
    A = 1.0/(2.0*a*a) - y/a
    B = -x/(2.0*a*a)
    xp, _r2, _r3 = np.roots([1.0, 0.0, A, B])  # problem can be reduced to solving this equation: x^3+Ax+B=0
    yp = a*pow(xp, 2)
    dx = xp - x
    dy = yp - y
    return dx, dy


class Parabola:
    def __init__(self, secondary, xrange, yrange, mask, max_t):
        """
        :param secondary: a secondary object
        :param xrange: range of x-values to consider (x_min, x_max) [mHz]
        :param yrange: range of y-values to consider (y_min, y_max) [us]
        :param mask: data points <= than this many pixels from origin/x-axis/y-axis will be masked
        :param max_t: max parabola thickness
        """
        functions.check_object_type(secondary, Secondary)

        self.sec = secondary
        self.sec_axes = self.sec.get_axes()
        self.sec_arr = self.sec.get_data()
        self.y_orig, self.x_orig = np.shape(self.sec_arr)  # get dimensions from the secondary array
        self.dy = abs(self.sec_axes[0][1]-self.sec_axes[0][0])  # sec_axes[0] is the delay/y-axis
        self.dx = abs(self.sec_axes[1][0][1]-self.sec_axes[1][0][0])
        # If specified x and y limits are outside of actual data, trim them to fit actual data
        self.min_x = xrange[0] if xrange[0] >= min(self.sec_axes[1]) else min(self.sec_axes[1])
        self.max_x = xrange[1] if xrange[1] <= max(self.sec_axes[1]) else max(self.sec_axes[1])
        self.min_y = yrange[0] if yrange[0] >= min(self.sec_axes[0]) else min(self.sec_axes[0])
        self.max_y = yrange[1] if yrange[1] <= max(self.sec_axes[0]) else max(self.sec_axes[0])
        self.mask_o, self.mask_x, self.mask_y = mask
        self.max_t = max_t

        # Check: min_y >= 0
        if self.min_y < 0.0:
            self.min_y = 0.0
            raise ValueError("Minimum y value is lower than 0")

        # Check that max's are more than min's
        if self.max_x < self.min_x or self.max_y < self.min_y:
            raise ValueError("Minimum values cannot be greater than maximum values")

        # Convert mask_o to an ellipse with correct units
        self.mask_ox = self.mask_o * self.dx
        self.mask_oy = self.mask_o * self.dy

        # Round mask to pixel boundary
        self.mask_x = np.floor(self.mask_x) + 0.5
        self.mask_y = np.floor(self.mask_y) + 0.5

        is_dB = True  # indicates power values are log (dB) values
