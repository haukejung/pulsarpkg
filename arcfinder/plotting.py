import matplotlib.pyplot as plt
from .multiprocessing_helper_functions import *
from . import computing
from matplotlib.backends.backend_pdf import PdfPages


cmap = 'viridis'  # set default colormap


class Pdf:
    def __init__(self, attr_dict):
        self.pdfp = PdfPages('query.pdf')
        self.info = self.pdfp.infodict()
        title = 'pulsarpkg result'
        if len(attr_dict) > 0:
            title += ' ('
        for attr, value in attr_dict.items():
            title += '{0}: {1}, '.format(attr, value)
        title = title[:-2]
        if len(attr_dict) > 0:
            title += ' )'
        self.info['Title'] = title
        self.info['Author'] = 'pulsarpkg'

    def save(self, fig):
        self.pdfp.savefig(fig)

    def __del__(self):
        self.pdfp.close()


def show_image(showme, axis_y=None, axis_x=None):
    """
    Shows an image with the given X and Y axes.
    :param showme: the 2D array to be shown
    :param axis_y: a 1D array containing the y axis to be used in the plot.
    :param axis_x: a 1D array containing the x axis to be used in the plot.
    :param cmap: matplotlib colormap
    :return:
    """
    if axis_x is None:
        axis_x = [i for i in range(len(showme[0]))]
    if axis_y is None:
        axis_y = [i for i in range(len(showme))]
    (x_min, x_max) = (min(axis_x), max(axis_x))
    (y_min, y_max) = (min(axis_y), max(axis_y))
    fig = plt.figure()
    plt.imshow(showme, aspect='auto', extent=[x_min, x_max, y_min, y_max], cmap=cmap)
    plt.colorbar()
    return fig


def show():
    plt.show()


def show_dyn(dyn_obj: computing.Dynamic, save=False, pdf: Pdf=None):
    fig = show_image(dyn_obj.dyn, dyn_obj.get_y_axis(), dyn_obj.get_x_axis())
    plt.title(dyn_obj.filename)
    plt.xlabel('Time (MJD - {0}) [s]'.format(dyn_obj.hdu_header['MJD']))
    plt.ylabel('Frequency [MHz]')
    if save:
        plt.savefig('{0}_dyn.pdf'.format(dyn_obj.filename), format='pdf')
    if pdf:
        pdf.save(fig)
    plt.close()
    return


def show_sec(sec_obj: computing.Secondary, save=False, pdf=None):
    """
    plots the secondary spectrum to the current figure in matplotlib
    :param sec_obj: "Secondary" object
    :return:
    :return:
    """
    fig = show_image(sec_obj.get_sec(), sec_obj.get_y_axis(), sec_obj.get_x_axis())
    if sec_obj.made_1D:
        overplot_parabolas(sec_obj, [min(sec_obj.etas), max(sec_obj.etas)])
    plt.title(sec_obj.observation_name)
    plt.ylabel('delay')
    plt.xlabel('fringe frequency')
    if save:
        plt.savefig('{0}_sec.pdf'.format(sec_obj.filename), format='pdf')
    if pdf:
        pdf.save(fig)
    return


def overplot_parabolas(sec_obj: computing.Secondary, etas, offsets=None):
    """
    plots parabolas over the secondary spectrum.
    :param sec_obj:
    :param etas: a list of the curvatures of parabolas desired
    :param offsets: a list of the y-offsets desired for the parabolas
    :return: nothing, but plots the parabolas to the current matplotlib figure.
    """
    if offsets is None:
        offsets = [0.]
    for eta in etas:
        for offset in offsets:
            eta = float(eta)
            axis_x = sec_obj.get_x_axis()
            plot_x = [x + offset for x in axis_x]
            axis_y = sec_obj.get_y_axis()
            parab = []
            for x in axis_x:
                y = eta * x ** 2 - eta * offset ** 2
                parab.append(y)
            plt.plot(plot_x, parab, 'b-')
            plt.xlim((min(axis_x), max(axis_x)))
            plt.ylim((min(axis_y), max(axis_y)))


def show_power_vs_eta(self, weird=False):
    """
    shows how much power is present in each eta value.
    Requires make_1D_by_quadratic to have been run, and simply plots
    the results from it to the current figure in matplotlib.
    :param weird:
    :return:
    """
    if not self.made_1D:
        print("make_1D_by_quadratic has not been run yet")
        return

    if not weird:
        plt.plot(self.etas, self.powers)
        plt.xlabel("eta")
        plt.ylabel("Power(dB), arbitrary scaling")
        plt.title("Power vs eta, " + self.observation_name)
        # self.overplot_parabolas([min(sec.etas),max(sec.etas)])
        return
    else:
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2)
        plt.plot([1 / eta ** 2 for eta in self.etas], self.powers)

        x_axis_points = np.linspace(1 / max(self.etas) ** 2, 1 / min(self.etas) ** 2, 10)
        x_axis = [round(1 / np.sqrt(x), 4) for x in x_axis_points]
        plt.xticks(x_axis_points, x_axis, rotation=90)
        plt.xlabel("eta")
        plt.ylabel("Power(dB), arbitrary scaling")
        plt.title("Power vs eta, " + self.observation_name)
        # self.overplot_parabolas([min(self.etas),max(self.etas)])
        return


def __give_eta_list(self, eta_range, num_etas, decimal_places=4):
    if num_etas is not 1:
        x_max = np.sqrt(1 / min(eta_range))
        x_min = np.sqrt(1 / max(eta_range))
        return [1 / x ** 2 for x in np.linspace(x_min, x_max, num_etas)]
    else:
        return [np.average(eta_range)]