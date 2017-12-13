import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


# convenience functions
def _weibull_ticks(y, pos):
    return "{:.0f}%".format(100 * (1 - np.exp(-np.exp(y))))


def _ftolnln(f):
    return np.log(-np.log(1.0 - np.asarray(f)))


class Analysis:
    """
    Based on life data, calculates a 2-parameter weibull _fit
    """
    def __init__(self, data):
        self._fits = {}

        dat = pd.DataFrame({'data': data})
        dat.index = np.arange(1, len(dat) + 1)

        # a suspension is when a unit is removed from test before it has failed
        dat['susp'] = [False if x else True for x in data]

        dat.sort_values('data', inplace=True)
        dat['rank'] = np.arange(1, len(dat) + 1)
        dat['f_rank'] = np.nan
        dat.loc[dat['susp'] == False, 'f_rank'] = np.arange(1,
                                                            len(dat[dat['susp'] == False]) + 1)
        di = dat['susp'] == False
        dat.loc[di, 'med_rank'] = self.med_ra(dat.loc[di, 'f_rank'])
        dat['rev_rank'] = dat['rank'].values[::-1]

        self.data = dat
        logger.debug('\n{}'.format(self.data))
        self._calc_adjrank()
        self._fit()

        fit = 'syx' if any(dat['susp']) else 'yx'
        logger.info('beta: {:.2f}, eta: {:.2f}'.format(
                    self._fits[fit]['beta'], self._fits[fit]['eta']))

        self.beta = self._fits[fit]['beta']
        self.eta = self._fits[fit]['eta']

    def _calc_adjrank(self):
        dat = self.data
        dat['adj_rank'] = np.nan
        fdat = dat[dat['susp'] == False]
        N = len(fdat)
        padj = [0]
        for i in range(N):
            n = fdat.index[i]
            pn = (fdat.loc[n, 'rev_rank'] * padj[-1] +
                  (len(dat) + 1.)) / (fdat.loc[n, 'rev_rank'] + 1)
            padj.append(pn)
            dat.loc[n, 'adj_rank'] = pn
        dat['adjm_rank'] = self.med_ra(dat['adj_rank'])

    def med_ra(self, i):
        """Calculate median rank using Bernard's approximation."""
        i = np.asarray(i)
        return (i - 0.3) / (len(i) + 0.4)

    def plot(self, file_name=None, **kwargs):
        dat = self.data

        susp = any(dat['susp'])
        fit = 'syx' if susp else 'yx'

        if susp:
            plt.semilogx(dat['data'], _ftolnln(dat['adjm_rank']), 'o')
        else:
            plt.semilogx(dat['data'], _ftolnln(dat['med_rank']), 'o')

        self.plot_fits(fit, **kwargs)

        ax = plt.gca()
        formatter = mpl.ticker.FuncFormatter(_weibull_ticks)
        ax.yaxis.set_major_formatter(formatter)
        yt_F = np.array([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                         0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
        yt_lnF = _ftolnln(yt_F)
        plt.yticks(yt_lnF)

        plt.ylim(_ftolnln([.01, .99]))

        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()

    def _fit(self):
        """
        Fit data.
        
        There are four fits.  X on Y and Y on X for data with no suspensions or
        with suspensions (prefixed by 's').
        """
        x0 = np.log(self.data.dropna()['data'].values)
        X = sm.add_constant(x0)
        Y = _ftolnln(self.data.dropna()['med_rank'])
        model = sm.OLS(Y, X)
        results = model.fit()

        xx = np.logspace(0, np.log(1000), 100, base=np.e)
        XX = sm.add_constant(np.log(xx))
        YY = results.predict(XX)
        eta = np.exp(-results.params[0] / results.params[1])

        self._fits['xy'] = {'results': results, 'model': model,
                           'line': np.row_stack([xx, YY]),
                           'beta': results.params[1],
                           'eta': eta}

        Yx = sm.add_constant(Y)
        model = sm.OLS(x0, Yx)
        results = model.fit()

        yy = _ftolnln(np.linspace(.001, .999, 100))
        YY = sm.add_constant(yy)
        XX = np.exp(results.predict(YY))
        eta = np.exp(results.predict([1, 0]))

        self._fits['yx'] = {'results': results, 'model': model,
                           'line': np.row_stack([XX, yy]),
                           'beta': 1 / results.params[1],
                           'eta': eta[0]}

        x0 = np.log(self.data.dropna()['data'].values)
        X = sm.add_constant(x0)
        Y = _ftolnln(self.data.dropna()['adjm_rank'])
        model = sm.OLS(Y, X)
        results = model.fit()

        xx = np.logspace(0, np.log(1000), 100, base=np.e)
        XX = sm.add_constant(np.log(xx))
        YY = results.predict(XX)
        eta = np.exp(-results.params[0] / results.params[1])

        self._fits['sxy'] = {'results': results, 'model': model,
                            'line': np.row_stack([xx, YY]),
                            'beta': results.params[1],
                            'eta': eta}

        Yx = sm.add_constant(Y)
        model = sm.OLS(x0, Yx)
        results = model.fit()

        YY = sm.add_constant(yy)
        XX = np.exp(results.predict(YY))
        eta = np.exp(results.predict([1, 0]))

        self._fits['syx'] = {'results': results, 'model': model,
                            'line': np.row_stack([XX, yy]),
                            'beta': 1 / results.params[1],
                            'eta': eta[0]}

    def plot_fits(self, fit='syx', **kwargs):
        dat = self._fits[fit]['line']

        plt.plot(dat[0], dat[1], **kwargs)


class Design:
    """
    Will determine the required test time required given the number of units
    under test and the target cycles OR it will determine the number of units
    given the test time and the target cycles.
    """

    def __init__(self, target_cycles,
                 reliability=0.9, confidence_level=0.9, expected_beta=2.0):
        """
        Initializes the Design class
        :param target_cycles: the target number of cycles
        :param reliability: the fraction of units still running after target_cycles
        :param confidence_level: the fractional level of confidence
        :param expected_beta: the anticipated level of beta (often worse-case)
        """
        if not 0.01 <= reliability <= 0.99:
            raise ValueError('The reliability must be between 0.01 and 0.99')
        if not 0.01 <= confidence_level <= 0.99:
            raise ValueError('The confidence level must be between 0.01 and 0.99')

        self.target_cycles = target_cycles
        self.reliability = reliability
        self.confidence_level = confidence_level
        self.beta = expected_beta

    def num_of_units(self, test_cycles):
        return self._calc_num_of_units(test_cycles)

    def num_of_cycles(self, num_of_units):
        return self._calc_test_cycles(num_of_units)

    def _calc_num_of_units(self, test_cycles):
        """
        Design a test, calculating the number of units
        required to run for the test duration / cycles

        :return: number of units required for the test
        """

        b = -np.log(self.reliability)
        c = b ** (1.0 / self.beta)

        ee = self.target_cycles / c

        units = np.log(1.0 - self.confidence_level) / (-(test_cycles / ee) ** self.beta)

        return units

    def _calc_test_cycles(self, number_of_units):
        """
        Design a test, calculating the test duration/cycles
        to prove the required reliability

        :return: the required duration or cycles
        """

        b = -np.log(self.reliability)
        c = b ** (1.0 / self.beta)

        ee = self.target_cycles / c

        cycles = (-np.log((1.0 - self.confidence_level) ** (1.0 / number_of_units))) ** (1.0 / self.beta) * ee
        return cycles


class Weibayes:

    def __init__(self, data, num_of_units=None, confidence_level=None, beta=2.0):
        if not 0.001 < confidence_level < 0.999:
            raise ValueError('confidence level must be between 0.01 and 0.99')

        if num_of_units:
            self.data = np.ones(num_of_units) * data
        else:
            self.data = np.asarray(data)

        self.beta = np.float(beta)
        self.confidence_level, self.r = None, None
        self.blife = None

        self.set_conf(confidence_level)

    def __str__(self):
        return f'weibayes: [eta: {self.eta:.02f}, beta: {self.beta:.02f}, cl: {self.confidence_level}]'

    def __repr__(self):
        return f"weibayes(beta={self.beta:.02f}, cl={self.confidence_level:.02f})"

    def run_calcs(self):
        self.calc()
        self.calc_icdf()
        self.calc_cdf()

    def set_conf(self, cl):
        confidence_levels = [0.5, cl]

        cl = np.asarray(confidence_levels)
        alpha = 1.0 - cl
        r = -np.log(alpha)

        self.confidence_level = cl
        self.r = r

        self.run_calcs()

    def calc(self, r=None):
        etaseries = np.empty((len(self.r), len(self.data)))
        for n, r in enumerate(self.r):
            etaseries[n, :] = ((self.data ** self.beta) / r)
        self.etaseries = etaseries
        self.eta = etaseries.sum(1) ** (1 / self.beta)

    def calc_cdf(self):
        """
        calculates the cumulative distribution function, saves within self.cdf
        :return: None
        """
        tmin = 10 ** (np.floor(np.log10(self.icdf.min())) - 1)
        tmax = 10 ** (np.floor(np.log10(self.icdf.max())) + 1)

        self.cdf_x = np.linspace(tmin, tmax, 1000)
        self.cdf = np.empty((len(self.eta), len(self.cdf_x)))

        for n, eta in enumerate(self.eta):
            self.cdf[n, :] = 1 - np.exp(- (self.cdf_x / eta) ** self.beta)

    def calc_icdf(self):
        """
        calculates the inverse cumulative distribution function
        :return: None
        """
        self.icdf_x = np.arange(.0001, .99, .0001)
        self.icdf = np.empty((len(self.eta), len(self.icdf_x)))

        tmp = pd.DataFrame(index=self.icdf_x)
        for n, eta in enumerate(self.eta):
            self.icdf[n, :] = eta * np.log(1. / (1 - self.icdf_x)) ** (1 / self.beta)
            tmp[self.confidence_level[n]] = self.icdf[n]

        self.blife = tmp.T  # transpose

        self.blife.index.name = 'B'

    def b_value(self, b):
        idxs = self.cdf <= 1
        bi = np.abs(self.cdf[idxs] - np.float(b) / 100.).argmin()
        return self.cdf_x[bi]

    def plot(self):
        for n, i in enumerate(self.confidence_level):
            plt.semilogx(self.cdf_x, _ftolnln(self.cdf[n]))
        ax = plt.gca()

        formatter = mpl.ticker.FuncFormatter(_weibull_ticks)
        ax.yaxis.set_major_formatter(formatter)
        yt_F = np.array([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                         0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
        yt_lnF = _ftolnln(yt_F)
        plt.yticks(yt_lnF)

        plt.ylim(yt_lnF[1], yt_lnF[-1])
        plt.xlim(self.cdf_x.min(), self.cdf_x.max())

        self.plot_annotate()
        plt.show()

        # plt.ylabel('failure rate')
        # plt.xlabel('time')

    def plot_annotate(self, b=None):
        ax = plt.gca()
        plt.text(.02, .95, 'beta: {:.0f}'.format(self.beta),
                 transform=ax.transAxes)

        ff = ["{:.5g}, ", ] * len(self.confidence_level)
        ff = "".join(ff).rstrip(", ")
        plt.text(.02, .85, 'eta: ' + ff.format(*self.eta),
                 transform=ax.transAxes)

        confidence_strings = [str(c) for c in self.confidence_level]
        confidence_string = ', '.join(confidence_strings)

        plt.text(.02, .90, f'cl: {confidence_string}',
                 transform=ax.transAxes)
        if b:
            plt.text(.02, .8, 'B{}: '.format(b) + ff.format(
                *self.blife[b].values.tolist()),
                     transform=ax.transAxes)

    def b_life(self, bs=(0.01, 0.02, 0.05, 0.10)):
        """
        Prints a table that contains the percent survival at confidence
        :param bs:
        :return:
        """
        string = str(self.blife[list(bs)].T)
        return string


# These functions need some work. I think they were supposed to calculation
# multiple confidence intervals at once.
def weibayes_calc(data, n=None, beta=2.0, cl=None):
    wcalcs = {}
    if cl:
        if type(cl) != type([]):
            cl = [cl, ]
        cl.append(50)
    else:
        cl = [50, ]
    for c in cl:
        wcalcs[c] = weibayes(data, N=n, beta=beta, cl=c)
    return wcalcs


def plot_weibayes(wcalcs):
    w = wcalcs.copy()
    w5 = w.pop(50)
    t = 'cdf'
    w5.plot(t, color=cc[0])
    for wi in w:
        w[wi].plot(t, color=cc[0])


def print_weibayes(wcalcs, bs=None):
    if not bs:
        bs = [1, 2, 5, 10]
    blife = pd.DataFrame(index=bs)
    blife.index.name = 'B'
    for w in sorted(wcalcs):
        x = wcalcs[w].blife[bs].T
        blife[w] = x
    print(blife)
