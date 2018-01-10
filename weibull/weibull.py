import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import scipy.stats
from scipy.special import gamma

rcParams.update({'figure.autolayout': True})

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN)


class ParameterError(Exception):
    def __init__(self, *args):
        default_str = 'Values for "beta" and "eta" not found; Run the "fit" method or assign values explicitly.'
        super().__init__(default_str, *args)


# convenience functions
def _weibull_ticks(y, _):
    return "{:.0f}%".format(100 * (1 - np.exp(-np.exp(y))))


def _ftolnln(f):
    return np.log(-np.log(1.0 - np.asarray(f)))


class Analysis:
    r"""
    Calculates and plots data points and curves for a standard 2-parameter Weibull for analyzing life data.

    :param data: A list or numpy array of life data, i.e. ``[127, 234, 329, 444]``
    :param suspended: A list or numpy array of suspensions as boolean values, i.e. ``[False, False, True, True]``. At any point which indicates ``True`` means that the test was stopped - or that the item was removed from the test - before the item failed.
    :param unit: The unit ('hour', 'minute', 'cycle', etc.).  This is used to add some useful information to the visualizations.  For instance, if the unit is ``hour``, then the x-axis will be labed in hours.

    :ivar beta: The current value of the shape parameter, :math:`\beta`.  This value is initially set to ``None``.  The proper value for ``beta`` will be calculated on call to the ``fit()`` method.  The user may also set this value directly.
    :ivar eta: The current value of the scale parameter, :math:`\eta`. This value is initially set to ``None``.  The proper value for ``beta`` will be calculated on call to the ``fit()`` method.  The user may also set this value directly.
    :ivar fit_test: Basic statistics regarding the results of ``fit()``, such as :math:`R^2` and P-value.
    """

    def __init__(self, data: list, suspended: list=None, unit: str='cycle'):

        self.x_unit = unit
        self.fit_test = None

        self.beta, self.eta = None, None

        dat = pd.DataFrame({'data': data})
        dat.index = np.arange(1, len(dat) + 1)

        # a suspension is when a unit is removed from test before it has failed
        if not suspended:
            dat['susp'] = [False if x else True for x in data]
            dat['data'].fillna(dat['data'].max(), inplace=True)
        else:
            dat['susp'] = suspended

        if dat['susp'].all():
            raise ValueError('data must contain at least one observed event')

        dat.sort_values('data', inplace=True)
        dat['rank'] = np.arange(1, len(dat) + 1)
        dat['f_rank'] = np.nan

        dat.loc[dat['susp'] == False, 'f_rank'] = np.arange(1,
                                                            len(dat[dat['susp'] == False]) + 1)
        di = dat['susp'] == False
        dat.loc[di, 'med_rank'] = self._med_ra(dat.loc[di, 'f_rank'])
        dat['reverse_rank'] = dat['rank'].values[::-1]

        self.data = dat
        logger.debug('\n{}'.format(self.data))

        self._calc_adjrank()

    def _calc_adjrank(self):
        dat = self.data
        dat['adj_rank'] = np.nan
        fdat = dat[dat['susp'] == False]
        N = len(fdat)
        padj = [0]
        for i in range(N):
            n = fdat.index[i]
            pn = (fdat.loc[n, 'reverse_rank'] * padj[-1] +
                  (len(dat) + 1.)) / (fdat.loc[n, 'reverse_rank'] + 1)
            padj.append(pn)
            dat.loc[n, 'adj_rank'] = pn

        dat['adjm_rank'] = self._med_ra(dat['adj_rank'])

    def _med_ra(self, i):
        """Calculate median rank using Bernard's approximation."""
        i = np.asarray(i)
        med_rank = (i - 0.3) / (len(i) + 0.4)

        return med_rank

    def _linear_regression(self):
        r"""
        Calculate :math:`\beta` and :math:`\eta` using a curve fit of the supplied data.

        :return: None
        """
        x0 = np.log(self.data.dropna()['data'].values)
        y = _ftolnln(self.data.dropna()['adjm_rank'])

        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y, x0)

        beta = 1.0/slope
        x_intercept = - intercept / beta
        eta = np.exp(-x_intercept/slope)

        self.beta = beta
        self.eta = eta

        logger.debug('beta: {:.2f}, eta: {:.2f}'.format(self.beta, self.eta))

        self.fit_test = pd.Series({'r_squared': r_value ** 2, 'p_value': p_value, 'fit method': 'linear regression'})

    def _maximum_likelihood_estimation(self):
        r"""
        Calculate :math:`\beta` and :math:`\eta` using the maximum likelihood estimation method.

        :return: None
        """
        data = self.data[['data', 'susp']].copy()
        df_failed = data[data.susp == False].copy()

        dtf_failed = df_failed["data"].values

        df_failed["ln_x_div_r"] = df_failed.apply(lambda s: np.log(s['data'])/len(df_failed), axis=1)

        dtf_all = self.data['data'].values

        # use Newton-Rhapson method for estimating the shape parameter

        # give initial value for the shape paramter:
        shape = (((6.0 / np.pi ** 2)
                 * (np.sum(np.log(dtf_all) ** 2)
                 - ((np.sum(np.log(dtf_all))) ** 2) / dtf_all.size))
                 / (dtf_all.size - 1)) ** -0.5

        # 10 iterations of the newton-rhapson method
        for i in range(1, 11):
            a = np.sum(np.log(dtf_failed) * 1.0) / dtf_failed.size
            b = np.sum(dtf_all ** shape)
            c = np.sum((dtf_all ** shape) * np.log(dtf_all))
            h = np.sum((dtf_all ** shape) * (np.log(dtf_all)) ** 2)

            shape = shape + (a + (1.0 / shape) - (c / b)) / ((1.0 / shape ** 2) + ((b * h) - c ** 2) / b ** 2)

        scale = (np.sum((dtf_all ** shape) / len(df_failed))) ** (1 / shape)

        self.beta = shape
        self.eta = scale

        self.fit_test = pd.Series({'fit method': 'maximum likelihood estimation'})

    def _confidence(self, confidence=0.95):
        r"""
        Calculate confidence intervals for :math:`\beta` and :math:`\eta` using the Fisher Matrix method.

        :return: None
        """
        # following the procedure as shown on page 54 of Weibull Analysis by Brian Dodson
        data = self.data[['data', 'susp']].copy().sort_values('susp')
        uncensored = data[data['susp'] == False]
        censored = data[data['susp'] == True]

        # step 3
        def calc(t):
            first_term = self.beta / self.eta ** 2
            second_term = ((t/self.eta) ** self.beta) * (self.beta / self.eta ** 2) * (self.beta + 1)

            return first_term - second_term

        data['step3'] = uncensored['data'].apply(func=calc)

        def calc(t):
            first_term = -1.0 / (self.beta ** 2)
            second_term = ((t / self.eta) ** self.beta) * (np.log(t / self.eta) ** 2)

            return first_term - second_term

        data['step4'] = uncensored['data'].apply(func=calc)

        def calc(t):
            first_term = -1.0 / self.eta
            second_term = ((t / self.eta) ** self.beta) * (1.0 / self.eta) * (self.beta * np.log(t / self.eta) + 1.0)

            return first_term + second_term

        data['step5'] = uncensored['data'].apply(func=calc)

        def calc(t):
            return -((t / self.eta) ** self.beta) * (self.beta / (self.eta ** 2)) * (self.beta + 1.0)

        data['step6'] = censored['data'].apply(func=calc)

        def calc(t):
            return -((t / self.eta) ** self.beta) * (np.log(t / self.eta) ** 2)

        data['step7'] = censored['data'].apply(func=calc)

        def calc(t):
            return ((t / self.eta) ** self.beta) * (1.0 / self.eta) * ((self.beta * np.log(t / self.eta)) + 1.0)

        data['step8'] = censored['data'].apply(func=calc)

        f11 = -np.sum(data['step3']) - np.sum(data['step6'])
        f12 = -np.sum(data['step5']) - np.sum(data['step8'])
        f22 = -np.sum(data['step4']) - np.sum(data['step7'])

        f = np.matrix([[f11, f12], [f12, f22]])
        fprime = np.linalg.inv(f)

        nd = scipy.stats.norm
        k_index = (1.0 - confidence)/2 + confidence
        k = nd.ppf(k_index)

        beta_lower = self.beta / (np.e ** (k * np.sqrt(fprime[1, 1]) / self.beta))
        beta_upper = self.beta * np.e ** (k * np.sqrt(fprime[1, 1]) / self.beta)

        eta_lower = self.eta / (np.e ** (k * np.sqrt(fprime[0, 0]) / self.eta))
        eta_upper = self.eta * np.e ** (k * np.sqrt(fprime[0, 0]) / self.eta)

        self.fit_test['confidence'] = confidence
        self.fit_test['beta lower limit'] = beta_lower
        self.fit_test['beta nominal'] = self.beta
        self.fit_test['beta upper limit'] = beta_upper
        self.fit_test['eta lower limit'] = eta_lower
        self.fit_test['eta nominal'] = self.eta
        self.fit_test['eta upper limit'] = eta_upper

    def fit(self, method: str='lr', confidence_level: float=0.95):
        r"""
        Calculate :math:`\beta` and :math:`\eta` using a linear regression
        or using the maximum likelihood method, depending on the 'method' value.

        :method: 'lr' for linear estimation or 'mle' for maximum likelihood estimation
        :return: None
        """
        if method not in ['lr', 'mle']:
            raise ValueError('The method specified must be '
                             'linear regression "lr" or maximum '
                             'likelihood estimation "mle"')

        if method is 'lr':
            if len(self.data) >= 15:
                logger.warning('the maximum likelihood method is likely '
                               'to yield better results with {} data points'.format(len(self.data)))

            self._linear_regression()
        elif method is 'mle':
            if len(self.data) < 15:
                logger.warning('the linear regression method is likely '
                               'to yield better results with {} data points'.format(len(self.data)))

            self._maximum_likelihood_estimation()

        self._confidence(confidence_level)

    def probplot(self, show: bool=True, file_name: str=None, **kwargs):
        r"""
        Generate a probability plot.  Use this to show the data points plotted with
        the beta and eta values.

        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :param kwargs: valid matplotlib options
        :return: None
        """
        if not self.eta or not self.beta:
            raise ParameterError

        susp = any(self.data['susp'])

        if susp:
            plt.semilogx(self.data['data'], _ftolnln(self.data['adjm_rank']), 'o')
        else:
            plt.semilogx(self.data['data'], _ftolnln(self.data['med_rank']), 'o')

        # calculate the ideal x and y values
        x_ideal = self.eta * np.random.weibull(self.beta, size=1000)
        x_ideal.sort()
        f = 1 - np.exp(-(x_ideal / self.eta) ** self.beta)
        x_ideal = x_ideal[f > 0.01]  # take f > 1%
        f = 1 - np.exp(-(x_ideal / self.eta) ** self.beta)
        x_ideal = x_ideal[f < 0.99]  # take f < 99%
        f = f[f < 0.99]
        y_ideal = np.log(-np.log(1 - f))

        plt.semilogx(x_ideal, y_ideal,
                     label="beta: {:.02f}\neta: {:.01f}".format(self.beta,
                                                                self.eta))
        plt.title("Weibull Probability Plot")
        plt.xlabel('{}s'.format(self.x_unit))
        plt.ylabel('Accumulated failures per {}'.format(self.x_unit))
        plt.legend(loc='lower right')

        # Generate ticks
        def weibull_CDF(y, _):
            return '{:.0f}%'.format((100 * (1 - np.exp(-np.exp(y)))))

        ax = plt.gca()
        formatter = mpl.ticker.FuncFormatter(weibull_CDF)
        ax.yaxis.set_major_formatter(formatter)

        yt_F = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                         0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
        yt_lnF = np.log(-np.log(1 - yt_F))
        plt.yticks(yt_lnF)
        ax.yaxis.grid()
        ax.xaxis.grid(which='both')

        if file_name:
            plt.savefig(file_name)

        if show:
            plt.show()

        return

    def pdf(self, show: bool=True, file_name: str=None):
        r"""
        Plot the probability density function

        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :return: None
        """
        if not self.eta or not self.beta:
            raise ParameterError

        x = np.linspace(0.01, self.eta*5, 100)
        y = scipy.stats.weibull_min.pdf(x, self.beta, 0, self.eta)

        self._plot_prob(x, y,
                        show=show, file_name=file_name,
                        title='Probability Density Function',
                        y_label='probability/{}'.format(self.x_unit))

    def sf(self, show: bool=True, file_name: str=None):
        r"""
        Plot the survival function

        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :return: None
        """
        if not self.eta or not self.beta:
            raise ParameterError

        x = np.linspace(0.01, self.eta * 5, 1000)
        y = scipy.stats.weibull_min.sf(x, self.beta, 0, self.eta)

        upper_upper = scipy.stats.weibull_min.sf(x,
                                                 self.fit_test['beta upper limit'],
                                                 0,
                                                 self.fit_test['eta upper limit'])
        lower_upper = scipy.stats.weibull_min.sf(x,
                                                 self.fit_test['beta lower limit'],
                                                 0,
                                                 self.fit_test['eta upper limit'])
        lower_lower = scipy.stats.weibull_min.sf(x,
                                                 self.fit_test['beta lower limit'],
                                                 0,
                                                 self.fit_test['eta lower limit'])
        upper_lower = scipy.stats.weibull_min.sf(x,
                                                 self.fit_test['beta upper limit'],
                                                 0,
                                                 self.fit_test['eta lower limit'])

        min_y = np.minimum(y, upper_upper)
        min_y = np.minimum(min_y, lower_upper)
        min_y = np.minimum(min_y, lower_lower)
        min_y = np.minimum(min_y, upper_lower)

        max_y = np.maximum(y, upper_upper)
        max_y = np.maximum(max_y, lower_upper)
        max_y = np.maximum(max_y, lower_lower)
        max_y = np.maximum(max_y, upper_lower)

        self._plot_prob(x, y, min_y, max_y,
                        show=show, file_name=file_name,
                        title='Survival Function',
                        y_label='probability of survival')

    def hazard(self, show: bool=True, file_name: str=None):
        r"""
        Plot the hazard (CDF) function

        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :return: None
        """
        self.cdf(show, file_name)

    def cdf(self, show: bool=True, file_name: str=None):
        r"""
        Plot the cumulative distribution function

        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :return: None
        """
        if not self.eta or not self.beta:
            raise ParameterError

        x = np.linspace(0.01, self.eta * 5, 1000)
        y = scipy.stats.weibull_min.cdf(x, self.beta, 0, self.eta)

        upper_upper = scipy.stats.weibull_min.cdf(x,
                                                  self.fit_test['beta upper limit'],
                                                  0,
                                                  self.fit_test['eta upper limit'])
        lower_upper = scipy.stats.weibull_min.cdf(x,
                                                  self.fit_test['beta lower limit'],
                                                  0,
                                                  self.fit_test['eta upper limit'])
        lower_lower = scipy.stats.weibull_min.cdf(x,
                                                  self.fit_test['beta lower limit'],
                                                  0,
                                                  self.fit_test['eta lower limit'])
        upper_lower = scipy.stats.weibull_min.cdf(x,
                                                  self.fit_test['beta upper limit'],
                                                  0,
                                                  self.fit_test['eta lower limit'])

        min_y = np.minimum(y, upper_upper)
        min_y = np.minimum(min_y, lower_upper)
        min_y = np.minimum(min_y, lower_lower)
        min_y = np.minimum(min_y, upper_lower)

        max_y = np.maximum(y, upper_upper)
        max_y = np.maximum(max_y, lower_upper)
        max_y = np.maximum(max_y, lower_lower)
        max_y = np.maximum(max_y, upper_lower)

        self._plot_prob(x, y, min_y, max_y,
                        show, file_name,
                        title='Hazard Function',
                        y_label='probability of failure')

    def fr(self, show: bool=True, file_name: str=None):
        r"""
        Plot failure rate as a function of cycles

        :param show: True if the item is to be shown now, False if other elements to be added later
        :param file_name: if file_name is stated, then the probplot will be saved as a PNG
        :return: None
        """
        if not self.eta or not self.beta:
            raise ParameterError

        x = np.linspace(0.01, self.eta * 5, 1000)
        y = (self.beta / self.eta) * (x / self.eta) ** (self.beta - 1)

        beta = self.fit_test['beta upper limit']
        eta = self.fit_test['eta upper limit']
        upper_upper = (beta / eta) * (x / eta) ** (beta - 1)

        beta = self.fit_test['beta lower limit']
        eta = self.fit_test['eta upper limit']
        lower_upper = (beta / eta) * (x / eta) ** (beta - 1)

        beta = self.fit_test['beta lower limit']
        eta = self.fit_test['eta lower limit']
        lower_lower = (beta / eta) * (x / eta) ** (beta - 1)

        beta = self.fit_test['beta upper limit']
        eta = self.fit_test['eta lower limit']
        upper_lower = (beta / eta) * (x / eta) ** (beta - 1)

        min_y = np.minimum(y, upper_upper)
        min_y = np.minimum(min_y, lower_upper)
        min_y = np.minimum(min_y, lower_lower)
        min_y = np.minimum(min_y, upper_lower)

        max_y = np.maximum(y, upper_upper)
        max_y = np.maximum(max_y, lower_upper)
        max_y = np.maximum(max_y, lower_lower)
        max_y = np.maximum(max_y, upper_lower)

        self._plot_prob(x, y, min_y, max_y,
                        show=show, file_name=file_name,
                        title='Failure Rate',
                        y_label='failures/{}'.format(self.x_unit))

    def _plot_prob(self, x: list, y: list,
                   min_y: list=None, max_y: list=None,
                   show: bool=True, file_name: str=None, title: str=None, y_label: str='probability'):
        r"""
        Base plot function used for the density function plotting

        :param x: the x values
        :param y: the y values
        :param min_y: the minimum y values (used to shade confidence limits)
        :param max_y: the maximum y values (used to shade confidence limits)
        :param show: True if the plot is to be shown, false if otherwise
        :param file_name: the file name to be passed to ``matplotlib.pyplot.savefig``
        :param title: the plot title
        :param y_label: the y-axis label
        :return: None
        """
        if min_y is not None and max_y is not None:
            if len(min_y) > 0 and len(max_y) > 0:
                plt.fill_between(x, min_y, max_y, alpha=0.25)
        plt.plot(x, y, label='beta: {:.02f}\neta: {:.01f}'.format(self.beta,
                                                                  self.eta))

        plt.legend()
        plt.xlim(0)
        plt.ylim(0)

        plt.xlabel('{}s'.format(self.x_unit))
        plt.ylabel(y_label)

        ax = plt.gca()
        ax.grid(True, which='both')

        if title:
            plt.title(title)

        if file_name:
            plt.savefig(file_name)

        if show:
            plt.show()

    def b(self, percent_failed: (float, str)=10.0):
        r"""
        Calculate the B-life value

        :param percent_failed: the number of elements that have failed as a percent (i.e. 10)
        :return: the life in cycles/hours/etc.
        """
        if not self.eta or not self.beta:
            raise ParameterError

        pf = float(percent_failed)

        if not 0.1 <= pf <= 99.0:
            raise ValueError('portion_failed must be between 0.001 and 0.999 (inclusive)')

        return scipy.stats.weibull_min.ppf(pf / 100, self.beta, 0, self.eta)

    @property
    def mean(self):
        r"""
        Calculates and returns mean life (aka, the MTTF) is the integral of the reliability function between 0 and inf,

        .. math::
            MTTF = \eta \Gamma(\frac{1}{\beta} + 1)

        where gamma function, :math:`\Gamma`, is evaluated at :math:`\frac{1}{\beta+1}`

        :return: the mean life of the product
        """
        if not self.eta or not self.beta:
            raise ParameterError

        return self.eta * gamma(1.0/self.beta + 1)

    @property
    def mttf(self):
        r"""
        Calculates and returns mean time between failures (MTTF)

        :return: the mean time to failure
        """
        if not self.eta or not self.beta:
            raise ParameterError

        return self.mean

    @property
    def median(self):
        r"""
        Calculates and returns median life of the product

        :return: The median life
        """
        if not self.eta or not self.beta:
            raise ParameterError

        return scipy.stats.weibull_min.ppf(0.5, self.beta, 0, self.eta)

    @property
    def characteristic_life(self):
        r"""
        Returns the current characteristic life of the product, aka :math:`\eta`

        :return: the characteristic life of the product
        """
        if not self.eta or not self.beta:
            raise ParameterError

        return self.eta


class Design:
    """
    Will determine the required test time required given the number of units
    under test and the target cycles OR it will determine the number of units
    given the test time and the target cycles.

    :param target_cycles: The target number of cycles/minutes/hours
    :param reliability: The fraction of units still running after target_cycles, 0.001 to 0.999
    :param confidence_level: The fractional level of confidence, 0.001 to 0.999
    :param expected_beta: The anticipated level of beta - often worse-case - based on historical data or other assumptions
    """

    def __init__(self, target_cycles: (int, float),
                 reliability: float=0.9, confidence_level: float=0.9,
                 expected_beta: float=2.0):
        if not 0.001 <= reliability <= 0.999:
            raise ValueError('The reliability must be between 0.01 and 0.99')
        if not 0.001 <= confidence_level <= 0.999:
            raise ValueError('The confidence level must be between 0.01 and 0.99')

        self.target_cycles = target_cycles
        self.reliability = reliability
        self.confidence_level = confidence_level
        self.beta = expected_beta

    def num_of_units(self, test_cycles: (int, float)):
        """
        Design a test, calculating the number of units required to run for the test duration / cycles in order to prove the reliability at target_cycles.

        :return: The number of units required
        """

        b = -np.log(self.reliability)
        c = b ** (1.0 / self.beta)

        ee = self.target_cycles / c

        units = np.log(1.0 - self.confidence_level) / (-(test_cycles / ee) ** self.beta)

        return units

    def num_of_cycles(self, number_of_units: int):
        """
        Design a test, calculating the test duration/cycles to prove the required reliability at target_cycles.

        :return: the required duration or cycles
        """

        b = -np.log(self.reliability)
        c = b ** (1.0 / self.beta)

        ee = self.target_cycles / c

        cycles = (-np.log((1.0 - self.confidence_level) ** (1.0 / number_of_units))) ** (1.0 / self.beta) * ee
        return cycles


class Weibayes:
    """
    Weibayes-style analysis of the data with a confidence level and beta.

    :param data: The data for each unit
    :param confidence_level: The fractional level of confidence, 0.001 to 0.999
    :param beta: The shape parameter
    """
    def __init__(self, data: list, confidence_level: float=None, beta: float=2.0):
        if not 0.001 < confidence_level < 0.999:
            raise ValueError('confidence level must be between 0.01 and 0.99')

        self.data = np.asarray(data)

        self.beta = np.float(beta)
        self.confidence_level, self.r = None, None
        self.blife = None

        self._set_confidence_level(confidence_level)

    def _set_confidence_level(self, confidence_level):
        cl = np.float(confidence_level)
        alpha = 1.0 - cl
        r = -np.log(alpha)

        self.confidence_level = cl
        self.r = r

        self._calc()
        self._calc_icdf()
        self._calc_cdf()

    def _calc(self):
        etaseries = np.empty((1, len(self.data)))

        etaseries[0, :] = ((self.data ** self.beta) / self.r)

        self.etaseries = etaseries
        self.eta = etaseries.sum(1) ** (1 / self.beta)

    def _calc_cdf(self):
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

    def _calc_icdf(self):
        """
        Calculates the inverse cumulative distribution function
        :return: None
        """
        self.icdf_x = np.arange(.0001, .99, .0001)
        self.icdf = np.empty((len(self.eta), len(self.icdf_x)))

        tmp = pd.DataFrame(index=self.icdf_x)
        self.icdf[0, :] = self.eta * np.log(1.0 / (1.0 - self.icdf_x)) ** (1.0 / self.beta)
        tmp[self.confidence_level] = self.icdf[0]

        self.blife = tmp.T  # transpose

        self.blife.index.name = 'B'

    def plot(self, confidence_level: float=None, show: bool=True, file_name: str=None):
        """
        Plot the linear plot line.

        :confidence_level: the desired confidence level
        :show: True if the plot is to be shown
        :file_name: Save the plot as "file_name"
        """
        if confidence_level:
            self._set_confidence_level(confidence_level)

        plt.semilogx(self.cdf_x, _ftolnln(self.cdf[0]))
        axis = plt.gca()

        axis.grid(True, which='both')

        formatter = mpl.ticker.FuncFormatter(_weibull_ticks)
        axis.yaxis.set_major_formatter(formatter)
        yt_F = np.array([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                         0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
        yt_lnF = _ftolnln(yt_F)
        plt.yticks(yt_lnF)

        plt.ylim(yt_lnF[1], yt_lnF[-1])
        plt.xlim(self.cdf_x.min(), self.cdf_x.max())

        self._plot_annotate()

        plt.ylabel('failure rate')
        plt.xlabel('cycles')

        if file_name:
            plt.savefig(file_name)
        if show:
            plt.show()

    def _plot_annotate(self):
        ax = plt.gca()
        ax.text(0.02, 0.95, 'beta: {:.0f}'.format(self.beta), transform=ax.transAxes)

        ax.text(.02, .90,
                'eta: {:.03g}'.format(self.eta[0]),
                transform=ax.transAxes)

        ax.text(.02, .85,
                'confidence level: {}'.format(self.confidence_level),
                transform=ax.transAxes)

    def b(self, b_spec: int=10, confidence_level: float=None):
        """
        Calculates the B-life

        :param b_spec: the B-specification (for instance, '10')
        :param confidence_level: the confidence level (usually between 0.01 and 0.99)
        :return: the B life
        """
        if not 1 <= b_spec <= 99:
            raise ValueError('b_spec must be between 1 and 99 (inclusive)')
        if confidence_level and not 0.001 < confidence_level < 0.999:
            raise ValueError('confidence level must be between 0.01 and 0.99')

        if confidence_level:
            self._set_confidence_level(confidence_level)

        b_spec_decimal = b_spec / 100.0
        return float(self.blife[b_spec_decimal].T)

