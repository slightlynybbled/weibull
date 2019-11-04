import warnings

from weibull.version import __version__
from weibull.weibull import Analysis, Design, Weibayes, ParameterError

__all__ = ['__version__',
           'Analysis', 'Design', 'Weibayes', 'ParameterError']

warnings.warn('The "weibull" module is being superceded by the "reliability" module and will no longer be maintained. '
              'Please refactor your code for this alternative model and thank you for using weibull!')
