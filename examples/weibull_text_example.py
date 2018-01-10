"""
Emulates the example in the book Weibull Analysis by Bryan Dodson

Data is sourced from page 22, table 2.5.  The author's estimation was based
on a probability plot and resulted in a beta of 2.974 and eta of 203.3 vs. this
package's beta of 3.47 with an eta of 187.9 using linear regression and 2.97 and
203.3 using maximum likelihood estimation.
"""

import weibull

fail_times = [
    42.1, 105.9, 151.3, 195.6,
    77.8, 117.0, 157.3, 207.0,
    83.3, 126.9, 163.8, 215.3,
    88.7, 138.7, 177.2, 217.4,
    101.8, 148.9, 194.3, 258.8
]

suspensions = [1, 0, 1, 1,
               0, 0, 0, 0,
               1, 0, 0, 1,
               1, 0, 1, 0,
               0, 0, 1, 1]

# this is where the actual analysis and curve fitting occur
analysis = weibull.Analysis(fail_times, suspensions, unit='hour')
analysis.fit(method='mle', confidence_level=0.6)

print(analysis.stats)

analysis.probplot(file_name='gallery-probplot.png')

analysis.pdf(file_name='gallery-pdf.png')
analysis.hazard(file_name='gallery-hazard.png')
analysis.sf(file_name='gallery-survival.png')
analysis.fr(file_name='gallery-fr.png')
