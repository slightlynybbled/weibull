"""
Emulates the example in the book Weibull Analysis by Bryan Dodson

Data is sourced from page 22, table 2.5.  The author's estimation was based
on a probability plot and resulted in a beta of 3.34 and eta of 191 vs. this
package's beta of 3.47 with an eta of 187.9.
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
analysis.fit()

analysis.probplot()

print(f'beta: {analysis.beta}\teta: {analysis.eta}')

"""
The reliability site at
http://reliabilityanalyticstoolkit.appspot.com/weibull_analysis_solution
has beta at 3.34 and eta at 190.3
"""
# a quick tool to format for the above website
for time, suspended in zip(fail_times, suspensions):
    print(f'{time} {"s" if suspended else "f"}')
