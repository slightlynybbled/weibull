import weibull

# the current run time of the test or the
# time that the test was suspended completely
current_run_time = 4200.0

fail_times = [77.8, 101.8, 105.9, 117.0, 126.9, 138.7,
              148.9, 157.3, 163.8, 207.0, 217.4]

suspended = [False] * len(fail_times)

analysis = weibull.Analysis(fail_times, suspended=suspended, unit='hour')
analysis.fit(method='mle')

analysis.probplot()
analysis.pdf()
analysis.sf()
analysis.hazard()
analysis.cdf()
analysis.fr()

print(f'beta: {analysis.beta}\teta: {analysis.eta}')
print(f'{analysis.stats}')
print(f'B2 life: {analysis.b(2):.02f}\nB10 life: {analysis.b(10):.02f}\nB50 life: {analysis.b(50):.02f}')
print(f'median: {analysis.median}')
print(f'mean: {analysis.mean}')
print(f'mttf: {analysis.mttf}')
