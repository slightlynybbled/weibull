import weibull

# the current run time of the test or the
# time that the test was suspended completely
current_run_time = 4200.0

fail_times = [77.8, 101.8, 105.9, 117.0, 126.9, 138.7,
              148.9, 157.3, 163.8, 207.0, 217.4]

suspended = [False] * len(fail_times)

analysis = weibull.Analysis(fail_times, suspended=suspended, unit='hour')
analysis.fit(method='mle', confidence_level=0.6)

analysis.probplot(watermark_text='blah blah')
analysis.pdf()
analysis.sf()
analysis.hazard()
analysis.cdf()
analysis.fr(watermark_text='hooray!')

print(analysis.stats)
