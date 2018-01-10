import weibull

# the current run time of the test or the
# time that the test was suspended completely
current_run_time = 4200.0

fail_times = [current_run_time] * 10
fail_times[7] = 1034.5
fail_times[8] = 2550.9
fail_times[6] = 3043.4

suspended = [True, True, True, True, True,
             False, False, False, True, True]

analysis = weibull.Analysis(fail_times, suspended=suspended, unit='hour')
analysis.fit(confidence_level=0.6)

analysis.probplot(file_name='weibull-fit.png')  # option to save as an image
analysis.pdf()
analysis.sf()
analysis.hazard()
analysis.cdf()
analysis.fr()

print(analysis.stats)
