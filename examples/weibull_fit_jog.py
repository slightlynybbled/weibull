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
analysis.fit()

analysis.probplot(show=False)

# beta and eta may be directly manipulated to try out different hypotheses
analysis.beta = 2.0
analysis.eta = 5000
analysis.probplot()


print(f'beta: {analysis.beta}\teta: {analysis.eta}')
print(f'B2 life: {analysis.b(2):.02f}\nB10 life: {analysis.b(10):.02f}\nB50 life: {analysis.b(50):.02f}')
print(f'median: {analysis.median}')
print(f'mean: {analysis.mean}')
print(f'mttf: {analysis.mttf}')
