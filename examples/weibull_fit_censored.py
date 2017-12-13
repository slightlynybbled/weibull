import weibull

# take real data and supply it for the failure times,
# leaving right-censored data as None
fail_times = [None] * 10
fail_times[5] = 461
fail_times[1] = 1444
fail_times[6] = 1444
fail_times[3] = 1783

analysis = weibull.Weibull(fail_times)
analysis.plot()
