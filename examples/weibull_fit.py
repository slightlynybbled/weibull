import weibull

# create blank table of failure times at test initialization
units_in_test = 9
fail_times = [None] * units_in_test  # when test is started, there are no failure times

# unit numbers and failure times
fail_times[8] = 6677
fail_times[0] = 8329
fail_times[1] = 8545

analysis = weibull.Weibull(fail_times)
analysis.plot()

