from scipy.stats import weibull_min
import numpy as np
import weibull

# simulate a complete set of failure times
N = 10

bins = 10
beta = 1.4
scale = 10000

x = scale * np.arange(1, N)
fail_times = weibull_min(beta, scale=scale).rvs(N)

analysis = weibull.Weibull(fail_times)
analysis.plot()

# take real data and supply it for the failure times
fail_times = [None] * 10
fail_times[5] = 461
fail_times[1] = 1444
fail_times[6] = 1444
fail_times[3] = 1783

analysis = weibull.Weibull(fail_times)
analysis.plot()
