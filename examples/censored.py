import pandas as pd
import weibull
import matplotlib.pyplot as plt

# create table of end times
units_in_test = 9
end_times = [0] * units_in_test

end_times[8] = 6677
end_times[0] = 8329
end_times[1] = 8545

suspensions = [False if x else True for x in end_times]
end_times = [x if x else max(end_times) for x in end_times]

analysis = weibull.Weibull(end_times, suspensions=suspensions)
analysis.fit()

analysis.plot(susp=1)
analysis.plot_fits('yx', linestyle='--')

plt.show()
