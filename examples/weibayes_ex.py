import weibull

import matplotlib.pyplot as plt

# 62 units run to 62 million time units
t = [6.2e7] * 62
ww = weibull.Weibayes(t, confidence_level=0.5, beta=2)

# display plot and lifetime block
print(f'B2 life: {ww.b(2)}')
print(f'B10 life: {ww.b(10)}')
ww.plot()

# annotate with B2 (98% survival) values
#ww.plot_annotate(2)

#plt.show()
