"""
This file being used to generate some of the plots for the documentation
"""
import weibull

plotter = weibull.Analysis([0], [False])

# Linearized CDF comparison
plotter.beta = 0.5
plotter.eta = 1.0
plotter.probplot(show=False)

plotter.beta = 1.0
plotter.probplot(show=False)

plotter.beta = 2.0
plotter.probplot(file_name='beta-effects-on-log-cdf.png')

# PDF comparison
plotter.beta = 0.5
plotter.pdf(show=False)

plotter.beta = 1.0
plotter.pdf(show=False)

plotter.beta = 2.0
plotter.pdf(file_name='beta-effects-on-pdf.png')

# failure rate comparison
plotter.beta = 0.5
plotter.fr(show=False)

plotter.beta = 1.0
plotter.fr(show=False)

plotter.beta = 2.0
plotter.fr(show=False)

plotter.beta = 2.4
plotter.fr(file_name='beta-effects-on-fr.png')

# eta effects on pdf
plotter.beta = 2.0
plotter.eta = 0.5
plotter.pdf(show=False)

plotter.eta = 1.0
plotter.pdf(show=False)

plotter.eta = 2.0
plotter.pdf(file_name='eta-effects-on-pdf.png')

# bathtub failure rate
# failure rate comparison
plotter.beta = 0.8
plotter.fr(show=False)

plotter.beta = 1.0
plotter.fr(show=False)

plotter.beta = 10
plotter.eta = 4.0

plotter.fr(file_name='bathtub-components.png')
