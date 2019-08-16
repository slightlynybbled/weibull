Analysis
========

The Analysis class is the primary class which will provide methods for analyzing your life data. This class is designed to take your data and calculate :math:`\beta` and :math:`\eta` values along with generating any appropriate plots for display of your data.

A typical use case of the ``Analysis`` case is as follows::

    import weibull

    # fail times include no censored data
    fail_times = [
        9402.7,
        6082.4,
        13367.2,
        10644.6,
        8632.0,
        3043.4,
        12860.2,
        1034.5,
        2550.9,
        3637.1
    ]

    # this is where the actual analysis and curve fitting occur
    analysis = weibull.Analysis(fail_times, unit='hour')

Fitting
-------

The ``fit()`` method is used to calculate appropriate :math:`\beta` and :math:`\eta` values, which are then stored into the class instance.  When ``fit()`` is called with no parameters, then the linear regression method of calculation is assumed::

    analysis.fit()

An alternative method is to use the Maximum Likelihood Estimation (MLE) method of fitting :math:`\beta` and :math:`\eta` to the data.  This may be done by specifying that the ``method='mle'``::

    analysis.fit(method='mle')

In many cases, the ``mle`` and ``lr`` methods will yield very similar values for :math:`\beta` and :math:`\eta`, but there are some cases in which one is preferred over the other.  It turns out that linear regression tends to work best for very low sample sizes, usually less than 15 while the maximum likelihood estimator works best for high sample sizes.  In both cases, the ``probplot()`` method should be used to verify that the data is a good fit.

To retrieve the :math:`\beta` and :math:`\eta` values, simply use the instance variables ``beta`` and ``eta``::

    print(f'beta: {analysis.beta: .02f}')
    print(f'eta: {analysis.eta: .02f}')

When using the ``fit()`` method, it is also possible to set the confidence levels

Use the ``stats()`` method to get a ``pandas.Series`` containing most internal estimates::

    $> analysis.stats
    fit method          maximum likelihood estimation
    confidence                                    0.6
    beta lower limit                          2.42828
    beta nominal                              2.97444
    beta upper limit                          3.64344
    eta lower limit                           186.483
    eta nominal                               203.295
    eta upper limit                           221.622
    mean life                                  181.47
    median life                               179.727
    b10 life                                   95.401
    dtype: object

Plotting
--------

One of the most often requested features of such a package is plotting the data, particularly in Jupyter Notebooks.  The ``weibull`` package comes with built-in methods to easily display and save standard plots with one-line methods.

Building on the ``analysis`` instance above, we will examine the probability plot::

    analysis.probplot()

.. image:: images/weibull-fit-10pt.png

There are some settings for the Weibull probability plots that will allow the user to add more information to the plot itself that are useful when the plot will be included in reports, on web pages, or other situations where the supporting information may not always accompany the plot.  These include being able to specify the Analyst's name and/or the Company or Organization name::

    analysis.analyst = 'John Q. Smith'
    analysis.company = 'Smith & Smith Engineering'

It should be noted that if either the Analyst or Company name is especially long, the annotation box (where they are displayed) will expand to fit the text which could obscure some portions of the Weibull plot.  To address this, you can place a line feed ``'\n'`` within either name which will force the text to wrap to a shorter length.  These must be defined before calling the ``probplot`` function.

Also, the user can now provide their own more descriptive title for the plot, use the default title, or suppress the title altogether, as follows::

    analysis.plot_title = 'Weibull Analysis of Golden Widget Failures'

Setting the title to a zero-length string (or not defining one at all) will cause the default title of 'Weibull Probability Plot' to be used.  Setting the title to ``None`` will suppress the title line altogether.  Any other value set for ``analysis.plot_title`` will be displayed on the Weibull plot.  This too must be defined before calling the ``probplot`` function.

By default, the probability plot will now be 12 inch x 8 inch (1200 x 800 pixels) to better display enhancements made to the plot format, but it is also possible to change the size of the plot to a different size or aspect ratio by specifying the ``figsize``::

    plt.figure(figsize=(8, 6))
    analysis.probplot()

where the values of 8 and 6 are in inches (at 100px/inch).  Be aware that generating plots smaller than the default size may result in the Y-Axis labels overlapping one another (especially when analyzing very large data sets), or the annotation box might covering some portion of the Weibull plot or the Characteristic Life line.  You may have to experiment with the size values to get a satisfactory result, and you may also need to use the image resizing features when inserting the image into a document or presentation.

We can also examine a number of other common function plots (only the hazard plot is shown, but the others are along the same line)::

    analysis.pdf()
    analysis.sf()
    analysis.hazard()
    analysis.cdf()

.. image:: images/weibull-hazard-10pt.png

Each of these functions will generate a plot that is suitable for publication or insertion into a Jupyter Notebook.  Again, note that some of these methods - such as ``hazard()`` and ``cdf()`` will produce the same plot with slightly different labeling.

Confidence Levels
^^^^^^^^^^^^^^^^^

Some plots will contain a shaded region which reflects the :ref:`confidence_levels`.  The confidence levels are calculate for :math:`\beta` and :math:`\eta` and the min/max values for :math:`\beta` and :math:`\eta` are explored rather than all possible values.  As a result, the visualizations shown are an approximation of the confidence limits.

Class Documentation
-------------------

.. autoclass:: weibull.Analysis
    :members:
