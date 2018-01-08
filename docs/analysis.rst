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
    analysis.fit()

When the ``fit()`` method is executed, the ``Analysis`` class will perform a curve-fit operation in which the Weibull distribution is fitted to the given data and the corresponding :math:`\beta` and :math:`\eta` are calculated.

To retrieve the :math:`\beta` and :math:`\eta` values, simply use the instance variables ``beta`` and ``eta``::

    print(f'beta: {analysis.beta: .02f}')
    print(f'eta: {analysis.eta: .02f}')

Useful Class Attributes
-----------------------

For user convenience, the ``mean``, ``median``, ``characteristic_life``, and ``mttf`` are defined as attributes of the class and may be called at any time after an initial curve fit.  Note that there is some overlap with other class variables.  For instance, the ``characteristic_life`` happens to be the same thing as ``eta``, but if a customer asks for the characteristic life, then having this available makes the code more readable and correspond more closely to the specification.

Life Calculations
-----------------

Perhaps one of the most often requested reliability attributes of an item is its B-life.  You can read the B-life as the life of a product until a certain percent of a product has failed.  For instance, it is common in bearings to refer to a B10 or L10, which is the life of a product until 10% of the items have failed.

To access the :ref:`b-life` of the ``Analysis`` class, simply call the ``b()`` method with the number specified as the sole parameter.  For instance::

    print(f'B10 life: {analysis.b(10):.0f}')
    print(f'B20 life: {analysis.b(20):.0f}')

Plotting
--------

One of the most often requested features of such a package is plotting the data, particularly in Jupyter Notebooks.  The ``weibull`` package comes with built-in methods to easily display and save standard plots with one-line methods.

Building on the ``analysis`` instance above, we will examine the probability plot::

    analysis.probplot()

.. image:: images/weibull-fit-10pt.png

We can also examine a number of other common function plots (only the hazard plot is shown, but the others are along the same line).::

    analysis.pdf()
    analysis.sf()
    analysis.hazard()
    analysis.cdf()

.. image:: images/weibull-hazard-10pt.png

Each of these functions will generate a plot that is suitable for publication or insertion into a Jupyter Notebook.  Again, note that some of these methods - such as ``hazard()`` and ``cdf()`` will produce the same plot with slightly different labeling.

Manipulating :math:`\beta` and :math:`\eta`
-------------------------------------------

It is possible to assign :math:`\beta` and :math:`\eta` using normal python commands::

    analysis.beta = 5.4
    analysis.eta = 6050

When these variables are manipulated, then any plotting functionality will assume the new values.  With that in mind, we can actually delay the plotting functionality using ``show=False`` in order to place one or more comparisons on the plot.::

    analysis.probplot(show=False)

    analysis.beta = 2.0  # assign new beta and eta values
    analysis.eta = 5000

    analysis.probplot()

.. image:: images/weibull-prob-plot-comparison.png

Class Documentation
-------------------

.. autoclass:: weibull.Analysis
    :members:
