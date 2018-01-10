Examples
========

Analysis
--------

Step 1: Determining :math:`\beta` and :math:`\eta` Values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before any suppositions may be gathered, it is appropriate to calculate :math:`\beta` and :math:`\eta` values.  Once we are satisfied that :math:`\beta` and :math:`\eta` match the raw data, we can move on to determining useful life characteristics for the product.

Example 1: Complete Test Data
*****************************

In this example, we will take a complete set of failure data that has no censorship and apply basic weibull analysis tool suite in order to achieve a simple, accurate, and useful analysis.::

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
    analysis.fit(method='mle')

    analysis.probplot()

By examining the probability plot, we can visually determine if the :math:`\beta` and :math:`\eta` are appropriately calculated.

By specifying a file name, the probability plot can be saved to a file ``analysis.probplot(file_name='prob.png')``.  This is optional, of course, and not required.

Example 2: Right-Censored Data
******************************

Often, it is necessary to use only the smallest amount of data in order to calculate the values for :math:`\beta` and :math:`\eta`.  For instance, a long-running test might have 10 units on the test bench, but only 3 of them have failed.  When the data is so small, the default linear regression fit method is probably going to yield better results than the maximum-likelihood estimation::

    current_run_time = 4200.0

    fail_times = [current_run_time] * 10
    fail_times[7] = 1034.5
    fail_times[8] = 2550.9
    fail_times[6] = 3043.4

    suspended = [True, True, True, True, True,
                 False, False, False, True, True]

    analysis = weibull.Analysis(fail_times, suspended=suspended, unit='hour')
    analysis.fit()

    analysis.probplot()

Again, we plot the raw data points against the calculated :math:`\beta` and :math:`\eta` in order to ensure that the linear regression is an appropriate fit for the data.  As more failures occur, more accurate curve fits may be run.
