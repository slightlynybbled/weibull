Examples
========

Determining :math:`\beta` and :math:`\eta` Values
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
    analysis = weibull.Analysis(fail_times,
                                unit='hour')
    analysis.fit(method='mle')

    analysis.probplot()

In this example, we chose to use the Maximum Likelihood Estimation method of estimating :math:`\beta` and :math:`\eta`, which is shown in the ``analysis.fit(method='mle)`` line.  If the ``fit()`` method were called with no parameters, it would - by default - have used linear regression.

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

    analysis = weibull.Analysis(fail_times,
                                suspended=suspended,
                                unit='hour')
    analysis.fit()

    analysis.probplot()

Again, we plot the raw data points against the calculated :math:`\beta` and :math:`\eta` in order to ensure that the linear regression is an appropriate fit for the data.  As more failures occur, more accurate curve fits may be run.

Life Calculations
^^^^^^^^^^^^^^^^^

Once :math:`\beta` and :math:`\eta` are determined, then they may be utilized to obtain the basic lifetime data that may be utilized for planning.  One common reliability metric is the :ref:`b-life`.  Obtaining a B10 life using the ``analysis`` object is trivial::

    print(f'B10 life: {analysis.b(10):.0f}')

As you can see, simply calling the ``b()`` function with the appropriate number as the parameter will return the B-life based on :math:`\beta` and :math:`\eta`.

Basic Life Statistics
^^^^^^^^^^^^^^^^^^^^^

For user convenience, the ``mean``, ``median``, ``characteristic_life``, and ``mttf`` are defined as attributes of the class and may be called at any time after an initial curve fit.  Note that there is some overlap with other class variables.  For instance, the ``characteristic_life`` happens to be the same thing as ``eta``, but if a customer asks for the characteristic life, then having this available makes the code more readable and correspond more closely to the specification.

Plotting
^^^^^^^^

We can also plot various functions of interest, such as the survival function and hazard functions, amongst others.::

    analysis.pdf()      # probability density function
    analysis.sf()       # survival function
    analysis.hazard()   # hazard function
    analysis.cdf()      # cumulative distribution function
    analysis.fr()       # failure rate

Each of these will generate a plot of the function.  For all plotting methods, if ``file_name`` is specified as a parameter, then the method will save to a file rather than display.  For instance::

    analysis.sf(file_name='survival_function.png')


Test Design
^^^^^^^^^^^

The Design class is to be utilized for two scenarios:

 - determine the required number of units to prove the target reliability given a test cycles/duration
 - determine the required number of cycles/duration to prove the target reliability given a number of units

To begin, first import and instantiate the Designer, which is the utility for the test designer. There are several parameters to consider and all of them are requirements or assumptions that must be entered as parameters for the Designer class:

 - target_cycles - the target to be proven in hours/days/weeks/cycles
 - reliability - defaults to 0.9
 - confidence_level - defaults to 0.95
 - expected_beta - an initial assumption for beta (defaults to 2)

Shown are two example calculations for a target lifetime of 10000 hours with a reliability of 0.9 at a confidence level of 0.5 and beta assumption of 1.5::

    import weibull

    designer = weibull.Design(
        target_cycles=10000,
        reliability=0.9,
        confidence_level=0.90,
        expected_beta=1.5
    )

    # The 'test_cycles' parameter can be in any units.
    # Days, weeks, hours, cycles, etc., so long
    #   as the target unit is consistent
    print(f'Minimum number of units for 10000 hour run:{designer.num_of_units(test_cycles=10000)}')
    print(f'Minimum hours for 20 units: {designer.num_of_cycles(num_of_units=20)}')

Weibayes Analysis
^^^^^^^^^^^^^^^^^

Use Weibayes analysis to assist with designing your test or evaluating reliability within a certain confidence interval based on historical data.

You have a product that needs to be tested to B2 life of 40 million time units with a confidence limit of 95%.  The product had an expected beta of 2 (lots of historical data there).  B2 life is the same as 98% survival.

Using the weibull test `Design` class, we need to run 62 units (the limit of our test rig) for 62 million time units with no failures::

    import weibull

    designer = weibull.Design(
        target_cycles=40e6,
        reliability=0.98,
        confidence_level=0.95,
        expected_beta=2
    )

    print(f'Minimum hours for 62 units: {designer.num_of_cycles(num_of_units=62)}')

Result::

    61860134.45191945

Weibayes analysis on the data would arrive at the same result.::

    import weibull

    # we want N units to run for H hours each
    N = 62
    H = 62.0e6

    run_times_desired = [H] * N
    weibayes = weibull.Weibayes(run_times_desired, confidence_level=0.95, beta=2)

    print(f'B2 life: {weibayes.b(2)}')

Results::

    B2 life: 40090439.86038491

Note that this `B2` matches very closely with `target_cycles` value found in the above iteration of the `Design` class.

We can further plot the data using `weibayes.plot()` resulting in:

.. image:: images/weibayes.png
