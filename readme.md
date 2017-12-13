# Temporary Disclaimer

I (slightlynybbled) forked this repository - including the readme - and I am currently modifying the library along with the readme to improve useability.  The below contains artifacts from the original readme that were written by the original author.  Where I have modified the library, I have modified the readme to reflect those changes.

As I advance in refactoring the library into a more complete and useable package, I will eventually replace all of the text of the readme with appropriate substitutions and remove this temporary disclaimer.

# Weibull Reliability Analysis

This is a rough collection of Weibull analysis routines.  I make no claim to the accuracy.

Routines are for low sample sizes.  I believe at large sample sizes, maximum likelihood methods are more accurate.  The Weibull fits here are done as Y on X and X on Y regressions - the equivalent to graphing on Weibull paper.  The class can handle right-censored data (data in which the test was stopped before all units have experienced failure).

A class for Weibayes analysis is also included.

Also included is are a few methods to determine test time or number of samples.

I wrote this while working at a manufacturing company.  Before I could polish it up, I left for a completely different job.  I doubt I will use this again, but I wanted to pull the work in progress out of the iPython notebook where it was and stick it in a file.  As a result, there are some duplicate functions in the file.  I couldn't be bothered to clean those up...

# Weibull Basic Fitting

The most fundamental weibull analysis is to calculate the values of beta and eta for a given list of life data.

Censored data can be fit as well as uncensored data.  Four fits are performed:

A basic example is shown here, but more complete examples may be found within the [examples](examples/) directory.

    import weibull
    
    # take real data and supply it for the failure times,
    # leaving right-censored data as None
    fail_times = [None] * 10
    fail_times[7] = 1034.5
    fail_times[8] = 2550.9
    fail_times[6] = 3043.4
    
    analysis = weibull.Analysis(fail_times)
    analysis.plot()
    
    print(f'beta: {analysis.beta}\teta: {analysis.eta}')

![weibull _fit](images/weibull-fit.png)

# Test design

The `Design` class is to be utilized for two scenarios:

 - determine the required number of units to prove the target reliability given a test cycles/duration
 - determine the required number of cycles/duration to prove the target reliability given a number of units
 
To begin, first import and instantiate the `Designer`, which is the utility for the test designer.  There are several parameters to consider and all of them are requirements or assumptions that must be entered as parameters for the `Designer` class:

 - `target_cycles` - the target to be proven in hours/days/weeks/cycles
 - `reliability` - defaults to 0.9
 - `confidence_level` - defaults to 0.95
 - `expected_beta` - an initial assumption for beta (defaults to 2)

Shown are two example calculations for a target lifetime of 10000 hours with a reliability of 0.9 at a confidence level of 0.5 and beta assumption of 1.5:

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
    print(f'Minimum number of units for 10000 hour run: {designer.num_of_units(test_cycles=10000)}')
    print(f'Minimum hours for 20 units: {designer.num_of_cycles(num_of_units=20)}')

# Weibayes

The company I worked for didn't like any failures in our testing.  They used canned software for Weibayes.  The software used a confidence limit of 50% as a default, so that's what a lot of our testing analysis ended up using.  The class here has 50% hardcoded in.  You can changed it, but the default is always there.

Confidence limits are expressed in the equation as `r`, where `r = -ln(cl)`.  A few values are shown in the table below.

 | Confidence   | r     | 
 | ------------ | ----- | 
 | 50%          | 0.693 | 
 | 63%          | 1.0   | 
 | 80%          | 1.61  | 
 | 90%          | 2.3   | 
 | 95%          | 3.0   | 
 | 99%          | 4.6   | 

 We had a product that needed to be tested to B2 life of 40 million time units with a confidence limit of 95%.  The product had an expected beta of 2 (lots of historical data there).  B2 life is the same as 98% survival.  The test design said we needed to run 62 units (the limit of our test rig) for 62 million time units with no failures:

    weibull.weib_t(62, 40e7, r=.98, cl=.95, beta=2)
    # 618601344.51919436

Weibayes analysis on the data would arrive at the same result.

    # 62 units run to 62 million time units
    t = [6.2e7]*62
    ww = weibull.weibayes(t, cl = 95, beta = 2)
    # display plot and lifetime block
    ww.display()
    # annotate with B2 (98% survival) values
    ww.plot_annotate(2)

    # output
    # B             50.0          95.0
    # 1.0   5.878481e+07  2.827654e+07
    # 2.0   8.334501e+07  4.009044e+07
    # 5.0   1.328021e+08  6.388021e+07
    # 10.0  1.903328e+08  9.155351e+07

![weibayes](images/weibayes.png)

# Contributions

Initial work on this repository was done by user [tgray](https://github.com/tgray).  You can still peruse the [original repository](https://github.com/tgray/weibull).
