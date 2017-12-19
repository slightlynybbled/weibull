.. _introduction-to-reliability-analysis:

Introduction to Reliability Analysis
====================================

Weibull Distribution
--------------------

In reliability analysis and, thus, in the ``weibull`` package, we are primarily concerned with the 2-parameter Weibull function defined herein as:

.. math::
  F(x) = \frac{\beta}{\eta} \left(\frac{x}{\eta}\right)^{\beta-1} e^{-\left(x/\eta\right)^\beta}

where:

 - :math:`\beta` or *beta* represents the **shape** parameter
 - :math:`\eta` or *eta* represents the **scale** parameter
 - :math:`x` represents the value at which the function is to be evaluated

Were one to plot the above :math:`F(x)` with given :math:`\beta` and :math:`\eta` values, one would get the probability density function, commonly shortened to PDF.  From the PDF alone, it is possible to derive the cumulative distribution function (a.k.a CDF and hazard functions), along wih the survival function which is very useful in reliability engineering.

Distribution Shape
******************

The **shape** parameter, :math:`\beta`, determines the overall shape of the distribution.  There are three primary regions in which :math:`\beta` may fall:

 - :math:`\beta < 1.0` Indicates infant mortality, or decreasing failures as time increases.  This is a distribution that may be observed when a phenomenon such as adhesive curing exists. As the adhesive cures, the product experiences fewer failures.
 - :math:`\beta = 1.0` Indicates 'random' or 'constant' failures.  This sort of distribution is most commonly applied to some electronic component categories, such as semiconductors.
 - :math:`\beta > 1.0` Indicates a wearout style of distribution.  This distribution is common 
