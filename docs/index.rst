.. weibull documentation master file, created by
   sphinx-quickstart on Tue Dec 19 10:10:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to weibull's documentation!
===================================

..  toctree::
    :maxdepth: 2
    :caption: Contents:

    Analysis <analysis.rst>

Introduction
------------

The weibull package is a package designed for easy reliability analysis using the weibull distribution.  This documentation will not make a high effort to explain Weibull analysis but will, instead, focus on the use of the package.

Installation
------------

To install weibull into a Python 3 environment, simply ``pip3 install weibull``.  This package does *not* support Python 2.

Dependencies
************

The ``weibull`` package is built on ``pandas``, ``numpy``, ``matplotlib``, and ``scipy`` libraries.  If you are having trouble installing these libraries, particularly within windows, then you may wish to use the Anaconda distribution of Python.

Test Coverage
-------------

Currently, the ``weibull`` package contains no testing.  This will be introduced soon.

Classes
-------

There are three primary classes available in ``weibull``:

 * ``Analysis``
 * ``Design``
 * ``Weibayes``

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
