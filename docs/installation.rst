Installation
============

Using ``pip``
-------------

Ideally, you should be able to pip install weibull and simply be finished. This package is a pure-python package, so it should work on any os. Unfortunately, this package utilizes other packages which may be more difficult to install. Please consult package documentation for more details.

This user had the most issues installing statsmodels in a Windows 10 environment. Other package dependencies - numpy, pandas, and matplotlib - installed without issue.

Using ``setup.py``
------------------

Additional dependencies are utilized by calling ``python setup.py install`` which ensure that the environment is appropriately set up for development.  For instance, the additional dependencies are ``flake8``, ``pytest``, and ``sphinx``.

Simply download a copy of the repository, ``cd`` into the directory, and ``python setup.py install``.

Using ``conda``
---------------
If you are having installation issues, perhaps try the Anaconda distribution! As I understand it, they have solved most of these installation problems for difficult packages!

Dependencies
------------

For most installations, you must have ``pandas``, ``numpy``, ``matplotlib``, and ``scipy`` properly installed into your environment.
