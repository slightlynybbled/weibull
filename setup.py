from setuptools import setup, find_packages
import os

# provide correct path for version
__version__ = None
here = os.path.dirname(os.path.dirname(__file__))
exec(open(os.path.join(here, 'weibull/version.py')).read())

requirements = [
    'pandas >= 0.20.0',
    'numpy >= 1.0',
    'matplotlib >= 2.0',
    'scipy >= 1.0.0'
]

setup_requirements = [
    'flake8 >= 3.5.0',
    'pytest >= 1.4.0'
]

setup(
    name='weibull',
    version=__version__,
    description='Weibull analysis and test design for reliability and life applications',
    author='Jason R. Jones',
    author_email='slightlynybbled@gmail.com',
    url='https://github.com/slightlynybbled/weibull',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    setup_requires=setup_requirements,
    zip_safe=True,
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English'
    ],
    keywords='weibull reliability'
)