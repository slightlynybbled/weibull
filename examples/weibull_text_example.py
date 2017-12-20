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

analysis.probplot()
analysis.hazard(file_name='hazard.png')

print(f'beta: {analysis.beta}\teta: {analysis.eta}')
print(f'B2 life: {analysis.b(2):.02f}\nB10 life: {analysis.b(10):.02f}\nB50 life: {analysis.b(50):.02f}')
print(f'median: {analysis.median}')
print(f'mean: {analysis.mean}')
print(f'mttf: {analysis.mttf}')
