import weibull

# take real data and supply it for the failure times,
# leaving right-censored data as None
fail_times = [None] * 10
fail_times[7] = 1034.5
fail_times[8] = 2550.9
fail_times[6] = 3043.4

analysis = weibull.Analysis(fail_times)

analysis.probplot(file_name='weibull-fit.png')  # option to save as an image
analysis.pdf()
analysis.sf()
analysis.hazard()
analysis.cdf()

print(f'beta: {analysis.beta}\teta: {analysis.eta}')
print(f'{analysis.fit_test}')
print(f'B2 life: {analysis.b(2):.02f}\nB10 life: {analysis.b(10):.02f}\nB50 life: {analysis.b(50):.02f}')
print(f'median: {analysis.median}')
print(f'mean: {analysis.mean}')
print(f'mttf: {analysis.mttf}')
