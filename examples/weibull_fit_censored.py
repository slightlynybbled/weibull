import weibull

# take real data and supply it for the failure times,
# leaving right-censored data as None
fail_times = [None] * 10
fail_times[7] = 1034.5
fail_times[8] = 2550.9
fail_times[6] = 3043.4

analysis = weibull.Weibull(fail_times)
analysis.plot(file_name='weibull-_fit.png')  # option to save as an image

print(f'beta: {analysis.beta}\teta: {analysis.eta}')
