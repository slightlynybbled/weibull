import weibull

# we want N units to run for H hours each
N = 62
H = 62.0e6

run_times_desired = [H] * N
weibayes = weibull.Weibayes(run_times_desired, confidence_level=0.95, beta=2)

# show the Bx life calculations
print(f'B2 life: {weibayes.b(2)}')
print(f'B10 life: {weibayes.b(10)}')

# probplot the confidences for each (note that the Bx
# values may be read from the Y axis)
weibayes.plot(file_name='weibayes.png')
weibayes.plot(0.5)

