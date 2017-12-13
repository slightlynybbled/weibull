import weibull

n = weibull.weib_n(48, t=100, r=0.9, cl=0.95, beta=1.5)
t = weibull.weib_t(20, t=100, r=0.9, cl=0.95, beta=1.5)

print(f'original n = {n}')
print(f'original t = {t}')

weibull.Design(number_of_units=20,
               test_cycles=48,
               target_cycles=100,
               reliability=0.9,
               confidence_level=0.95,
               expected_beta=1.5)

