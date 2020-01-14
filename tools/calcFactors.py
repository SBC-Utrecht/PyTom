from functools import reduce

def calc_factors(n):    
    return[e for e in list(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))) if e!=n and e!=1]
    


def calc_prime_factors(n):
    N = n
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)

    return [e for e in factors if e != N]

def calc_fast_gpu_dimensions(start=32, end=1000, allowed_factors=[1,2,3,5,7,11], verbose=False):
    list_fast_dimensions = []
    primes = []
    
    for i in range(2, end):
        factors = calc_prime_factors(i)

        if not factors:
            continue
        
        fast = True

        # Multiple factors of 11 will not be fast so we will exclude these numbers
        if factors.count(11) > 1: 
            fast = False
            if verbose: print(f'factor 11 occured {factors.count(11)} times in {i}, and therefore not selected.')
            continue

        # If one factor does not occur in the list of allowed factors, the number i will be slow.
        for factor in factors:
            if factor == i: continue
            if not factor in allowed_factors:
                fast =False
                break

        # If i can be factorized into small primes, the corresponding fft on GPU will have a high
        # probability to be fast. So i is added to a list of fast z-dimensions.
        if fast:
            if i > start-1:
                list_fast_dimensions += [i]
            if not i in allowed_factors:
                allowed_factors += [i]
            if verbose: print(i,factors)

    return list_fast_dimensions

if __name__=='__main__':
    import sys
    sys.argv += [32,1000] 
    print(calc_fast_gpu_dimensions(int(sys.argv[1]), int(sys.argv[2])+1))
