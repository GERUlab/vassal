import time

def mytimer(func):
    def f(*args, **kwargs):
        before = time.time()
        result = func(*args, **kwargs)
        after = time.time()
        print 'elapsed time: {0:.3} s'.format(after - before)
        return result
    return f

