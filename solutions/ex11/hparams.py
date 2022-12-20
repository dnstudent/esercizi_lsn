
def range_domain(hparam, step=1):
    """Utility function to extract the domain from a Discrete or IntInterval hparam"""
    try:
        return hparam.domain.values
    except:
        return range(hparam.domain.min_value, hparam.domain.max_value + 1, step)


def hppow(base, hparam, dtype):
    """Elevates base to an exponent in the domain of hparam sampled uniformly"""
    return dtype(pow(base, hparam.domain.sample_uniform()))
