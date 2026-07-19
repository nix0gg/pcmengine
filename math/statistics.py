import math

def mean(data):
    return sum(data) / len(data)

def mean_grouped(midpoints, frequencies):
    n = sum(frequencies)
    return sum(m * f for m, f in zip(midpoints, frequencies)) / n

def median(data):
    s = sorted(data)
    n = len(s)
    if n % 2 == 0:
        return (s[n//2 - 1] + s[n//2]) / 2
    return s[n//2]

def median_grouped(l, h, f, cf, n):
    # l = lower boundary, h = class width, f = freq of median class
    # cf = cumulative freq before median class, n = total freq
    return l + ((n/2 - cf) / f) * h

def mode(data):
    from collections import Counter
    c = Counter(data)
    return c.most_common(1)[0][0]

def mode_grouped(l, h, f1, f0, f2):
    # l = lower boundary of modal class, f1 = freq of modal class
    # f0 = freq before, f2 = freq after
    return l + (f1 - f0) / (2*f1 - f0 - f2) * h

def range_stat(data):
    return max(data) - min(data)

def mean_deviation_mean(data):
    m = mean(data)
    return sum(abs(x - m) for x in data) / len(data)

def mean_deviation_median(data):
    med = median(data)
    return sum(abs(x - med) for x in data) / len(data)

def mean_deviation_grouped(midpoints, frequencies, about='mean'):
    n = sum(frequencies)
    if about == 'mean':
        a = mean_grouped(midpoints, frequencies)
    else:
        a = median(midpoints)
    return sum(f * abs(m - a) for m, f in zip(midpoints, frequencies)) / n

def variance(data):
    m = mean(data)
    return sum((x - m) ** 2 for x in data) / len(data)

def variance_grouped(midpoints, frequencies):
    n = sum(frequencies)
    m = mean_grouped(midpoints, frequencies)
    return sum(f * (x - m) ** 2 for x, f in zip(midpoints, frequencies)) / n

def standard_deviation(data):
    return math.sqrt(variance(data))

def standard_deviation_grouped(midpoints, frequencies):
    return math.sqrt(variance_grouped(midpoints, frequencies))

def coefficient_of_variation(data):
    return (standard_deviation(data) / mean(data)) * 100

def coefficient_of_variation_grouped(midpoints, frequencies):
    sd = standard_deviation_grouped(midpoints, frequencies)
    m = mean_grouped(midpoints, frequencies)
    return (sd / m) * 100

def variance_shortcut(data):
    n = len(data)
    return sum(x**2 for x in data)/n - (sum(data)/n)**2

def combined_mean(n1, x1, n2, x2):
    return (n1 * x1 + n2 * x2) / (n1 + n2)

def combined_variance(n1, x1, v1, n2, x2, v2):
    xc = combined_mean(n1, x1, n2, x2)
    d1 = x1 - xc
    d2 = x2 - xc
    return (n1 * (v1 + d1**2) + n2 * (v2 + d2**2)) / (n1 + n2)

__all__ = [name for name in globals() if not name.startswith("_")]
