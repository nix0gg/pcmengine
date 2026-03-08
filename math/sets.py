def custom_union(set1, set2):
    result = set1.copy()
    for element in set2:
        result.add(element)
    return result

def custom_intersection(set1, set2):
    result = set()
    for element in set1:
        if element in set2:
            result.add(element)
    return result

def custom_difference(set1, set2):
    result = set()
    for element in set1:
        if element not in set2:
            result.add(element)
    return result
