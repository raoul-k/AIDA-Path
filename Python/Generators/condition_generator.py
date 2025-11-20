
from itertools import chain, combinations, product

# cond_gen0: [{...}, k]
# cond_gen1: [{...}, {...}, ..., k]
# cond_gen2: [[{.},...], [{.},...], (l,m,n)]
# cond_gen3: [[{.},...], [{.},...]]

# Condition Generator 0
#######################
#  Descr: generate the elements of the powerset with at least k elements
#         (alternative function name: cond_gen0)
#  Input: list or set of symptoms, condition k
# Output: itertools.chain (needs to be iterated over and transformed to list with sets)
def powerset(s, k=0):
    s = list(s)
    return chain.from_iterable(combinations(s,r) for r in range(k,len(s)+1))

# Condition Generator 1
#######################
#  Descr: generate a list with sets with at least k elements from k different sets
#  Input: list with sets e.g. [{...},{...},...,{...}], condition k      
# Output: list with sets
def cond_gen1(sets, k):
    # generate the powerset
    pset = []
    pset = [{*s} for s in powerset(set.union(*sets))]
    
    # check if the condition is met for each element
    fin = []
    for ps in pset:
        c = len([1 for s in sets if len(s&ps)])
        if c >= k:
            fin.append(ps)
    return fin

# Condition Generator 2
#######################
#  Descr: generate a list with sets with at least k elements from first, l elements from second and m elements from both lists
#  Input: list with 2 lists containing sets e.g. [[{.},...],[{.},...]], conditions k in the form (l,m,n)
# Output: list with sets
def cond_gen2(lists, k):
    # generate the powersets for each list
    # make the product of each element in each set and unite them
    pset0, pset1, prod_pset, prod_uset = [], [], [], []
    pset0 = [{*s} for s in powerset(set.union(*lists[0]))]
    pset1 = [{*s} for s in powerset(set.union(*lists[1]))]
    prod_pset = [set.union(*p) for p in product(pset0, pset1)]
    
    # check if the conditions are met for each element
    fin = []
    for ps in prod_pset:
        c0 = len([1 for s in lists[0] if len(s&ps)])
        c1 = len([1 for s in lists[1] if len(s&ps)])
        
        if c0>=k[0] and c1>=k[1] and c0+c1>=k[2]:
            fin.append(ps)
    return fin

# Condition Generator 3
#######################
#  Descr: just cartesian product of 2 lists
#  Input: list with 2 lists containing sets e.g. [[{.},...],[{.},...]], conditions k in the form (l,m,-n)
# Output: list with sets
def cond_gen3(lists):
    return [set.union(*p) for p in product(lists[0], lists[1])]