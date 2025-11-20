
import pandas as pd
import numpy as np
import copy

from combine_generators import gen_all_comb

#  Descr: get domain of the generators
#  Input: specific generator list (list with lists)
# Output: list with all symptoms (domain)
def getdom(gens):
    all_lists = []

    def recurse(nested):
        for element in nested:
            if isinstance(element, set):
                all_lists.append(list(element))
            elif isinstance(element, list):
                recurse(element)
    recurse(gens)
    
    return [item for sublist in all_lists for item in sublist]

#  Descr: single integer to binary combination based on a specific domain
#  Input: single integer, domain (list of symptoms), bin = as binary or symptom list
# Output: symptom set or binary representation
def int2comb(integer, domain, bin=True):
    temp = format(integer,'0'+str(len(domain))+'b')
    if bin:
        return temp
    else:
        return set([domain[idx] for idx in [i for i,s in enumerate(temp) if s=='1']])

#  Descr: integer list to binary combination based on a specific domain
#  Input: integer list, domain (list of symptoms), bin = as binary or symptom list
# Output: list of symptom sets or binary representation list
def int2comb_list(integers, domain, bin=True):
    temp = [format(i,'0'+str(len(domain))+'b') for i in integers]
    if bin:
        return temp
    else:
        return [set([domain[idx] for idx in [i for i,s in enumerate(t) if s=='1']]) for t in temp]

# only for small gens
def gen2CSV(generators, gen_name):
    domain = getdom(generators)
    combis = gen_all_comb(generators, ret=True)
    combis_df = [[gen_name]+[1 if element in comb else 0 for element in domain] for comb in combis]
    combis_df.insert(0,['name']+domain)
    combis_pd = pd.DataFrame(combis_df)
    combis_pd.to_csv(gen_name+'.csv', header=False, index=False)


def collect_sets(nested):
    all_sets = set()
    for item in nested:
        if isinstance(item, list):
            all_sets.update(collect_sets(item))
        elif isinstance(item, set):
            all_sets.update(item)
    return all_sets

def generators2matrix(generators, domain):
    list_combinations = gen_all_comb(generators, ret=True)
    return np.array([[1 if domain_element in combination else 0 for domain_element in domain] for combination in list_combinations])

def union_list(list1, list2):
    dcopy = copy.deepcopy(list1)
    [dcopy.append(l2) for l2 in list2 if l2 not in dcopy]
    return dcopy

def collect_excl_sets(nested):
    all_sets = set()
    for gen in nested:
        if isinstance(gen[0], list) and len(gen)==2:
            all_sets.update(*gen[0])
            all_sets.update(*gen[1])
    return all_sets

#set '{}' to 'set()'
def replace_dicts(nested):
    for i, elem in enumerate(nested):
        if isinstance(elem, dict):
            nested[i] = set()
        elif isinstance(elem, list):
            nested[i] = replace_dicts(elem)
    return nested