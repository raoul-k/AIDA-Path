
# Conditional Generators

import copy

from general_methods import *

#1: [[],[],(,,)]
#2: [[],[]]
#3: [{},{},...,k]
#4: [{},k]
#5: [{}] or [[{}]]
def collection_generators(generators):
    coll_of_gens = []
    #[generator number, generator, domain, numbers]
    for gen in generators:
        collected = collect_sets(gen)
        if isinstance(gen[0], list):
            if len(gen)>2: #[[],[],(,,)]
                coll_of_gens.append([1, gen, collected, gen[2]])
            elif len(gen)==2: #[[],[]]
                coll_of_gens.append([2, gen, collected, -1])
            else: #[]
                coll_of_gens.append([5, gen[0], collected, -1])   
        else:
            if len(gen)>2: #[{},{},...,k]
                coll_of_gens.append([3, gen, collected, gen[-1]])  
            elif len(gen)==2: #[{},k]
                coll_of_gens.append([4, gen, collected, gen[-1]]) 
            else: #[{}]
                coll_of_gens.append([5, gen, collected, -1])
    return coll_of_gens

def gen_simplifier_by_ones(generator, ones):
    cog = collection_generators(generator)
    ones_removed = []
    
    for symptom in ones:
        for gen in cog:
            if symptom in gen[2] and gen[0]!=2:
                gen[2].remove(symptom)
                ones_removed.append(symptom)
    
    new_generator = []
    new_generator.append([set(ones_removed)])

    #1: [[],[],(,,)] - !
    #2: [[],[]]
    #3: [{},{},...,k] - !
    #4: [{},k] - !
    #5: [{}] or [[{}]] - !
    for gen in cog:
        match gen[0]:
            case 1:
                new_tupl = list(gen[3])
                new_set1 = []
                new_set2 = []
                
                for h in gen[1][0]:
                    h_inter = set.intersection(h, gen[2])
                    if len(h) - len(h_inter) > 0:
                        new_tupl[0] = max(0, new_tupl[0]-1)
                        new_tupl[2] = max(0, new_tupl[2]-1)
                        h_inter = set()
                    if h_inter:
                        new_set1.append(h_inter)
                
                for h in gen[1][1]:
                    h_inter = set.intersection(h, gen[2])
                    if len(h) - len(h_inter) > 0:
                        new_tupl[1] = max(0, new_tupl[1]-1)
                        new_tupl[2] = max(0, new_tupl[2]-1)
                        h_inter = set()
                    if h_inter:
                        new_set2.append(h_inter)

                if new_set1 and new_set2:
                    new_generator.append([new_set1, new_set2, tuple(new_tupl)])
                elif new_set1 or new_set2:
                    new_generator.append([*(new_set1+new_set2), tuple(new_tupl)[2]])
            case 2:
                new_generator.append(gen[1])
            case 3:
                new_set = []
                new_k = gen[3]
                for h in gen[1][:-1]:
                    h_inter = h.intersection(gen[2])
                    if len(h)-len(h_inter) > 0:
                        new_k = new_k - 1
                    if h_inter:
                        new_set.append(h_inter)

                if new_set:
                    new_set.append(max(0, new_k))
                    new_generator.append(new_set)
            case 4: 
                if len(gen[2]) > 0:
                    new_generator.append([gen[2], max(0,gen[3]-(len(gen[1][0])-len(gen[2])))])
            case 5: 
                if len(gen[2]) > 0:
                    new_generator.append([gen[2]])
            case _:
                print('Should not be here!')

    return new_generator

def gen_simplifier_by_zeroes(generator, zeroes):
    cog = collection_generators(generator)
    new_generator = []

    #1: [[],[],(,,)] - !
    #2: [[],[]]
    #3: [{},{},...,k] - !
    #4: [{},k] - !
    #5: [{}] or [[{}]]    
    for gen in cog:
        match gen[0]:
            case 1:
                new_gen1 = copy.deepcopy(gen[1][0])
                new_gen2 = copy.deepcopy(gen[1][1])
                intersection = set.intersection(zeroes,set.union(*(new_gen1+new_gen2)))
                
                while len(intersection)>0:
                    x = intersection.pop()
                    temp1 = [s-set(x) for s in new_gen1 if len(s-set(x))>0]
                    temp2 = [s-set(x) for s in new_gen2 if len(s-set(x))>0]
                
                    if len(temp1)+len(temp2) >= gen[1][-1][2]:
                        if len(temp1) >= gen[1][-1][0]:
                            new_gen1 = temp1
                    
                        if len(temp2) >= gen[1][-1][1]:
                            new_gen2 = temp2
                
                new_generator.append([new_gen1, new_gen2, gen[1][-1]])
            case 2:
                new_generator.append(gen[1])
            case 3:
                new_gen = copy.deepcopy(gen[1])[:-1]
                intersection = set.intersection(zeroes,set.union(*new_gen))
                
                while len(intersection) > 0:
                    x = intersection.pop()
                    temp = [s-set(x) for s in new_gen if len(s-set(x))>0]
                    if len(temp) >= gen[1][-1]:
                        new_gen = temp
                
                new_gen.append(gen[1][-1])
                new_generator.append(new_gen)
            case 4: 
                new_gen = copy.deepcopy(gen[1])
                intersection = set.intersection(new_gen[0], zeroes)
                
                while len(intersection)>0 and len(new_gen[0])>new_gen[1]:
                    new_gen[0].remove(intersection.pop())
                
                new_generator.append(new_gen)
            case 5: 
                new_generator.append(gen[1])
            case _:
                print('Should not be here!')

    return new_generator

def gen_simplifier(generator1, generator2):
    gen1 = copy.deepcopy(generator1)
    gen2 = copy.deepcopy(generator2)
    gen1 = replace_dicts(gen1)
    gen2 = replace_dicts(gen2)
    domain1 = getdom(gen1)
    domain2 = getdom(gen2)
    
    gen1_excl = collect_excl_sets(gen1)
    gen2_excl = collect_excl_sets(gen2)
    gen1_ones = [ele for ele in domain1 if ele not in gen1_excl]
    gen2_ones = [ele for ele in domain2 if ele not in gen2_excl]
    gen1_zero = [ele for ele in union_list(domain1, domain2) if ele not in domain1]
    gen2_zero = [ele for ele in union_list(domain1, domain2) if ele not in domain2]

    union_excl = union_list(list(gen1_excl), list(gen2_excl))
    gen_inter_ones = [ele for ele in gen1_ones if ele in gen2_ones] #gen_inter_ones = [ele for ele in gen2_ones if ele in gen1_ones]
    gen1_inter_zero = [ele for ele in gen1_ones if ele in gen2_zero]
    gen2_inter_zero = [ele for ele in gen2_ones if ele in gen1_zero]
    
    gen1 = gen_simplifier_by_ones(gen1, set(gen_inter_ones))
    gen1 = gen_simplifier_by_zeroes(gen1, set(gen1_inter_zero))
    gen2 = gen_simplifier_by_ones(gen2, set(gen_inter_ones))
    gen2 = gen_simplifier_by_zeroes(gen2, set(gen2_inter_zero))

    return gen1, gen2