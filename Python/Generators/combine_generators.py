
import itertools
import numpy as np

from condition_generator import *
from general_methods import *

#  Descr: yield results from itertools product function
#  Input: list of generated lists
# Output: cartesian product of all elements (sets) between all lists
def gen_comb_yield(*list):
    for combi in itertools.product(*list):
        yield combi

#  Descr: writes the combined (cartesian product) sets as bytes into a binary file,
#         the first byte is reserved for domain length so the rest can be read back with the correct byte size,
#         it is dumped in chunks into the file to lessen the number of write operations
#  Input: file path, domain, list of generated lists
# Output: - (file)
def gen_comb_save(file_path, domain, *list):
    with open(file_path,'wb') as f:
        chunk = 5_000
        chunk_list = []
        f.write(len(domain).to_bytes(1, byteorder='big'))
        
        for combi in itertools.product(*list):
            yield combi
            combi_union = [set.union(*combi)]
            s = ["".join(map(str,[1 if element in cset else 0 for element in domain])) for cset in combi_union]
            chunk_list.append(*s)
            
            if len(chunk_list) >= chunk:
                chunk_list_bytes = [int(a,2).to_bytes((len(a)+7) // 8, byteorder='big') for a in chunk_list]
                f.write(b''.join(chunk_list_bytes))
                chunk_list = []
                
        if chunk_list:
            chunk_list_bytes = [int(a,2).to_bytes((len(a)+7) // 8, byteorder='big') for a in chunk_list]
            f.write(b''.join(chunk_list_bytes))
        
        print("Last 3 Elements:", chunk_list[-3:])
        print()

#  Descr: reads a binary file and converts the data back to list of symptom combinations (as an integer list)
#  Input: file path
# Output: numpy array in unsigned integer 64 format
def read_all_comb(file_path):
    chunk = 5_000 #chunk size to read
    chunk_np = 1_000_000 #chunk size to write to numpy array
    chunk_list = []
    chunk_np_list = []
    combis = np.empty([0,1],dtype=np.uint64)
    counter = 0
    
    with open(file_path,"rb") as f:
        bytes2read = f.read(1) # first byte is always the domain length
        bytes2read = (int.from_bytes(bytes2read, byteorder='big')+7) // 8
        #print("Domain Length: ", int.from_bytes(bytes2read, byteorder='big'))
        #print("Bytes to Read: ", bytes2read)
        
        while True:
            entry = f.read(bytes2read)
            if not entry:
                break
            chunk_list.append(entry)
            
            counter = counter + 1
            if counter >= chunk:
                chunk_np_list += [np.uint64(int.from_bytes(a, byteorder='big')) for a in chunk_list]
                chunk_list = []
                counter = 0
                
            if len(chunk_np_list) >= chunk_np:
                combis = np.append(combis, np.array(chunk_np_list))
                chunk_np_list = []
                
        if chunk_list:
            chunk_np_list += [np.uint64(int.from_bytes(a, byteorder='big')) for a in chunk_list]
            combis = np.append(combis, np.array(chunk_np_list))
    
    return combis

#  Descr: expands all the generators into lists, combines them (cartesion product) and saves them in a binary file
#  Input: file path, list of the generators
# Output: - (file)
def gen_all_comb_save(file_path, list_of_gens):
    list_of_combis = []
    for gen in list_of_gens:
        # first generator element is a list -> #2/#3
        if isinstance(gen[0], list):
            # generator len is over 2 -> #2
            if len(gen)>2:
                generated = cond_gen2(gen[:-1], gen[-1])
                list_of_combis.append(generated)
            elif len(gen)==2:
                generated = cond_gen3(gen)
                list_of_combis.append(generated) 
            else:
                list_of_combis.append(gen[0])
        # first generator element is a set -> #0/#1        
        else:
            # generator len is over 2 -> #1
            if len(gen)>2:
                generated = cond_gen1(gen[:-1], gen[-1])
                list_of_combis.append(generated)
            elif len(gen)==2:
                generated = [{*x} for x in powerset(gen[0], gen[1])]
                list_of_combis.append(generated)
            else:
                list_of_combis.append(gen)
    
    domain = getdom(list_of_gens)
    combi_generator = gen_comb_save(file_path, domain, *list_of_combis)    
    counter = 0
    for c in combi_generator:
        counter = counter + 1
        
    #print(counter,"profiles generated!")
    print(f'{counter:,}',"profiles generated!")

#  Descr: expands all the generators into lists, combines them (cartesion product) and prints them
#  Input: file path, list of the generators
# Output: - (cartesian products per print)
def gen_all_comb(list_of_gens, ret=False):
    list_of_combis = []
    for gen in list_of_gens:
        # first generator element is a list -> #2/#3
        if isinstance(gen[0], list):
            # generator len is over 2 -> #2
            if len(gen)>2:
                generated = cond_gen2(gen[:-1], gen[-1])
                list_of_combis.append(generated)
            elif len(gen)==2:
                generated = cond_gen3(gen)
                list_of_combis.append(generated) 
            else:
                list_of_combis.append(gen[0])
        # first generator element is a set -> #0/#1        
        else:
            # generator len is over 2 -> #1
            if len(gen)>2:
                generated = cond_gen1(gen[:-1], gen[-1])
                list_of_combis.append(generated)
            elif len(gen)==2:
                generated = [{*x} for x in powerset(gen[0], gen[1])]
                list_of_combis.append(generated)
            else:
                list_of_combis.append(gen)
    
    combi_generator = gen_comb_yield(*list_of_combis)    
    
    if ret:
        return_list = []
        for c in combi_generator:
            return_list.append(set.union(*c))
            counter  = counter + 1
        print()
        print(f'{counter:,}',"profiles generated!")
        
        return return_list
    else:
        counter = 0
        for c in combi_generator:
            print(set.union(*c))
            counter  = counter + 1
        print()
        print(f'{counter:,}',"profiles generated!")
