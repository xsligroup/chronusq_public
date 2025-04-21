import os
import pickle

import sys
import pdaggerq

def load_pq_output(pq_filepath):
    pq_file = open(pq_filepath, "rb")
    pq_data = pickle.load(pq_file)
    pq_file.close()

    return pq_data[0], pq_data[1]

def main():
    all_contracted_strings = []
    pq = pdaggerq.pq_helper("fermi")
   
    # < L2|H|R2 >
    # set right and left-hand operators
    pq.set_right_operators_type('DIP')
    pq.set_left_operators([['a*(i)', 'a*(j)']])
    pq.set_right_operators([['a(l)','a(k)']])
    
    print('')
    print('    H(i,j;k,l) = <0|i* j* e(-T) H e(T) l k |0>')
    print('')
    
    pq.add_st_operator(1.0,['f'],['t1','t2'])
    pq.add_st_operator(1.0,['v'],['t1','t2'])
    
    pq.simplify()
    # grab list of fully-contracted strings, then print
    all_contracted_strings.append(pq.fully_contracted_strings())
    for my_term in pq.fully_contracted_strings():
        print(my_term)
    pq.clear()

    # < L3|H|R3 >
    # set right and left-hand operators
    pq.set_right_operators_type('DIP')
    pq.set_left_operators([['a*(i)', 'a*(j)', 'a*(k)', 'a(a)']])
    pq.set_right_operators([['a*(b)','a(n)','a(m)','a(l)',]])
    
    print('')
    print('    H(ijka;lmnb) = <0|i* j* k* a e(-T) H e(T) b* n m l|0>')
    print('')
    
    pq.add_st_operator(1.0,['f'],['t1','t2'])
    pq.add_st_operator(1.0,['v'],['t1','t2'])
    
    pq.simplify()
    # grab list of fully-contracted strings, then print
    all_contracted_strings.append(pq.fully_contracted_strings())
    for my_term in pq.fully_contracted_strings():
        print(my_term)
    pq.clear()


    # < L2|H|R3 >
    # set right and left-hand operators
    pq.set_right_operators_type('DIP')
    pq.set_left_operators([['a*(i)', 'a*(j)']])
    pq.set_right_operators([['a*(b)','a(n)','a(m)','a(l)']])
    
    print('')
    print('    H(i,j;lmnb) = <0|i* j* e(-T) H e(T) b* n m l|0>')
    print('')
    
    pq.add_st_operator(1.0,['f'],['t1','t2'])
    pq.add_st_operator(1.0,['v'],['t1','t2'])
    
    pq.simplify()
    # grab list of fully-contracted strings, then print
    all_contracted_strings.append(pq.fully_contracted_strings())
    for my_term in pq.fully_contracted_strings():
        print(my_term)
    pq.clear()

    # < L3|H|R3 >
    # set right and left-hand operators
    pq.set_right_operators_type('DIP')
    pq.set_left_operators([['a*(i)', 'a*(j)', 'a*(k)', 'a(a)']])
    pq.set_right_operators([['a(m)','a(l)',]])
    
    print('')
    print('    H(ijka;lm) = <0|i* j* k* a e(-T) H e(T) m l|0>')
    print('')
    
    pq.add_st_operator(1.0,['f'],['t1','t2'])
    pq.add_st_operator(1.0,['v'],['t1','t2'])
    
    pq.simplify()
    # grab list of fully-contracted strings, then print
    all_contracted_strings.append(pq.fully_contracted_strings())
    for my_term in pq.fully_contracted_strings():
        print(my_term)
    pq.clear()




    # remove '"' from eqNames
    # before (), you can name them whatever as long as it contains its featureted such as t3,u0,u1,r0,r1,etc;
    # within (), the information is related with dimention of tensor
    # you can use ab-g for vir, ij-o for occ
    eqTensorTypes = ['H(i,j,k,l)', 'H(i,j,k,a,l,m,n,b)', 'H(i,j,l,m,n,b)','H(i,j,k,a,l,m)' ]

#eqNames = [x.replace('"', '') for x in eqNames]

    ## replace ', ' with ',' from eqNames
    #eqTensorTypes = [x.replace(', ', ',') for x in eqTensorTypes]
    
    ## polariton stuff
    ## replace all 'd+' and 'd-' with 'dp'
    #for eq, i in enumerate(all_eq):
    #    for term, j in enumerate(i):
    #        for tensor, k in enumerate(j):
    #           all_eq[eq][term][tensor] = k.replace('d+', 'dp').replace('d-', 'dp')

    tabuild = pdaggerq.tabuilder()

    tabuild.set_options({
        't1_transform': False,
        'verbose': False,
        #'t1_transform': False,
        #'make_scalars' : True,
        #'batched' : False,
        #'iterative_merge' : True,
        #'permuted_merge' : True,
        'reuse_permutations' : True,
        'sigma_vectors' : {'r0','r1','r2','r3','s0','s1','s2','s3','l0','l1','l2','l3','m0','m1','m2','m3'}, #determines name of trial vectors,
        #'conditions' : {'t3','u0','u1','u2','u3','r3','s0','s1','s2','s3','l3','m0','m1','m2','m3'}, # if include t3, put in conditions.
        #'max_contraction_ops' : -1,
        #'depth' : -1,
        #'max_temps' : -1,
        #'num_threads' : 12,
    })

    tabuild.build(eqTensorTypes, all_contracted_strings)
    #tabuild.optimize()
    tabuild.reorder()
    tabuild.substitute()
    tabuild.merge_terms()
    #tabuild.merge_permutations()
    tabuild.reuse()
    print_str = tabuild.str()

    #formatting
    import re

    eri_match = re.compile(r"eri_([ab]+?_[ov]+?)\(")
    f_match = re.compile(r"f_([ab]+?_[ov]+?)\(")
    Id_match = re.compile(r"Id_([ab]+?_[ov]+?)\(")
    d_match = re.compile(r"d_([ab]+?_[ov]+?)\(")
    dp_match = re.compile(r"d[+-]_([ab]+?_[ov]+?)\(")

    print_str = eri_match.sub(r'V_blks_["\1"](', print_str)
    print_str = f_match.sub(r'F_blks_["\1"](', print_str)
    print_str = Id_match.sub(r'Id_blks_["\1"](', print_str)
    print_str = d_match.sub(r'Id_blks_["\1"](', print_str)
    print_str = dp_match.sub(r'dp_\1(', print_str)

    print(print_str, flush=True)

    tabuild.analysis()

if __name__ == '__main__':
    main()








