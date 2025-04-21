import re

# change
#   conj(this->antiSymMoints["vvoo"]("l,k,a,b"))
# to 
#   conj(this->antiSymMoints["vvoo"]("a,b,l,k"))

def reorder_conj(file):

    with open(filename, 'r') as file:
        lines = file.readlines()

    #conj_pat = re.compile(r'conj\((this->antiSymMoints\["[vo][vo][vo][vo]"\])\(([a-o]),([a-o]),([a-o]),([a-o])\)(.*)')
    pattern = r'(conj\(this->antiSymMoints\["[vo][vo][vo][vo]"\])\("([a-o]),([a-o]),([a-o]),([a-o])"\)(.*)'
    replacement = r'\1("\4,\5,\2,\3")\6'

    with open(filename, 'w') as file:
        for line in lines:
            match = re.search(pattern, line)
            if match:
                new_line = re.sub(pattern,replacement, line)
                file.write(new_line)
            else:
                file.write(line)




filename = 'eom_dip_full_hamiltonian.code'
reorder_conj(filename)
