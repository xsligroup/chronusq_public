import re

####
# change 
#    tempOps["vooo_46"] = TArrayD(); 
# to 
#    TAManager::get().free("vooo", std::move(tempOps["vooo_46"]));

import re

def replace_free_lines(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    pattern = r'tempOps\["([a-zA-Z]+)_([0-9]+)"\] = TArrayD\(\);'
    replacement = r'TAmanager.free("\1", std::move(tempOps["\1_\2"]));'

    with open(filename, 'w') as file:
        for line in lines:
            new_line = re.sub(pattern, replacement, line)
            file.write(new_line)

def add_malloc_lines(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

               #tempOps["vooo_0"]("a,j,k,i") = 
    pattern = r'tempOps\["([a-zA-Z]+)_([0-9]+)"\]\(".+"\) = .*'
    replacement = r'tempOps["\1_\2"] = TAmanager.malloc<MatsT>("\1");'
#    replacement = r'hello"\1"_"\2"'

    with open(filename, 'w') as file:
        for line in lines:
            match = re.search(pattern, line)
            if match:
                new_line = re.sub(pattern, replacement, line)
                file.write(new_line)
            file.write(line)


filename = 'eom_dip_full_hamiltonian.code'
replace_free_lines(filename)
add_malloc_lines(filename)

