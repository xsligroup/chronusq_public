
cp eom_dip_full_hamiltonian.out eom_dip_full_hamiltonian.code

replace 'f_oo' 'this->fockMatrix_ta["oo"]'             eom_dip_full_hamiltonian.code 
replace 'f_ov' 'this->fockMatrix_ta["ov"]'             eom_dip_full_hamiltonian.code
replace 'f_vv' 'this->fockMatrix_ta["vv"]'             eom_dip_full_hamiltonian.code
replace 'f_vo' 'this->fockMatrix_ta["vo"]'             eom_dip_full_hamiltonian.code
replace 't1_vo' 'this->T1_'                            eom_dip_full_hamiltonian.code
replace 't2_vvoo' 'this->T2_'                          eom_dip_full_hamiltonian.code  
replace 'r2_oo' 'R2'                                   eom_dip_full_hamiltonian.code
replace 'r3_vooo' 'R3'                                 eom_dip_full_hamiltonian.code  
replace 'l2_oo' 'L2'                                   eom_dip_full_hamiltonian.code
replace 'l3_vooo' 'L3'                                 eom_dip_full_hamiltonian.code  
replace 'eri_oooo' 'this->antiSymMoints["oooo"]'       eom_dip_full_hamiltonian.code
replace 'eri_vooo' 'this->antiSymMoints["vooo"]'       eom_dip_full_hamiltonian.code
replace 'eri_vovo' 'this->antiSymMoints["vovo"]'       eom_dip_full_hamiltonian.code
replace 'eri_oovv' 'conj(this->antiSymMoints["vvoo"]'  eom_dip_full_hamiltonian.code
replace 'eri_oovo' 'conj(this->antiSymMoints["vooo"]'  eom_dip_full_hamiltonian.code
replace 'eri_vovv' 'conj(this->antiSymMoints["vvvo"]'  eom_dip_full_hamiltonian.code
sed -i 's/conj\([^)]*"\)/conj\1)/'                     eom_dip_full_hamiltonian.code
replace ' d_oo' ' Id_oo' eom_dip_full_hamiltonian.code
replace ' d_vv' ' Id_vv' eom_dip_full_hamiltonian.code
python3 conj_reorder.py
python3 ta_free.py
