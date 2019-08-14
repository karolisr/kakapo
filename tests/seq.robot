*** Settings ***
Library   kakapo.seq

*** Variables ***
${nt1}  GGNGCNGTNCTNTTRATHCCNTTYTAYTGYATGCAYAARCGNAGRTGGTCNAGYACNGAYGARAAYCARRAYSARTARTGA
${aa1}  GAVLLIPFYCMHKRRWSSTDENQXX**

${nt2}  attgcagatctgactcagaagatctttgaccttcgaggcaagtttaagcggcccaccctgcggagagtgaggatctctgcagatgccatgatgcaggcgctgctgggggcccgggctaaggagtccctgacctgcgggcccacctcaagcaggtga
${aa2}  IADLTQKIFDLRGKFKRPTLRRVRISADAMMQALLGARAKESLTCGPTSSR*

*** Test Cases ***
seq.reverse_complement
    ${x}  Reverse Complement          TGCAYRKMWSVHDB
    Should Be Equal As Strings  ${x}  VHDBSWKMYRTGCA

seq.translate 1
    ${x}  Translate  ${nt1}  ${tt1}  ${startc1}
    Should Be Equal  ${x}    ${aa1}

seq.translate 2
    ${x}  Translate  ${nt2}  ${tt1}  ${startc1}
    Should Be Equal  ${x}    ${aa2}
