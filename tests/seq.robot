*** Settings ***
Library   kakapo.seq

*** Variables ***
${nt1}  GGNGCNGTNCTNTTRATHCCNTTYTAYTGYATGCAYAARCGNAGRTGGTCNAGYACNGAYGARAAYCARRAYSARTARTGA
${aa1}  GAVLLIPFYCMHKRRWSSTDENQXX**

*** Test Cases ***
seq.reverse_complement
    ${x}  Reverse Complement          TGCAYRKMWSVHDB
    Should Be Equal As Strings  ${x}  VHDBSWKMYRTGCA

seq.translate
    ${x}  Translate  ${nt1}  ${tt1}
    Should Be Equal  ${x}    ${aa1}
