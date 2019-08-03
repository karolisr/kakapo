*** Settings ***
Library  kakapo.seq

*** Test Cases ***
Reverse Complement
    ${result}=  Reverse Complement
  ...  seq=AAAA
    Should be equal  ${result}  TTTT
