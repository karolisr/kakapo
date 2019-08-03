*** Settings ***
Suite Setup  Setup before all tests

*** Keywords ***
Setup before all tests
    evaluate  sys.path.append('${EXECDIR}')  modules=os, sys
