*** Settings ***
Library   OperatingSystem
Library   kakapo.workflow
Resource  resources/resource_directories.robot
Resource  resources/resource_taxonomy.robot
Suite Setup      Setup before all tests
Suite Teardown   Cleanup after all tests

*** Keywords ***
Setup before all tests

    Evaluate  os.chdir('tests')  modules=os, sys

    Prepare Output Directories
        ...    dir_out=${dir_out}
        ...    prj_name=${prj_name}

    ${tax}  Init Taxonomy
    Set Global Variable  ${tax}

Cleanup after all tests

    Remove Directory
    ...    path=${dir_out}
    ...    recursive=True

    Remove Directory
    ...    path=${dir_tax}
    ...    recursive=True
