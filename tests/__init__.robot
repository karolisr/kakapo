*** Settings ***
Library   OperatingSystem
Library   kakapo.workflow
Library   ncbi_taxonomy_local.Taxonomy  WITH NAME  Tax
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

    &{tt1}  Tax.trans_table_for_genetic_code_id  1
    Set Global Variable  ${tt1}

Cleanup after all tests

    Remove Directory
    ...    path=${dir_out}
    ...    recursive=True

    Remove Directory
    ...    path=${dir_tax}
    ...    recursive=True
