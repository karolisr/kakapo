*** Settings ***
Library   OperatingSystem
Library   kakapo.workflow
Library   kakapo.translation_tables
Resource  resources/resource_directories.robot
Resource  resources/resource_taxonomy.robot
Suite Setup      Setup before all tests
Suite Teardown   Cleanup after all tests

*** Keywords ***
Setup before all tests

    Evaluate  os.chdir('tests')  modules=os

    Prepare Output Directories
        ...    dir_out=${dir_out}
        ...    prj_name=${prj_name}

    ${tax}  Init Taxonomy
    Set Global Variable  ${tax}

    ${TT}  Get Library Instance  kakapo.translation_tables
    &{tt1}  Set Variable  ${TT.TranslationTable(1).table_ambiguous}
    Set Global Variable  ${tt1}
    @{startc1}  Set Variable  ${TT.TranslationTable(1).start_codons_ambiguous}
    Set Global Variable  ${startc1}

Cleanup after all tests

    Remove Directory
    ...    path=${dir_out}
    ...    recursive=True

    Remove Directory
    ...    path=${dir_tax}
    ...    recursive=True
