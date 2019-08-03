*** Settings ***
Library   OperatingSystem
Library   kakapo.workflow
Resource  resources/resource_directories.robot

Suite Setup      Setup before all tests
Suite Teardown   Cleanup after all tests

*** Keywords ***
Setup before all tests
    Prepare Output Directories
        ...    dir_out=${dir_out}
        ...    prj_name=${prj_name}

Cleanup after all tests
    Remove Directory
    ...    path=${dir_out}
    ...    recursive=True
