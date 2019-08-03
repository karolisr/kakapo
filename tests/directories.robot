*** Settings ***
Library   OperatingSystem
Library   kakapo.workflow

Resource  resources/resource_directories.robot

*** Test Cases ***
workflow.prepare_output_directories

    ${result}  Prepare Output Directories
    ...    dir_out=${dir_out}
    ...    prj_name=${prj_name}

    Length Should Be  ${result}  19

Output Directory

    Directory Should Exist  ${dir_out}
    ${dir_count}  Count Directories In Directory  ${dir_out}
    Should Be Equal As Integers  ${dir_count}  ${dir_out_count}

Cache Directory

    Directory Should Exist  ${dir_cache}
    ${dir_count}  Count Directories In Directory  ${dir_cache}
    Should Be Equal As Integers  ${dir_count}  ${dir_cache_count}

Cache Subdirectories

    Directory Should Exist  ${dir_cache_fq_minlen}
    Directory Should Exist  ${dir_cache_pfam_acc}
    Directory Should Exist  ${dir_cache_prj}

Temp Directory

    Directory Should Exist  ${dir_temp}
    ${dir_count}  Count Directories In Directory  ${dir_temp}
    Should Be Equal As Integers  ${dir_count}  ${dir_temp_count}

Global Directory

    Directory Should Exist  ${dir_global}
    ${dir_count}  Count Directories In Directory  ${dir_global}
    Should Be Equal As Integers  ${dir_count}  ${dir_global_count}

Project Directory

    Directory Should Exist  ${dir_prj}
    ${dir_count}  Count Directories In Directory  ${dir_prj}
    Should Be Equal As Integers  ${dir_count}  ${dir_prj_count}

Project Subdirectories

    Directory Should Exist  ${dir_prj_queries}
    Directory Should Exist  ${dir_prj_blast_results_fa_trim}
    Directory Should Exist  ${dir_prj_vsearch_results_fa_trim}
    Directory Should Exist  ${dir_prj_spades_assemblies}
    Directory Should Exist  ${dir_prj_blast_assmbl}
    Directory Should Exist  ${dir_prj_assmbl_blast_results}
    Directory Should Exist  ${dir_prj_transcripts}
    Directory Should Exist  ${dir_prj_transcripts_combined}
