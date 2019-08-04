*** Settings ***
Library  ncbi_taxonomy_local

*** Variables ***
${dir_tax}                temp/ncbi-taxonomy-db

*** Keywords ***
Init Taxonomy
    ${value}  taxonomy  ${dir_tax}
    [Return]  ${value}
