*** Settings ***
Library   kakapo.workflow
Resource  resources/resource_directories.robot
Resource  resources/resource_taxonomy.robot

*** Test Cases ***
workflow.descending_tax_ids
    @{taxids}  Evaluate  [3701, 3702]
    @{desc_taxids_expected}  Evaluate  [3702, 29726, 38785, 45249, 45251, 59689, 59690, 59691, 63677, 63680, 81970, 81971, 81972, 97979, 97980, 97982, 302551, 347883, 347884, 347885, 378006, 395843, 412662, 675858, 675859, 864766, 864768, 869750, 869751, 1219863, 1240361, 1328956, 1532269, 1547868, 1547871, 1547872, 1547873, 1547874, 1746102, 1837063, 2203885, 2203886, 2203887, 2203888, 2203889, 2203890, 2486701]
    @{desc_taxids}  Descending Tax Ids
    ...    tax_ids_user=${taxids}
    ...    taxonomy=${tax}
    @{desc_taxids}  Evaluate  sorted(${desc_taxids})
    Should Be Equal  ${desc_taxids_expected}  ${desc_taxids}

workflow.pfam_uniprot_accessions
    @{pfacc}  Evaluate  ['PF18245']
    @{taxids}  Evaluate  [7215]
    @{uniprot_acc_expected}  Evaluate  ['A0A1W4VG33', 'B3MPT4', 'B4I753', 'B4MA44', 'B4NPI5', 'E1JJR3', 'Q29IV9', 'Q9VWI1']
    @{desc_taxids}  Descending Tax Ids
    ...    tax_ids_user=${taxids}
    ...    taxonomy=${tax}
    @{uniprot_acc}  Pfam Uniprot Accessions
    ...    pfam_acc=${pfacc}
    ...    tax_ids=${desc_taxids}
    ...    dir_cache_pfam_acc=${dir_cache_pfam_acc}
    @{uniprot_acc}  Evaluate  sorted(${uniprot_acc})
    Should Be Equal  ${uniprot_acc_expected}  ${uniprot_acc}
