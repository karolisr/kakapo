*** Variables ***
${prj_name}               robot

${dir_out}                temp/robot
${dir_out_count}          4

${dir_cache}              ${dir_out}/00-cache
${dir_cache_count}        3

${dir_cache_fq_minlen}    ${dir_cache}/min-acceptable-read-lengths
${dir_cache_pfam_acc}     ${dir_cache}/pfam-uniprot-accessions
${dir_cache_prj}          ${dir_cache}/projects

${dir_temp}               ${dir_out}/00-temp
${dir_temp_count}         0

${dir_global}             ${dir_out}/01-global
${dir_global_count}       4

${dir_prj}                ${dir_out}/02-project-specific/${prj_name}
${dir_prj_count}          8

${dir_prj_queries}                     ${dir_prj}/01-queries
${dir_prj_blast_results_fa_trim}       ${dir_prj}/02-trimmed-fa-blast-results
${dir_prj_vsearch_results_fa_trim}     ${dir_prj}/03-trimmed-fa-vsearch-results
${dir_prj_spades_assemblies}           ${dir_prj}/04-spades-assemblies
${dir_prj_blast_assmbl}                ${dir_prj}/05-assemblies-blast-db-data
${dir_prj_assmbl_blast_results}        ${dir_prj}/06-assemblies-blast-results
${dir_prj_transcripts}                 ${dir_prj}/07-transcripts
${dir_prj_transcripts_combined}        ${dir_prj}/08-transcripts-combined
