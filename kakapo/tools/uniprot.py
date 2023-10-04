"""UniProtKB Proteins."""
import re
from collections.abc import Collection, Mapping
from typing import Any, Generator, Union

from kakapo.tools.ebi_domain_search import pfam_entry
from kakapo.utils.http import Response, get
from kakapo.utils.logging import Log

# https://www.uniprot.org/help/api_queries
# https://www.uniprot.org/help/query-fields
# https://www.uniprot.org/help/return_fields

# Possible values for the Accept request-header field:
UNIPROT_ACC_HEAD: dict[str, dict[str, str]] = {
    'fasta': {'Accept': 'text/plain; format=fasta'},
    'gff': {'Accept': 'text/plain; format=gff'},
    'json': {'Accept': 'application/json'},
    'list': {'Accept': 'text/plain; format=list'},
    'obo': {'Accept': 'text/plain; format=obo'},
    'rdf': {'Accept': 'application/rdf+xml'},
    'tsv': {'Accept': 'text/plain; format=tsv'},
    'txt': {'Accept': 'text/plain; format=flatfile'},
    'xlsx': {'Accept': 'application/vnd.ms-excel'},
    'xml': {'Accept': 'application/xml'},
}

UNIPROT_FIELD_IO_MAP: dict[str, tuple[str, ...]] = {
    'accession': ('primaryAccession',),
    'organism_id': ('organism', 'taxonId',),
    'id': ('uniProtkbId',),
}


UNIPROT_API_URL = 'https://rest.uniprot.org/uniprotkb/search'
REGEX_NEXT_LINK = re.compile(r'<(.+)>; rel="next"')


def _valid_response_formats() -> tuple[str, ...]:
    return tuple(UNIPROT_ACC_HEAD.keys())


def _make_query(pfam_acc: str, taxids: Collection[int],
                only_reviewed: bool = False) -> str:
    taxn_term: str = f'(taxonomy_id: {" OR taxonomy_id:".join(map(str, taxids))})'
    pfam_term: str = f'(xref:pfam-{pfam_acc})'
    revd_term: str = '(reviewed:true)'
    terms: list[str] = [taxn_term, pfam_term]

    if only_reviewed is True:
        terms.append(revd_term)

    query: str = ' AND '.join(terms)
    return query


def _link_to_next_page(headers: Mapping[str, str]) -> Union[str, None]:
    if 'Link' in headers:
        match = REGEX_NEXT_LINK.match(headers['Link'])
        if match is not None:
            return match.group(1)


def _get_page(params: Union[dict, None], headers: dict
              ) -> Generator[tuple[Response, int], Any, None]:
    url: Union[str, None] = UNIPROT_API_URL
    while url is not None:
        response = get(url, params, headers)
        assert response is not None
        total = int(response.headers['X-Total-Results'])
        yield response, total
        url = _link_to_next_page(response.headers)
        params = None


def _uniprot_search_paged_get(params: dict, page_size: int,
                              response_format: str,
                              fields: list[str], logger=Log) -> Any:

    params['size'] = page_size

    if response_format in ('fasta',):
        fields = []

    params['fields'] = ','.join(fields)

    assert response_format in _valid_response_formats()
    headers = UNIPROT_ACC_HEAD[response_format]

    results: list = list()
    for response, total in _get_page(params, headers):
        if response_format == 'json':
            results += response.json()['results']
        elif response_format == 'fasta':
            results += response.text.strip()[1:].split('\n>')
        # _ = len(str(total))
        # logger.msg(f'\t{len(results):>{_}}/{total}'.rstrip(), '')

    results_parsed: Any = None

    if response_format == 'json':
        results_parsed = list()
        for r in results:
            for f in fields:
                assert f in UNIPROT_FIELD_IO_MAP
                value = UNIPROT_FIELD_IO_MAP[f]
                if len(value) == 1:
                    results_parsed.append(r[value[0]])
                elif len(value) == 2:
                    results_parsed.append(r[value[0]][value[1]])

    elif response_format == 'fasta':
        results_parsed = '>' + '\n>'.join(results)

    return results_parsed


def uniprot_entries_for_pfam_term(
        pfam_term: str,
        taxids: Collection[int],
        response_format: str = 'json',
        fields: list[str] = ['accession']):

    pfam_accs = pfam_entry(pfam_term)
    if len(pfam_accs) > 0:
        assert len(pfam_accs) == 1
        pfam_term = pfam_accs[0]['acc']

    query: str = _make_query(pfam_term, taxids)
    params = {'query': query}

    results_parsed: Any = _uniprot_search_paged_get(
        params, 500, response_format, fields)

    return results_parsed
