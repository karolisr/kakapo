"""
Wrap NCBI's Entrez Programming Utilities (E-utilities).

More information on E-utilities at:
    http://www.ncbi.nlm.nih.gov/books/NBK25497

Entrez Unique Identifiers (UIDs) for selected databases:
    http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T
    ._entrez_unique_identifiers_ui/?report=objectonly

Valid values of &retmode and &rettype for EFetch (null = empty string):
    https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T
    ._valid_values_of__retmode_and/?report=objectonly
"""

import os
from collections.abc import Iterable
from functools import wraps
from typing import Callable as CallableT
from typing import Iterable as IterableT
from xml.etree import ElementTree

from multipledispatch import dispatch
from requests.models import Response

from kakapo.tools.parsers import eutils_loc_str
from kakapo.tools.parsers import parse_efetch_sra_csv_text
from kakapo.tools.parsers import seq_records_gb
from kakapo.utils.http import post
from kakapo.utils.logging import Log

EUTILS_NS = dict()
ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def find_api_key() -> str:
    global_variables = globals()
    if 'ENTREZ_KEY' in global_variables:
        api_key = global_variables['ENTREZ_KEY']
    elif 'ENTREZ_KEY' in os.environ:
        api_key = os.environ['ENTREZ_KEY']
    else:
        Log.wrn('Warning:', 'ENTREZ_KEY is not defined.')
        api_key = None

    return api_key


def ncbi_api_key(f: CallableT) -> CallableT:
    @wraps(f)
    def wrapper(api_key: str = None, *args, **kwargs):
        if api_key is None:
            api_key = find_api_key()
        return f(api_key=api_key, *args, **kwargs)

    return wrapper


def history(f: CallableT) -> CallableT:
    @wraps(f)
    def wrapper(usehistory: bool = False, web_env: str = None,
                query_key: str = None, *args,
                **kwargs):
        if usehistory is True or web_env is not None:
            usehistory = 'y'
        else:
            usehistory = None

        return f(usehistory=usehistory, web_env=web_env, query_key=query_key,
                 *args, **kwargs)

    return wrapper


def with_history(f: CallableT) -> CallableT:
    @wraps(f)
    def wrapper(data: dict):

        db = data['db']
        query_keys = data['query_keys']
        web_env = data['web_env']

        new_data = {'db': db, 'query_key': None, 'web_env': web_env}

        if len(query_keys) == 1:
            new_data['query_key'] = query_keys[0]
            return f(new_data)

        elif len(query_keys) > 1:
            ret_list = list()
            for k in query_keys:
                new_data['query_key'] = k
                ret_list += f(new_data)
            return ret_list

    return wrapper


@ncbi_api_key
@history
def eutil(util, db, api_key=None, bdata=None, cmd=None, complexity=None,
          datetype=None, dbfrom=None, email=None, field=None, holding=None,
          ids=None, ids_separate=None, idtype='acc', linkname=None,
          location=None, maxdate=None, mindate=None, query_key=None,
          reldate=None, retmax=None, retmode=None, retstart=None, rettype=None,
          seq_start=None, seq_stop=None, sort=None, strand=None, term=None,
          tool=None, usehistory=None, version=None, web_env=None) -> Response:
    url = ENTREZ_BASE_URL + util

    if ids is not None:
        assert isinstance(ids, IterableT)
        ids = ','.join(ids)

    params = {'api_key': api_key, 'bdata': bdata, 'cmd': cmd,
              'complexity': complexity, 'datetype': datetype, 'db': db,
              'dbfrom': dbfrom, 'email': email, 'field': field,
              'holding': holding, 'id': ids, 'idtype': idtype,
              'linkname': linkname, 'location': location, 'maxdate': maxdate,
              'mindate': mindate, 'query_key': query_key, 'reldate': reldate,
              'retmax': retmax, 'retmode': retmode, 'retstart': retstart,
              'rettype': rettype, 'seq_start': seq_start, 'seq_stop': seq_stop,
              'sort': sort, 'strand': strand, 'term': term, 'tool': tool,
              'usehistory': usehistory, 'version': version, 'WebEnv': web_env}

    params = [(k, params[k]) for k in params]

    if ids_separate is not None:
        params += ids_separate

    return post(url, params, retmode)


# E-utilities ----------------------------------------------------------------
def einfo(db: str) -> dict:
    r = eutil(util='einfo.fcgi', db=db, version='2.0', retmode='json')
    if r is not None:
        return r.json()


def esearch(db: str, term: str, web_env: str = None,
            query_key: int = None, **kwargs) -> dict:
    r = eutil(util='esearch.fcgi', rettype='uilist', retmode='json', db=db,
              term=term, usehistory=True, web_env=web_env,
              query_key=query_key, **kwargs)

    if r is not None:
        esearchresult = r.json()['esearchresult']
        web_env = esearchresult['webenv']
        query_key = int(esearchresult['querykey'])
        return epost(db=db, web_env=web_env, query_key=query_key)


def epost(db: str, ids: IterableT[str] = None, web_env: str = None,
          **kwargs) -> dict:
    r = eutil(util='epost.fcgi', db=db, ids=ids, web_env=web_env,
              retmode='xml', **kwargs)

    query_keys = list()
    web_env = None

    if r is not None:
        root = ElementTree.fromstring(r.text)
        for child in root:
            if child.tag == 'QueryKey':
                query_keys.append(int(child.text))
            if child.tag == 'WebEnv':
                web_env = child.text

    return {'db': db, 'query_keys': query_keys, 'web_env': web_env}


def esummary(db: str, ids: IterableT[str] = None, **kwargs) -> list:
    return_list = []

    r = eutil(util='esummary.fcgi', db=db, ids=ids, version='2.0',
              retmode='json', **kwargs)
    if r is not None:
        parsed = r.json()['result']
        keys = parsed['uids']

        for k in keys:
            return_list.append(parsed[k])

        return return_list


def efetch(db: str, ids: IterableT[str] = None, rettype: str = None,
           retmode: str = None, **kwargs) -> str:
    r = eutil(util='efetch.fcgi', db=db, ids=ids, rettype=rettype,
              retmode=retmode, **kwargs)
    if r is not None:
        return r.text


def elink(db: str, dbfrom: str, cmd: str, ids: IterableT[str] = None,
          **kwargs) -> list:
    if ids is not None:
        ids = [('id', x) for x in ids]
    r = eutil(util='elink.fcgi', db=db, dbfrom=dbfrom, cmd=cmd,
              ids_separate=ids, retmode='json', **kwargs)

    if r is not None:
        r_dict = r.json()
        linksets = r_dict['linksets']
        return linksets


def espell(db: str, term: str, **kwargs) -> str:
    r = eutil(util='espell.fcgi', db=db, term=term, retmode='xml', **kwargs)

    corrected_query = None
    if r is not None:
        root = ElementTree.fromstring(r.text)
        for child in root:
            if child.tag == 'CorrectedQuery':
                corrected_query = child.text
                break

    return corrected_query


# Secondary functions (private) ----------------------------------------------
def _strip_split(txt: str) -> list:
    return txt.strip().split('\n')


def _elink_parse(linksets, dbfrom, db) -> tuple:
    ids_from = []
    ids_to = []
    for ls in linksets:
        if dbfrom == ls['dbfrom']:
            ids_from += ls['ids']
            linksetdbs = ls['linksetdbs']
            for lsdb in linksetdbs:
                if db == lsdb['dbto']:
                    ids_to += lsdb['links']
    return ids_from, ids_to


def _elink_runner(dbfrom: str, db: str, ids: IterableT[str]) -> tuple:
    linkname = dbfrom + '_' + db
    linksets = elink(cmd='neighbor', dbfrom=dbfrom, db=db, linkname=linkname,
                     ids=ids)
    ids_from, ids_to = _elink_parse(linksets, dbfrom, db)
    return ids_from, ids_to


# Secondary functions --------------------------------------------------------
@dispatch(dict, namespace=EUTILS_NS)
@with_history
def accs(data: dict) -> list:
    txt = efetch(rettype='uilist', retmode='text', **data)
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(str, Iterable, namespace=EUTILS_NS)
def accs(db: str, ids: IterableT[str]) -> list:
    txt = efetch(db=db, ids=ids, rettype='uilist', retmode='text')
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(str, Iterable, namespace=EUTILS_NS)
def taxids(db: str, ids: IterableT[str]) -> dict:
    dbfrom = db
    db = 'taxonomy'
    ids_from, ids_to = _elink_runner(dbfrom, db, ids)
    ids_to = [int(x) for x in ids_to]
    return dict(zip(ids_from, ids_to))


@dispatch(dict, namespace=EUTILS_NS)
def taxids(data: dict) -> dict:
    ids = accs(data)
    return taxids(data['db'], ids)


@dispatch(dict, namespace=EUTILS_NS)
def cds_accs(data_protein: dict) -> dict:
    assert data_protein['db'] == 'protein'
    ids_protein = accs(data_protein)
    return cds_accs(ids_protein)


@dispatch(Iterable, namespace=EUTILS_NS)
def cds_accs(ids_protein: IterableT[str]) -> dict:
    dbfrom = 'protein'
    db = 'nuccore'
    ids_from, ids_to = _elink_runner(dbfrom, db, ids_protein)
    return dict(zip(ids_from, ids_to))


def seqs(db: str, ids: IterableT[str], rettype: str = 'gb',
         retmode: str = 'xml') -> list:
    gb_txt = efetch(db=db, ids=ids, rettype=rettype, retmode=retmode,
                    usehistory=True)
    return seq_records_gb(gb_txt)


def cds(ids_protein: IterableT[str]) -> list:
    rec_list = seqs('protein', ids_protein)
    ret_list = list()
    for rec in rec_list:
        cds_info = rec.coded_by
        if cds_info is not None:
            cds_loc_str = eutils_loc_str(cds_info)
            cds_acc = cds_info['external_ref']
            cds_gb = efetch(db='nuccore', ids=[cds_acc], location=cds_loc_str,
                            rettype='gb', retmode='xml')
            cds_seq_record = seq_records_gb(cds_gb)[0]
            ret_list.append(cds_seq_record)
    return ret_list


def sra_run_info(ids_srr: IterableT[str]) -> list:
    efetch_csv_txt = efetch(db='sra', ids=ids_srr, rettype='runinfo',
                            retmode='csv')
    ret_list = parse_efetch_sra_csv_text(efetch_csv_txt)
    return ret_list


# test_data = epost(db='protein', ids=['GER25982.1', 'NP_001098858'])

# accs('protein', ['GER25982.1', 'NP_001098858'])
# accs(test_data)

# taxids('protein', ['GER25982.1', 'NP_001098858'])
# taxids(test_data)

# cds_accs(['GER25982.1', 'NP_001098858'])
# cds_accs(test_data)

# gb_txt = efetch(db='protein', ids=['NP_001098858'],
#                 rettype='gb', retmode='xml')

# gb_txt = efetch(db='nuccore', ids=['NM_001105388.1'],
#                 rettype='gb', retmode='xml')

# seqs('nuccore', ['KF764990.1'])
# seqs('protein', ['GER25982.1', 'NP_001098858'])

# cds(['GER25982.1', 'NP_001098858'])

# sra_run_info(['SRR000060'])
