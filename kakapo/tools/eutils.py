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
from io import StringIO
from time import sleep
from typing import Callable as CallableT
from typing import Iterable as IterableT
from xml.etree import ElementTree

from multipledispatch import dispatch
from requests.models import Response

from kakapo.tools.bioio import read_fasta
from kakapo.tools.parsers import eutils_loc_str
from kakapo.tools.parsers import parse_efetch_sra_csv_text
from kakapo.tools.parsers import parse_esummary_xml_text
from kakapo.tools.parsers import seq_records_gb
from kakapo.tools.seq import SEQ_TYPE_NT, SEQ_TYPE_AA
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
    def wrapper(data: dict, **kwargs):
        # def wrapper(data: dict, rettype: str = 'gb', retmode: str = 'xml'):

        db = data['db']
        query_keys = data['query_keys']
        web_env = data['web_env']

        new_data = {'db': db, 'query_key': None, 'web_env': web_env}

        if len(query_keys) == 0:
            return f(new_data)

        elif len(query_keys) == 1:
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


def esearch(db: str, term: str, usehistory: bool = False, web_env: str = None,
            query_key: int = None, **kwargs) -> dict:
    r = eutil(util='esearch.fcgi', rettype='uilist', retmode='json', db=db,
              term=term, usehistory=usehistory, web_env=web_env,
              query_key=query_key, **kwargs)
    if r is not None:
        return r.json()


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


def esummary(db: str, ids: IterableT[str] = None, **kwargs) -> dict:
    r = eutil(util='esummary.fcgi', db=db, ids=ids, version='1.0',
              retmode='xml', **kwargs)
    if r is not None:
        return r.text


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


def _elink_ok(dbfrom: str, ids: IterableT[str],
              linkname: str) -> tuple:
    ids = gis(dbfrom, ids)
    links_available = elink(db='', dbfrom=dbfrom, cmd='acheck', ids=ids,
                            linkname=linkname)
    ret_val = list()
    for ok in links_available:
        idchecklist = ok['idchecklist']
        idlinksets = idchecklist['idlinksets']
        for idlinkset in idlinksets:
            acc = idlinkset['id']
            for linkinfo in idlinkset['linkinfos']:
                if linkname == linkinfo['linkname']:
                    ret_val.append(acc)
                    break

    return tuple(ret_val)


def _elink_parse(db: str, dbfrom: str, linksets: IterableT[dict]) -> tuple:
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


def _elink_runner(db: str, dbfrom: str,
                  ids: IterableT[str], linkname: str = None) -> tuple:
    if type(db) not in (str, ) or type(dbfrom) not in (str, ):
        raise Exception('Parameters db and dbfrom must be strings.')
    if linkname is None:
        linkname = dbfrom + '_' + db
    ids_ok = _elink_ok(dbfrom=dbfrom, ids=ids, linkname=linkname)
    linksets = elink(db=db, dbfrom=dbfrom, cmd='neighbor', linkname=linkname,
                     ids=ids_ok)
    ids_from, ids_to = _elink_parse(db, dbfrom, linksets)
    return ids_from, ids_to


def _history_server_data_ok(data):
    if data['query_key'] is None or data['web_env'] is None:
        return False
    else:
        return True


# Secondary functions --------------------------------------------------------
def search(db: str, term: str) -> dict:
    ret_dict = {'db': db, 'query_keys': [], 'web_env': None, 'count': 0}
    esearchresult = None
    for i in range(5):
        json = esearch(db=db, term=term, usehistory=True)
        if json is not None:
            esearchresult = json['esearchresult']
            if 'count' not in esearchresult:
                esearchresult = None
                continue
            else:
                break

    if esearchresult is not None:
        count = int(esearchresult['count'])
        if count == 0:
            return ret_dict
        ret_dict['query_keys'] = [int(esearchresult['querykey'])]
        ret_dict['web_env'] = esearchresult['webenv']

    return ret_dict


@dispatch(dict, namespace=EUTILS_NS)
@with_history
def summary(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = esummary(**data)
    return parse_esummary_xml_text(txt)


@dispatch(str, Iterable, namespace=EUTILS_NS)
def summary(db: str, ids: IterableT[str]) -> list:
    txt = esummary(db=db, ids=ids)
    return parse_esummary_xml_text(txt)


@dispatch(dict, namespace=EUTILS_NS)
@with_history
def accs(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = efetch(rettype='uilist', retmode='text', **data)
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(str, Iterable, namespace=EUTILS_NS)
def accs(db: str, ids: IterableT[str]) -> list:
    txt = efetch(db=db, ids=ids, rettype='uilist', retmode='text')
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(dict, namespace=EUTILS_NS)
@with_history
def gis(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = efetch(rettype='uilist', retmode='text', idtype=None, **data)
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(str, Iterable, namespace=EUTILS_NS)
def gis(db: str, ids: IterableT[str]) -> list:
    txt = efetch(db=db, ids=ids, rettype='uilist', retmode='text',
                 idtype=None)
    ret_list = _strip_split(txt)
    return ret_list


@dispatch(str, Iterable, namespace=EUTILS_NS)
def taxids(db: str, ids: IterableT[str]) -> dict:
    dbfrom = db
    db = 'taxonomy'
    ids_from, ids_to = _elink_runner(db=db, dbfrom=dbfrom, ids=ids)
    ids_to = [int(x) for x in ids_to]
    return dict(zip(ids_from, ids_to))


@dispatch(dict, namespace=EUTILS_NS)
def taxids(data: dict) -> dict:
    ids = accs(data)
    return taxids(data['db'], ids)


def _process_downloaded_seq_data(efetch_txt: str, db: str, rettype: str,
                                 retmode: str):
    rec_list = list()
    if rettype == 'gb' and retmode == 'xml':
        rec_list = seq_records_gb(efetch_txt)
    elif rettype == 'fasta' and retmode == 'text':
        seq_type = None
        if db == 'nuccore':
            seq_type = SEQ_TYPE_NT
        elif db == 'protein':
            seq_type = SEQ_TYPE_AA
        rec_list = read_fasta(StringIO(efetch_txt), seq_type, parse_def=True)
    return rec_list


@dispatch(dict, namespace=EUTILS_NS)
@with_history
def seqs(data: dict, rettype: str = 'gb', retmode: str = 'xml') -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = efetch(rettype=rettype, retmode=retmode, **data)
    db = data['db']
    return _process_downloaded_seq_data(txt, db, rettype, retmode)


@dispatch(str, Iterable, namespace=EUTILS_NS)
def seqs(db: str, ids: IterableT[str], rettype: str = 'gb',
         retmode: str = 'xml', location: dict = None) -> list:
    if location is None:
        ids_requested = set(ids)
        r_temp = list()
        epost_r = epost(db, ids)
        r = seqs(epost_r, rettype=rettype, retmode=retmode)
        r_temp += r
        if len(r) != len(ids_requested):
            # DEBUG
            # print('DEBUG : eutils.seqs : len(r) != len(ids_requested) ::',
            #       len(r), len(ids_requested))
            # END DEBUG
            ids_dnld = set([x.accession_version for x in r_temp])
            ids_remaining = ids_requested - ids_dnld
            r = seqs(db, list(ids_remaining), rettype=rettype, retmode=retmode)
            r_temp += r
            return r_temp
        else:
            return r_temp
    else:
        txt = efetch(db=db, ids=ids, rettype=rettype, retmode=retmode,
                     location=location)
        return _process_downloaded_seq_data(txt, db, rettype, retmode)


@dispatch(Iterable, namespace=EUTILS_NS)
def cds(ids_protein: IterableT[str]) -> list:
    prot_recs = seqs('protein', ids_protein)
    ret_list = list()
    for prot_rec in prot_recs:
        prot_def = prot_rec.definition
        prot_acc_ver = prot_rec.accession_version
        if prot_rec.coded_by is not None:
            cds_loc_str = eutils_loc_str(prot_rec.coded_by)
            cds_acc = prot_rec.coded_by['external_ref']

            cds_rec = seqs('nuccore', [cds_acc], rettype='fasta',
                           retmode='text', location=cds_loc_str)

            if len(cds_rec) > 0:
                cds_rec = cds_rec[0]
            else:
                continue

            cds_def = cds_rec.definition
            cds_acc_ver = cds_rec.accession_version
            cds_def_new = 'CDS for: aa_acc[' + prot_acc_ver + '] ' + \
                          prot_def[0].title() + prot_def[1:] + \
                          '. Extracted from: nt_acc[' + cds_acc_ver + '] ' + \
                          cds_def[0].title() + cds_def[1:] + '.'
            cds_rec.definition = cds_def_new
            cds_rec.accession = prot_rec.accession
            cds_rec.version = prot_rec.version
            ret_list.append(cds_rec)
            sleep(0.5)
    return ret_list


@dispatch(dict, namespace=EUTILS_NS)
def cds(data_protein: dict) -> dict:
    ok = _history_server_data_ok(data_protein)
    if ok is False:
        return list()
    assert data_protein['db'] == 'protein'
    ids_protein = accs(data_protein)
    return cds(ids_protein)


@dispatch(Iterable, namespace=EUTILS_NS)
def sra_run_info(ids_srr: IterableT[str]) -> list:
    efetch_csv_txt = efetch(db='sra', ids=ids_srr, rettype='runinfo',
                            retmode='csv')
    ret_list = parse_efetch_sra_csv_text(efetch_csv_txt)
    ret_list = [x for x in ret_list if x['Run'] in ids_srr]
    return ret_list


@dispatch(str, namespace=EUTILS_NS)
def sra_run_info(id_srr: str) -> list:
    return sra_run_info([id_srr])
