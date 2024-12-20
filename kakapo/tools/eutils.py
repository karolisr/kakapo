"""
Wrap NCBI's Entrez Programming Utilities (E-utilities).

More information on E-utilities at:
    http://www.ncbi.nlm.nih.gov/books/NBK25497

Entrez Unique Identifiers (UIDs) for selected databases:
    http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

Valid values of &retmode and &rettype for EFetch (null = empty string):
    https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
"""

import os
from collections.abc import Iterable
from copy import deepcopy
from functools import wraps
from io import StringIO
from sys import exit
from time import sleep
from typing import Any, Callable, Union
from xml.etree import ElementTree

from requests import Response

from kakapo.tools.bioio import read_fasta
from kakapo.tools.parsers import (eutils_loc_str, parse_efetch_sra_xml_text,
                                  parse_esummary_xml_text, seq_records_gb)
from kakapo.tools.seq import SEQ_TYPE_AA, SEQ_TYPE_NT, SeqRecord
from kakapo.utils.http import post
from kakapo.utils.logging import Log

EUTILS_NS = dict()
ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def find_api_key() -> Union[str, None]:
    api_key: Union[str, None] = None
    global_variables: dict[str, str] = globals()

    if 'ENTREZ_KEY' in global_variables:
        api_key = global_variables['ENTREZ_KEY']
        # if api_key != '':
        #     Log.msg('ENTREZ_KEY found in global_variables', api_key)

    elif 'ENTREZ_KEY' in os.environ:
        api_key = os.environ['ENTREZ_KEY']
        # if api_key != '':
        #     Log.msg('ENTREZ_KEY found in os.environ', api_key)

    # ToDo: make this an exception?
    if api_key == '' or api_key is None:
        api_key = None
        Log.err('Warning:', 'ENTREZ_KEY is not defined.')
    else:
        global_variables['ENTREZ_KEY'] = api_key

    return api_key


def ncbi_api_key(f: Callable) -> Callable:
    @wraps(f)
    def wrapper(api_key: Union[str, None] = None, *args, **kwargs) -> Callable:
        if api_key is None:
            api_key = find_api_key()
        if api_key == '':
            api_key = None
        return f(api_key=api_key, *args, **kwargs)
    return wrapper


def history(f: Callable) -> Callable:
    @wraps(f)
    def wrapper(usehistory: bool = False, web_env: Union[str, None] = None,
                query_key: Union[str, None] = None, *args, **kwargs
                ) -> Callable:
        uh: Union[str, None] = None
        if usehistory is True or web_env is not None:
            uh = 'y'
        return f(usehistory=uh, web_env=web_env, query_key=query_key, *args,
                 **kwargs)
    return wrapper


def with_history(f: Callable) -> Callable:
    @wraps(f)
    def wrapper(data: dict, *args, **kwargs) -> Any:

        db: str = data['db']
        query_keys: list = data['query_keys']
        web_env: str = data['web_env']

        new_data = {
            'db': db, 'query_key': None, 'web_env': web_env}

        if len(query_keys) == 0:
            return f(new_data, *args, **kwargs)

        elif len(query_keys) == 1:
            new_data['query_key'] = query_keys[0]
            return f(new_data, *args, **kwargs)

        elif len(query_keys) > 1:
            ret_list = list()
            for k in query_keys:
                new_data['query_key'] = k
                ret_list += f(new_data, *args, **kwargs)
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

    url: str = ENTREZ_BASE_URL + util

    ids_joined: Union[str, None] = None
    if ids is not None:
        assert isinstance(ids, Iterable)
        ids_joined = ','.join(ids)

    params = {'api_key': api_key, 'bdata': bdata, 'cmd': cmd,
              'complexity': complexity, 'datetype': datetype, 'db': db,
              'dbfrom': dbfrom, 'email': email, 'field': field,
              'holding': holding, 'id': ids_joined, 'idtype': idtype,
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
def einfo(db: str) -> Union[dict, None]:
    r: Response = eutil(util='einfo.fcgi', db=db, version='2.0',
                        retmode='json')
    if r is not None:
        return r.json()


def esearch(db: str, term: str, usehistory: bool = False,
            web_env: Union[str, None] = None,
            query_key: Union[int, None] = None, **kwargs) -> Union[dict, None]:
    r = eutil(util='esearch.fcgi', rettype='uilist', retmode='json',
              db=db, term=term, usehistory=usehistory, web_env=web_env,
              query_key=query_key, **kwargs)
    if r is not None:
        return r.json()


def epost(db: str, ids: Union[Iterable[str], None] = None,
          web_env: Union[str, None] = None,
          **kwargs) -> dict:
    r: Response = eutil(util='epost.fcgi', db=db, ids=ids, web_env=web_env,
                        retmode='xml', **kwargs)

    query_keys = list()
    web_env = None

    if r is not None:
        root = ElementTree.fromstring(r.text)
        for child in root:
            if child.tag == 'QueryKey':
                assert child.text is not None
                query_keys.append(int(child.text))
            if child.tag == 'WebEnv':
                web_env = child.text

    return {'db': db, 'query_keys': query_keys, 'web_env': web_env}


def esummary(db: str, ids: Union[Iterable[Union[str, int]], None] = None,
             **kwargs) -> Union[str, None]:
    if ids is None:
        return None
    ids = [str(id) for id in ids]
    r: Response = eutil(util='esummary.fcgi', db=db, ids=ids, version='1.0',
                        retmode='xml', **kwargs)
    if r is not None:
        return r.text


def efetch(db: str, ids: Union[Iterable[str], None] = None,
           rettype: Union[str, None] = None,
           retmode: Union[str, None] = None, **kwargs) -> Union[str, None]:
    r: Response = eutil(util='efetch.fcgi', db=db, ids=ids, rettype=rettype,
                        retmode=retmode, **kwargs)
    if r is not None:
        return r.text


def elink(db: str, dbfrom: str, cmd: str,
          ids: Union[Iterable[str], None] = None, **kwargs
          ) -> Union[list, None]:
    ids_separate = None
    if ids is not None:
        ids_separate = [('id', x) for x in ids]
    r: Response = eutil(util='elink.fcgi', db=db, dbfrom=dbfrom, cmd=cmd,
                        ids_separate=ids_separate, retmode='json', **kwargs)

    if r is not None:
        r_dict = r.json()
        linksets = r_dict['linksets']
        return linksets


def espell(db: str, term: str, **kwargs) -> Union[str, None]:
    r: Response = eutil(util='espell.fcgi', db=db, term=term, retmode='xml',
                        **kwargs)

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


def _elink_ok(dbfrom: str, ids: Iterable[str],
              linkname: str) -> tuple:
    ids = gis(dbfrom, ids)
    links_available: Union[list, None] = elink(
        db='', dbfrom=dbfrom, cmd='acheck', ids=ids, linkname=linkname)
    ret_val = list()
    if links_available is not None:
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


def _elink_parse(db: str, dbfrom: str, linksets: Iterable[dict]) -> tuple:
    ids_from = []
    ids_to = []
    for ls in linksets:
        # if 'linksetdbs' not in ls:
        #     continue
        if dbfrom == ls['dbfrom']:
            ids_from += ls['ids']
            linksetdbs = ls['linksetdbs']
            for lsdb in linksetdbs:
                if db == lsdb['dbto']:
                    ids_to += lsdb['links']
    return ids_from, ids_to


def _elink_runner(db: str, dbfrom: str, ids: Iterable[str],
                  linkname: Union[str, None] = None) -> tuple:
    if type(db) not in (str, ) or type(dbfrom) not in (str, ):
        raise Exception('Parameters db and dbfrom must be strings.')
    if linkname is None:
        linkname = dbfrom + '_' + db
    ids_ok: tuple = _elink_ok(dbfrom=dbfrom, ids=ids, linkname=linkname)
    ids_from = []
    ids_to = []
    linksets: Union[list, None] = elink(db=db, dbfrom=dbfrom, cmd='neighbor',
                                        linkname=linkname, ids=ids_ok)
    if linksets is not None:
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
        if json is not None and 'esearchresult' in json:
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


@with_history
def summary_with_data(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = esummary(**data)
    if txt is None:
        Log.err('kakapo experienced a problem downloading data from NCBI.', '')
        exit(0)
    return parse_esummary_xml_text(txt)


def summary(db: str, ids: Iterable[Union[str, int]]) -> list:
    txt = esummary(db=db, ids=ids)
    if txt is None:
        ids = [str(id) for id in ids]
        Log.err(f'kakapo experienced a problem downloading data from an NCBI'
                f' database {db}:',
                f' {",".join(ids)}')
        exit(0)
    return parse_esummary_xml_text(txt)


@with_history
def accs_with_data(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = efetch(rettype='uilist', retmode='text', **data)
    ret_list = []
    if txt is not None:
        ret_list = _strip_split(txt)
    return ret_list


def accs(db: str, ids: Iterable[str]) -> list:
    txt = efetch(db=db, ids=ids, rettype='uilist', retmode='text')
    ret_list = []
    if txt is not None:
        ret_list = _strip_split(txt)
    return ret_list


@with_history
def gis_with_data(data: dict) -> list:
    ok = _history_server_data_ok(data)
    if ok is False:
        return list()
    txt = efetch(rettype='uilist', retmode='text', idtype=None, **data)
    ret_list = []
    if txt is not None:
        ret_list = _strip_split(txt)
    return ret_list


def gis(db: str, ids: Iterable[str]) -> list:
    txt = efetch(db=db, ids=ids, rettype='uilist', retmode='text',
                 idtype=None)
    ret_list = []
    if txt is not None:
        ret_list = _strip_split(txt)
    return ret_list


def taxids(db: str, ids: Iterable[str]) -> dict:
    dbfrom = db
    db = 'taxonomy'
    ids_from, ids_to = _elink_runner(db=db, dbfrom=dbfrom, ids=ids)
    ids_to = [x for x in ids_to]
    return dict(zip(ids_from, ids_to))


def taxids_with_data(data: dict) -> dict:
    ids = accs_with_data(data)
    return taxids(data['db'], ids)


def _process_downloaded_seq_data(efetch_txt: str, db: str, rettype: str,
                                 retmode: str):
    rec_list: list[SeqRecord] = list()
    if rettype == 'gb' and retmode == 'xml':
        rec_list = seq_records_gb(efetch_txt)
    elif rettype == 'fasta' and retmode == 'text':
        seq_type = None
        if db == 'nuccore':
            seq_type = SEQ_TYPE_NT
        elif db == 'protein':
            seq_type = SEQ_TYPE_AA
        _ = read_fasta(StringIO(efetch_txt), seq_type, parse_def=True)
        for r in _:
            assert isinstance(r, SeqRecord)
            rec_list.append(r)
    return rec_list


@with_history
def seqs_with_data(data: dict, rettype: str = 'gb', retmode: str = 'xml'
                   ) -> list[SeqRecord]:
    ret_list: list[SeqRecord] = []
    ok = _history_server_data_ok(data)
    if ok is False:
        return ret_list
    txt = efetch(rettype=rettype, retmode=retmode, **data)
    db = data['db']
    if txt is not None:
        ret_list = _process_downloaded_seq_data(txt, db, rettype, retmode)
    return ret_list


def seqs(db: str, ids: Iterable[str], rettype: str = 'gb',
         retmode: str = 'xml', location: Union[str, None] = None
         ) -> list[SeqRecord]:
    if location is None:
        ids_requested = set(ids)
        r_temp = list()
        epost_r = epost(db, ids)
        r = seqs_with_data(data=epost_r, rettype=rettype, retmode=retmode)
        r_temp += r
        if len(r) != len(ids_requested):
            # DEBUG
            # print('DEBUG : eutils.seqs : len(r) != len(ids_requested) ::',
            #       len(r), len(ids_requested))
            # END DEBUG
            ids_dnld = set([x.accession_version for x in r_temp])
            ids_remaining = ids_requested - ids_dnld

            # ----------------------------------------------------------------
            # Check if accession is actually in the specified database,
            #   otherwise this would recurse to infinity!
            # Examples:
            # eutils.seqs('protein', ['QYF06680.1'])  # OK
            # eutils.seqs('nuccore', ['QYF06680.1'])  # KO
            for id in deepcopy(ids_remaining):
                check_id = len(
                    summary_with_data(search(db=db, term=id))) == 1
                if check_id is False:
                    print(f'Record {id} was not found in the database "{db}".')
                    ids_remaining.remove(id)
            # ----------------------------------------------------------------
            r = seqs(db, list(ids_remaining), rettype=rettype,
                     retmode=retmode)
            r_temp += r
            return r_temp
        else:
            return r_temp
    else:
        txt = efetch(db=db, ids=ids, rettype=rettype,
                     retmode=retmode, location=location)
        if txt is not None:
            return _process_downloaded_seq_data(txt, db, rettype, retmode)
        return []


def cds(ids_protein: Iterable[str]) -> list:
    prot_recs = seqs('protein', ids_protein)
    ret_list = list()
    for prot_rec in prot_recs:
        prot_def = prot_rec.definition
        prot_acc_ver: Union[str, None] = prot_rec.accession_version
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

            assert prot_acc_ver is not None
            assert cds_acc_ver is not None

            cds_def_new = 'CDS for: aa_acc[' + prot_acc_ver + '] ' \
                + prot_def[0].title() + prot_def[1:] \
                + '. Extracted from: nt_acc[' + cds_acc_ver + '] ' \
                + cds_def[0].title() + cds_def[1:] + '.'

            cds_rec.definition = cds_def_new
            cds_rec.accession = prot_rec.accession
            cds_rec.version = prot_rec.version
            ret_list.append(cds_rec)
            sleep(0.5)

    return ret_list


def cds_with_data(data_protein: dict) -> list:
    ok = _history_server_data_ok(data_protein)
    if ok is False:
        return list()
    assert data_protein['db'] == 'protein'
    ids_protein = accs_with_data(data_protein)
    return cds(ids_protein)


def sra_run_info(srr: Union[Iterable[str], str]) -> list:
    if isinstance(srr, str):
        srr = [srr]
    # CSV output example: https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&acc=SRR13805638,SRR13805642
    efetch_xml_txt = efetch(db='sra', ids=srr, rettype='runinfo',
                            retmode='xml')
    ret_list = parse_efetch_sra_xml_text(efetch_xml_txt)
    ret_list = [x for x in ret_list if x['Run'] in srr]
    return ret_list
