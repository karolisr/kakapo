from kakapo.tools.eutils import efetch
from kakapo.tools.eutils import einfo
from kakapo.tools.eutils import elink
from kakapo.tools.eutils import epost
from kakapo.tools.eutils import esearch
from kakapo.tools.eutils import espell
from kakapo.tools.eutils import esummary


def test_einfo():
    r = einfo(db='protein')
    assert 'dbinfo' in r['einforesult']


def test_esearch():
    r = esearch(db='protein',
                term='("ribonuclease t2"[title] AND "RNase"[title])')
    assert 'XP_001384105.2' in r['esearchresult']['idlist']


def test_epost():
    r = epost(db='protein', ids=['1746590493'])
    assert r['query_keys'][0] == 1


def test_esummary():
    r = esummary(db='protein', ids=['15718680', 'NP_001098858.1', '119703751'])
    assert r['result']['119703751']['accessionversion'] == 'NP_034713.2'


def test_efetch():
    r = efetch(db='protein', ids=['15718680', 'NP_001098858.1', '119703751'],
               rettype='acc', retmode='text')
    assert r == 'NP_005537.3\nNP_001098858.1\nNP_034713.2\n'


def test_elink():
    r = elink(db='nuccore', dbfrom='protein', linkname='protein_nuccore',
              cmd='neighbor', ids=['GER25982.1'])
    assert r[0]['linksetdbs'][0]['links'] == ['BKCP01000669.1']


def test_espell():
    r = espell(db='taxonomy', term='Arabidois OR Solanm')
    assert r == 'arabidopsis or solanum'
