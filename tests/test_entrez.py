from kakapo.tools.eutils import search


def test_esearch():
    esearch_results = search(db='protein', term='AAB87127.1')

    assert esearch_results['db'] == 'protein'
    assert isinstance(esearch_results['query_keys'], list)
    assert isinstance(esearch_results['web_env'], str)
