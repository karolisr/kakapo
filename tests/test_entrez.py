from kakapo.tools.eutils import esearch


def test_esearch():
    esearch_results = esearch(term='AAB87127.1',
                              db='protein')

    assert esearch_results['db'] == 'protein'
    assert isinstance(esearch_results['query_keys'], list)
    assert isinstance(esearch_results['web_env'], str)
