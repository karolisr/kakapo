from kakapo.tools.eutils import esearch


def test_esearch():
    esearch_results = esearch(term='AAB87127.1',
                              db='protein',
                              api_key=None,
                              ret_type='uilist')

    assert esearch_results['db'] == 'protein'
    assert isinstance(esearch_results['counts'], list)
    assert isinstance(esearch_results['query_keys'], list)
    assert isinstance(esearch_results['web_env'], str)
    assert esearch_results['counts'][0] == 1
