from kakapo.utils.http import get
from requests.exceptions import RetryError


class TestGet:
    def test_get_retry_fail(self):
        try:
            r = get(url='http://httpbin.org/status/500',
                    response_format='text')
        except RetryError as e:
            pass

    def test_get_success(self):
        r = get(url='http://httpbin.org/status/200',
                response_format='text')
        assert r.status_code == 200

    def test_get_eventual_success(self):
        r = get(url='http://httpbin.org/status/200,500',
                response_format='text')
        assert r.status_code == 200


def test_post():
    assert True
