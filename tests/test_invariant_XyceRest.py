import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'utils', 'XyceCInterface'))
from utils.XyceCInterface.XyceRest import app


@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client


@pytest.mark.parametrize("method,endpoint,headers", [
    ("POST", "/xyce_open", {}),
    ("POST", "/xyce_initialize", {"Authorization": "Bearer expired.token.value"}),
    ("POST", "/xyce_getcircuitvalue", {"Authorization": "Bearer malformed"}),
    ("GET", "/status", {"Authorization": ""}),
    ("POST", "/xyce_getdacnames", {"Authorization": "InvalidScheme abc123"}),
])
def test_unauthenticated_requests_are_rejected(client, method, endpoint, headers):
    """Invariant: Protected endpoints must reject unauthenticated requests with 401 or 403."""
    if method == "POST":
        response = client.post(endpoint, headers=headers, json={})
    else:
        response = client.get(endpoint, headers=headers)

    assert response.status_code in (401, 403), (
        f"Endpoint {endpoint} returned {response.status_code} without valid auth — "
        f"expected 401 or 403. Unauthenticated access must be denied."
    )