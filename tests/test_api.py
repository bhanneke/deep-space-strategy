"""The comparison service: one query with an impulsive and an electric ship, validation of bad input."""
import pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'flyby'))
import pytest
pytest.importorskip('fastapi')
pytestmark = pytest.mark.skipif(not (ROOT / 'data' / 'alternatives.npz').exists(), reason='run make data first')

@pytest.fixture(scope='module')
def client():
    from fastapi.testclient import TestClient
    from engine.api import app
    return TestClient(app)

def test_meta(client):
    m = client.get('/api/meta').json()
    assert {'earth', 'ceres', 'eml1'} <= {n['id'] for n in m['nodes']}
    assert 'solar-electric cargo' in m['presets']

def test_compare_impulsive_and_electric(client):
    from engine.ship import PRESETS
    q = {'origin': 'earth', 'destination': 'mars', 'depart_from': '2029-01-01', 'depart_until': '2029-11-30', 'max_legs': 2,
         'ships': [PRESETS['methalox tug'].to_dict(), PRESETS['solar-electric cargo'].to_dict()]}
    r = client.post('/api/compare', json=q)
    assert r.status_code == 200
    imp, lt = r.json()['ships']
    assert imp['results'] and imp['results'][0]['model'] == 'impulsive'
    for x in imp['results']:
        assert x['max_payload_kg'] >= x['payload_delivered_kg'] - 1
        assert x['depart'] >= '2029-01-01'
    rows = lt['results']
    for x in rows:
        assert x['model'] == 'continuous thrust' and [p['kind'] for p in x['phases']] == ['planet-side', 'heliocentric', 'planet-side']
        assert x['check']['position_error_km'] < 50
        assert x['depart'] >= '2029-01-01'               # the window bounds the departure from the parking orbit
    for a, b in zip(rows, rows[1:]):                # alternatives: each slower option uses less fuel
        assert a['days'] <= b['days'] and a['fuel_left_kg'] < b['fuel_left_kg']

def test_bad_input_is_rejected(client):
    q = {'origin': 'earth', 'destination': 'earth', 'depart_from': '2029-01-01', 'depart_until': '2029-02-01', 'ships': []}
    assert client.post('/api/compare', json=q).status_code == 400
    q.update(destination='mars', depart_until='2028-01-01')
    assert client.post('/api/compare', json=q).status_code == 400

def test_page_is_served(client):
    assert 'Ship Comparison' in client.get('/').text
