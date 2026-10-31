"""Comparison service: journeys for several ships between two sites, impulsive and low thrust side by side.

Run (in the solmaps-env container, from /work): uvicorn engine.api:app --host 0.0.0.0 --port 8077
Then open http://localhost:8077/
"""
import datetime as dt, pathlib, sys
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
from typing import Optional
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field

from engine.conventions import SITES, NAMES, FUEL_SOURCES, end_state, export as conventions
from engine.data import Alternatives
from engine.journey import Network, plan, max_payload
from engine.ship import Ship, PRESETS

ROOT = pathlib.Path(__file__).resolve().parent.parent
app = FastAPI(title='solmaps comparison')
ALTS = Alternatives(ROOT / 'data' / 'alternatives.npz')
NET = Network(ALTS)
JD_UNIX = 2440587.5

def jd_of(date):
    d = dt.date.fromisoformat(date)
    return JD_UNIX + (d - dt.date(1970, 1, 1)).days

def iso(jd):
    return (dt.date(1970, 1, 1) + dt.timedelta(days=jd - JD_UNIX)).isoformat()

class ShipIn(BaseModel):
    name: str
    dry_kg: float = Field(gt=0)
    payload_kg: float = Field(ge=0)
    fuel_capacity_kg: float = Field(gt=0)
    thrust_n: float = Field(gt=0)
    isp_s: float = Field(gt=0)
    propulsion: str = 'chemical'

class Query(BaseModel):
    origin: str
    destination: str
    depart_from: str
    depart_until: str
    arrive_by: Optional[str] = None
    ships: list[ShipIn]
    refuel_at: Optional[list[str]] = None
    max_legs: int = Field(default=3, ge=1, le=5)

@app.get('/api/meta')
def meta():
    nodes = SITES + ['eml1']
    return {'nodes': [{'id': n, 'name': NAMES[n], 'end_state': end_state(n)['label'], 'fuel_source': n in FUEL_SOURCES} for n in nodes],
            'presets': {k: v.to_dict() for k, v in PRESETS.items()},
            'dates': {'first': iso(ALTS.jd(0)), 'last': iso(ALTS.jd(ALTS.n_dates - 1))}, 'conventions': conventions()}

def impulsive_rows(q, ship):
    res = plan(ALTS, ship, q.origin, q.destination, jd_of(q.depart_from), jd_of(q.depart_until),
               arrive_by_jd=jd_of(q.arrive_by) if q.arrive_by else None, refuel_at=q.refuel_at, max_legs=q.max_legs, net=NET)
    rows = []
    for j in res:
        rows.append({'model': 'impulsive', 'depart': iso(j['depart_jd']), 'arrive': iso(j['arrive_jd']), 'days': j['days'],
                     'dv_total': j['dv_total'], 'fuel_loaded_kg': j['fuel_loaded_kg'], 'fuel_left_kg': j['fuel_left_kg'],
                     'payload_delivered_kg': j['payload_delivered_kg'], 'max_payload_kg': round(max_payload(ship, j)),
                     'legs': [{**l, 'depart': iso(l['depart_jd']), 'arrive': iso(l['arrive_jd'])} for l in j['legs']],
                     'warnings': sorted({w for l in j['legs'] for w in l['warnings']})})
    return rows

def lowthrust_rows(q, ship):
    from engine import refine
    t0 = (jd_of(q.depart_from) - 2451544.5, jd_of(q.depart_until) - 2451544.5)
    last = jd_of(q.arrive_by) - jd_of(q.depart_until) if q.arrive_by else 1500
    if last < 60:
        return [{'model': 'low thrust', 'feasible': False, 'reason': 'deadline leaves less than 60 days after the window'}]
    r = refine.lowthrust_transfer(q.origin, q.destination, ship, t0, (60, last))
    if not r.get('leg'):
        return [{'model': 'low thrust', 'feasible': False, 'reason': r.get('reason', 'no feasible transfer')}]
    leg = r['leg']
    return [{'model': 'low thrust', 'feasible': r['feasible'], 'reason': r.get('reason', ''),
             'depart': iso(r['depart_mjd2000'] + 2451544.5), 'arrive': iso(r['arrive_mjd2000'] + 2451544.5),
             'days': round(r['arrive_mjd2000'] - r['depart_mjd2000'], 1), 'dv_total': round(sum(p['dv_kms'] for p in r['phases']), 3),
             'fuel_loaded_kg': ship.fuel_capacity_kg, 'fuel_left_kg': round(ship.fuel_capacity_kg - r['fuel_kg'], 1),
             'payload_delivered_kg': r['payload_delivered_kg'], 'phases': r['phases'], 'cached': r['cached'],
             'check': leg['check'], 'vinf_slack_ms': [round(leg['vinf_dep_ms'], 1), round(leg['vinf_arr_ms'], 1)]}]

@app.post('/api/compare')
def compare(q: Query):
    nodes = SITES + ['eml1']
    if q.origin not in nodes or q.destination not in nodes or q.origin == q.destination:
        raise HTTPException(400, 'choose two different sites')
    if jd_of(q.depart_until) < jd_of(q.depart_from):
        raise HTTPException(400, 'the departure window ends before it starts')
    out = []
    for s in q.ships:
        ship = Ship(**s.model_dump())
        if ship.low_thrust and (q.origin == 'eml1' or q.destination == 'eml1'):
            rows = [{'model': 'low thrust', 'feasible': False, 'reason': 'low-thrust transfers start and end at sites, not at L1'}]
        else:
            rows = lowthrust_rows(q, ship) if ship.low_thrust else impulsive_rows(q, ship)
        out.append({'ship': ship.to_dict(), 'dv_capacity_kms': round(ship.dv_capacity_kms(), 3), 'results': rows})
    return {'query': q.model_dump(), 'ships': out}

@app.get('/', response_class=HTMLResponse)
def page():
    return (ROOT / 'engine' / 'ui' / 'index.html').read_text()
