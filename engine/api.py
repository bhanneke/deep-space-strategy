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
    return (dt.date(1970, 1, 1) + dt.timedelta(days=round(jd - JD_UNIX, 6))).isoformat()   # rounding keeps a date at a bound

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

CAPS = (1.0, 0.6, 0.4, 0.25, 0.15, 0.1)   # continuous thrust: flight-time caps as shares of the time to the deadline (at least 60 days)

def lowthrust_rows(q, ship):
    """continuous-thrust alternatives: the fuel-optimal transfer under each flight-time cap, kept when no faster
    alternative uses less fuel (Pareto on days and fuel), fastest first"""
    from engine import refine
    t0 = (jd_of(q.depart_from) - 2451544.5, jd_of(q.depart_until) - 2451544.5)
    last = jd_of(q.arrive_by) - jd_of(q.depart_until) if q.arrive_by else 1500
    if last < 60:
        return [{'model': 'continuous thrust', 'feasible': False, 'reason': 'deadline leaves less than 60 days after the window'}]
    found, reason = [], 'no feasible transfer'
    for cap in sorted({max(60.0, round(last * c)) for c in CAPS}):
        r = refine.lowthrust_transfer(q.origin, q.destination, ship, t0, (min(60.0, 0.5 * cap), cap))
        if r.get('leg') and r['feasible']:
            found.append(r)
        else:
            reason = r.get('reason') or reason
    if not found:
        return [{'model': 'continuous thrust', 'feasible': False, 'reason': reason}]
    days = lambda r: r['arrive_mjd2000'] - r['depart_mjd2000']
    found.sort(key=lambda r: (days(r), r['fuel_kg']))
    front = []
    for r in found:
        if not front or r['fuel_kg'] < 0.995 * front[-1]['fuel_kg']:     # slower only when it saves at least 0.5% fuel
            front.append(r)
    rows = []
    for r in front:
        leg = r['leg']
        rows.append({'model': 'continuous thrust', 'feasible': True, 'reason': '',
                     'depart': iso(r['depart_mjd2000'] + 2451544.5), 'arrive': iso(r['arrive_mjd2000'] + 2451544.5),
                     'days': round(days(r), 1), 'dv_total': round(sum(p['dv_kms'] for p in r['phases']), 3),
                     'fuel_loaded_kg': ship.fuel_capacity_kg, 'fuel_left_kg': round(ship.fuel_capacity_kg - r['fuel_kg'], 1),
                     'payload_delivered_kg': r['payload_delivered_kg'], 'phases': r['phases'], 'cached': r['cached'],
                     'check': leg['check'], 'vinf_slack_ms': [round(leg['vinf_dep_ms'], 1), round(leg['vinf_arr_ms'], 1)]})
    return rows

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
        rows = impulsive_rows(q, ship)               # short burns: valid legs only (journey.IMPULSIVE_SHARE)
        if ship.low_thrust:
            if q.origin == 'eml1' or q.destination == 'eml1':
                rows += [{'model': 'continuous thrust', 'feasible': False, 'reason': 'continuous-thrust transfers start and end at sites, not at L1'}]
            else:
                rows += lowthrust_rows(q, ship)
        out.append({'ship': ship.to_dict(), 'dv_capacity_kms': round(ship.dv_capacity_kms(), 3), 'results': rows})
    return {'query': q.model_dump(), 'ships': out}

@app.get('/', response_class=HTMLResponse)
def page():
    return (ROOT / 'engine' / 'ui' / 'index.html').read_text()
