# Checks for the Horizons small bodies: interpolation error between samples, agreement with the SBDB two-body orbit
# near its epoch, and how far two-body motion drifts by 2050, 2075 and 2100.
import sys, json
sys.path.insert(0, '/work/flyby')
import numpy as np, pykep as pk
import bodies

atlas = bodies.load_atlas()
for name in ['ceres', 'vesta', 'psyche', '433', '101955', '162173']:
    H = bodies.HorizonsBody(name); K = bodies.small_body(name, atlas)
    # interpolation: midpoint between two samples against a Horizons sample treated as unknown (skip one row)
    rows = H.rows; errs = []
    for i in range(100, 3000, 300):
        a, m, b = rows[i], rows[i + 1], rows[i + 2]
        sub = bodies.HorizonsBody(name); sub.rows = [a, b]; sub.t0 = a[0] - 2451544.5; sub.h = b[0] - a[0]
        r = np.array(sub.eph(m[0] - 2451544.5)[0]) / pk.AU
        errs.append(np.linalg.norm(r - np.array(m[1:4])) * pk.AU / 1000)
    el = {b['id']: b for b in atlas['bodies']}[name]['orbit']['elements']
    ep = el['epoch_jd_tdb'] - 2451544.5
    drift = []
    for jd in (ep + 2451544.5, 2469807.5, 2478938.5, 2488069.5):        # epoch, 2050, 2075, 2100
        t = jd - 2451544.5
        drift.append(np.linalg.norm(np.array(H.eph(t)[0]) - np.array(K.eph(t)[0])) / pk.AU)
    print(f'{name:7s} interpolation error over a 20-day gap: max {max(errs):8.1f} km | two-body minus Horizons, AU: '
          f'epoch {drift[0]:.5f}, 2050 {drift[1]:.4f}, 2075 {drift[2]:.4f}, 2100 {drift[3]:.4f}')
