# Check DE440 bodies: positions match jpl_lp before 2050, a flyby problem builds, Lambert speed.
import time, sys
sys.path.insert(0, '/work/flyby')
import pykep as pk
import bodies

def main():
    atlas = bodies.load_atlas()
    mars, lp = bodies.body('mars', atlas), pk.planet(pk.udpla.jpl_lp('mars'))
    r1, _ = mars.eph(pk.epoch('2040-06-01')); r2, _ = lp.eph(pk.epoch('2040-06-01'))
    print('DE440 vs jpl_lp Mars 2040, AU:', round(max(abs(a - b) for a, b in zip(r1, r2)) / pk.AU, 6))
    print('mars mu_self', mars.get_mu_self(), 'radius', mars.get_radius())
    seq = [bodies.body(n, atlas) for n in ('earth', 'mars', 'jupiter')]
    udp = pk.trajopt.mga_1dsm(seq=seq, t0=[pk.epoch('2080-01-01').mjd2000, pk.epoch('2082-01-01').mjd2000],
                              tof=[[100, 500], [300, 2500]], vinf=[0.5, 9], add_vinf_dep=False, add_vinf_arr=True)
    lb, ub = udp.get_bounds(); x = [(a + b) / 2 for a, b in zip(lb, ub)]
    print('flyby problem in 2080 evaluates:', udp.fitness(x)[0] > 0)
    earth, jup = bodies.body('earth', atlas), bodies.body('jupiter', atlas)
    t = time.time(); n = 0
    for k in range(2000):
        t0 = 10000 + k * 10
        ra, va = earth.eph(t0); rb, vb = jup.eph(t0 + 1000)
        pk.lambert_problem(ra, rb, 1000 * pk.DAY2SEC, pk.MU_SUN); n += 1
    print('Lambert + 2 ephemeris calls:', round((time.time() - t) / n * 1e6, 1), 'microseconds each')

if __name__ == '__main__':
    main()
