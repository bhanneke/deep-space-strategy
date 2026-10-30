# Quick end-to-end check of routes.py with a tiny optimiser budget.
import sys
sys.path.insert(0, '/work/flyby')
import routes, pykep as pk

def main():
    t0 = [pk.epoch('2031-01-01').mjd2000, pk.epoch('2033-01-01').mjd2000]
    for seq in [['earth', 'jupiter'], ['earth', 'venus', 'earth', 'jupiter'], ['earth', 'ceres']]:
        print(routes.solve(seq, t0, islands=4, gens=150, pop=20, evolves=1, seed=1), flush=True)

if __name__ == '__main__':      # pygmo starts worker processes that re-import this module
    main()
