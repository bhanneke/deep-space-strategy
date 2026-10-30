# First learned surrogate: an MLP that predicts the cheapest rendezvous cost (vinf_dep + vinf_arr) and its
# flight time from the two heliocentric states at departure. Trained on make_dataset.py samples.
import sys, time, json
import numpy as np
import torch, torch.nn as nn

K2 = 2.9591220828e-4                          # AU^3/day^2

def features(X):
    r1, v1, r2, v2 = X[:, 0:3], X[:, 3:6], X[:, 6:9], X[:, 9:12]
    ang = np.arctan2(r1[:, 1], r1[:, 0])      # rotate about the ecliptic pole so body 1 sits on the +x axis
    c, s = np.cos(-ang), np.sin(-ang)
    def rot(v):
        return np.stack([c * v[:, 0] - s * v[:, 1], s * v[:, 0] + c * v[:, 1], v[:, 2]], axis=1)
    R1, V1, R2, V2 = rot(r1), rot(v1), rot(r2), rot(v2)
    def sma(r, v):                            # vis-viva semi-major axis, AU
        return 1.0 / (2.0 / np.linalg.norm(r, axis=1) - (v ** 2).sum(1) / K2)
    a1, a2 = sma(r1, v1), sma(r2, v2)
    V1, V2 = V1 * 58.13, V2 * 58.13           # AU/day -> units near 1 (Earth's speed is about 0.0172 AU/day)
    return np.column_stack([R1[:, 0], R1[:, 2], V1, R2, V2, np.log(a1), np.log(a2)]).astype(np.float32)

def main():
    path, epochs = sys.argv[1], int(sys.argv[2])
    torch.set_num_threads(int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    d = np.load(path)
    X = d['X']
    X = X[(X[:, 14] < 60) & np.isfinite(X).all(1)]            # drop extreme transfers above 60 km/s
    F = features(X)
    Y = np.column_stack([np.log(X[:, 14]), np.log(X[:, 15])]).astype(np.float32)
    rng = np.random.default_rng(0); idx = rng.permutation(len(F))
    ntr, nva = int(0.8 * len(F)), int(0.1 * len(F))
    tr, va, te = idx[:ntr], idx[ntr:ntr + nva], idx[ntr + nva:]
    mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
    Fn = (F - mu) / sd
    net = nn.Sequential(nn.Linear(F.shape[1], 256), nn.SiLU(), nn.Linear(256, 256), nn.SiLU(),
                        nn.Linear(256, 256), nn.SiLU(), nn.Linear(256, 256), nn.SiLU(), nn.Linear(256, 2))
    opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    Xt, Yt = torch.tensor(Fn[tr]), torch.tensor(Y[tr])
    Xv, Yv = torch.tensor(Fn[va]), torch.tensor(Y[va])
    t0 = time.time()
    for ep in range(epochs):
        net.train(); perm = torch.randperm(len(Xt))
        for k in range(0, len(Xt), 4096):
            b = perm[k:k + 4096]
            loss = nn.functional.smooth_l1_loss(net(Xt[b]), Yt[b], beta=0.05)
            opt.zero_grad(); loss.backward(); opt.step()
        sched.step()
        if ep % 5 == 4 or ep == epochs - 1:
            net.eval()
            with torch.no_grad():
                vl = nn.functional.l1_loss(net(Xv)[:, 0], Yv[:, 0]).item()
            print(f'epoch {ep + 1:3d}  val mean |log error| of cost {vl:.4f}  ({time.time() - t0:.0f} s)', flush=True)
    net.eval()
    with torch.no_grad():
        P = net(torch.tensor(Fn[te])).numpy()
    true, pred = X[te, 14], np.exp(P[:, 0])
    err = np.abs(pred - true)
    from scipy.stats import spearmanr
    res = {'test_samples': int(len(te)), 'mae_kms': float(err.mean()), 'median_abs_err_kms': float(np.median(err)),
           'p95_abs_err_kms': float(np.percentile(err, 95)), 'share_within_0_5_kms': float((err < 0.5).mean()),
           'median_rel_err': float(np.median(err / true)), 'spearman_rank': float(spearmanr(pred, true).correlation),
           'tof_median_rel_err': float(np.median(np.abs(np.exp(P[:, 1]) - X[te, 15]) / X[te, 15]))}
    print(json.dumps(res, indent=1))
    torch.save({'state': net.state_dict(), 'mu': mu, 'sd': sd}, '/work/surrogate/mlp_v1.pt')
    json.dump(res, open('/work/surrogate/mlp_v1_metrics.json', 'w'), indent=1)

if __name__ == '__main__':
    main()
