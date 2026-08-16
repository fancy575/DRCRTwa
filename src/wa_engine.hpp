#ifndef DRCRTWA_WA_ENGINE_HPP
#define DRCRTWA_WA_ENGINE_HPP

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace drcrtwa {

inline bool finite_number(double x) { return R_finite(x); }

inline double clamp_exp_arg(double x) {
  if (x > 25.0) return 25.0;
  if (x < -25.0) return -25.0;
  return x;
}

inline double safe_exp(double x) { return std::exp(clamp_exp_arg(x)); }
inline int idx2(int i, int g, int G) { return i * G + g; }
inline int idx3(int k, int i, int g, int N, int G) {
  return (k * N + i) * G + g;
}

inline double dot_row(const std::vector<double>& X, int i, int p,
                      const std::vector<double>& beta) {
  double out = 0.0;
  const int off = i * p;
  for (int k = 0; k < p; ++k) out += X[off + k] * beta[k];
  return out;
}

inline bool solve_linear(std::vector<double> A,
                         std::vector<double> b,
                         int p,
                         std::vector<double>& x,
                         double ridge = 0.0) {
  x.assign(p, 0.0);
  if (p == 0) return true;
  for (int i = 0; i < p; ++i) A[i * p + i] += ridge;
  for (int col = 0; col < p; ++col) {
    int pivot = col;
    double best = std::fabs(A[col * p + col]);
    for (int r = col + 1; r < p; ++r) {
      const double v = std::fabs(A[r * p + col]);
      if (v > best) {
        best = v;
        pivot = r;
      }
    }
    if (!finite_number(best) || best < 1e-13) return false;
    if (pivot != col) {
      for (int c = col; c < p; ++c)
        std::swap(A[col * p + c], A[pivot * p + c]);
      std::swap(b[col], b[pivot]);
    }
    const double diag = A[col * p + col];
    for (int c = col; c < p; ++c) A[col * p + c] /= diag;
    b[col] /= diag;
    for (int r = 0; r < p; ++r) {
      if (r == col) continue;
      const double fac = A[r * p + col];
      if (fac == 0.0) continue;
      for (int c = col; c < p; ++c)
        A[r * p + c] -= fac * A[col * p + c];
      b[r] -= fac * b[col];
    }
  }
  x = b;
  return true;
}

inline bool invert_matrix(const std::vector<double>& A,
                          int p,
                          std::vector<double>& inv,
                          double ridge = 0.0) {
  inv.assign(p * p, 0.0);
  if (p == 0) return true;
  for (int j = 0; j < p; ++j) {
    std::vector<double> e(p, 0.0), sol;
    e[j] = 1.0;
    if (!solve_linear(A, e, p, sol, ridge)) return false;
    for (int i = 0; i < p; ++i) inv[i * p + j] = sol[i];
  }
  return true;
}

inline int grid_bin(double t, double tau, int G) {
  if (t <= 0.0) return 0;
  int g = static_cast<int>(std::floor(t * static_cast<double>(G) / tau));
  if (g < 0) g = 0;
  if (g >= G) g = G - 1;
  return g;
}

struct DataSet {
  int n_subjects;
  int n_units;
  int G;
  int n_event_types;
  double tau;
  std::vector<int> unit;
  std::vector<int> A;
  std::vector<int> cluster_size;
  std::vector<double> prob0;
  std::vector<double> prob1;
  std::vector<double> death_time;
  std::vector<double> censor_time;
  std::vector<double> observed_time;
  std::vector<unsigned char> death_event;
  std::vector<unsigned char> censor_event;
  // recurrent_time[k][i] contains exact observed type-k event times for subject i.
  std::vector< std::vector< std::vector<double> > > recurrent_time;
  std::vector<unsigned char> Y;
  std::vector<unsigned char> Ydag;
  std::vector<unsigned char> dD;
  std::vector<unsigned char> dC;
  // Type-specific recurrent increments, indexed by (k, i, g).
  std::vector<double> dR;
  std::vector<double> dt;
};

// kind = -2 for censoring, -1 for terminal event, and 0,...,K-1 for
// recurrent event types.
struct EventGroup {
  double time;
  std::vector<int> subject;
};

inline std::vector<EventGroup> collect_event_groups(const DataSet& dat,
                                                     int arm,
                                                     int kind) {
  std::vector< std::pair<double, int> > raw;
  raw.reserve(dat.n_subjects);
  for (int i = 0; i < dat.n_subjects; ++i) {
    if (dat.A[i] != arm) continue;
    if (kind == -1) {
      if (dat.death_event[i])
        raw.push_back(std::make_pair(dat.death_time[i], i));
    } else if (kind == -2) {
      if (dat.censor_event[i])
        raw.push_back(std::make_pair(dat.censor_time[i], i));
    } else {
      const std::vector<double>& ev = dat.recurrent_time[kind][i];
      for (std::size_t q = 0; q < ev.size(); ++q) {
        if (ev[q] <= dat.observed_time[i] + 1e-12 &&
            ev[q] <= dat.tau + 1e-12) {
          raw.push_back(std::make_pair(ev[q], i));
        }
      }
    }
  }
  std::sort(raw.begin(), raw.end(),
            [](const std::pair<double, int>& a,
               const std::pair<double, int>& b) {
              return (a.first < b.first) ||
                (a.first == b.first && a.second < b.second);
            });
  std::vector<EventGroup> groups;
  for (std::size_t q = 0; q < raw.size(); ++q) {
    if (groups.empty() ||
        std::fabs(raw[q].first - groups.back().time) > 1e-11) {
      EventGroup eg;
      eg.time = raw[q].first;
      eg.subject.push_back(raw[q].second);
      groups.push_back(eg);
    } else {
      groups.back().subject.push_back(raw[q].second);
    }
  }
  return groups;
}

struct CoxEval {
  std::vector<double> score;
  std::vector<double> info;
  double loglik;
  bool valid;
};

struct CoxFit {
  int p;
  int arm;
  int kind;
  int G;
  std::vector<double> X;
  std::vector<double> beta;
  std::vector<double> base_inc;
  std::vector<double> if_beta;
  std::vector<double> if_base;
  bool converged;
  bool valid;
  int iterations;
  int n_events;
};

inline CoxEval evaluate_cox(const DataSet& dat,
                            const std::vector<double>& X,
                            int p,
                            const std::vector<double>& beta,
                            const std::vector<EventGroup>& groups,
                            const std::vector<int>& risk_subject,
                            const std::vector<double>& risk_stop) {
  (void) dat;
  CoxEval out;
  out.score.assign(p, 0.0);
  out.info.assign(p * p, 0.0);
  out.loglik = 0.0;
  out.valid = true;
  const int nr = static_cast<int>(risk_subject.size());
  if (nr == 0 || groups.empty()) {
    out.valid = false;
    return out;
  }

  std::vector<double> rev0(nr + 1, 0.0);
  std::vector<double> rev1((nr + 1) * p, 0.0);
  std::vector<double> rev2((nr + 1) * p * p, 0.0);
  for (int pos = nr - 1; pos >= 0; --pos) {
    const int i = risk_subject[pos];
    const double rr = safe_exp(dot_row(X, i, p, beta));
    rev0[pos] = rev0[pos + 1] + rr;
    for (int k = 0; k < p; ++k) {
      const double xk = X[i * p + k];
      rev1[pos * p + k] = rev1[(pos + 1) * p + k] + rr * xk;
      for (int l = 0; l < p; ++l) {
        rev2[(pos * p + k) * p + l] =
          rev2[((pos + 1) * p + k) * p + l] +
          rr * xk * X[i * p + l];
      }
    }
  }

  for (std::size_t j = 0; j < groups.size(); ++j) {
    const double t = groups[j].time;
    const int pos = static_cast<int>(
      std::lower_bound(risk_stop.begin(), risk_stop.end(), t - 1e-12) -
        risk_stop.begin());
    if (pos >= nr) {
      out.valid = false;
      return out;
    }
    const double s0 = rev0[pos];
    if (!finite_number(s0) || s0 <= 1e-12) {
      out.valid = false;
      return out;
    }
    const double d = static_cast<double>(groups[j].subject.size());
    std::vector<double> bar(p, 0.0);
    for (int k = 0; k < p; ++k) bar[k] = rev1[pos * p + k] / s0;

    out.loglik -= d * std::log(s0);
    for (std::size_t e = 0; e < groups[j].subject.size(); ++e) {
      const int i = groups[j].subject[e];
      out.loglik += dot_row(X, i, p, beta);
      for (int k = 0; k < p; ++k)
        out.score[k] += X[i * p + k] - bar[k];
    }
    for (int k = 0; k < p; ++k) {
      for (int l = 0; l < p; ++l) {
        const double v = rev2[(pos * p + k) * p + l] / s0 -
          bar[k] * bar[l];
        out.info[k * p + l] += d * v;
      }
    }
  }
  return out;
}

inline CoxFit fit_cox_lwyy(const DataSet& dat,
                           const std::vector<double>& X,
                           int p,
                           int arm,
                           int kind,
                           int max_iter = 40,
                           double tol = 1e-8) {
  CoxFit fit;
  fit.p = p;
  fit.arm = arm;
  fit.kind = kind;
  fit.G = dat.G;
  fit.X = X;
  fit.beta.assign(p, 0.0);
  fit.base_inc.assign(dat.G, 0.0);
  fit.if_beta.assign(dat.n_units * p, 0.0);
  fit.if_base.assign(dat.n_units * dat.G, 0.0);
  fit.converged = false;
  fit.valid = true;
  fit.iterations = 0;

  std::vector<EventGroup> groups = collect_event_groups(dat, arm, kind);
  fit.n_events = 0;
  for (std::size_t j = 0; j < groups.size(); ++j)
    fit.n_events += static_cast<int>(groups[j].subject.size());

  // A zero-event nuisance process has a zero baseline increment. This is a
  // valid boundary fit and is preferable to failing an otherwise estimable
  // while-alive analysis.
  if (fit.n_events == 0) {
    fit.converged = true;
    return fit;
  }

  std::vector<int> risk_subject;
  std::vector<double> risk_stop;
  for (int i = 0; i < dat.n_subjects; ++i) {
    if (dat.A[i] == arm) {
      risk_subject.push_back(i);
      risk_stop.push_back(dat.observed_time[i]);
    }
  }
  std::vector<int> ord(risk_subject.size());
  for (std::size_t q = 0; q < ord.size(); ++q)
    ord[q] = static_cast<int>(q);
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return risk_stop[a] < risk_stop[b]; });
  std::vector<int> rs2(ord.size());
  std::vector<double> st2(ord.size());
  for (std::size_t q = 0; q < ord.size(); ++q) {
    rs2[q] = risk_subject[ord[q]];
    st2[q] = risk_stop[ord[q]];
  }
  risk_subject.swap(rs2);
  risk_stop.swap(st2);

  CoxEval cur = evaluate_cox(dat, fit.X, fit.p, fit.beta, groups,
                             risk_subject, risk_stop);
  if (!cur.valid) {
    fit.valid = false;
    return fit;
  }

  if (p == 0) {
    fit.converged = true;
  } else {
    for (int it = 0; it < max_iter; ++it) {
      fit.iterations = it + 1;
      std::vector<double> step;
      bool ok = false;
      double ridge = 1e-10;
      for (int attempt = 0; attempt < 10 && !ok; ++attempt) {
        ok = solve_linear(cur.info, cur.score, fit.p, step, ridge);
        ridge *= 10.0;
      }
      if (!ok) break;
      double maxstep = 0.0;
      for (int k = 0; k < fit.p; ++k)
        maxstep = std::max(maxstep, std::fabs(step[k]));
      double fac = 1.0;
      bool accepted = false;
      std::vector<double> cand(fit.p, 0.0);
      CoxEval ce;
      for (int half = 0; half < 24; ++half) {
        for (int k = 0; k < fit.p; ++k) {
          cand[k] = std::max(-8.0,
                             std::min(8.0, fit.beta[k] + fac * step[k]));
        }
        ce = evaluate_cox(dat, fit.X, fit.p, cand, groups,
                          risk_subject, risk_stop);
        if (ce.valid && finite_number(ce.loglik) &&
            ce.loglik >= cur.loglik - 1e-8) {
          accepted = true;
          break;
        }
        fac *= 0.5;
      }
      if (!accepted) break;
      fit.beta = cand;
      cur = ce;
      if (fac * maxstep < tol) {
        fit.converged = true;
        break;
      }
    }
    cur = evaluate_cox(dat, fit.X, fit.p, fit.beta, groups,
                       risk_subject, risk_stop);
    if (!cur.valid) {
      fit.valid = false;
      return fit;
    }
    if (!fit.converged) {
      double maxscore = 0.0;
      for (int k = 0; k < fit.p; ++k)
        maxscore = std::max(maxscore, std::fabs(cur.score[k]));
      if (maxscore < 1e-5 * std::max(1, fit.n_events))
        fit.converged = true;
    }
  }

  std::vector<double> inv_info;
  if (p > 0) {
    bool inv_ok = false;
    double ridge = 1e-10;
    for (int attempt = 0; attempt < 12 && !inv_ok; ++attempt) {
      inv_ok = invert_matrix(cur.info, fit.p, inv_info, ridge);
      ridge *= 10.0;
    }
    if (!inv_ok) {
      fit.valid = false;
      return fit;
    }
  }

  const int U = dat.n_units;
  std::vector<double> unit_score(U * fit.p, 0.0);
  std::vector<double> rr(dat.n_subjects, 0.0);
  for (int i = 0; i < dat.n_subjects; ++i) {
    if (dat.A[i] == arm)
      rr[i] = safe_exp(dot_row(fit.X, i, fit.p, fit.beta));
  }

  struct JumpInfo {
    double time;
    double s0;
    double dlam;
    std::vector<double> bar;
    std::vector<int> event_subject;
  };
  std::vector<JumpInfo> jump;
  jump.reserve(groups.size());

  const int nr = static_cast<int>(risk_subject.size());
  std::vector<double> rev0(nr + 1, 0.0);
  std::vector<double> rev1((nr + 1) * fit.p, 0.0);
  for (int pos = nr - 1; pos >= 0; --pos) {
    const int i = risk_subject[pos];
    rev0[pos] = rev0[pos + 1] + rr[i];
    for (int k = 0; k < fit.p; ++k) {
      rev1[pos * fit.p + k] =
        rev1[(pos + 1) * fit.p + k] + rr[i] * fit.X[i * fit.p + k];
    }
  }
  for (std::size_t j = 0; j < groups.size(); ++j) {
    const int pos = static_cast<int>(
      std::lower_bound(risk_stop.begin(), risk_stop.end(),
                       groups[j].time - 1e-12) - risk_stop.begin());
    if (pos >= nr || rev0[pos] <= 1e-12) {
      fit.valid = false;
      return fit;
    }
    JumpInfo ji;
    ji.time = groups[j].time;
    ji.s0 = rev0[pos];
    ji.dlam = static_cast<double>(groups[j].subject.size()) / ji.s0;
    ji.bar.assign(fit.p, 0.0);
    ji.event_subject = groups[j].subject;
    for (int k = 0; k < fit.p; ++k)
      ji.bar[k] = rev1[pos * fit.p + k] / ji.s0;
    for (std::size_t e = 0; e < groups[j].subject.size(); ++e) {
      const int i = groups[j].subject[e];
      const int u = dat.unit[i];
      for (int k = 0; k < fit.p; ++k)
        unit_score[u * fit.p + k] += fit.X[i * fit.p + k] - ji.bar[k];
    }
    fit.base_inc[grid_bin(ji.time, dat.tau, dat.G)] += ji.dlam;
    jump.push_back(ji);
  }

  // Add the risk-set compensator to the robust unit score. This is the
  // martingale score residual for Cox fits and the mean-zero LWYY rate-score
  // residual for recurrent-event fits.
  std::vector<double> jump_time(jump.size(), 0.0);
  std::vector<double> cum_lam(jump.size() + 1, 0.0);
  std::vector<double> cum_bar_lam((jump.size() + 1) * fit.p, 0.0);
  for (std::size_t j = 0; j < jump.size(); ++j) {
    jump_time[j] = jump[j].time;
    cum_lam[j + 1] = cum_lam[j] + jump[j].dlam;
    for (int k = 0; k < fit.p; ++k) {
      cum_bar_lam[(j + 1) * fit.p + k] =
        cum_bar_lam[j * fit.p + k] + jump[j].dlam * jump[j].bar[k];
    }
  }
  for (int i = 0; i < dat.n_subjects; ++i) {
    if (dat.A[i] != arm) continue;
    const std::size_t nj = static_cast<std::size_t>(
      std::upper_bound(jump_time.begin(), jump_time.end(),
                       dat.observed_time[i] + 1e-12) - jump_time.begin());
    if (nj == 0) continue;
    const int u = dat.unit[i];
    const double lam = cum_lam[nj];
    for (int k = 0; k < fit.p; ++k) {
      const double comp = rr[i] *
        (fit.X[i * fit.p + k] * lam -
         cum_bar_lam[nj * fit.p + k]);
      unit_score[u * fit.p + k] -= comp;
    }
  }

  for (int u = 0; u < U; ++u) {
    for (int k = 0; k < fit.p; ++k) {
      double val = 0.0;
      for (int l = 0; l < fit.p; ++l)
        val += inv_info[k * fit.p + l] * unit_score[u * fit.p + l];
      fit.if_beta[u * fit.p + k] = static_cast<double>(U) * val;
    }
  }

  // Breslow baseline-increment influence contribution, including the
  // regression-coefficient estimation effect.
  std::vector<double> risk_u(U, 0.0), event_u(U, 0.0);
  for (std::size_t j = 0; j < jump.size(); ++j) {
    std::fill(risk_u.begin(), risk_u.end(), 0.0);
    std::fill(event_u.begin(), event_u.end(), 0.0);
    for (int i = 0; i < dat.n_subjects; ++i) {
      if (dat.A[i] == arm &&
          dat.observed_time[i] + 1e-12 >= jump[j].time) {
        risk_u[dat.unit[i]] += rr[i];
      }
    }
    for (std::size_t e = 0; e < jump[j].event_subject.size(); ++e)
      event_u[dat.unit[jump[j].event_subject[e]]] += 1.0;
    const int g = grid_bin(jump[j].time, dat.tau, dat.G);
    for (int u = 0; u < U; ++u) {
      double beta_part = 0.0;
      for (int k = 0; k < fit.p; ++k)
        beta_part += jump[j].bar[k] * fit.if_beta[u * fit.p + k];
      const double residual_part = static_cast<double>(U) *
        (event_u[u] - jump[j].dlam * risk_u[u]) / jump[j].s0;
      fit.if_base[idx2(u, g, dat.G)] +=
        residual_part - jump[j].dlam * beta_part;
    }
  }
  return fit;
}

struct MethodResult {
  double nu;
  double mu;
  double rate;
  double se;
  bool valid;
  std::vector<double> if_rate;
};

struct ArmMethodResults {
  MethodResult method[3];
};

enum MethodIndex { DR = 0, IPCW = 1, OR = 2 };

struct FunctionalAdjoint {
  bool valid;
  double nu;
  double mu;
  double rate;
  std::vector<double> Bmean;
  std::vector<double> ADmean;
  std::vector<double> ARmean;
  std::vector<double> incD;
  std::vector<double> incR;
  std::vector<double> S;
  std::vector<double> adjB;
  std::vector<double> adjAD;
  std::vector<double> adjAR;
};

inline FunctionalAdjoint evaluate_functional(
    const DataSet& dat,
    const std::vector<double>& Bunit,
    const std::vector<double>& ADunit,
    const std::vector<double>& ARunit) {
  FunctionalAdjoint f;
  const int U = dat.n_units;
  const int G = dat.G;
  f.valid = true;
  f.nu = 0.0;
  f.mu = 0.0;
  f.rate = NA_REAL;
  f.Bmean.assign(G, 0.0);
  f.ADmean.assign(G, 0.0);
  f.ARmean.assign(G, 0.0);
  f.incD.assign(G, 0.0);
  f.incR.assign(G, 0.0);
  f.S.assign(G + 1, 1.0);
  for (int g = 0; g < G; ++g) {
    for (int u = 0; u < U; ++u) {
      const int ug = idx2(u, g, G);
      f.Bmean[g] += Bunit[ug];
      f.ADmean[g] += ADunit[ug];
      f.ARmean[g] += ARunit[ug];
    }
    f.Bmean[g] /= static_cast<double>(U);
    f.ADmean[g] /= static_cast<double>(U);
    f.ARmean[g] /= static_cast<double>(U);
    if (!finite_number(f.Bmean[g]) || f.Bmean[g] <= 1e-12) {
      f.valid = false;
      return f;
    }
    f.incD[g] = f.ADmean[g] / f.Bmean[g];
    f.incR[g] = f.ARmean[g] / f.Bmean[g];
    if (!finite_number(f.incD[g]) || !finite_number(f.incR[g]) ||
        f.incD[g] >= 0.999999) {
      f.valid = false;
      return f;
    }
    f.nu += f.S[g] * dat.dt[g];
    f.mu += f.S[g] * f.incR[g];
    f.S[g + 1] = f.S[g] * (1.0 - f.incD[g]);
  }
  if (!finite_number(f.nu) || f.nu <= 1e-12 || !finite_number(f.mu)) {
    f.valid = false;
    return f;
  }
  f.rate = f.mu / f.nu;

  f.adjB.assign(G, 0.0);
  f.adjAD.assign(G, 0.0);
  f.adjAR.assign(G, 0.0);
  const double a_mu = 1.0 / f.nu;
  const double a_nu = -f.mu / (f.nu * f.nu);
  double adj_S_next = 0.0;
  for (int g = G - 1; g >= 0; --g) {
    const double adj_incR = a_mu * f.S[g];
    const double adj_incD = -adj_S_next * f.S[g];
    const double adj_S = a_nu * dat.dt[g] + a_mu * f.incR[g] +
      adj_S_next * (1.0 - f.incD[g]);
    adj_S_next = adj_S;
    f.adjAD[g] = adj_incD / f.Bmean[g];
    f.adjAR[g] = adj_incR / f.Bmean[g];
    f.adjB[g] =
      -(adj_incD * f.incD[g] + adj_incR * f.incR[g]) / f.Bmean[g];
  }
  return f;
}

inline double model_if_dot(const CoxFit& fit,
                           const std::vector<double>& grad_beta,
                           const std::vector<double>& grad_base,
                           int u) {
  double out = 0.0;
  for (int k = 0; k < fit.p; ++k)
    out += grad_beta[k] * fit.if_beta[u * fit.p + k];
  for (int g = 0; g < fit.G; ++g)
    out += grad_base[g] * fit.if_base[idx2(u, g, fit.G)];
  return out;
}

inline ArmMethodResults compute_arm_methods(
    const DataSet& dat,
    int arm,
    const CoxFit& fitC,
    const CoxFit& fitD,
    const std::vector<CoxFit>& fitR,
    const std::vector<double>& event_weights,
    const std::vector<double>& target_weight) {
  const int N = dat.n_subjects;
  const int U = dat.n_units;
  const int G = dat.G;
  const int K = dat.n_event_types;
  const int NG = N * G;
  const int UG = U * G;
  ArmMethodResults out;
  for (int m = 0; m < 3; ++m) {
    out.method[m].valid = true;
    out.method[m].nu = NA_REAL;
    out.method[m].mu = NA_REAL;
    out.method[m].rate = NA_REAL;
    out.method[m].se = NA_REAL;
    out.method[m].if_rate.assign(U, 0.0);
  }
  bool nuisance_valid = fitC.valid && fitD.valid &&
    static_cast<int>(fitR.size()) == K &&
    static_cast<int>(event_weights.size()) == K;
  for (int k = 0; k < K; ++k) nuisance_valid = nuisance_valid && fitR[k].valid;
  if (!nuisance_valid) {
    for (int m = 0; m < 3; ++m) out.method[m].valid = false;
    return out;
  }

  std::vector<double> cumC(NG, 0.0), cumD(NG, 0.0),
    Uaug(NG, 1.0), dLC(NG, 0.0), dLD(NG, 0.0), mC(NG, 0.0);
  std::vector< std::vector<double> > dLR(
    K, std::vector<double>(NG, 0.0));
  std::vector< std::vector<double> > B(
    3, std::vector<double>(UG, 0.0));
  std::vector< std::vector<double> > AD(
    3, std::vector<double>(UG, 0.0));
  std::vector< std::vector<double> > AR(
    3, std::vector<double>(UG, 0.0));

  for (int i = 0; i < N; ++i) {
    const double rrC = safe_exp(dot_row(fitC.X, i, fitC.p, fitC.beta));
    const double rrD = safe_exp(dot_row(fitD.X, i, fitD.p, fitD.beta));
    std::vector<double> rrR(K, 1.0);
    for (int k = 0; k < K; ++k)
      rrR[k] = safe_exp(dot_row(fitR[k].X, i, fitR[k].p, fitR[k].beta));
    double cC = 0.0, cD = 0.0, ua = 1.0;
    const double pa = (arm == 1) ? dat.prob1[i] : dat.prob0[i];
    const double xi = (dat.A[i] == arm) ? 1.0 / pa : 0.0;
    const int u = dat.unit[i];
    const double tw = target_weight[i];
    for (int g = 0; g < G; ++g) {
      const int id = idx2(i, g, G);
      const int ug = idx2(u, g, G);
      cumC[id] = cC;
      cumD[id] = cD;
      Uaug[id] = ua;
      dLC[id] = rrC * fitC.base_inc[g];
      dLD[id] = rrD * fitD.base_inc[g];
      double model_weighted_recur = 0.0;
      double observed_weighted_recur = 0.0;
      for (int k = 0; k < K; ++k) {
        dLR[k][id] = rrR[k] * fitR[k].base_inc[g];
        model_weighted_recur += event_weights[k] * dLR[k][id];
        observed_weighted_recur += event_weights[k] *
          dat.dR[idx3(k, i, g, N, G)];
      }
      const double Ksurv = std::exp(-cC);
      const double Hsurv = std::exp(-cD);
      const double inv = xi / std::max(Ksurv, 1e-12);
      const double aug = 1.0 - xi * ua;
      const double qD = Hsurv * dLD[id];
      const double qR = Hsurv * model_weighted_recur;
      const double y = static_cast<double>(dat.Y[id]);
      const double dd = static_cast<double>(dat.dD[id]);
      const double dc = static_cast<double>(dat.dC[id]);
      const double dr = observed_weighted_recur;

      B[DR][ug] += tw * (inv * y + aug * Hsurv);
      AD[DR][ug] += tw * (inv * dd + aug * qD);
      AR[DR][ug] += tw * (inv * dr + aug * qR);
      B[IPCW][ug] += tw * inv * y;
      AD[IPCW][ug] += tw * inv * dd;
      AR[IPCW][ug] += tw * inv * dr;
      B[OR][ug] += tw * Hsurv;
      AD[OR][ug] += tw * qD;
      AR[OR][ug] += tw * qR;

      mC[id] = dc - static_cast<double>(dat.Ydag[id]) * dLC[id];
      ua -= mC[id] * std::exp(cC + cD);
      cC += dLC[id];
      cD += dLD[id];
    }
  }

  for (int m = 0; m < 3; ++m) {
    FunctionalAdjoint f = evaluate_functional(dat, B[m], AD[m], AR[m]);
    MethodResult mr;
    mr.valid = f.valid;
    mr.nu = f.nu;
    mr.mu = f.mu;
    mr.rate = f.rate;
    mr.se = NA_REAL;
    mr.if_rate.assign(U, 0.0);
    if (!f.valid) {
      out.method[m] = mr;
      continue;
    }

    // Empirical contribution with nuisance functions held fixed.
    for (int u = 0; u < U; ++u) {
      double val = 0.0;
      for (int g = 0; g < G; ++g) {
        const int ug = idx2(u, g, G);
        val += f.adjB[g] * (B[m][ug] - f.Bmean[g]) +
          f.adjAD[g] * (AD[m][ug] - f.ADmean[g]) +
          f.adjAR[g] * (AR[m][ug] - f.ARmean[g]);
      }
      mr.if_rate[u] = val;
    }

    std::vector<double> gBC(fitC.p, 0.0), gBD(fitD.p, 0.0);
    std::vector<double> gLC(G, 0.0), gLD(G, 0.0);
    std::vector< std::vector<double> > gBR(K), gLR(K);
    for (int k = 0; k < K; ++k) {
      gBR[k].assign(fitR[k].p, 0.0);
      gLR[k].assign(G, 0.0);
    }

    // Reverse-mode derivative of the scalar arm rate with respect to all
    // fitted Cox/LWYY coefficients and Breslow increments.
    for (int i = 0; i < N; ++i) {
      std::vector<double> aCumC(G, 0.0), aCumD(G, 0.0),
        aU(G, 0.0), aDLC(G, 0.0), aDLD(G, 0.0);
      std::vector< std::vector<double> > aDLR(
        K, std::vector<double>(G, 0.0));
      const double pa = (arm == 1) ? dat.prob1[i] : dat.prob0[i];
      const double xi = (dat.A[i] == arm) ? 1.0 / pa : 0.0;
      const double coeff = target_weight[i] / static_cast<double>(U);

      for (int g = 0; g < G; ++g) {
        const int id = idx2(i, g, G);
        const double Hsurv = std::exp(-cumD[id]);
        const double inv = xi * std::exp(cumC[id]);
        const double aug = 1.0 - xi * Uaug[id];
        double model_weighted_recur = 0.0;
        double observed_weighted_recur = 0.0;
        for (int k = 0; k < K; ++k) {
          model_weighted_recur += event_weights[k] * dLR[k][id];
          observed_weighted_recur += event_weights[k] *
            dat.dR[idx3(k, i, g, N, G)];
        }
        const double qD = Hsurv * dLD[id];
        const double qR = Hsurv * model_weighted_recur;
        const double aB = coeff * f.adjB[g];
        const double aAD = coeff * f.adjAD[g];
        const double aAR = coeff * f.adjAR[g];
        const double y = static_cast<double>(dat.Y[id]);
        const double dd = static_cast<double>(dat.dD[id]);
        const double dr = observed_weighted_recur;

        if (m == DR) {
          const double aInv = aB * y + aAD * dd + aAR * dr;
          const double aAug = aB * Hsurv + aAD * qD + aAR * qR;
          aU[g] += -xi * aAug;
          double aH = aB * aug;
          const double aqD = aAD * aug;
          const double aqR = aAR * aug;
          aH += aqD * dLD[id] + aqR * model_weighted_recur;
          aDLD[g] += aqD * Hsurv;
          for (int k = 0; k < K; ++k)
            aDLR[k][g] += aqR * Hsurv * event_weights[k];
          aCumC[g] += aInv * inv;
          aCumD[g] += -aH * Hsurv;
        } else if (m == IPCW) {
          const double aInv = aB * y + aAD * dd + aAR * dr;
          aCumC[g] += aInv * inv;
        } else {
          const double aH = aB + aAD * dLD[id] +
            aAR * model_weighted_recur;
          aDLD[g] += aAD * Hsurv;
          for (int k = 0; k < K; ++k)
            aDLR[k][g] += aAR * Hsurv * event_weights[k];
          aCumD[g] += -aH * Hsurv;
        }
      }

      if (m == DR) {
        double adj_next = 0.0;
        for (int g = G - 1; g >= 0; --g) {
          const int id = idx2(i, g, G);
          if (adj_next != 0.0) {
            const double invKH = std::exp(cumC[id] + cumD[id]);
            const double aM = -adj_next * invKH;
            const double aInvKH = -adj_next * mC[id];
            aCumC[g] += aInvKH * invKH;
            aCumD[g] += aInvKH * invKH;
            aDLC[g] += -static_cast<double>(dat.Ydag[id]) * aM;
          }
          adj_next = aU[g] + adj_next;
        }
      }

      double carryC = 0.0, carryD = 0.0;
      for (int g = G - 1; g >= 0; --g) {
        aDLC[g] += carryC;
        aDLD[g] += carryD;
        carryC += aCumC[g];
        carryD += aCumD[g];
      }

      const int offC = i * fitC.p;
      const int offD = i * fitD.p;
      const double rrC = safe_exp(dot_row(fitC.X, i, fitC.p, fitC.beta));
      const double rrD = safe_exp(dot_row(fitD.X, i, fitD.p, fitD.beta));
      std::vector<double> rrR(K, 1.0);
      for (int k = 0; k < K; ++k)
        rrR[k] = safe_exp(dot_row(fitR[k].X, i, fitR[k].p,
                                  fitR[k].beta));

      for (int g = 0; g < G; ++g) {
        gLC[g] += aDLC[g] * rrC;
        gLD[g] += aDLD[g] * rrD;
        for (int k = 0; k < fitC.p; ++k) {
          gBC[k] += aDLC[g] * dLC[idx2(i, g, G)] *
            fitC.X[offC + k];
        }
        for (int k = 0; k < fitD.p; ++k) {
          gBD[k] += aDLD[g] * dLD[idx2(i, g, G)] *
            fitD.X[offD + k];
        }
        for (int e = 0; e < K; ++e) {
          gLR[e][g] += aDLR[e][g] * rrR[e];
          const int offR = i * fitR[e].p;
          for (int k = 0; k < fitR[e].p; ++k) {
            gBR[e][k] += aDLR[e][g] * dLR[e][idx2(i, g, G)] *
              fitR[e].X[offR + k];
          }
        }
      }
    }

    for (int u = 0; u < U; ++u) {
      mr.if_rate[u] += model_if_dot(fitC, gBC, gLC, u) +
        model_if_dot(fitD, gBD, gLD, u);
      for (int k = 0; k < K; ++k)
        mr.if_rate[u] += model_if_dot(fitR[k], gBR[k], gLR[k], u);
    }

    const double mn = std::accumulate(mr.if_rate.begin(),
                                      mr.if_rate.end(), 0.0) / U;
    double ss = 0.0;
    for (int u = 0; u < U; ++u)
      ss += (mr.if_rate[u] - mn) * (mr.if_rate[u] - mn);
    mr.se = (U > 1) ?
      std::sqrt(ss / (static_cast<double>(U) * (U - 1.0))) : NA_REAL;
    out.method[m] = mr;
  }
  return out;
}

inline Rcpp::List fit_to_list(const CoxFit& fit) {
  return Rcpp::List::create(
    Rcpp::_["beta"] = Rcpp::wrap(fit.beta),
    Rcpp::_["base_increment"] = Rcpp::wrap(fit.base_inc),
    Rcpp::_["converged"] = fit.converged,
    Rcpp::_["valid"] = fit.valid,
    Rcpp::_["iterations"] = fit.iterations,
    Rcpp::_["n_events"] = fit.n_events
  );
}

}  // namespace drcrtwa

#endif
