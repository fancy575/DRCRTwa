// Compiled semiparametric engine for DRCRTwa.
//
// This implementation is adapted from the C++ simulation and applied-data
// engines accompanying the doubly robust while-alive methodology. It retains
// the Cox/LWYY nuisance fitting, censoring-martingale augmentation, and full
// independent-unit influence-function calculation, generalized here to an
// arbitrary number of recurrent event types and prespecified event weights.

#include <Rcpp.h>
#include "wa_engine.hpp"

using namespace Rcpp;
using namespace drcrtwa;

namespace {

std::vector<double> matrix_to_row_major(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  std::vector<double> out(n * p, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      out[i * p + j] = X(i, j);
  return out;
}

DataSet build_dataset(
    const IntegerVector& arm,
    const IntegerVector& unit,
    const IntegerVector& cluster_size,
    const NumericVector& prob0,
    const NumericVector& prob1,
    const NumericVector& final_time,
    const IntegerVector& final_status,
    const std::vector< std::vector< std::vector<double> > >& all_recurrent,
    int terminal_code,
    double tau,
    int n_grid) {
  const int N = arm.size();
  const int K = static_cast<int>(all_recurrent.size());
  DataSet dat;
  dat.n_subjects = N;
  dat.n_units = 0;
  dat.G = n_grid;
  dat.n_event_types = K;
  dat.tau = tau;
  dat.unit.resize(N);
  dat.A.resize(N);
  dat.cluster_size.resize(N);
  dat.prob0.resize(N);
  dat.prob1.resize(N);
  dat.death_time.assign(N, std::numeric_limits<double>::infinity());
  dat.censor_time.assign(N, std::numeric_limits<double>::infinity());
  dat.observed_time.assign(N, tau);
  dat.death_event.assign(N, 0);
  dat.censor_event.assign(N, 0);
  dat.recurrent_time.assign(
    K, std::vector< std::vector<double> >(N));
  dat.Y.assign(N * n_grid, 0);
  dat.Ydag.assign(N * n_grid, 0);
  dat.dD.assign(N * n_grid, 0);
  dat.dC.assign(N * n_grid, 0);
  dat.dR.assign(K * N * n_grid, 0.0);
  dat.dt.assign(n_grid, tau / static_cast<double>(n_grid));

  int max_unit = -1;
  for (int i = 0; i < N; ++i) {
    dat.A[i] = arm[i];
    dat.unit[i] = unit[i];
    dat.cluster_size[i] = cluster_size[i];
    dat.prob0[i] = prob0[i];
    dat.prob1[i] = prob1[i];
    max_unit = std::max(max_unit, dat.unit[i]);
    dat.observed_time[i] = std::min(tau, final_time[i]);

    if (final_status[i] == terminal_code && final_time[i] <= tau + 1e-12) {
      dat.death_event[i] = 1;
      dat.death_time[i] = final_time[i];
    }
    if (final_status[i] == 0 && final_time[i] < tau - 1e-12) {
      dat.censor_event[i] = 1;
      dat.censor_time[i] = final_time[i];
    }

    for (int k = 0; k < K; ++k) {
      const std::vector<double>& ev = all_recurrent[k][i];
      for (std::size_t q = 0; q < ev.size(); ++q) {
        if (finite_number(ev[q]) && ev[q] > 0.0 &&
            ev[q] <= dat.observed_time[i] + 1e-12 &&
            ev[q] <= tau + 1e-12) {
          dat.recurrent_time[k][i].push_back(ev[q]);
        }
      }
      std::sort(dat.recurrent_time[k][i].begin(),
                dat.recurrent_time[k][i].end());
    }
  }
  dat.n_units = max_unit + 1;

  for (int i = 0; i < N; ++i) {
    for (int g = 0; g < n_grid; ++g) {
      const double left = tau * static_cast<double>(g) /
        static_cast<double>(n_grid);
      const double right = tau * static_cast<double>(g + 1) /
        static_cast<double>(n_grid);
      const int id = idx2(i, g, n_grid);
      const bool risk = dat.observed_time[i] + 1e-12 >= left;
      dat.Y[id] = static_cast<unsigned char>(risk);
      dat.Ydag[id] = static_cast<unsigned char>(risk);
      if (dat.death_event[i] && dat.death_time[i] > left &&
          dat.death_time[i] <= right + 1e-12) {
        dat.dD[id] = 1;
      }
      if (dat.censor_event[i] && dat.censor_time[i] > left &&
          dat.censor_time[i] <= right + 1e-12) {
        dat.dC[id] = 1;
      }
      for (int k = 0; k < K; ++k) {
        double nr = 0.0;
        const std::vector<double>& ev = dat.recurrent_time[k][i];
        for (std::size_t q = 0; q < ev.size(); ++q) {
          if (ev[q] > left && ev[q] <= right + 1e-12) nr += 1.0;
        }
        dat.dR[idx3(k, i, g, N, n_grid)] = nr;
      }
    }
  }
  return dat;
}

DataFrame arm_table(const ArmMethodResults& x,
                    int arm,
                    const std::string& estimand,
                    double horizon) {
  CharacterVector methods = CharacterVector::create("dr", "ipcw", "or");
  CharacterVector estimands(3);
  NumericVector horizons(3, horizon);
  IntegerVector arms(3, arm);
  NumericVector estimate(3), std_error(3), burden(3), rmst(3);
  LogicalVector valid(3);
  for (int m = 0; m < 3; ++m) {
    estimands[m] = estimand;
    estimate[m] = x.method[m].rate;
    std_error[m] = x.method[m].se;
    burden[m] = x.method[m].mu;
    rmst[m] = x.method[m].nu;
    valid[m] = x.method[m].valid;
  }
  return DataFrame::create(
    _["time"] = horizons,
    _["estimand"] = estimands,
    _["arm"] = arms,
    _["method"] = methods,
    _["estimate"] = estimate,
    _["std_error"] = std_error,
    _["weighted_burden"] = burden,
    _["rmst"] = rmst,
    _["valid"] = valid
  );
}

DataFrame contrast_table(const ArmMethodResults& arm0,
                         const ArmMethodResults& arm1,
                         const std::string& estimand,
                         double horizon) {
  CharacterVector methods = CharacterVector::create("dr", "ipcw", "or");
  CharacterVector estimands(3);
  NumericVector horizons(3, horizon);
  NumericVector estimate(3), std_error(3), rate0(3), rate1(3),
    burden0(3), burden1(3), rmst0(3), rmst1(3),
    if_mean(3), if_sd(3), if_max_abs(3);
  LogicalVector valid(3);
  for (int m = 0; m < 3; ++m) {
    estimands[m] = estimand;
    rate0[m] = arm0.method[m].rate;
    rate1[m] = arm1.method[m].rate;
    burden0[m] = arm0.method[m].mu;
    burden1[m] = arm1.method[m].mu;
    rmst0[m] = arm0.method[m].nu;
    rmst1[m] = arm1.method[m].nu;
    valid[m] = arm0.method[m].valid && arm1.method[m].valid;
    estimate[m] = rate1[m] - rate0[m];
    const int U = static_cast<int>(arm0.method[m].if_rate.size());
    if (valid[m] && U > 1) {
      std::vector<double> IF(U, 0.0);
      for (int u = 0; u < U; ++u)
        IF[u] = arm1.method[m].if_rate[u] - arm0.method[m].if_rate[u];
      const double mn = std::accumulate(IF.begin(), IF.end(), 0.0) /
        static_cast<double>(U);
      double ss = 0.0, mx = 0.0;
      for (int u = 0; u < U; ++u) {
        ss += (IF[u] - mn) * (IF[u] - mn);
        mx = std::max(mx, std::fabs(IF[u]));
      }
      if_mean[m] = mn;
      if_sd[m] = std::sqrt(ss / static_cast<double>(U - 1));
      if_max_abs[m] = mx;
      std_error[m] = if_sd[m] / std::sqrt(static_cast<double>(U));
    } else {
      if_mean[m] = NA_REAL;
      if_sd[m] = NA_REAL;
      if_max_abs[m] = NA_REAL;
      std_error[m] = NA_REAL;
    }
  }
  return DataFrame::create(
    _["time"] = horizons,
    _["estimand"] = estimands,
    _["method"] = methods,
    _["estimate"] = estimate,
    _["std_error"] = std_error,
    _["rate_control"] = rate0,
    _["rate_treatment"] = rate1,
    _["weighted_burden_control"] = burden0,
    _["weighted_burden_treatment"] = burden1,
    _["rmst_control"] = rmst0,
    _["rmst_treatment"] = rmst1,
    _["if_mean"] = if_mean,
    _["if_sd"] = if_sd,
    _["if_max_abs"] = if_max_abs,
    _["valid"] = valid
  );
}

List influence_table(const ArmMethodResults& arm0,
                     const ArmMethodResults& arm1) {
  const char* method_name[3] = {"dr", "ipcw", "or"};
  List out(3);
  CharacterVector nm(3);
  for (int m = 0; m < 3; ++m) {
    const int U = static_cast<int>(arm0.method[m].if_rate.size());
    NumericVector if0(U), if1(U), contrast(U);
    for (int u = 0; u < U; ++u) {
      if0[u] = arm0.method[m].if_rate[u];
      if1[u] = arm1.method[m].if_rate[u];
      contrast[u] = if1[u] - if0[u];
    }
    out[m] = List::create(
      _["control"] = if0,
      _["treatment"] = if1,
      _["contrast"] = contrast
    );
    nm[m] = method_name[m];
  }
  out.attr("names") = nm;
  return out;
}

List nuisance_list(const CoxFit fitC[2],
                   const CoxFit fitD[2],
                   const std::vector<CoxFit> fitR[2],
                   const IntegerVector& recurrent_codes) {
  List recurrent0(fitR[0].size()), recurrent1(fitR[1].size());
  CharacterVector rn(recurrent_codes.size());
  for (int k = 0; k < recurrent_codes.size(); ++k) {
    recurrent0[k] = fit_to_list(fitR[0][k]);
    recurrent1[k] = fit_to_list(fitR[1][k]);
    rn[k] = std::to_string(recurrent_codes[k]);
  }
  recurrent0.attr("names") = rn;
  recurrent1.attr("names") = rn;
  return List::create(
    _["censoring"] = List::create(
      _["control"] = fit_to_list(fitC[0]),
      _["treatment"] = fit_to_list(fitC[1])
    ),
    _["terminal"] = List::create(
      _["control"] = fit_to_list(fitD[0]),
      _["treatment"] = fit_to_list(fitD[1])
    ),
    _["recurrent"] = List::create(
      _["control"] = recurrent0,
      _["treatment"] = recurrent1
    )
  );
}

DataFrame diagnostic_table(const CoxFit fitC[2],
                           const CoxFit fitD[2],
                           const std::vector<CoxFit> fitR[2],
                           const IntegerVector& recurrent_codes,
                           double horizon) {
  const int K = recurrent_codes.size();
  const int nrow = 4 + 2 * K;
  NumericVector time(nrow, horizon);
  IntegerVector arm(nrow), event_code(nrow), n_events(nrow), iterations(nrow);
  CharacterVector model(nrow);
  LogicalVector converged(nrow), valid(nrow);
  int r = 0;
  for (int a = 0; a < 2; ++a) {
    model[r] = "censoring";
    arm[r] = a;
    event_code[r] = 0;
    n_events[r] = fitC[a].n_events;
    iterations[r] = fitC[a].iterations;
    converged[r] = fitC[a].converged;
    valid[r] = fitC[a].valid;
    ++r;
    model[r] = "terminal";
    arm[r] = a;
    event_code[r] = NA_INTEGER;
    n_events[r] = fitD[a].n_events;
    iterations[r] = fitD[a].iterations;
    converged[r] = fitD[a].converged;
    valid[r] = fitD[a].valid;
    ++r;
    for (int k = 0; k < K; ++k) {
      model[r] = "recurrent";
      arm[r] = a;
      event_code[r] = recurrent_codes[k];
      n_events[r] = fitR[a][k].n_events;
      iterations[r] = fitR[a][k].iterations;
      converged[r] = fitR[a][k].converged;
      valid[r] = fitR[a][k].valid;
      ++r;
    }
  }
  return DataFrame::create(
    _["time"] = time,
    _["arm"] = arm,
    _["model"] = model,
    _["event_code"] = event_code,
    _["n_events"] = n_events,
    _["iterations"] = iterations,
    _["converged"] = converged,
    _["valid"] = valid
  );
}

}  // anonymous namespace

// [[Rcpp::export]]
List drcrtwa_fit_cpp(
    IntegerVector arm,
    IntegerVector unit,
    IntegerVector cluster_size,
    NumericMatrix X_censoring,
    NumericMatrix X_terminal,
    NumericMatrix X_recurrent,
    NumericVector final_time,
    IntegerVector final_status,
    List recurrent_times,
    IntegerVector recurrent_codes,
    NumericVector event_weights,
    int terminal_code,
    NumericVector horizons,
    IntegerVector n_grid,
    NumericVector prob0,
    NumericVector prob1,
    List target_weights,
    CharacterVector estimand_names,
    bool keep_influence,
    bool keep_nuisance) {
  const int N = arm.size();
  const int K = recurrent_codes.size();
  const int H = horizons.size();
  const int L = target_weights.size();
  if (N < 2) stop("At least two participants are required.");
  if (unit.size() != N || cluster_size.size() != N ||
      X_censoring.nrow() != N || X_terminal.nrow() != N ||
      X_recurrent.nrow() != N || final_time.size() != N ||
      final_status.size() != N || prob0.size() != N || prob1.size() != N) {
    stop("Subject-level inputs have inconsistent dimensions.");
  }
  if (K < 1 || recurrent_times.size() != K || event_weights.size() != K)
    stop("Recurrent-event codes, histories, and weights have inconsistent dimensions.");
  if (H < 1 || n_grid.size() != H)
    stop("The horizon and integration-grid vectors must have equal positive length.");
  if (L < 1 || estimand_names.size() != L)
    stop("Target weights and estimand names have inconsistent dimensions.");

  std::vector< std::vector< std::vector<double> > > all_recurrent(
    K, std::vector< std::vector<double> >(N));
  for (int k = 0; k < K; ++k) {
    List by_subject = recurrent_times[k];
    if (by_subject.size() != N)
      stop("Each recurrent-event history list must have one element per participant.");
    for (int i = 0; i < N; ++i) {
      NumericVector ev = by_subject[i];
      all_recurrent[k][i].reserve(ev.size());
      for (int q = 0; q < ev.size(); ++q)
        if (R_finite(ev[q])) all_recurrent[k][i].push_back(ev[q]);
    }
  }

  std::vector<double> Xc = matrix_to_row_major(X_censoring);
  std::vector<double> Xd = matrix_to_row_major(X_terminal);
  std::vector<double> Xr = matrix_to_row_major(X_recurrent);
  std::vector<double> weights(K, 0.0);
  for (int k = 0; k < K; ++k) weights[k] = event_weights[k];

  std::vector< std::vector<double> > target(L, std::vector<double>(N, 0.0));
  for (int l = 0; l < L; ++l) {
    NumericVector w = target_weights[l];
    if (w.size() != N)
      stop("Each target-weight vector must have one element per participant.");
    for (int i = 0; i < N; ++i) target[l][i] = w[i];
  }

  List horizon_results(H);
  List final_nuisance;
  for (int h = 0; h < H; ++h) {
    if (!R_finite(horizons[h]) || horizons[h] <= 0.0 || n_grid[h] < 2)
      stop("All analysis times and integration-grid sizes must be positive.");
    DataSet dat = build_dataset(
      arm, unit, cluster_size, prob0, prob1, final_time, final_status,
      all_recurrent, terminal_code, horizons[h], n_grid[h]);

    CoxFit fitC[2], fitD[2];
    std::vector<CoxFit> fitR[2];
    for (int a = 0; a < 2; ++a) {
      fitC[a] = fit_cox_lwyy(dat, Xc, X_censoring.ncol(), a, -2);
      fitD[a] = fit_cox_lwyy(dat, Xd, X_terminal.ncol(), a, -1);
      fitR[a].reserve(K);
      for (int k = 0; k < K; ++k)
        fitR[a].push_back(fit_cox_lwyy(
          dat, Xr, X_recurrent.ncol(), a, k));
    }

    List contrast_by_estimand(L), arms_by_estimand(L), influence_by_estimand;
    CharacterVector level_names(L);
    if (keep_influence) influence_by_estimand = List(L);
    for (int l = 0; l < L; ++l) {
      ArmMethodResults a0 = compute_arm_methods(
        dat, 0, fitC[0], fitD[0], fitR[0], weights, target[l]);
      ArmMethodResults a1 = compute_arm_methods(
        dat, 1, fitC[1], fitD[1], fitR[1], weights, target[l]);
      const std::string level = as<std::string>(estimand_names[l]);
      level_names[l] = level;
      contrast_by_estimand[l] = contrast_table(a0, a1, level, horizons[h]);
      Function rbind("rbind");
      arms_by_estimand[l] = rbind(
        arm_table(a0, 0, level, horizons[h]),
        arm_table(a1, 1, level, horizons[h]));
      if (keep_influence)
        influence_by_estimand[l] = influence_table(a0, a1);
    }
    contrast_by_estimand.attr("names") = level_names;
    arms_by_estimand.attr("names") = level_names;
    if (keep_influence) influence_by_estimand.attr("names") = level_names;

    List one = List::create(
      _["contrasts"] = contrast_by_estimand,
      _["arms"] = arms_by_estimand,
      _["diagnostics"] = diagnostic_table(
        fitC, fitD, fitR, recurrent_codes, horizons[h])
    );
    if (keep_influence) one["influence"] = influence_by_estimand;
    horizon_results[h] = one;

    if (keep_nuisance && h == H - 1)
      final_nuisance = nuisance_list(fitC, fitD, fitR, recurrent_codes);
  }

  int n_units_out = 0;
  for (int i = 0; i < N; ++i)
    n_units_out = std::max(n_units_out, unit[i] + 1);

  List out = List::create(
    _["horizons"] = horizon_results,
    _["n_subjects"] = N,
    _["n_units"] = n_units_out,
    _["terminal_code"] = terminal_code,
    _["recurrent_codes"] = recurrent_codes,
    _["event_weights"] = event_weights
  );
  if (keep_nuisance) out["nuisance"] = final_nuisance;
  return out;
}
