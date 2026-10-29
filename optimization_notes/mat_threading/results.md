## Runs

| Run | Variable sites | min/M states | Wall (s) | STL step (us) | propose | evaluate | store+accept/restore+adapt |
|---|---|---|---|---|---|---|---|
| big-2-local | 18217 | 1.48 |  | 80.47 | 31.87 | 44.05 | 4.55 |
| big-2-path | 18217 | 4.42 |  | 171.28 | 120.58 | 45.12 | 5.58 |
| big-3-gamma-local | 2405 | 3.80 |  | 99.96 | 45.48 | 50.38 | 4.10 |
| big-3-local | 2405 | 1.18 |  | 62.00 | 29.75 | 28.72 | 3.53 |
| big-3-path | 2405 | 3.12 |  | 140.65 | 106.29 | 30.10 | 4.27 |
| sim-100k-r005 | 13907 | 1.08 | 33 | 46.68 | 13.37 | 29.14 | 4.17 |
| sim-100k-r02 | 44866 | 1.68 | 50 | 53.75 | 18.82 | 30.77 | 4.18 |
| sim-20k-r005 | 2769 | 0.86 | 26 | 42.24 | 11.08 | 27.37 | 3.80 |
| sim-20k-r02 | 8988 | 0.99 | 30 | 44.67 | 12.55 | 28.43 | 3.69 |
| sim-300k-r02 | 134581 | 3.44 |  | 67.73 | 27.90 | 35.36 | 4.47 |
| small-2-local | 294 | 0.37 | 22 | 17.89 | 9.50 | 7.05 | 1.34 |
| small-2-path | 294 | 0.76 | 23 | 37.20 | 28.05 | 7.48 | 1.67 |
| wide-1-local | 31181 | 0.96 |  | 25.89 | 20.55 | 4.30 | 1.05 |
| wide-1-path | 31181 | 4.57 |  | 72.02 | 64.77 | 5.59 | 1.65 |

## Step profile of every operator (microseconds per step)

### big-2-local

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 820524 | 3.61 | 31.87 | 44.05 | 0.78 | 0.16 | 80.47 |
| subtreeJump(treeModel) | 82527 | 3.69 | 93.60 | 42.72 | 1.13 | 0.02 | 141.16 |
| mutationHistoryGibbs(treeModel) | 82214 | 3.71 | 39.41 | 7.21 | 0.39 | 0.02 | 50.74 |
| scale(clock.rate) | 1621 | 4.63 | 5.88 | 1824.90 | 4.47 | 1.57 | 1841.44 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1584 | 5.39 | 304.10 | 1365.21 | 8.77 | 1.81 | 1685.28 |
| scale(constant.popSize) | 1530 | 4.59 | 6.54 | 15.38 | 0.80 | 1.36 | 28.67 |
| (all operators) | 990000 | 3.63 | 37.99 | 45.86 | 0.80 | 0.15 | 88.43 |

### big-2-path

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 410591 | 4.02 | 120.58 | 45.12 | 1.37 | 0.19 | 171.28 |
| subtreeJump(treeModel) | 40978 | 4.19 | 1133.46 | 50.72 | 3.10 | 0.03 | 1191.50 |
| mutationHistoryGibbs(treeModel) | 41157 | 4.25 | 41.99 | 7.83 | 0.71 | 0.02 | 54.81 |
| scale(constant.popSize) | 768 | 5.32 | 10.25 | 15.68 | 1.16 | 2.10 | 34.50 |
| scale(clock.rate) | 773 | 5.29 | 8.34 | 1894.34 | 7.13 | 2.27 | 1917.35 |
| up(nodeHeights(treeModel))/down(clock.rate) | 733 | 6.22 | 329.18 | 1387.55 | 10.57 | 2.84 | 1736.36 |
| (all operators) | 495000 | 4.06 | 197.86 | 47.32 | 1.48 | 0.17 | 250.89 |

### big-3-gamma-local

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 808021 | 2.67 | 45.48 | 50.38 | 1.32 | 0.11 | 99.96 |
| mutationHistoryGibbs(treeModel) | 80951 | 3.00 | 58.94 | 25.87 | 0.80 | 0.03 | 88.62 |
| subtreeJump(treeModel) | 80494 | 2.99 | 89.60 | 49.03 | 1.70 | 0.03 | 143.36 |
| scale(alpha) | 2695 | 3.90 | 5.88 | 3490.58 | 12.02 | 2.12 | 3514.50 |
| siteRateCategoryGibbs(treeModel) | 7974 | 3.41 | 8429.43 | 3469.65 | 5.44 | 0.12 | 11908.04 |
| scale(kappa) | 891 | 4.00 | 6.96 | 3381.18 | 23.80 | 2.10 | 3418.04 |
| scale(constant.popSize) | 2680 | 3.88 | 5.02 | 9.23 | 1.00 | 1.63 | 20.75 |
| scale(clock.rate) | 2734 | 3.82 | 5.68 | 3377.03 | 8.46 | 2.11 | 3397.09 |
| up(nodeHeights(treeModel))/down(clock.rate) | 2685 | 4.74 | 179.46 | 2149.93 | 10.39 | 2.35 | 2346.87 |
| deltaExchange(frequencies) | 875 | 4.00 | 11.75 | 3354.44 | 10.49 | 2.82 | 3383.50 |
| (all operators) | 990000 | 2.75 | 117.67 | 105.86 | 1.44 | 0.12 | 227.84 |

### big-3-local

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeJump(treeModel) | 81251 | 2.63 | 67.25 | 27.89 | 1.29 | 0.03 | 99.08 |
| mutationHistoryGibbs(treeModel) | 81183 | 2.63 | 39.45 | 6.56 | 0.59 | 0.02 | 49.25 |
| subtreeLeap(treeModel) | 817454 | 2.43 | 29.75 | 28.72 | 1.00 | 0.10 | 62.00 |
| up(nodeHeights(treeModel))/down(clock.rate) | 2713 | 4.19 | 169.05 | 721.03 | 7.18 | 2.08 | 903.53 |
| scale(clock.rate) | 2819 | 3.44 | 5.15 | 1013.34 | 5.06 | 1.79 | 1028.78 |
| scale(constant.popSize) | 2756 | 3.43 | 5.31 | 9.02 | 0.86 | 1.57 | 20.19 |
| deltaExchange(frequencies) | 924 | 3.83 | 11.27 | 984.42 | 4.86 | 2.55 | 1006.93 |
| scale(kappa) | 900 | 3.58 | 6.74 | 976.60 | 5.85 | 1.73 | 994.50 |
| (all operators) | 990000 | 2.48 | 33.83 | 33.24 | 1.02 | 0.11 | 70.67 |

### big-3-path

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 408479 | 2.75 | 106.29 | 30.10 | 1.39 | 0.13 | 140.65 |
| subtreeJump(treeModel) | 40897 | 3.02 | 610.49 | 34.27 | 3.17 | 0.05 | 651.01 |
| mutationHistoryGibbs(treeModel) | 40700 | 3.10 | 40.84 | 7.16 | 0.78 | 0.03 | 51.91 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1343 | 4.82 | 171.39 | 693.28 | 6.66 | 2.46 | 878.61 |
| scale(constant.popSize) | 1337 | 3.98 | 7.22 | 9.13 | 1.03 | 1.87 | 23.23 |
| scale(kappa) | 430 | 4.15 | 10.39 | 1005.66 | 4.39 | 2.17 | 1026.76 |
| scale(clock.rate) | 1360 | 3.90 | 7.13 | 1041.45 | 4.98 | 2.27 | 1059.73 |
| deltaExchange(frequencies) | 454 | 4.23 | 14.51 | 1022.93 | 4.65 | 3.15 | 1049.48 |
| (all operators) | 495000 | 2.82 | 142.04 | 34.84 | 1.51 | 0.13 | 181.34 |

### sim-100k-r005

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 408202 | 2.45 | 13.37 | 29.14 | 1.61 | 0.11 | 46.68 |
| mutationHistoryGibbs(treeModel) | 40895 | 2.65 | 50.12 | 7.36 | 1.01 | 0.03 | 61.17 |
| subtreeJump(treeModel) | 40939 | 2.65 | 111.86 | 29.24 | 2.27 | 0.03 | 146.05 |
| scale(clock.rate) | 1371 | 3.36 | 5.88 | 1097.61 | 7.38 | 2.14 | 1116.37 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1355 | 4.40 | 191.84 | 791.38 | 8.73 | 2.59 | 998.94 |
| deltaExchange(frequencies) | 462 | 3.94 | 13.77 | 1093.15 | 24.23 | 3.46 | 1138.57 |
| scale(constant.popSize) | 1301 | 3.35 | 6.52 | 9.17 | 1.26 | 1.82 | 22.13 |
| scale(kappa) | 475 | 3.56 | 8.84 | 1068.44 | 6.91 | 1.96 | 1089.71 |
| (all operators) | 495000 | 2.50 | 25.00 | 34.33 | 1.67 | 0.12 | 63.62 |

### sim-100k-r02

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 408624 | 2.78 | 18.82 | 30.77 | 1.28 | 0.12 | 53.75 |
| mutationHistoryGibbs(treeModel) | 40677 | 3.05 | 129.63 | 9.16 | 1.05 | 0.03 | 142.92 |
| subtreeJump(treeModel) | 40814 | 3.01 | 359.71 | 36.25 | 2.35 | 0.04 | 401.35 |
| scale(clock.rate) | 1356 | 3.86 | 6.99 | 1269.02 | 8.11 | 2.35 | 1290.33 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1319 | 4.77 | 192.77 | 959.71 | 9.92 | 2.68 | 1169.86 |
| deltaExchange(frequencies) | 461 | 4.10 | 14.87 | 1248.18 | 10.78 | 3.48 | 1281.40 |
| scale(constant.popSize) | 1300 | 4.02 | 7.58 | 9.40 | 1.06 | 1.93 | 23.99 |
| scale(kappa) | 449 | 3.98 | 9.97 | 1249.53 | 7.90 | 2.35 | 1273.72 |
| (all operators) | 495000 | 2.83 | 56.42 | 37.49 | 1.40 | 0.13 | 98.27 |

### sim-20k-r005

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeJump(treeModel) | 40770 | 2.61 | 42.18 | 26.46 | 1.57 | 0.03 | 72.85 |
| subtreeLeap(treeModel) | 408769 | 2.43 | 11.08 | 27.37 | 1.26 | 0.11 | 42.24 |
| mutationHistoryGibbs(treeModel) | 40512 | 2.62 | 18.67 | 6.45 | 0.73 | 0.04 | 28.51 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1327 | 4.10 | 180.57 | 720.50 | 8.50 | 2.37 | 916.04 |
| deltaExchange(frequencies) | 481 | 3.80 | 13.38 | 967.55 | 7.48 | 2.72 | 994.92 |
| scale(constant.popSize) | 1311 | 3.36 | 6.04 | 8.90 | 1.33 | 1.44 | 21.07 |
| scale(clock.rate) | 1371 | 3.33 | 5.79 | 990.64 | 5.55 | 1.78 | 1007.09 |
| scale(kappa) | 459 | 3.49 | 8.09 | 950.06 | 6.69 | 1.80 | 970.12 |
| (all operators) | 495000 | 2.47 | 14.69 | 31.83 | 1.28 | 0.12 | 50.38 |

### sim-20k-r02

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 408314 | 2.64 | 12.55 | 28.43 | 0.94 | 0.11 | 44.67 |
| subtreeJump(treeModel) | 40985 | 2.78 | 85.92 | 27.17 | 1.29 | 0.03 | 117.19 |
| mutationHistoryGibbs(treeModel) | 40829 | 2.80 | 39.83 | 6.53 | 0.57 | 0.03 | 49.76 |
| deltaExchange(frequencies) | 449 | 3.82 | 14.01 | 1002.25 | 4.28 | 2.66 | 1027.02 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1380 | 4.43 | 179.56 | 726.54 | 6.42 | 2.17 | 919.12 |
| scale(kappa) | 445 | 3.61 | 8.33 | 1004.48 | 6.36 | 1.90 | 1024.68 |
| scale(clock.rate) | 1288 | 3.49 | 5.66 | 1015.35 | 4.44 | 1.65 | 1030.59 |
| scale(constant.popSize) | 1310 | 3.46 | 6.23 | 9.30 | 0.85 | 1.24 | 21.09 |
| (all operators) | 495000 | 2.68 | 21.31 | 32.74 | 0.97 | 0.11 | 57.81 |

### sim-300k-r02

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 407888 | 3.10 | 27.90 | 35.36 | 1.25 | 0.12 | 67.73 |
| subtreeJump(treeModel) | 41324 | 3.45 | 1182.19 | 56.50 | 3.72 | 0.08 | 1245.94 |
| mutationHistoryGibbs(treeModel) | 40843 | 3.51 | 327.25 | 13.62 | 1.56 | 0.04 | 345.98 |
| up(nodeHeights(treeModel))/down(clock.rate) | 1339 | 5.56 | 208.71 | 1337.98 | 10.60 | 3.06 | 1565.90 |
| scale(clock.rate) | 1356 | 4.36 | 7.22 | 1658.89 | 7.62 | 2.63 | 1680.73 |
| scale(constant.popSize) | 1342 | 4.19 | 7.93 | 9.55 | 1.30 | 2.12 | 25.10 |
| deltaExchange(frequencies) | 429 | 4.40 | 16.25 | 1647.53 | 6.39 | 3.88 | 1678.45 |
| scale(kappa) | 479 | 4.30 | 10.24 | 1646.58 | 8.69 | 2.60 | 1672.43 |
| (all operators) | 495000 | 3.18 | 149.31 | 46.19 | 1.54 | 0.14 | 200.36 |

### small-2-local

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 789540 | 0.60 | 9.50 | 7.05 | 0.58 | 0.16 | 17.89 |
| up(nodeHeights(treeModel))/down(clock.rate) | 14059 | 1.20 | 29.03 | 118.09 | 1.77 | 0.48 | 150.56 |
| mutationHistoryGibbs(treeModel) | 78778 | 0.74 | 11.87 | 2.77 | 0.31 | 0.02 | 15.70 |
| scale(constant.popSize) | 14156 | 0.82 | 0.74 | 1.92 | 0.50 | 0.34 | 4.32 |
| subtreeJump(treeModel) | 79055 | 0.73 | 15.38 | 6.50 | 0.68 | 0.02 | 23.31 |
| scale(clock.rate) | 14412 | 0.83 | 0.91 | 166.08 | 1.23 | 0.40 | 169.45 |
| (all operators) | 990000 | 0.64 | 10.18 | 10.48 | 0.59 | 0.15 | 22.05 |

### small-2-path

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 394920 | 0.72 | 28.05 | 7.48 | 0.78 | 0.17 | 37.20 |
| mutationHistoryGibbs(treeModel) | 39457 | 0.87 | 11.65 | 2.99 | 0.40 | 0.03 | 15.94 |
| up(nodeHeights(treeModel))/down(clock.rate) | 7112 | 1.62 | 31.61 | 125.13 | 1.86 | 0.72 | 160.94 |
| subtreeJump(treeModel) | 39226 | 0.87 | 99.70 | 7.64 | 1.26 | 0.03 | 109.50 |
| scale(clock.rate) | 7056 | 1.08 | 1.49 | 174.90 | 1.39 | 0.63 | 179.50 |
| scale(constant.popSize) | 7229 | 1.09 | 1.46 | 1.96 | 0.59 | 0.53 | 5.64 |
| (all operators) | 495000 | 0.77 | 31.71 | 11.13 | 0.81 | 0.17 | 44.58 |

### wide-1-local

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 793810 | 0.41 | 20.55 | 4.30 | 0.47 | 0.17 | 25.89 |
| subtreeJump(treeModel) | 79863 | 0.48 | 243.00 | 5.93 | 0.80 | 0.02 | 250.23 |
| mutationHistoryGibbs(treeModel) | 79283 | 0.49 | 199.67 | 3.27 | 0.37 | 0.02 | 203.81 |
| scale(constant.popSize) | 37044 | 0.66 | 0.37 | 0.90 | 0.45 | 0.21 | 2.58 |
| (all operators) | 990000 | 0.43 | 52.09 | 4.22 | 0.48 | 0.14 | 57.36 |

### wide-1-path

| Operator | Steps | store | propose | evaluate | accept/restore | adapt | total |
|---|---|---|---|---|---|---|---|
| subtreeLeap(treeModel) | 396735 | 0.66 | 64.77 | 5.59 | 0.80 | 0.19 | 72.02 |
| mutationHistoryGibbs(treeModel) | 39717 | 0.92 | 215.57 | 4.43 | 0.88 | 0.04 | 221.83 |
| scale(constant.popSize) | 18709 | 1.26 | 1.61 | 0.98 | 0.68 | 0.38 | 4.91 |
| subtreeJump(treeModel) | 39839 | 0.81 | 2167.99 | 13.75 | 2.89 | 0.07 | 2185.52 |
| (all operators) | 495000 | 0.72 | 243.76 | 5.98 | 0.97 | 0.18 | 251.60 |

## MAT phases

Microseconds per call, and per subtree move (local.before + local.after, or path.*) or per evaluation (lik.*), from the whole run.

### big-2-local

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 712195 | 10.621 | 14.91 | 14.91 |
| local.runs | 1089452 | 1.465 | 1.35 | 2.06 |
| local.chains | 712195 | 2.042 | 2.87 | 2.87 |
| local.scoreOld | 712195 | 8.427 | 11.83 | 11.83 |
| local.after | 712195 | 19.408 | 27.25 | 27.25 |
| local.newChain | 712195 | 1.507 | 2.12 | 2.12 |
| local.sampleNew | 712195 | 12.060 | 16.93 | 16.93 |
| gibbs | 82214 | 3.226 | 39.24 | 4.53 |
| branches | 794409 | 5.905 | 7.43 | 8.29 |
| store.updateBranch | 2083088 | 3.097 | 1.49 | 4.35 |
| store.updateBranchSites | 3035452 | 2.142 | 0.71 | 3.01 |
| store.journalClone | 1840203 | 2.437 | 1.32 | 3.42 |
| store.commit | 1702492 | 0.077 | 0.05 | 0.11 |
| store.revert | 712492 | 0.215 | 0.30 | 0.30 |
| lik.updateMAT | 1782765 | 0.028 | 0.02 | 0.03 |
| lik.ensureRates | 3752412 | 0.060 | 0.02 | 0.06 |
| lik.calculate | 988356 | 12.121 | 12.26 | 12.26 |
| lik.branchTerms | 988356 | 8.871 | 8.98 | 8.98 |
| lik.sum | 988356 | 3.025 | 3.06 | 3.06 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 712195 |  |
| path_moves | 0 |  |
| gibbs_updates | 82214 |  |
| candidate_sites | 7927919 | 11.13 |
| changed_sites | 38604543 | 54.21 |
| hot_sites | 13129051 | 18.43 |
| branch_scan_patterns | 37944230345 | 53277.87 |
| branch_site_updates | 229798500 | 322.66 |
| journal_clone_bytes | 33530338863 | 47080.28 |
| runs_recomputed | 1089452 | 1.53 |
| terms_recomputed | 24144932 | 33.90 |
| evaluations | 988356 |  |
| transition_cache_misses | 13528643 | 19.00 |
| transition_cache_clears | 826 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### big-2-path

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.runs | 73718 | 0.132 | 1.78 | 0.37 |
| path.before | 358309 | 43.096 | 120.28 | 120.28 |
| path.candidates | 358309 | 3.869 | 10.80 | 10.80 |
| path.scoreOld | 358309 | 39.141 | 109.24 | 109.24 |
| path.after | 358309 | 50.291 | 140.36 | 140.36 |
| path.sampleNew | 358309 | 43.998 | 122.79 | 122.79 |
| gibbs | 41157 | 1.715 | 41.68 | 4.79 |
| branches | 399466 | 6.269 | 15.69 | 17.50 |
| store.updateBranch | 1048745 | 1.490 | 1.42 | 4.16 |
| store.updateBranchSites | 5577633 | 4.107 | 0.74 | 11.46 |
| store.journalClone | 3113456 | 4.057 | 1.30 | 11.32 |
| store.commit | 850875 | 0.080 | 0.09 | 0.22 |
| store.revert | 355875 | 0.217 | 0.61 | 0.61 |
| lik.updateMAT | 895328 | 0.024 | 0.03 | 0.05 |
| lik.ensureRates | 1653103 | 0.028 | 0.02 | 0.06 |
| lik.calculate | 495862 | 6.720 | 13.55 | 13.55 |
| lik.branchTerms | 495862 | 5.057 | 10.20 | 10.20 |
| lik.sum | 495862 | 1.517 | 3.06 | 3.06 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 0 |  |
| path_moves | 358309 |  |
| gibbs_updates | 41157 |  |
| candidate_sites | 209822983 | 585.59 |
| changed_sites | 24310623 | 67.85 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 19104809605 | 53319.37 |
| branch_site_updates | 633202534 | 1767.20 |
| journal_clone_bytes | 56730281776 | 158327.82 |
| runs_recomputed | 73718 | 0.21 |
| terms_recomputed | 17004944 | 47.46 |
| evaluations | 495862 |  |
| transition_cache_misses | 6549952 | 18.28 |
| transition_cache_clears | 400 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### big-3-gamma-local

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 708044 | 15.307 | 21.62 | 21.62 |
| local.runs | 924752 | 0.343 | 0.37 | 0.48 |
| local.chains | 708044 | 1.405 | 1.98 | 1.98 |
| local.scoreOld | 708044 | 13.660 | 19.29 | 19.29 |
| local.after | 708044 | 24.850 | 35.10 | 35.10 |
| local.newChain | 708044 | 1.372 | 1.94 | 1.94 |
| local.sampleNew | 708044 | 17.747 | 25.06 | 25.06 |
| gibbs | 80951 | 4.740 | 58.55 | 6.69 |
| branches | 788995 | 5.831 | 7.39 | 8.24 |
| store.updateBranch | 2074448 | 1.640 | 0.79 | 2.32 |
| store.updateBranchSites | 2114556 | 3.154 | 1.49 | 4.45 |
| store.journalClone | 1367933 | 0.480 | 0.35 | 0.68 |
| store.setCategory | 19177470 | 40.010 | 2.09 | 56.51 |
| store.commit | 1696281 | 0.067 | 0.04 | 0.09 |
| store.revert | 706281 | 0.242 | 0.34 | 0.34 |
| lik.updateMAT | 1786879 | 0.041 | 0.02 | 0.04 |
| lik.ensureRates | 3557529 | 0.517 | 0.15 | 0.52 |
| lik.calculate | 989910 | 84.877 | 85.74 | 85.74 |
| lik.branchTerms | 989910 | 82.188 | 83.03 | 83.03 |
| lik.sum | 989910 | 2.036 | 2.06 | 2.06 |
| categories.resample | 7974 | 67.208 | 8428.38 | 94.92 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 708044 |  |
| path_moves | 0 |  |
| gibbs_updates | 80951 |  |
| candidate_sites | 47554832 | 67.16 |
| changed_sites | 188713095 | 266.53 |
| hot_sites | 71934757 | 101.60 |
| branch_scan_patterns | 4994808555 | 7054.38 |
| branch_site_updates | 538817899 | 760.99 |
| journal_clone_bytes | 3295350597 | 4654.16 |
| runs_recomputed | 924752 | 1.31 |
| terms_recomputed | 49171361 | 69.45 |
| evaluations | 989910 |  |
| transition_cache_misses | 160199102 | 226.26 |
| transition_cache_clears | 17304 | 0.02 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 7853802 | 11.09 |

### big-3-local

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 715706 | 9.393 | 13.12 | 13.12 |
| local.runs | 935571 | 0.337 | 0.36 | 0.47 |
| local.chains | 715706 | 1.347 | 1.88 | 1.88 |
| local.scoreOld | 715706 | 7.858 | 10.98 | 10.98 |
| local.after | 715706 | 16.871 | 23.57 | 23.57 |
| local.newChain | 715706 | 1.171 | 1.64 | 1.64 |
| local.sampleNew | 715706 | 9.993 | 13.96 | 13.96 |
| gibbs | 81183 | 3.182 | 39.19 | 4.45 |
| branches | 796889 | 5.821 | 7.31 | 8.13 |
| store.updateBranch | 2097567 | 1.648 | 0.79 | 2.30 |
| store.updateBranchSites | 2128785 | 3.200 | 1.50 | 4.47 |
| store.journalClone | 1378694 | 0.556 | 0.40 | 0.78 |
| store.commit | 1702823 | 0.064 | 0.04 | 0.09 |
| store.revert | 712823 | 0.198 | 0.28 | 0.28 |
| lik.updateMAT | 1784740 | 0.034 | 0.02 | 0.03 |
| lik.ensureRates | 3572618 | 0.140 | 0.04 | 0.14 |
| lik.calculate | 987851 | 12.969 | 13.13 | 13.13 |
| lik.branchTerms | 987851 | 10.967 | 11.10 | 11.10 |
| lik.sum | 987851 | 1.710 | 1.73 | 1.73 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 715706 |  |
| path_moves | 0 |  |
| gibbs_updates | 81183 |  |
| candidate_sites | 48010147 | 67.08 |
| changed_sites | 187254635 | 261.64 |
| hot_sites | 72110211 | 100.75 |
| branch_scan_patterns | 5050111968 | 7056.13 |
| branch_site_updates | 540545834 | 755.26 |
| journal_clone_bytes | 3321273846 | 4640.56 |
| runs_recomputed | 935571 | 1.31 |
| terms_recomputed | 27186502 | 37.99 |
| evaluations | 987851 |  |
| transition_cache_misses | 16970618 | 23.71 |
| transition_cache_clears | 3445 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### big-3-path

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.runs | 54224 | 0.025 | 0.47 | 0.07 |
| path.before | 361735 | 27.684 | 76.53 | 76.53 |
| path.candidates | 361735 | 5.191 | 14.35 | 14.35 |
| path.scoreOld | 361735 | 22.384 | 61.88 | 61.88 |
| path.after | 361735 | 38.563 | 106.61 | 106.61 |
| path.sampleNew | 361735 | 29.482 | 81.50 | 81.50 |
| gibbs | 40700 | 1.646 | 40.45 | 4.55 |
| branches | 402435 | 9.060 | 22.51 | 25.05 |
| store.updateBranch | 1060528 | 0.732 | 0.69 | 2.02 |
| store.updateBranchSites | 5087949 | 7.514 | 1.48 | 20.77 |
| store.journalClone | 2873021 | 0.869 | 0.30 | 2.40 |
| store.commit | 851231 | 0.066 | 0.08 | 0.18 |
| store.revert | 356231 | 0.217 | 0.61 | 0.60 |
| lik.updateMAT | 896087 | 0.017 | 0.02 | 0.03 |
| lik.ensureRates | 1660257 | 0.070 | 0.04 | 0.14 |
| lik.calculate | 493652 | 6.829 | 13.83 | 13.83 |
| lik.branchTerms | 493652 | 5.791 | 11.73 | 11.73 |
| lik.sum | 493652 | 0.858 | 1.74 | 1.74 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 0 |  |
| path_moves | 361735 |  |
| gibbs_updates | 40700 |  |
| candidate_sites | 395316589 | 1092.83 |
| changed_sites | 97540112 | 269.65 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 2553607452 | 7059.33 |
| branch_site_updates | 1495971727 | 4135.55 |
| journal_clone_bytes | 6921107589 | 19133.09 |
| runs_recomputed | 54224 | 0.15 |
| terms_recomputed | 17560012 | 48.54 |
| evaluations | 493652 |  |
| transition_cache_misses | 8236311 | 22.77 |
| transition_cache_clears | 1656 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### sim-100k-r005

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 207281 | 2.818 | 13.60 | 13.60 |
| local.runs | 360498 | 0.263 | 0.73 | 1.27 |
| local.chains | 207281 | 0.822 | 3.97 | 3.97 |
| local.scoreOld | 207281 | 1.931 | 9.32 | 9.32 |
| local.after | 207281 | 6.126 | 29.55 | 29.55 |
| local.newChain | 207281 | 0.087 | 0.42 | 0.42 |
| local.sampleNew | 207281 | 4.899 | 23.63 | 23.63 |
| gibbs | 40895 | 2.038 | 49.82 | 9.83 |
| branches | 248176 | 1.233 | 4.97 | 5.95 |
| store.updateBranch | 587333 | 0.652 | 1.11 | 3.15 |
| store.updateBranchSites | 527720 | 0.284 | 0.54 | 1.37 |
| store.journalClone | 329798 | 0.464 | 1.41 | 2.24 |
| store.commit | 852340 | 0.021 | 0.03 | 0.10 |
| store.revert | 357340 | 0.080 | 0.22 | 0.39 |
| lik.updateMAT | 749307 | 0.041 | 0.05 | 0.08 |
| lik.ensureRates | 1204764 | 0.074 | 0.06 | 0.15 |
| lik.calculate | 501131 | 7.309 | 14.58 | 14.59 |
| lik.branchTerms | 501131 | 6.246 | 12.46 | 12.46 |
| lik.sum | 501131 | 0.871 | 1.74 | 1.74 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 207281 |  |
| path_moves | 0 |  |
| gibbs_updates | 40895 |  |
| candidate_sites | 3776124 | 18.22 |
| changed_sites | 33394300 | 161.11 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 8170389363 | 39416.97 |
| branch_site_updates | 50856275 | 245.35 |
| journal_clone_bytes | 4587819978 | 22133.34 |
| runs_recomputed | 360498 | 1.74 |
| terms_recomputed | 12703770 | 61.29 |
| evaluations | 501131 |  |
| transition_cache_misses | 8335017 | 40.21 |
| transition_cache_clears | 1744 | 0.01 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### sim-100k-r02

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 154944 | 5.177 | 33.41 | 33.41 |
| local.runs | 294616 | 0.494 | 1.68 | 3.19 |
| local.chains | 154944 | 1.965 | 12.68 | 12.68 |
| local.scoreOld | 154944 | 3.154 | 20.35 | 20.36 |
| local.after | 154944 | 15.996 | 103.23 | 103.24 |
| local.newChain | 154944 | 0.098 | 0.63 | 0.63 |
| local.sampleNew | 154944 | 13.363 | 86.24 | 86.24 |
| gibbs | 40677 | 5.255 | 129.19 | 33.92 |
| branches | 195621 | 2.888 | 14.76 | 18.64 |
| store.updateBranch | 433837 | 1.273 | 2.93 | 8.22 |
| store.updateBranchSites | 470320 | 0.895 | 1.90 | 5.78 |
| store.journalClone | 276796 | 0.766 | 2.77 | 4.94 |
| store.commit | 852545 | 0.024 | 0.03 | 0.15 |
| store.revert | 357545 | 0.077 | 0.21 | 0.50 |
| lik.updateMAT | 692682 | 0.026 | 0.04 | 0.05 |
| lik.ensureRates | 1043247 | 0.094 | 0.09 | 0.19 |
| lik.calculate | 497061 | 8.417 | 16.93 | 16.93 |
| lik.branchTerms | 497061 | 7.350 | 14.79 | 14.79 |
| lik.sum | 497061 | 0.867 | 1.74 | 1.74 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 154944 |  |
| path_moves | 0 |  |
| gibbs_updates | 40677 |  |
| candidate_sites | 11355678 | 73.29 |
| changed_sites | 126895355 | 818.98 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 19466266190 | 125634.20 |
| branch_site_updates | 195198414 | 1259.80 |
| journal_clone_bytes | 12419836520 | 80156.94 |
| runs_recomputed | 294616 | 1.90 |
| terms_recomputed | 12531157 | 80.88 |
| evaluations | 497061 |  |
| transition_cache_misses | 8173983 | 52.75 |
| transition_cache_clears | 1690 | 0.01 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### sim-20k-r005

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 294857 | 1.908 | 6.47 | 6.47 |
| local.runs | 418124 | 0.107 | 0.25 | 0.36 |
| local.chains | 294857 | 0.463 | 1.57 | 1.57 |
| local.scoreOld | 294857 | 1.356 | 4.60 | 4.60 |
| local.after | 294857 | 3.177 | 10.78 | 10.77 |
| local.newChain | 294857 | 0.083 | 0.28 | 0.28 |
| local.sampleNew | 294857 | 2.472 | 8.39 | 8.38 |
| gibbs | 40512 | 0.745 | 18.38 | 2.53 |
| branches | 335369 | 0.632 | 1.88 | 2.14 |
| store.updateBranch | 853191 | 0.352 | 0.41 | 1.19 |
| store.updateBranchSites | 609445 | 0.108 | 0.18 | 0.37 |
| store.journalClone | 416205 | 0.211 | 0.51 | 0.72 |
| store.commit | 852317 | 0.022 | 0.03 | 0.07 |
| store.revert | 357317 | 0.085 | 0.24 | 0.29 |
| lik.updateMAT | 832867 | 0.028 | 0.03 | 0.06 |
| lik.ensureRates | 1463093 | 0.077 | 0.05 | 0.15 |
| lik.calculate | 497498 | 6.391 | 12.85 | 12.85 |
| lik.branchTerms | 497498 | 5.357 | 10.77 | 10.77 |
| lik.sum | 497498 | 0.870 | 1.75 | 1.75 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 294857 |  |
| path_moves | 0 |  |
| gibbs_updates | 40512 |  |
| candidate_sites | 1114120 | 3.78 |
| changed_sites | 7833703 | 26.57 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 2365898643 | 8023.88 |
| branch_site_updates | 11405724 | 38.68 |
| journal_clone_bytes | 1154136465 | 3914.22 |
| runs_recomputed | 418124 | 1.42 |
| terms_recomputed | 12752750 | 43.25 |
| evaluations | 497498 |  |
| transition_cache_misses | 8400705 | 28.49 |
| transition_cache_clears | 1741 | 0.01 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### sim-20k-r02

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 228853 | 2.519 | 11.01 | 11.01 |
| local.runs | 397689 | 0.213 | 0.53 | 0.93 |
| local.chains | 228853 | 0.660 | 2.89 | 2.88 |
| local.scoreOld | 228853 | 1.791 | 7.83 | 7.83 |
| local.after | 228853 | 5.073 | 22.17 | 22.17 |
| local.newChain | 228853 | 0.080 | 0.35 | 0.35 |
| local.sampleNew | 228853 | 4.032 | 17.62 | 17.62 |
| gibbs | 40829 | 1.616 | 39.59 | 7.06 |
| branches | 269682 | 1.033 | 3.83 | 4.51 |
| store.updateBranch | 649079 | 0.531 | 0.82 | 2.32 |
| store.updateBranchSites | 551814 | 0.240 | 0.44 | 1.05 |
| store.journalClone | 351179 | 0.357 | 1.02 | 1.56 |
| store.commit | 852404 | 0.021 | 0.02 | 0.09 |
| store.revert | 357404 | 0.071 | 0.20 | 0.31 |
| lik.updateMAT | 764643 | 0.018 | 0.02 | 0.04 |
| lik.ensureRates | 1263178 | 0.069 | 0.05 | 0.14 |
| lik.calculate | 494961 | 6.434 | 13.00 | 13.00 |
| lik.branchTerms | 494961 | 5.419 | 10.95 | 10.95 |
| lik.sum | 494961 | 0.860 | 1.74 | 1.74 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 228853 |  |
| path_moves | 0 |  |
| gibbs_updates | 40829 |  |
| candidate_sites | 3365189 | 14.70 |
| changed_sites | 28190278 | 123.18 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 5836518368 | 25503.35 |
| branch_site_updates | 42218231 | 184.48 |
| journal_clone_bytes | 3157801568 | 13798.38 |
| runs_recomputed | 397689 | 1.74 |
| terms_recomputed | 12290823 | 53.71 |
| evaluations | 494961 |  |
| transition_cache_misses | 8078508 | 35.30 |
| transition_cache_clears | 1669 | 0.01 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### sim-300k-r02

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 117676 | 10.286 | 87.41 | 87.41 |
| local.runs | 237819 | 0.978 | 4.11 | 8.31 |
| local.chains | 117676 | 5.108 | 43.41 | 43.41 |
| local.scoreOld | 117676 | 5.116 | 43.48 | 43.48 |
| local.after | 117676 | 48.513 | 412.25 | 412.26 |
| local.newChain | 117676 | 0.133 | 1.13 | 1.13 |
| local.sampleNew | 117676 | 42.034 | 357.20 | 357.20 |
| gibbs | 40843 | 13.343 | 326.68 | 113.39 |
| branches | 158519 | 7.571 | 47.76 | 64.34 |
| store.updateBranch | 329780 | 2.950 | 8.94 | 25.07 |
| store.updateBranchSites | 426441 | 2.737 | 6.42 | 23.26 |
| store.journalClone | 240013 | 1.413 | 5.89 | 12.01 |
| store.commit | 852400 | 0.028 | 0.03 | 0.24 |
| store.revert | 357400 | 0.086 | 0.24 | 0.73 |
| lik.updateMAT | 653322 | 0.020 | 0.03 | 0.04 |
| lik.ensureRates | 929517 | 0.159 | 0.17 | 0.32 |
| lik.calculate | 494803 | 11.839 | 23.93 | 23.93 |
| lik.branchTerms | 494803 | 10.704 | 21.63 | 21.63 |
| lik.sum | 494803 | 0.868 | 1.75 | 1.75 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 117676 |  |
| path_moves | 0 |  |
| gibbs_updates | 40843 |  |
| candidate_sites | 26548163 | 225.60 |
| changed_sites | 379722799 | 3226.85 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 44383441300 | 377166.47 |
| branch_site_updates | 586248150 | 4981.88 |
| journal_clone_bytes | 32302149605 | 274500.74 |
| runs_recomputed | 237819 | 2.02 |
| terms_recomputed | 12765190 | 108.48 |
| evaluations | 494803 |  |
| transition_cache_misses | 8301886 | 70.55 |
| transition_cache_clears | 1703 | 0.01 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### small-2-local

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 668535 | 2.860 | 4.28 | 4.28 |
| local.runs | 831221 | 0.061 | 0.07 | 0.09 |
| local.chains | 668535 | 0.415 | 0.62 | 0.62 |
| local.scoreOld | 668535 | 2.287 | 3.42 | 3.42 |
| local.after | 668535 | 4.277 | 6.40 | 6.40 |
| local.newChain | 668535 | 0.094 | 0.14 | 0.14 |
| local.sampleNew | 668535 | 3.501 | 5.24 | 5.24 |
| gibbs | 78778 | 0.918 | 11.66 | 1.37 |
| branches | 747313 | 0.663 | 0.89 | 0.99 |
| store.updateBranch | 1947242 | 0.257 | 0.13 | 0.38 |
| store.updateBranchSites | 1274584 | 0.127 | 0.10 | 0.19 |
| store.journalClone | 902989 | 0.098 | 0.11 | 0.15 |
| store.commit | 1701029 | 0.041 | 0.02 | 0.06 |
| store.revert | 711029 | 0.105 | 0.15 | 0.16 |
| lik.updateMAT | 1717993 | 0.027 | 0.02 | 0.03 |
| lik.ensureRates | 3133841 | 0.051 | 0.02 | 0.05 |
| lik.calculate | 970680 | 6.480 | 6.68 | 6.68 |
| lik.branchTerms | 970680 | 5.949 | 6.13 | 6.13 |
| lik.sum | 970680 | 0.313 | 0.32 | 0.32 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 668535 |  |
| path_moves | 0 |  |
| gibbs_updates | 78778 |  |
| candidate_sites | 1726129 | 2.58 |
| changed_sites | 4065887 | 6.08 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 579806382 | 867.28 |
| branch_site_updates | 8349200 | 12.49 |
| journal_clone_bytes | 269090722 | 402.51 |
| runs_recomputed | 831221 | 1.24 |
| terms_recomputed | 18517400 | 27.70 |
| evaluations | 970680 |  |
| transition_cache_misses | 10366357 | 15.51 |
| transition_cache_clears | 633 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### small-2-path

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.runs | 1214 | 0.001 | 0.79 | 0.00 |
| path.before | 334702 | 6.558 | 19.59 | 19.59 |
| path.candidates | 334702 | 0.472 | 1.41 | 1.41 |
| path.scoreOld | 334702 | 6.009 | 17.95 | 17.95 |
| path.after | 334702 | 7.560 | 22.59 | 22.59 |
| path.sampleNew | 334702 | 6.710 | 20.05 | 20.05 |
| gibbs | 39457 | 0.450 | 11.41 | 1.34 |
| branches | 374159 | 0.820 | 2.19 | 2.45 |
| store.updateBranch | 975141 | 0.124 | 0.13 | 0.37 |
| store.updateBranchSites | 3377290 | 0.349 | 0.10 | 1.04 |
| store.journalClone | 1989136 | 0.142 | 0.07 | 0.42 |
| store.commit | 850553 | 0.035 | 0.04 | 0.10 |
| store.revert | 355553 | 0.106 | 0.30 | 0.32 |
| lik.updateMAT | 859344 | 0.014 | 0.02 | 0.03 |
| lik.ensureRates | 1568205 | 0.024 | 0.02 | 0.05 |
| lik.calculate | 485185 | 3.439 | 7.09 | 7.09 |
| lik.branchTerms | 485185 | 3.161 | 6.51 | 6.52 |
| lik.sum | 485185 | 0.157 | 0.32 | 0.32 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 0 |  |
| path_moves | 334702 |  |
| gibbs_updates | 39457 |  |
| candidate_sites | 4212462 | 12.59 |
| changed_sites | 4531147 | 13.54 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 290387888 | 867.60 |
| branch_site_updates | 53601320 | 160.15 |
| journal_clone_bytes | 592762528 | 1771.02 |
| runs_recomputed | 1214 | 0.00 |
| terms_recomputed | 11618953 | 34.71 |
| evaluations | 485185 |  |
| transition_cache_misses | 5176250 | 15.47 |
| transition_cache_clears | 316 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### wide-1-local

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.before | 405999 | 7.235 | 17.82 | 17.82 |
| local.runs | 844549 | 1.776 | 2.10 | 4.37 |
| local.chains | 405999 | 1.383 | 3.41 | 3.41 |
| local.scoreOld | 405999 | 5.751 | 14.17 | 14.17 |
| local.after | 405999 | 27.560 | 67.88 | 67.88 |
| local.newChain | 405999 | 0.697 | 1.72 | 1.72 |
| local.sampleNew | 405999 | 17.117 | 42.16 | 42.16 |
| gibbs | 79283 | 15.819 | 199.53 | 38.96 |
| branches | 485282 | 12.389 | 25.53 | 30.51 |
| store.updateBranch | 1159648 | 4.944 | 4.26 | 12.18 |
| store.updateBranchSites | 1620490 | 6.113 | 3.77 | 15.06 |
| store.journalClone | 975205 | 1.412 | 1.45 | 3.48 |
| store.commit | 1706756 | 0.041 | 0.02 | 0.10 |
| store.revert | 716756 | 0.091 | 0.13 | 0.22 |
| lik.updateMAT | 1438263 | 0.024 | 0.02 | 0.03 |
| lik.ensureRates | 2589405 | 0.046 | 0.02 | 0.05 |
| lik.calculate | 952981 | 2.551 | 2.68 | 2.68 |
| lik.branchTerms | 952981 | 2.248 | 2.36 | 2.36 |
| lik.sum | 952981 | 0.125 | 0.13 | 0.13 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 405999 |  |
| path_moves | 0 |  |
| gibbs_updates | 79283 |  |
| candidate_sites | 28175102 | 69.40 |
| changed_sites | 227460776 | 560.25 |
| hot_sites | 29080277 | 71.63 |
| branch_scan_patterns | 36162032445 | 89069.26 |
| branch_site_updates | 1152741934 | 2839.27 |
| journal_clone_bytes | 30411767925 | 74906.02 |
| runs_recomputed | 844549 | 2.08 |
| terms_recomputed | 7079504 | 17.44 |
| evaluations | 952981 |  |
| transition_cache_misses | 2398977 | 5.91 |
| transition_cache_clears | 146 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

### wide-1-path

| Phase | Calls | Total (s) | us/call | us per move or evaluation |
|---|---|---|---|---|
| local.runs | 22815 | 0.066 | 2.90 | 0.33 |
| path.before | 202072 | 45.491 | 225.12 | 225.12 |
| path.candidates | 202072 | 3.137 | 15.52 | 15.52 |
| path.scoreOld | 202072 | 42.279 | 209.23 | 209.23 |
| path.after | 202072 | 65.756 | 325.41 | 325.41 |
| path.sampleNew | 202072 | 54.702 | 270.71 | 270.71 |
| gibbs | 39717 | 8.547 | 215.20 | 42.30 |
| branches | 241789 | 12.159 | 50.29 | 60.17 |
| store.updateBranch | 576448 | 2.679 | 4.65 | 13.26 |
| store.updateBranchSites | 1833153 | 8.461 | 4.62 | 41.87 |
| store.journalClone | 1083906 | 1.771 | 1.63 | 8.76 |
| store.commit | 853219 | 0.033 | 0.04 | 0.16 |
| store.revert | 358219 | 0.085 | 0.24 | 0.42 |
| lik.updateMAT | 718108 | 0.015 | 0.02 | 0.03 |
| lik.ensureRates | 1161969 | 0.026 | 0.02 | 0.05 |
| lik.calculate | 476319 | 1.952 | 4.10 | 4.10 |
| lik.branchTerms | 476319 | 1.763 | 3.70 | 3.70 |
| lik.sum | 476319 | 0.068 | 0.14 | 0.14 |

| Counter | Total | per move |
|---|---|---|
| local_moves | 0 |  |
| path_moves | 202072 |  |
| gibbs_updates | 39717 |  |
| candidate_sites | 291938586 | 1444.73 |
| changed_sites | 166843974 | 825.67 |
| hot_sites | 0 | 0.00 |
| branch_scan_patterns | 17975408220 | 88955.46 |
| branch_site_updates | 2377596672 | 11766.09 |
| journal_clone_bytes | 33801608610 | 167275.07 |
| runs_recomputed | 22815 | 0.11 |
| terms_recomputed | 5341465 | 26.43 |
| evaluations | 476319 |  |
| transition_cache_misses | 1200862 | 5.94 |
| transition_cache_clears | 73 | 0.00 |
| parsimony_full | 0 | 0.00 |
| parsimony_incremental | 0 | 0.00 |
| category_changes | 0 | 0.00 |

## MAT phases per step of each tree operator (microseconds per step)

### big-2-local

| Phase | subtreeLeap (820524 steps) | subtreeJump (82527 steps) | mutationHistoryGibbs (82214 steps) |
|---|---|---|---|
| local.before | 11.19 | 17.47 |  |
| local.runs | 1.44 | 2.80 | 0.66 |
| local.chains | 2.11 | 3.78 |  |
| local.scoreOld | 8.92 | 13.44 |  |
| local.after | 17.94 | 56.79 |  |
| local.newChain | 1.37 | 4.63 |  |
| local.sampleNew | 11.23 | 34.48 |  |
| branches | 5.17 | 17.29 | 2.90 |
| store.updateBranch | 3.28 | 4.89 |  |
| store.updateBranchSites | 1.32 | 10.79 | 2.08 |
| store.journalClone | 1.89 | 6.48 | 4.25 |
| store.commit | 0.07 | 0.16 | 0.03 |
| store.revert | 0.21 | 0.47 |  |
| lik.updateMAT | 0.03 | 0.04 | 0.04 |
| lik.ensureRates | 0.06 | 0.10 | 0.05 |
| lik.calculate | 7.27 | 8.14 | 6.69 |
| lik.branchTerms | 3.99 | 4.81 | 3.37 |
| lik.sum | 3.06 | 3.06 | 3.06 |
| gibbs |  |  | 39.24 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| local_moves | 0.77 | 1.00 |  |
| candidate_sites | 8.39 | 12.67 |  |
| changed_sites | 32.50 | 137.91 | 6.81 |
| hot_sites | 9.80 | 61.69 |  |
| branch_scan_patterns | 40757.84 | 54545.10 |  |
| branch_site_updates | 97.29 | 1769.68 | 47.68 |
| journal_clone_bytes | 26337.47 | 90024.68 | 54617.79 |
| runs_recomputed | 1.07 | 2.07 | 0.45 |
| terms_recomputed | 15.46 | 22.52 | 15.51 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 5.76 | 5.29 | 2.43 |
| transition_cache_clears | 0.00 | 0.00 |  |
| gibbs_updates |  |  | 1.00 |

### big-2-path

| Phase | subtreeLeap (410591 steps) | subtreeJump (40978 steps) | mutationHistoryGibbs (41157 steps) |
|---|---|---|---|
| local.runs | 0.23 | 0.67 | 0.24 |
| path.before | 54.28 | 507.83 |  |
| path.candidates | 7.01 | 24.21 |  |
| path.scoreOld | 47.09 | 483.32 |  |
| path.after | 62.22 | 603.80 |  |
| path.sampleNew | 53.49 | 537.76 |  |
| branches | 8.44 | 65.23 | 3.16 |
| store.updateBranch | 3.15 | 4.83 |  |
| store.updateBranchSites | 4.22 | 55.67 | 2.28 |
| store.journalClone | 6.25 | 32.25 | 4.15 |
| store.commit | 0.14 | 0.51 | 0.04 |
| store.revert | 0.38 | 1.52 |  |
| lik.updateMAT | 0.05 | 0.06 | 0.05 |
| lik.ensureRates | 0.06 | 0.07 | 0.05 |
| lik.calculate | 8.41 | 13.80 | 7.35 |
| lik.branchTerms | 5.07 | 10.22 | 4.01 |
| lik.sum | 3.07 | 3.13 | 3.06 |
| gibbs |  |  | 41.68 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| path_moves | 0.77 | 1.00 |  |
| candidate_sites | 380.59 | 1306.97 |  |
| changed_sites | 41.71 | 167.37 | 7.89 |
| branch_scan_patterns | 41085.35 | 54554.50 |  |
| branch_site_updates | 473.28 | 10654.65 | 55.23 |
| journal_clone_bytes | 85681.36 | 471032.28 | 54628.03 |
| runs_recomputed | 0.13 | 0.39 | 0.13 |
| terms_recomputed | 24.33 | 55.22 | 20.20 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 5.78 | 6.25 | 1.90 |
| transition_cache_clears | 0.00 | 0.00 |  |
| gibbs_updates |  |  | 1.00 |

### big-3-gamma-local

| Phase | subtreeLeap (808021 steps) | mutationHistoryGibbs (80951 steps) | subtreeJump (80494 steps) |
|---|---|---|---|
| local.before | 16.58 |  | 23.73 |
| local.runs | 0.35 | 0.15 | 0.61 |
| local.chains | 1.50 |  | 2.45 |
| local.scoreOld | 14.82 |  | 20.91 |
| local.after | 25.83 |  | 49.45 |
| local.newChain | 1.42 |  | 2.75 |
| local.sampleNew | 18.48 |  | 34.95 |
| branches | 5.57 | 5.30 | 11.17 |
| store.updateBranch | 1.77 |  | 2.59 |
| store.updateBranchSites | 2.77 | 4.55 | 6.84 |
| store.journalClone | 0.39 | 1.29 | 0.73 |
| store.commit | 0.07 | 0.03 | 0.11 |
| store.revert | 0.25 |  | 0.45 |
| lik.updateMAT | 0.04 | 0.05 | 0.05 |
| lik.ensureRates | 0.15 | 0.14 | 0.17 |
| lik.calculate | 28.17 | 24.98 | 27.99 |
| lik.branchTerms | 25.96 | 22.76 | 25.74 |
| lik.sum | 1.94 | 1.93 | 1.94 |
| gibbs |  | 58.55 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| local_moves | 0.78 |  | 1.00 |
| candidate_sites | 51.65 |  | 72.35 |
| changed_sites | 201.22 | 5.08 | 319.45 |
| hot_sites | 73.14 |  | 159.44 |
| branch_scan_patterns | 5463.23 |  | 7210.48 |
| branch_site_updates | 522.41 | 35.47 | 1414.10 |
| journal_clone_bytes | 2815.54 | 7216.76 | 5418.20 |
| runs_recomputed | 0.93 | 0.40 | 1.76 |
| terms_recomputed | 19.58 | 19.37 | 21.52 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 44.59 | 30.39 | 41.55 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  | 1.00 |  |

### big-3-local

| Phase | subtreeJump (81251 steps) | mutationHistoryGibbs (81183 steps) | subtreeLeap (817454 steps) |
|---|---|---|---|
| local.before | 14.76 |  | 10.02 |
| local.runs | 0.59 | 0.15 | 0.34 |
| local.chains | 2.37 |  | 1.41 |
| local.scoreOld | 12.09 |  | 8.41 |
| local.after | 36.90 |  | 16.97 |
| local.newChain | 2.39 |  | 1.20 |
| local.sampleNew | 22.86 |  | 9.95 |
| branches | 11.13 | 5.13 | 5.51 |
| store.updateBranch | 2.61 |  | 1.76 |
| store.updateBranchSites | 6.90 | 4.49 | 2.78 |
| store.journalClone | 0.86 | 1.31 | 0.47 |
| store.commit | 0.10 | 0.03 | 0.06 |
| store.revert | 0.35 |  | 0.21 |
| lik.updateMAT | 0.04 | 0.04 | 0.03 |
| lik.ensureRates | 0.12 | 0.09 | 0.09 |
| lik.calculate | 7.04 | 5.77 | 6.72 |
| lik.branchTerms | 5.01 | 3.77 | 4.75 |
| lik.sum | 1.74 | 1.73 | 1.73 |
| gibbs |  | 39.19 |  |

| Counter per step | subtreeJump | mutationHistoryGibbs | subtreeLeap |
|---|---|---|---|
| local_moves | 1.00 |  | 0.78 |
| candidate_sites | 72.34 |  | 51.54 |
| changed_sites | 316.14 | 4.94 | 197.16 |
| hot_sites | 159.61 |  | 72.35 |
| branch_scan_patterns | 7211.79 |  | 5461.04 |
| branch_site_updates | 1422.46 | 34.53 | 516.44 |
| journal_clone_bytes | 5409.58 | 7215.93 | 2808.63 |
| runs_recomputed | 1.76 | 0.39 | 0.93 |
| terms_recomputed | 18.27 | 14.83 | 16.32 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 6.04 | 2.60 | 6.50 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  | 1.00 |  |

### big-3-path

| Phase | subtreeLeap (408479 steps) | subtreeJump (40897 steps) | mutationHistoryGibbs (40700 steps) |
|---|---|---|---|
| local.runs | 0.04 | 0.13 | 0.05 |
| path.before | 43.10 | 246.45 |  |
| path.candidates | 10.20 | 25.04 |  |
| path.scoreOld | 32.67 | 221.02 |  |
| path.after | 59.70 | 346.63 |  |
| path.sampleNew | 44.66 | 274.84 |  |
| branches | 14.59 | 70.84 | 4.99 |
| store.updateBranch | 1.54 | 2.55 |  |
| store.updateBranchSites | 11.56 | 63.99 | 4.31 |
| store.journalClone | 1.41 | 5.99 | 1.15 |
| store.commit | 0.12 | 0.36 | 0.03 |
| store.revert | 0.38 | 1.47 |  |
| lik.updateMAT | 0.03 | 0.05 | 0.05 |
| lik.ensureRates | 0.09 | 0.10 | 0.10 |
| lik.calculate | 7.27 | 10.72 | 6.35 |
| lik.branchTerms | 5.23 | 8.47 | 4.35 |
| lik.sum | 1.74 | 1.77 | 1.73 |
| gibbs |  |  | 40.45 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| path_moves | 0.79 | 1.00 |  |
| candidate_sites | 826.09 | 1415.15 |  |
| changed_sites | 203.85 | 343.62 | 5.33 |
| branch_scan_patterns | 5529.52 | 7211.15 |  |
| branch_site_updates | 2204.66 | 14521.86 | 37.22 |
| journal_clone_bytes | 11090.89 | 51275.35 | 7216.35 |
| runs_recomputed | 0.09 | 0.29 | 0.09 |
| terms_recomputed | 23.60 | 42.73 | 19.69 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 6.41 | 5.45 | 2.77 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  |  | 1.00 |

### sim-100k-r005

| Phase | subtreeLeap (408202 steps) | mutationHistoryGibbs (40895 steps) | subtreeJump (40939 steps) |
|---|---|---|---|
| local.before | 4.83 |  | 20.66 |
| local.runs | 0.51 |  | 1.34 |
| local.chains | 1.03 |  | 9.85 |
| local.scoreOld | 3.68 |  | 10.46 |
| local.after | 6.93 |  | 80.51 |
| local.newChain | 0.15 |  | 0.59 |
| local.sampleNew | 5.02 |  | 69.64 |
| branches | 1.66 | 4.22 | 9.30 |
| store.updateBranch | 1.26 |  | 3.36 |
| store.updateBranchSites | 0.14 | 2.83 | 2.73 |
| store.journalClone | 0.56 | 4.39 | 1.35 |
| store.commit | 0.04 | 0.03 | 0.05 |
| store.revert | 0.16 |  | 0.36 |
| lik.updateMAT | 0.08 | 0.07 | 0.11 |
| lik.ensureRates | 0.08 | 0.10 | 0.14 |
| lik.calculate | 7.69 | 6.67 | 9.02 |
| lik.branchTerms | 5.61 | 4.59 | 6.89 |
| lik.sum | 1.76 | 1.77 | 1.78 |
| gibbs |  | 49.82 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| local_moves | 0.41 |  | 1.00 |
| candidate_sites | 6.72 |  | 25.28 |
| changed_sites | 8.32 | 69.01 | 663.84 |
| branch_scan_patterns | 15844.97 |  | 41584.85 |
| branch_site_updates | 9.93 | 478.57 | 665.15 |
| journal_clone_bytes | 5668.72 | 41675.85 | 13911.00 |
| runs_recomputed | 0.69 |  | 1.94 |
| terms_recomputed | 14.13 | 17.08 | 15.32 |
| evaluations | 1.02 | 1.01 | 1.02 |
| transition_cache_misses | 6.08 | 2.61 | 6.38 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  | 1.00 |  |

### sim-100k-r02

| Phase | subtreeLeap (408624 steps) | mutationHistoryGibbs (40677 steps) | subtreeJump (40814 steps) |
|---|---|---|---|
| local.before | 7.06 |  | 56.15 |
| local.runs | 0.87 |  | 3.40 |
| local.chains | 1.57 |  | 32.41 |
| local.scoreOld | 5.39 |  | 23.32 |
| local.after | 10.07 |  | 291.12 |
| local.newChain | 0.14 |  | 0.96 |
| local.sampleNew | 7.02 |  | 257.15 |
| branches | 2.79 | 13.17 | 29.70 |
| store.updateBranch | 2.14 |  | 9.79 |
| store.updateBranchSites | 0.25 | 9.92 | 9.53 |
| store.journalClone | 0.73 | 8.88 | 2.63 |
| store.commit | 0.05 | 0.03 | 0.06 |
| store.revert | 0.14 |  | 0.46 |
| lik.updateMAT | 0.05 | 0.06 | 0.08 |
| lik.ensureRates | 0.10 | 0.13 | 0.15 |
| lik.calculate | 8.61 | 8.19 | 12.54 |
| lik.branchTerms | 6.55 | 6.09 | 10.41 |
| lik.sum | 1.75 | 1.76 | 1.78 |
| gibbs |  | 129.19 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| local_moves | 0.28 |  | 1.00 |
| candidate_sites | 17.72 |  | 100.80 |
| changed_sites | 19.38 | 274.45 | 2641.58 |
| branch_scan_patterns | 34233.27 |  | 134212.03 |
| branch_site_updates | 24.18 | 1899.15 | 2647.80 |
| journal_clone_bytes | 12532.34 | 134412.55 | 44870.00 |
| runs_recomputed | 0.52 |  | 2.05 |
| terms_recomputed | 13.93 | 17.88 | 14.22 |
| evaluations | 1.01 | 1.00 | 1.01 |
| transition_cache_misses | 5.93 | 2.94 | 5.46 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  | 1.00 |  |

### sim-20k-r005

| Phase | subtreeJump (40770 steps) | subtreeLeap (408769 steps) | mutationHistoryGibbs (40512 steps) |
|---|---|---|---|
| local.before | 8.43 | 3.83 |  |
| local.runs | 0.42 | 0.22 |  |
| local.chains | 3.19 | 0.81 |  |
| local.scoreOld | 4.91 | 2.83 |  |
| local.after | 23.25 | 5.45 |  |
| local.newChain | 0.36 | 0.17 |  |
| local.sampleNew | 19.67 | 4.09 |  |
| branches | 2.86 | 1.09 | 1.74 |
| store.updateBranch | 1.07 | 0.75 |  |
| store.updateBranchSites | 0.75 | 0.09 | 0.96 |
| store.journalClone | 0.49 | 0.33 | 1.42 |
| store.commit | 0.05 | 0.05 | 0.03 |
| store.revert | 0.31 | 0.18 |  |
| lik.updateMAT | 0.07 | 0.06 | 0.05 |
| lik.ensureRates | 0.11 | 0.10 | 0.09 |
| lik.calculate | 6.97 | 6.68 | 5.75 |
| lik.branchTerms | 4.92 | 4.65 | 3.73 |
| lik.sum | 1.77 | 1.76 | 1.75 |
| gibbs |  |  | 18.38 |

| Counter per step | subtreeJump | subtreeLeap | mutationHistoryGibbs |
|---|---|---|---|
| local_moves | 1.00 | 0.62 |  |
| candidate_sites | 5.01 | 2.23 |  |
| changed_sites | 143.50 | 3.44 | 14.24 |
| branch_scan_patterns | 8289.55 | 4961.08 |  |
| branch_site_updates | 143.83 | 3.80 | 98.41 |
| journal_clone_bytes | 2773.00 | 1723.67 | 8306.13 |
| runs_recomputed | 1.75 | 0.85 |  |
| terms_recomputed | 15.57 | 14.51 | 14.74 |
| evaluations | 1.01 | 1.01 | 1.00 |
| transition_cache_misses | 5.99 | 6.30 | 2.78 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  |  | 1.00 |

### sim-20k-r02

| Phase | subtreeLeap (408314 steps) | subtreeJump (40985 steps) | mutationHistoryGibbs (40829 steps) |
|---|---|---|---|
| local.before | 4.54 | 16.21 |  |
| local.runs | 0.42 | 0.98 |  |
| local.chains | 0.89 | 7.21 |  |
| local.scoreOld | 3.52 | 8.68 |  |
| local.after | 6.51 | 58.92 |  |
| local.newChain | 0.15 | 0.49 |  |
| local.sampleNew | 4.81 | 50.42 |  |
| branches | 1.45 | 7.27 | 3.47 |
| store.updateBranch | 1.06 | 2.39 |  |
| store.updateBranchSites | 0.13 | 2.30 | 2.29 |
| store.journalClone | 0.46 | 0.97 | 3.21 |
| store.commit | 0.04 | 0.05 | 0.02 |
| store.revert | 0.14 | 0.30 |  |
| lik.updateMAT | 0.04 | 0.05 | 0.04 |
| lik.ensureRates | 0.08 | 0.11 | 0.09 |
| lik.calculate | 6.78 | 7.26 | 5.75 |
| lik.branchTerms | 4.79 | 5.25 | 3.75 |
| lik.sum | 1.74 | 1.75 | 1.74 |
| gibbs |  |  | 39.59 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| local_moves | 0.46 | 1.00 |  |
| candidate_sites | 6.21 | 20.22 |  |
| changed_sites | 7.97 | 553.21 | 55.42 |
| branch_scan_patterns | 11594.75 | 26893.29 |  |
| branch_site_updates | 9.37 | 554.47 | 383.75 |
| journal_clone_bytes | 4137.28 | 8992.00 | 26940.54 |
| runs_recomputed | 0.78 | 1.94 |  |
| terms_recomputed | 14.03 | 13.96 | 15.53 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 6.12 | 5.37 | 2.58 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  |  | 1.00 |

### sim-300k-r02

| Phase | subtreeLeap (407888 steps) | subtreeJump (41324 steps) | mutationHistoryGibbs (40843 steps) |
|---|---|---|---|
| local.before | 10.15 | 148.69 |  |
| local.runs | 1.50 | 8.84 |  |
| local.chains | 2.67 | 97.24 |  |
| local.scoreOld | 7.39 | 50.89 |  |
| local.after | 15.70 | 1018.96 |  |
| local.newChain | 0.15 | 1.76 |  |
| local.sampleNew | 10.01 | 918.39 |  |
| branches | 5.36 | 87.67 | 43.12 |
| store.updateBranch | 4.24 | 29.50 |  |
| store.updateBranchSites | 0.47 | 27.16 | 34.79 |
| store.journalClone | 1.02 | 5.72 | 18.64 |
| store.commit | 0.05 | 0.09 | 0.03 |
| store.revert | 0.14 | 0.67 |  |
| lik.updateMAT | 0.04 | 0.06 | 0.06 |
| lik.ensureRates | 0.15 | 0.23 | 0.19 |
| lik.calculate | 11.78 | 26.38 | 12.09 |
| lik.branchTerms | 9.68 | 24.05 | 9.92 |
| lik.sum | 1.75 | 1.82 | 1.78 |
| gibbs |  |  | 326.68 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| local_moves | 0.19 | 1.00 |  |
| candidate_sites | 34.54 | 301.53 |  |
| changed_sites | 35.70 | 8012.09 | 834.19 |
| branch_scan_patterns | 68022.69 | 402618.37 |  |
| branch_site_updates | 45.64 | 8030.68 | 5772.61 |
| journal_clone_bytes | 25192.78 | 134585.00 | 403122.33 |
| runs_recomputed | 0.37 | 2.14 |  |
| terms_recomputed | 14.16 | 14.72 | 17.29 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 5.93 | 5.63 | 2.50 |
| transition_cache_clears | 0.00 | 0.00 | 0.00 |
| gibbs_updates |  |  | 1.00 |

### small-2-local

| Phase | subtreeLeap (789540 steps) | mutationHistoryGibbs (78778 steps) | subtreeJump (79055 steps) |
|---|---|---|---|
| local.before | 3.16 |  | 4.62 |
| local.runs | 0.07 | 0.00 | 0.11 |
| local.chains | 0.45 |  | 0.75 |
| local.scoreOld | 2.54 |  | 3.58 |
| local.after | 4.66 |  | 7.57 |
| local.newChain | 0.10 |  | 0.16 |
| local.sampleNew | 3.82 |  | 6.14 |
| branches | 0.61 | 1.28 | 1.06 |
| store.updateBranch | 0.28 |  | 0.45 |
| store.updateBranchSites | 0.08 | 0.64 | 0.18 |
| store.journalClone | 0.08 | 0.29 | 0.11 |
| store.commit | 0.04 | 0.02 | 0.05 |
| store.revert | 0.11 |  | 0.14 |
| lik.updateMAT | 0.03 | 0.03 | 0.03 |
| lik.ensureRates | 0.05 | 0.05 | 0.07 |
| lik.calculate | 2.75 | 2.25 | 2.63 |
| lik.branchTerms | 2.22 | 1.64 | 2.07 |
| lik.sum | 0.32 | 0.33 | 0.32 |
| gibbs |  | 11.66 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| local_moves | 0.75 |  | 1.00 |
| candidate_sites | 1.86 |  | 3.25 |
| changed_sites | 3.00 | 8.80 | 12.73 |
| branch_scan_patterns | 646.45 |  | 877.96 |
| branch_site_updates | 3.25 | 60.21 | 13.16 |
| journal_clone_bytes | 222.49 | 886.89 | 298.00 |
| runs_recomputed | 0.90 | 0.01 | 1.56 |
| terms_recomputed | 11.73 | 11.37 | 12.11 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 3.72 | 0.23 | 2.85 |
| transition_cache_clears | 0.00 |  | 0.00 |
| gibbs_updates |  | 1.00 |  |

### small-2-path

| Phase | subtreeLeap (394920 steps) | mutationHistoryGibbs (39457 steps) | subtreeJump (39226 steps) |
|---|---|---|---|
| local.runs | 0.00 | 0.00 | 0.01 |
| path.before | 12.16 |  | 44.79 |
| path.candidates | 0.98 |  | 2.20 |
| path.scoreOld | 11.02 |  | 42.28 |
| path.after | 14.05 |  | 51.30 |
| path.sampleNew | 12.39 |  | 46.33 |
| branches | 1.49 | 1.26 | 4.65 |
| store.updateBranch | 0.27 |  | 0.48 |
| store.updateBranchSites | 0.57 | 0.58 | 2.59 |
| store.journalClone | 0.27 | 0.27 | 0.67 |
| store.commit | 0.07 | 0.03 | 0.12 |
| store.revert | 0.22 |  | 0.48 |
| lik.updateMAT | 0.03 | 0.04 | 0.04 |
| lik.ensureRates | 0.05 | 0.05 | 0.07 |
| lik.calculate | 2.94 | 2.45 | 3.38 |
| lik.branchTerms | 2.38 | 1.86 | 2.72 |
| lik.sum | 0.32 | 0.33 | 0.35 |
| gibbs |  | 11.41 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| path_moves | 0.75 |  | 1.00 |
| candidate_sites | 7.84 |  | 28.44 |
| changed_sites | 7.78 | 8.78 | 28.36 |
| branch_scan_patterns | 648.09 |  | 878.14 |
| branch_site_updates | 73.85 | 60.00 | 562.64 |
| journal_clone_bytes | 1129.02 | 886.70 | 2852.82 |
| runs_recomputed | 0.00 | 0.00 | 0.01 |
| terms_recomputed | 16.43 | 13.93 | 22.56 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 3.72 | 0.16 | 2.90 |
| transition_cache_clears | 0.00 |  | 0.00 |
| gibbs_updates |  | 1.00 |  |

### wide-1-local

| Phase | subtreeLeap (793810 steps) | subtreeJump (79863 steps) | mutationHistoryGibbs (79283 steps) |
|---|---|---|---|
| local.before | 6.11 | 29.89 |  |
| local.runs | 1.60 | 5.15 | 1.24 |
| local.chains | 0.95 | 7.91 |  |
| local.scoreOld | 5.06 | 21.70 |  |
| local.after | 13.53 | 210.65 |  |
| local.newChain | 0.49 | 3.91 |  |
| local.sampleNew | 7.95 | 135.27 |  |
| branches | 4.98 | 69.76 | 36.17 |
| store.updateBranch | 3.75 | 24.63 |  |
| store.updateBranchSites | 0.75 | 37.66 | 31.70 |
| store.journalClone | 0.88 | 4.55 | 4.47 |
| store.commit | 0.04 | 0.07 | 0.02 |
| store.revert | 0.09 | 0.22 |  |
| lik.updateMAT | 0.02 | 0.04 | 0.04 |
| lik.ensureRates | 0.04 | 0.10 | 0.06 |
| lik.calculate | 2.52 | 4.20 | 2.76 |
| lik.branchTerms | 2.21 | 3.84 | 2.39 |
| lik.sum | 0.13 | 0.13 | 0.15 |
| gibbs |  |  | 199.53 |

| Counter per step | subtreeLeap | subtreeJump | mutationHistoryGibbs |
|---|---|---|---|
| local_moves | 0.41 | 1.00 |  |
| candidate_sites | 16.72 | 186.63 |  |
| changed_sites | 32.10 | 1667.33 | 868.09 |
| hot_sites | 3.29 | 331.47 |  |
| branch_scan_patterns | 36325.16 | 91741.61 |  |
| branch_site_updates | 69.13 | 8296.23 | 5490.54 |
| journal_clone_bytes | 19214.40 | 98423.91 | 92059.53 |
| runs_recomputed | 0.78 | 2.25 | 0.57 |
| terms_recomputed | 6.84 | 11.27 | 9.50 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 2.85 | 1.67 | 0.03 |
| transition_cache_clears | 0.00 | 0.00 |  |
| gibbs_updates |  |  | 1.00 |

### wide-1-path

| Phase | subtreeLeap (396735 steps) | mutationHistoryGibbs (39717 steps) | subtreeJump (39839 steps) |
|---|---|---|---|
| local.runs | 0.06 | 0.36 | 0.74 |
| path.before | 27.17 |  | 871.29 |
| path.candidates | 2.82 |  | 50.61 |
| path.scoreOld | 24.22 |  | 820.10 |
| path.after | 36.00 |  | 1292.01 |
| path.sampleNew | 28.42 |  | 1090.03 |
| branches | 7.29 | 37.40 | 195.36 |
| store.updateBranch | 4.05 |  | 26.86 |
| store.updateBranchSites | 2.35 | 32.72 | 156.32 |
| store.journalClone | 2.29 | 4.91 | 16.74 |
| store.commit | 0.06 | 0.03 | 0.17 |
| store.revert | 0.14 |  | 0.65 |
| lik.updateMAT | 0.03 | 0.05 | 0.06 |
| lik.ensureRates | 0.05 | 0.07 | 0.11 |
| lik.calculate | 3.53 | 3.70 | 10.14 |
| lik.branchTerms | 3.17 | 3.24 | 9.48 |
| lik.sum | 0.13 | 0.17 | 0.21 |
| gibbs |  | 215.20 |  |

| Counter per step | subtreeLeap | mutationHistoryGibbs | subtreeJump |
|---|---|---|---|
| path_moves | 0.41 |  | 1.00 |
| candidate_sites | 178.51 |  | 5550.28 |
| changed_sites | 64.33 | 876.36 | 2673.61 |
| branch_scan_patterns | 36096.18 |  | 91738.96 |
| branch_site_updates | 416.25 | 5523.40 | 50028.40 |
| journal_clone_bytes | 43774.31 | 92036.46 | 320775.99 |
| runs_recomputed | 0.02 | 0.12 | 0.24 |
| terms_recomputed | 9.85 | 11.33 | 24.71 |
| evaluations | 1.00 | 1.00 | 1.00 |
| transition_cache_misses | 2.85 | 0.02 | 1.70 |
| transition_cache_clears | 0.00 |  | 0.00 |
| gibbs_updates |  | 1.00 |  |

## Overhead of the profiles

| Run | min/M (profiled) | min/M (unprofiled) | Overhead % | STL us (profiled) | STL us (unprofiled) |
|---|---|---|---|---|---|
| big-2-local | 1.48 | 1.60 | -7.5 | 80.47 | 87.76 |
| big-3-local | 1.18 | 1.16 | 1.7 | 62.00 | 60.98 |
| sim-100k-r02 | 1.68 | 1.54 | 9.1 | 53.75 | 49.84 |
| small-2-local | 0.37 | 0.33 | 12.1 | 17.89 | 16.28 |
