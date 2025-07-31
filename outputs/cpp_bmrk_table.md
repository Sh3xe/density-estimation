### bfgs
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|1.00|0.00|5.06e-14|7.68e-14|8.26e-27|2.47e-26|2\textmu s|4\textmu s|
|sphere_10d|1.00|0.00|6.60e-13|3.26e-13|5.38e-25|5.19e-25|6\textmu s|3\textmu s|
|sphere_30d|1.00|0.00|3.93e-12|1.95e-12|1.91e-23|1.86e-23|20\textmu s|3\textmu s|
|schwefel_2d|8.33|1.32|6.17e+02|1.14e+02|6.52e+02|8.60e+01|18\textmu s|9\textmu s|
|schwefel_10d|20.43|3.60|1.38e+03|1.13e+02|3.27e+03|2.89e+02|314\textmu s|96\textmu s|
|schwefel_30d|24.23|3.11|2.38e+03|1.09e+02|1.02e+04|3.68e+02|1ms|444\textmu s|
|rastrigin_2d|59.43|136.63|2.87e+00|1.05e+00|9.35e+00|6.34e+00|112\textmu s|275\textmu s|
|rastrigin_10d|72.83|121.15|7.40e+00|1.89e+00|5.85e+01|2.86e+01|1ms|2ms|
|rastrigin_30d|92.13|72.97|1.31e+01|2.24e+00|1.78e+02|5.75e+01|6ms|7ms|
|schaffer_f6|3.87|0.51|6.49e+00|3.59e+00|4.74e-02|4.18e-02|2\textmu s|1\textmu s|
|rosenbrock|297.13|235.97|4.25e-04|3.49e-05|3.64e-08|5.36e-09|199\textmu s|169\textmu s|
### lbfgs30
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|1.00|0.00|5.06e-14|7.68e-14|1.91e+00|1.93e+00|1\textmu s|4\textmu s|
|sphere_10d|1.00|0.00|6.60e-13|3.26e-13|1.01e+01|4.21e+00|2\textmu s|0\textmu s|
|sphere_30d|1.00|0.00|3.93e-12|1.95e-12|3.02e+01|8.94e+00|6\textmu s|0\textmu s|
|schwefel_2d|7.13|1.74|6.08e+02|1.05e+02|6.50e+02|8.36e+01|7\textmu s|4\textmu s|
|schwefel_10d|12.40|2.75|1.38e+03|1.12e+02|3.29e+03|2.90e+02|66\textmu s|17\textmu s|
|schwefel_30d|14.93|2.36|2.38e+03|1.13e+02|1.02e+04|3.69e+02|751\textmu s|211\textmu s|
|rastrigin_2d|17.57|58.66|2.72e+00|8.31e-01|8.09e+00|4.49e+00|70\textmu s|313\textmu s|
|rastrigin_10d|23.83|37.35|6.55e+00|1.51e+00|4.53e+01|2.06e+01|325\textmu s|764\textmu s|
|rastrigin_30d|19.13|22.24|1.10e+01|1.90e+00|1.24e+02|4.46e+01|1ms|2ms|
|schaffer_f6|3.87|0.51|6.49e+00|3.59e+00|4.74e-02|4.18e-02|3\textmu s|1\textmu s|
|rosenbrock|196.47|215.30|4.39e-04|1.90e-05|4.04e-08|4.88e-09|341\textmu s|408\textmu s|
### cg_fr
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|1.00|0.00|5.06e-14|7.68e-14|8.26e-27|2.47e-26|0\textmu s|0\textmu s|
|sphere_10d|1.00|0.00|6.60e-13|3.26e-13|5.38e-25|5.19e-25|1\textmu s|0\textmu s|
|sphere_30d|1.00|0.00|3.93e-12|1.95e-12|1.91e-23|1.86e-23|5\textmu s|1\textmu s|
|schwefel_2d|488.23|64.45|1.24e+13|6.81e+13|9.29e+25|5.09e+26|469\textmu s|119\textmu s|
|schwefel_10d|488.07|65.36|1.35e+13|7.39e+13|1.09e+26|5.98e+26|5ms|1ms|
|schwefel_30d|494.67|29.21|5.74e+11|3.15e+12|1.98e+23|1.08e+24|53ms|16ms|
|rastrigin_2d|500.00|0.00|8.01e+03|3.37e+04|1.16e+10|6.16e+10|590\textmu s|86\textmu s|
|rastrigin_10d|500.00|0.00|4.21e+01|8.70e+01|9.04e+04|3.02e+05|9ms|2ms|
|rastrigin_30d|500.00|0.00|3.23e+02|1.51e+03|2.31e+07|1.25e+08|58ms|13ms|
|schaffer_f6|8.50|3.07|1.30e+06|3.72e+06|5.00e-01|9.58e-07|12\textmu s|4\textmu s|
|rosenbrock|29.53|26.42|2.51e+39|1.03e+40|2.73e+163|inf|16\textmu s|14\textmu s|
### cg_pr
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|1.00|0.00|5.06e-14|7.68e-14|8.26e-27|2.47e-26|0\textmu s|0\textmu s|
|sphere_10d|1.00|0.00|6.60e-13|3.26e-13|5.38e-25|5.19e-25|1\textmu s|0\textmu s|
|sphere_30d|1.00|0.00|3.93e-12|1.95e-12|1.91e-23|1.86e-23|5\textmu s|0\textmu s|
|schwefel_2d|43.73|3.87|6.09e+02|1.01e+02|6.63e+02|7.48e+01|15\textmu s|2\textmu s|
|schwefel_10d|50.83|2.72|1.38e+03|1.07e+02|3.32e+03|2.38e+02|215\textmu s|44\textmu s|
|schwefel_30d|54.47|2.22|2.38e+03|1.09e+02|1.02e+04|3.66e+02|2ms|519\textmu s|
|rastrigin_2d|32.90|7.71|5.36e+24|2.93e+25|8.60e+51|4.71e+52|43\textmu s|10\textmu s|
|rastrigin_10d|40.93|8.38|6.86e+00|1.61e+00|4.99e+01|2.26e+01|695\textmu s|146\textmu s|
|rastrigin_30d|42.03|6.59|1.11e+01|2.01e+00|1.28e+02|4.78e+01|5ms|883\textmu s|
|schaffer_f6|49.53|49.36|6.38e+00|3.73e+00|4.71e-02|4.21e-02|50\textmu s|48\textmu s|
|rosenbrock|39.13|97.15|1.18e+80|6.41e+80|inf|-nan|21\textmu s|57\textmu s|
### cg_prp
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|1.00|0.00|5.06e-14|7.68e-14|8.26e-27|2.47e-26|0\textmu s|0\textmu s|
|sphere_10d|1.00|0.00|6.60e-13|3.26e-13|5.38e-25|5.19e-25|1\textmu s|0\textmu s|
|sphere_30d|1.00|0.00|3.93e-12|1.95e-12|1.91e-23|1.86e-23|4\textmu s|0\textmu s|
|schwefel_2d|43.73|3.87|6.09e+02|1.01e+02|6.63e+02|7.48e+01|16\textmu s|3\textmu s|
|schwefel_10d|50.83|2.72|1.38e+03|1.07e+02|3.32e+03|2.38e+02|216\textmu s|30\textmu s|
|schwefel_30d|54.47|2.22|2.38e+03|1.09e+02|1.02e+04|3.66e+02|2ms|529\textmu s|
|rastrigin_2d|32.90|7.71|5.36e+24|2.93e+25|8.60e+51|4.71e+52|48\textmu s|13\textmu s|
|rastrigin_10d|40.93|8.38|6.86e+00|1.61e+00|4.99e+01|2.26e+01|666\textmu s|154\textmu s|
|rastrigin_30d|42.03|6.59|1.11e+01|2.01e+00|1.28e+02|4.78e+01|5ms|895\textmu s|
|schaffer_f6|49.53|49.36|6.38e+00|3.73e+00|4.71e-02|4.21e-02|48\textmu s|49\textmu s|
|rosenbrock|39.13|97.15|1.18e+80|6.41e+80|inf|-nan|19\textmu s|49\textmu s|
### nelder_mead
| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std
|-|-|-|-|-|-|-|-|-|
|sphere_2d|16.77|2.64|2.99e-03|1.89e-03|1.24e-05|1.39e-05|1\textmu s|4\textmu s|
|sphere_10d|169.07|20.54|5.97e-03|1.13e-03|3.68e-05|1.54e-05|23\textmu s|6\textmu s|
|sphere_30d|500.00|0.00|1.51e-01|7.15e-02|2.79e-02|3.10e-02|219\textmu s|38\textmu s|
|schwefel_2d|40.33|4.10|6.07e+02|9.57e+01|6.62e+02|7.28e+01|4\textmu s|0\textmu s|
|schwefel_10d|475.00|23.32|1.37e+03|1.09e+02|3.33e+03|2.29e+02|115\textmu s|15\textmu s|
|schwefel_30d|500.00|0.00|2.38e+03|1.07e+02|1.05e+04|3.94e+02|380\textmu s|51\textmu s|
|rastrigin_2d|42.23|9.78|1.69e+00|9.44e-01|3.75e+00|3.33e+00|6\textmu s|1\textmu s|
|rastrigin_10d|497.23|10.88|4.14e+00|9.20e-01|1.90e+01|8.68e+00|260\textmu s|40\textmu s|
|rastrigin_30d|500.00|0.00|8.53e+00|1.24e+00|1.18e+02|3.10e+01|932\textmu s|292\textmu s|
|schaffer_f6|14.60|3.69|6.38e+00|4.08e+00|4.76e-02|5.05e-02|2\textmu s|0\textmu s|
|rosenbrock|80.60|23.91|4.90e-03|4.66e-03|2.51e-05|4.35e-05|5\textmu s|3\textmu s|
