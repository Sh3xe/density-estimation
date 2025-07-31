### LBFGS30
Method|Iters|x\_diff|f\_diff|duration
-|-|-|-|-
sphere\_2d (fdaPDE)|1.0 \textpm 0.0 |5.06e-14 \textpm 7.68e-14|1.91e+00\textpm 1.93e+00|1ns\textpm 4ns
sphere\_2d (scipy)|1.9 \textpm 0.3 |1.29e-08 \textpm 1.90e-08|5.16e-16\textpm 1.70e-15|388\textmu s\textpm 233\textmu s
sphere\_10d (fdaPDE)|1.0 \textpm 0.0 |6.60e-13 \textpm 3.26e-13|1.01e+01\textpm 4.21e+00|2ns\textpm 0ns
sphere\_10d (scipy)|2.0 \textpm 0.0 |3.64e-07 \textpm 2.90e-07|2.13e-13\textpm 3.20e-13|425\textmu s\textpm 45\textmu s
sphere\_30d (fdaPDE)|1.0 \textpm 0.0 |3.93e-12 \textpm 1.95e-12|3.02e+01\textpm 8.94e+00|6ns\textpm 0ns
sphere\_30d (scipy)|2.0 \textpm 0.0 |4.14e-06 \textpm 2.41e-06|2.27e-11\textpm 2.40e-11|716\textmu s\textpm 15\textmu s
rastrigin\_2d (fdaPDE)|17.6 \textpm 58.7 |2.72e+00 \textpm 8.31e-01|8.09e+00\textpm 4.49e+00|70ns\textpm 313ns
rastrigin\_2d (scipy)|5.5 \textpm 1.5 |2.70e+00 \textpm 8.82e-01|8.09e+00\textpm 4.75e+00|836\textmu s\textpm 237\textmu s
rastrigin\_10d (fdaPDE)|23.8 \textpm 37.3 |6.55e+00 \textpm 1.51e+00|4.53e+01\textpm 2.06e+01|325ns\textpm 764ns
rastrigin\_10d (scipy)|7.5 \textpm 2.2 |6.37e+00 \textpm 1.19e+00|4.21e+01\textpm 1.51e+01|1ms\textpm 514\textmu s
rastrigin\_30d (fdaPDE)|19.1 \textpm 22.2 |1.10e+01 \textpm 1.90e+00|1.24e+02\textpm 4.46e+01|1\textmu s\textpm 2\textmu s
rastrigin\_30d (scipy)|9.8 \textpm 1.7 |1.06e+01 \textpm 1.43e+00|1.16e+02\textpm 3.05e+01|10ms\textpm 2ms
schwefel\_2d (fdaPDE)|7.1 \textpm 1.7 |6.08e+02 \textpm 1.05e+02|6.50e+02\textpm 8.36e+01|7ns\textpm 4ns
schwefel\_2d (scipy)|5.5 \textpm 1.7 |6.12e+02 \textpm 1.19e+02|6.46e+02\textpm 8.51e+01|840\textmu s\textpm 291\textmu s
schwefel\_10d (fdaPDE)|12.4 \textpm 2.7 |1.38e+03 \textpm 1.12e+02|3.29e+03\textpm 2.90e+02|66ns\textpm 17ns
schwefel\_10d (scipy)|9.1 \textpm 3.3 |1.38e+03 \textpm 1.07e+02|3.31e+03\textpm 2.26e+02|3ms\textpm 1ms
schwefel\_30d (fdaPDE)|14.9 \textpm 2.4 |2.38e+03 \textpm 1.13e+02|1.02e+04\textpm 3.69e+02|751ns\textpm 211ns
schwefel\_30d (scipy)|9.5 \textpm 1.7 |2.38e+03 \textpm 1.11e+02|1.02e+04\textpm 3.65e+02|16ms\textpm 2ms
rosenbrock (fdaPDE)|196.5 \textpm 215.3 |4.39e-04 \textpm 1.90e-05|4.04e-08\textpm 4.88e-09|341ns\textpm 408ns
rosenbrock (scipy)|33.8 \textpm 3.6 |4.61e-04 \textpm 4.21e-04|1.31e-07\textpm 1.70e-07|3ms\textpm 378\textmu s
schaffer\_f6 (fdaPDE)|3.9 \textpm 0.5 |6.49e+00 \textpm 3.59e+00|4.74e-02\textpm 4.18e-02|3ns\textpm 1ns
schaffer\_f6 (scipy)|3.5 \textpm 0.6 |6.70e+00 \textpm 3.57e+00|4.94e-02\textpm 4.36e-02|507\textmu s\textpm 98\textmu s

### Nelder-Mead
Method|Iters|x\_diff|f\_diff|duration
-|-|-|-|-
sphere\_2d (fdaPDE)|16.8 \textpm 2.6 |2.99e-03 \textpm 1.89e-03|1.24e-05\textpm 1.39e-05|1ns\textpm 4ns
sphere\_2d (scipy)|50.1 \textpm 3.8 |3.32e-06 \textpm 1.05e-06|1.21e-11\textpm 8.12e-12|846\textmu s\textpm 223\textmu s
sphere\_10d (fdaPDE)|169.1 \textpm 20.5 |5.97e-03 \textpm 1.13e-03|3.68e-05\textpm 1.54e-05|23ns\textpm 6ns
sphere\_10d (scipy)|500.0 \textpm 0.0 |5.94e-01 \textpm 3.29e-01|4.57e-01\textpm 5.40e-01|6ms\textpm 157\textmu s
sphere\_30d (fdaPDE)|500.0 \textpm 0.0 |1.51e-01 \textpm 7.15e-02|2.79e-02\textpm 3.10e-02|219ns\textpm 38ns
sphere\_30d (scipy)|500.0 \textpm 0.0 |3.71e+00 \textpm 7.61e-01|1.43e+01\textpm 5.82e+00|6ms\textpm 154\textmu s
rastrigin\_2d (fdaPDE)|42.2 \textpm 9.8 |1.69e+00 \textpm 9.44e-01|3.75e+00\textpm 3.33e+00|6ns\textpm 1ns
rastrigin\_2d (scipy)|39.1 \textpm 3.9 |2.92e+00 \textpm 1.14e+00|9.82e+00\textpm 7.24e+00|716\textmu s\textpm 56\textmu s
rastrigin\_10d (fdaPDE)|497.2 \textpm 10.9 |4.14e+00 \textpm 9.20e-01|1.90e+01\textpm 8.68e+00|260ns\textpm 40ns
rastrigin\_10d (scipy)|500.0 \textpm 0.0 |6.43e+00 \textpm 1.19e+00|4.67e+01\textpm 1.73e+01|11ms\textpm 332\textmu s
rastrigin\_30d (fdaPDE)|500.0 \textpm 0.0 |8.53e+00 \textpm 1.24e+00|1.18e+02\textpm 3.10e+01|932ns\textpm 292ns
rastrigin\_30d (scipy)|500.0 \textpm 0.0 |1.07e+01 \textpm 1.43e+00|2.24e+02\textpm 3.76e+01|19ms\textpm 464\textmu s
schwefel\_2d (fdaPDE)|40.3 \textpm 4.1 |6.07e+02 \textpm 9.57e+01|6.62e+02\textpm 7.28e+01|4ns\textpm 0ns
schwefel\_2d (scipy)|55.1 \textpm 5.4 |6.10e+02 \textpm 9.59e+01|6.65e+02\textpm 7.42e+01|1ms\textpm 121\textmu s
schwefel\_10d (fdaPDE)|475.0 \textpm 23.3 |1.37e+03 \textpm 1.09e+02|3.33e+03\textpm 2.29e+02|115ns\textpm 15ns
schwefel\_10d (scipy)|500.0 \textpm 0.0 |1.38e+03 \textpm 1.11e+02|3.33e+03\textpm 2.26e+02|13ms\textpm 576\textmu s
schwefel\_30d (fdaPDE)|500.0 \textpm 0.0 |2.38e+03 \textpm 1.07e+02|1.05e+04\textpm 3.94e+02|380ns\textpm 51ns
schwefel\_30d (scipy)|500.0 \textpm 0.0 |2.38e+03 \textpm 1.12e+02|1.08e+04\textpm 4.60e+02|23ms\textpm 480\textmu s
rosenbrock (fdaPDE)|80.6 \textpm 23.9 |4.90e-03 \textpm 4.66e-03|2.51e-05\textpm 4.35e-05|5ns\textpm 3ns
rosenbrock (scipy)|105.8 \textpm 10.9 |2.09e-06 \textpm 1.12e-06|6.48e-12\textpm 1.05e-11|1ms\textpm 176\textmu s
schaffer\_f6 (fdaPDE)|14.6 \textpm 3.7 |6.38e+00 \textpm 4.08e+00|4.76e-02\textpm 5.05e-02|2ns\textpm 0ns
schaffer\_f6 (scipy)|53.9 \textpm 4.0 |6.49e+00 \textpm 3.59e+00|4.74e-02\textpm 4.18e-02|954\textmu s\textpm 95\textmu s

### CGFR
Method|Iters|x\_diff|f\_diff|duration
-|-|-|-|-
sphere\_2d (fdaPDE)|1.0 \textpm 0.0 |5.06e-14 \textpm 7.68e-14|8.26e-27\textpm 2.47e-26|0ns\textpm 0ns
sphere\_2d (scipy)|1.7 \textpm 0.7 |1.83e-08 \textpm 2.80e-08|1.09e-15\textpm 3.21e-15|422\textmu s\textpm 95\textmu s
sphere\_10d (fdaPDE)|1.0 \textpm 0.0 |6.60e-13 \textpm 3.26e-13|5.38e-25\textpm 5.19e-25|1ns\textpm 0ns
sphere\_10d (scipy)|1.1 \textpm 0.4 |3.45e-07 \textpm 5.36e-07|3.97e-13\textpm 1.53e-12|473\textmu s\textpm 70\textmu s
sphere\_30d (fdaPDE)|1.0 \textpm 0.0 |3.93e-12 \textpm 1.95e-12|1.91e-23\textpm 1.86e-23|5ns\textpm 1ns
sphere\_30d (scipy)|1.9 \textpm 0.7 |9.64e-07 \textpm 2.18e-06|5.53e-12\textpm 1.76e-11|1ms\textpm 353\textmu s
rastrigin\_2d (fdaPDE)|500.0 \textpm 0.0 |8.01e+03 \textpm 3.37e+04|1.16e+10\textpm 6.16e+10|590ns\textpm 86ns
rastrigin\_2d (scipy)|5.9 \textpm 1.5 |2.61e+00 \textpm 9.84e-01|7.79e+00\textpm 5.23e+00|1ms\textpm 500\textmu s
rastrigin\_10d (fdaPDE)|500.0 \textpm 0.0 |4.21e+01 \textpm 8.70e+01|9.04e+04\textpm 3.02e+05|9\textmu s\textpm 2\textmu s
rastrigin\_10d (scipy)|7.9 \textpm 1.9 |6.37e+00 \textpm 1.19e+00|4.21e+01\textpm 1.51e+01|3ms\textpm 815\textmu s
rastrigin\_30d (fdaPDE)|500.0 \textpm 0.0 |3.23e+02 \textpm 1.51e+03|2.31e+07\textpm 1.25e+08|58\textmu s\textpm 13\textmu s
rastrigin\_30d (scipy)|9.4 \textpm 1.3 |1.06e+01 \textpm 1.39e+00|1.15e+02\textpm 2.94e+01|17ms\textpm 2ms
schwefel\_2d (fdaPDE)|488.2 \textpm 64.4 |1.24e+13 \textpm 6.81e+13|9.29e+25\textpm 5.09e+26|469ns\textpm 119ns
schwefel\_2d (scipy)|6.6 \textpm 1.9 |6.09e+02 \textpm 1.01e+02|6.63e+02\textpm 7.48e+01|1ms\textpm 414\textmu s
schwefel\_10d (fdaPDE)|488.1 \textpm 65.4 |1.35e+13 \textpm 7.39e+13|1.09e+26\textpm 5.98e+26|5\textmu s\textpm 1\textmu s
schwefel\_10d (scipy)|12.5 \textpm 2.5 |1.38e+03 \textpm 1.13e+02|3.29e+03\textpm 2.96e+02|7ms\textpm 2ms
schwefel\_30d (fdaPDE)|494.7 \textpm 29.2 |5.74e+11 \textpm 3.15e+12|1.98e+23\textpm 1.08e+24|53\textmu s\textpm 16\textmu s
schwefel\_30d (scipy)|15.0 \textpm 2.3 |2.39e+03 \textpm 1.10e+02|1.02e+04\textpm 3.66e+02|41ms\textpm 11ms
rosenbrock (fdaPDE)|29.5 \textpm 26.4 |2.51e+39 \textpm 1.03e+40|2.73e+163\textpm inf|16ns\textpm 14ns
rosenbrock (scipy)|27.2 \textpm 5.2 |1.12e-05 \textpm 1.66e-05|1.83e-10\textpm 9.00e-10|5ms\textpm 1ms
schaffer\_f6 (fdaPDE)|8.5 \textpm 3.1 |1.30e+06 \textpm 3.72e+06|5.00e-01\textpm 9.13e-07|12ns\textpm 4ns
schaffer\_f6 (scipy)|2.1 \textpm 0.4 |6.49e+00 \textpm 3.59e+00|4.74e-02\textpm 4.18e-02|656\textmu s\textpm 187\textmu s

### CGPR
Method|Iters|x\_diff|f\_diff|duration
-|-|-|-|-
sphere\_2d (fdaPDE)|1.0 \textpm 0.0 |5.06e-14 \textpm 7.68e-14|8.26e-27\textpm 2.47e-26|0ns\textpm 0ns
sphere\_2d (scipy)|1.7 \textpm 0.7 |1.83e-08 \textpm 2.80e-08|1.09e-15\textpm 3.21e-15|422\textmu s\textpm 95\textmu s
sphere\_10d (fdaPDE)|1.0 \textpm 0.0 |6.60e-13 \textpm 3.26e-13|5.38e-25\textpm 5.19e-25|1ns\textpm 0ns
sphere\_10d (scipy)|1.1 \textpm 0.4 |3.45e-07 \textpm 5.36e-07|3.97e-13\textpm 1.53e-12|473\textmu s\textpm 70\textmu s
sphere\_30d (fdaPDE)|1.0 \textpm 0.0 |3.93e-12 \textpm 1.95e-12|1.91e-23\textpm 1.86e-23|5ns\textpm 0ns
sphere\_30d (scipy)|1.9 \textpm 0.7 |9.64e-07 \textpm 2.18e-06|5.53e-12\textpm 1.76e-11|1ms\textpm 353\textmu s
rastrigin\_2d (fdaPDE)|32.9 \textpm 7.7 |5.36e+24 \textpm 2.93e+25|8.60e+51\textpm 4.71e+52|43ns\textpm 10ns
rastrigin\_2d (scipy)|5.9 \textpm 1.5 |2.61e+00 \textpm 9.84e-01|7.79e+00\textpm 5.23e+00|1ms\textpm 500\textmu s
rastrigin\_10d (fdaPDE)|40.9 \textpm 8.4 |6.86e+00 \textpm 1.61e+00|4.99e+01\textpm 2.26e+01|695ns\textpm 146ns
rastrigin\_10d (scipy)|7.9 \textpm 1.9 |6.37e+00 \textpm 1.19e+00|4.21e+01\textpm 1.51e+01|3ms\textpm 815\textmu s
rastrigin\_30d (fdaPDE)|42.0 \textpm 6.6 |1.11e+01 \textpm 2.01e+00|1.28e+02\textpm 4.78e+01|5\textmu s\textpm 883ns
rastrigin\_30d (scipy)|9.4 \textpm 1.3 |1.06e+01 \textpm 1.39e+00|1.15e+02\textpm 2.94e+01|17ms\textpm 2ms
schwefel\_2d (fdaPDE)|43.7 \textpm 3.9 |6.09e+02 \textpm 1.01e+02|6.63e+02\textpm 7.48e+01|15ns\textpm 2ns
schwefel\_2d (scipy)|6.6 \textpm 1.9 |6.09e+02 \textpm 1.01e+02|6.63e+02\textpm 7.48e+01|1ms\textpm 414\textmu s
schwefel\_10d (fdaPDE)|50.8 \textpm 2.7 |1.38e+03 \textpm 1.07e+02|3.32e+03\textpm 2.38e+02|215ns\textpm 44ns
schwefel\_10d (scipy)|12.5 \textpm 2.5 |1.38e+03 \textpm 1.13e+02|3.29e+03\textpm 2.96e+02|7ms\textpm 2ms
schwefel\_30d (fdaPDE)|54.5 \textpm 2.2 |2.38e+03 \textpm 1.09e+02|1.02e+04\textpm 3.66e+02|2\textmu s\textpm 519ns
schwefel\_30d (scipy)|15.0 \textpm 2.3 |2.39e+03 \textpm 1.10e+02|1.02e+04\textpm 3.66e+02|41ms\textpm 11ms
rosenbrock (fdaPDE)|39.1 \textpm 97.1 |1.18e+80 \textpm 6.41e+80|inf\textpm nan|21ns\textpm 57ns
rosenbrock (scipy)|27.2 \textpm 5.2 |1.12e-05 \textpm 1.66e-05|1.83e-10\textpm 9.00e-10|5ms\textpm 1ms
schaffer\_f6 (fdaPDE)|49.5 \textpm 49.4 |6.38e+00 \textpm 3.73e+00|4.71e-02\textpm 4.21e-02|50ns\textpm 48ns
schaffer\_f6 (scipy)|2.1 \textpm 0.4 |6.49e+00 \textpm 3.59e+00|4.74e-02\textpm 4.18e-02|656\textmu s\textpm 187\textmu s

### CGPRP
Method|Iters|x\_diff|f\_diff|duration
-|-|-|-|-
sphere\_2d (fdaPDE)|1.0 \textpm 0.0 |5.06e-14 \textpm 7.68e-14|8.26e-27\textpm 2.47e-26|0ns\textpm 0ns
sphere\_2d (scipy)|1.7 \textpm 0.7 |1.83e-08 \textpm 2.80e-08|1.09e-15\textpm 3.21e-15|422\textmu s\textpm 95\textmu s
sphere\_10d (fdaPDE)|1.0 \textpm 0.0 |6.60e-13 \textpm 3.26e-13|5.38e-25\textpm 5.19e-25|1ns\textpm 0ns
sphere\_10d (scipy)|1.1 \textpm 0.4 |3.45e-07 \textpm 5.36e-07|3.97e-13\textpm 1.53e-12|473\textmu s\textpm 70\textmu s
sphere\_30d (fdaPDE)|1.0 \textpm 0.0 |3.93e-12 \textpm 1.95e-12|1.91e-23\textpm 1.86e-23|4ns\textpm 0ns
sphere\_30d (scipy)|1.9 \textpm 0.7 |9.64e-07 \textpm 2.18e-06|5.53e-12\textpm 1.76e-11|1ms\textpm 353\textmu s
rastrigin\_2d (fdaPDE)|32.9 \textpm 7.7 |5.36e+24 \textpm 2.93e+25|8.60e+51\textpm 4.71e+52|48ns\textpm 13ns
rastrigin\_2d (scipy)|5.9 \textpm 1.5 |2.61e+00 \textpm 9.84e-01|7.79e+00\textpm 5.23e+00|1ms\textpm 500\textmu s
rastrigin\_10d (fdaPDE)|40.9 \textpm 8.4 |6.86e+00 \textpm 1.61e+00|4.99e+01\textpm 2.26e+01|666ns\textpm 154ns
rastrigin\_10d (scipy)|7.9 \textpm 1.9 |6.37e+00 \textpm 1.19e+00|4.21e+01\textpm 1.51e+01|3ms\textpm 815\textmu s
rastrigin\_30d (fdaPDE)|42.0 \textpm 6.6 |1.11e+01 \textpm 2.01e+00|1.28e+02\textpm 4.78e+01|5\textmu s\textpm 895ns
rastrigin\_30d (scipy)|9.4 \textpm 1.3 |1.06e+01 \textpm 1.39e+00|1.15e+02\textpm 2.94e+01|17ms\textpm 2ms
schwefel\_2d (fdaPDE)|43.7 \textpm 3.9 |6.09e+02 \textpm 1.01e+02|6.63e+02\textpm 7.48e+01|16ns\textpm 3ns
schwefel\_2d (scipy)|6.6 \textpm 1.9 |6.09e+02 \textpm 1.01e+02|6.63e+02\textpm 7.48e+01|1ms\textpm 414\textmu s
schwefel\_10d (fdaPDE)|50.8 \textpm 2.7 |1.38e+03 \textpm 1.07e+02|3.32e+03\textpm 2.38e+02|216ns\textpm 30ns
schwefel\_10d (scipy)|12.5 \textpm 2.5 |1.38e+03 \textpm 1.13e+02|3.29e+03\textpm 2.96e+02|7ms\textpm 2ms
schwefel\_30d (fdaPDE)|54.5 \textpm 2.2 |2.38e+03 \textpm 1.09e+02|1.02e+04\textpm 3.66e+02|2\textmu s\textpm 529ns
schwefel\_30d (scipy)|15.0 \textpm 2.3 |2.39e+03 \textpm 1.10e+02|1.02e+04\textpm 3.66e+02|41ms\textpm 11ms
rosenbrock (fdaPDE)|39.1 \textpm 97.1 |1.18e+80 \textpm 6.41e+80|inf\textpm nan|19ns\textpm 49ns
rosenbrock (scipy)|27.2 \textpm 5.2 |1.12e-05 \textpm 1.66e-05|1.83e-10\textpm 9.00e-10|5ms\textpm 1ms
schaffer\_f6 (fdaPDE)|49.5 \textpm 49.4 |6.38e+00 \textpm 3.73e+00|4.71e-02\textpm 4.21e-02|48ns\textpm 49ns
schaffer\_f6 (scipy)|2.1 \textpm 0.4 |6.49e+00 \textpm 3.59e+00|4.74e-02\textpm 4.18e-02|656\textmu s\textpm 187\textmu s

