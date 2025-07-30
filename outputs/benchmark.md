### horseshoe
Optimizer|Iterations|Loss|CV error|lambda|cv_duration|duration|
-|-|-|-|-|-|-
lbfgs30|83|2797.24|1.85504|0.01|1s|58ms|
grad_descent|500|2758.25|1.88091|0.001|9s|1s|
cg_pr|90|2758.26|1.88271|0.001|2s|205ms|
cg_fr|221|-nan|1.79769e+308|1e-05|7s|439ms|

### gaussian_square
Optimizer|Iterations|Loss|CV error|lambda|cv_duration|duration|
-|-|-|-|-|-|-
lbfgs30|105|195.627|-0.759535|0.000316228|6s|99ms|
grad_descent|390|197.021|-0.760926|0.000316228|39s|1s|
cg_pr|131|195.831|-0.760193|0.000316228|21s|410ms|
cg_fr|500|321.189|-0.733355|0.00316228|28s|1s|

### infections_southampton
Optimizer|Iterations|Loss|CV error|lambda|cv_duration|duration|
-|-|-|-|-|-|-
lbfgs30|133|81036.5|7.20417|0.01|16s|403ms|
grad_descent|500|81037.8|7.19859|0.01|101s|3s|
cg_pr|80|81040.4|7.20035|0.01|49s|699ms|
cg_fr|350|81036.7|7.20262|0.01|86s|2s|

