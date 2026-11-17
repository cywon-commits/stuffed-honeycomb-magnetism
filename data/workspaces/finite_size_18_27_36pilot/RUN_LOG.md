# Measured run log

The two 18-spin shapes ran concurrently under a two-CPU quota. Per-point walls overlap; their sum is not a single elapsed turnaround. Target-H counts exclude auxiliary S² actions and some validation work. No unconverged eigenvalues are used.

## Final verification

{
  "status": "PASS",
  "number_of_18spin_couplings_and_shapes": 122,
  "number_of_18spin_momentum_minima": 732,
  "number_of_18spin_returned_eigenpairs": 3660,
  "max_ground_H_residual": 4.691651745880047e-14,
  "max_all_eigenpair_H_residual": 2.2006568501437907e-12,
  "max_ground_S2_residual": 4.907653602423959e-10,
  "max_all_eigenpair_S2_residual": 1.1321672950292963e-09,
  "max_18spin_total_spin_sum_rule_error": 1.857541898075965e-14,
  "max_27spin_total_spin_sum_rule_error": 5.320854867818525e-12,
  "max_18spin_energy_from_G_error": 1.4210854715202004e-14,
  "note": "Global ground-state coverage is at all sampled points. Numerical scans cannot rule out arbitrarily narrow re-entrant plateaus between unchanged endpoints."
}

## Phase controller records

[
  {
    "phase": "finish18.py",
    "returncode": 0,
    "wall_seconds": 55.10984520799684
  },
  {
    "phase": "audit27.py",
    "returncode": 0,
    "wall_seconds": 32.710356812996906
  },
  {
    "phase": "pilot36.py",
    "returncode": 0,
    "wall_seconds": 1.8180341240004054
  }
]

## Environment

{
  "CPU_model": "INTEL(R) XEON(R) PLATINUM 8573C",
  "CPU_quota": "200000 100000",
  "RAM_limit_bytes": 8589934592,
  "Python": "3.12.14",
  "NumPy": "2.3.5",
  "OpenMP_threads": 2,
  "BLAS_threads": 1
}

## Near-crossing spin resolution

{
  "status": "PASS",
  "rows": [
    {
      "shape": "K18",
      "r": 0.501706589,
      "old_max_S2_residual": 1.210452139149408e-07,
      "new_max_S2_residual": 1.343823093515504e-11,
      "ground_energy_difference": 1.7763568394002505e-15,
      "method": "Joint H/S2 subspace diagonalization, then H diagonalization separately for each integer S."
    },
    {
      "shape": "C18",
      "r": 0.3346873744,
      "old_max_S2_residual": 2.3908413376657875e-07,
      "new_max_S2_residual": 2.4145621422594323e-11,
      "ground_energy_difference": 5.329070518200751e-15,
      "method": "Joint H/S2 subspace diagonalization, then H diagonalization separately for each integer S."
    },
    {
      "shape": "C18",
      "r": 0.5042126561,
      "old_max_S2_residual": 1.4528769678056387e-07,
      "new_max_S2_residual": 3.5958663698417085e-11,
      "ground_energy_difference": 4.440892098500626e-15,
      "method": "Joint H/S2 subspace diagonalization, then H diagonalization separately for each integer S."
    }
  ]
}

## Original scan progress


### K18_scan.log

```text

{'number_of_index6_sublattices': 12, 'valid_K_compatible_shapes': 3, 'all_K_compatible_shapes_pointgroup_equivalent': True, 'control_T': [[1, 2], [0, 6]], 'geometry_metric': 'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.', 'caution': 'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'}
{"shape": "K18", "r": 0.1, "E": -7.02418677324936, "S": 0, "q": [0.0, 0.5], "pH": -0.5984507135247523, "wall_seconds": 4.811423379993357}
{"shape": "K18", "r": 0.15, "E": -7.038674131286463, "S": 0, "q": [0.0, 0.5], "pH": -0.5962685872930532, "wall_seconds": 4.577936490997672}
{"shape": "K18", "r": 0.2, "E": -7.0601021049226915, "S": 0, "q": [0.0, 0.5], "pH": -0.5928338541133272, "wall_seconds": 4.639649983000709}
{"shape": "K18", "r": 0.225, "E": -7.073705087618323, "S": 0, "q": [0.0, 0.5], "pH": -0.5905353756418664, "wall_seconds": 4.72260077700048}
{"shape": "K18", "r": 0.25, "E": -7.089417994728179, "S": 0, "q": [0.0, 0.5], "pH": -0.587768799196314, "wall_seconds": 4.519339779006259}
{"shape": "K18", "r": 0.275, "E": -7.107402300639545, "S": 0, "q": [0.0, 0.5], "pH": -0.5844543902572349, "wall_seconds": 4.152845333002915}
{"shape": "K18", "r": 0.3, "E": -7.12784169897646, "S": 0, "q": [0.0, 0.5], "pH": -0.5804907273352342, "wall_seconds": 4.47672387700004}
{"shape": "K18", "r": 0.35, "E": -7.176952883060874, "S": 0, "q": [0.0, 0.5], "pH": -0.57004984397884, "wall_seconds": 4.2145277009985875}
{"shape": "K18", "r": 0.4, "E": -7.277648294880409, "S": 2, "q": [0.0, 0.0], "pH": 0.955917336583381, "wall_seconds": 3.9291653129985207}
{"shape": "K18", "r": 0.45, "E": -7.397752855351421, "S": 2, "q": [0.0, 0.0], "pH": 0.9483033821878221, "wall_seconds": 5.493029771998408}
{"shape": "K18", "r": 0.5, "E": -7.519200723337365, "S": 2, "q": [0.0, 0.0], "pH": 0.9406387756140762, "wall_seconds": 5.169513649001601}
{"shape": "K18", "r": 0.55, "E": -7.745559598948251, "S": 1, "q": [0.0, 0.0], "pH": 0.7965539717118484, "wall_seconds": 5.165116226999089}
{"shape": "K18", "r": 0.6, "E": -7.979475647652568, "S": 1, "q": [0.0, 0.0], "pH": 0.8069567624238424, "wall_seconds": 10.012048542004777}
{"shape": "K18", "r": 0.65, "E": -8.21533032180254, "S": 1, "q": [0.0, 0.0], "pH": 0.8087959145175281, "wall_seconds": 8.661535327002639}
{"shape": "K18", "r": 0.675, "E": -8.333759173488389, "S": 1, "q": [0.0, 0.0], "pH": 0.8082951495292233, "wall_seconds": 14.199665584994364}
{"shape": "K18", "r": 0.68, "E": -8.357478830563611, "S": 1, "q": [0.0, 0.0], "pH": 0.8081227886024381, "wall_seconds": 7.651057636998303}
{"shape": "K18", "r": 0.7, "E": -8.4524615579071, "S": 1, "q": [0.0, 0.0], "pH": 0.8072446030561296, "wall_seconds": 8.774585249993834}
{"shape": "K18", "r": 0.375, "E": -7.218135778815232, "S": 2, "q": [0.0, 0.0], "pH": 0.9596903914383685, "wall_seconds": 7.362030042997503}
{"shape": "K18", "r": 0.525, "E": -7.629835057789947, "S": 1, "q": [0.0, 0.0], "pH": 0.7832827197771877, "wall_seconds": 8.26014167899848}
{"shape": "K18", "r": 0.3625, "E": -7.191129378099118, "S": 0, "q": [0.0, 0.5], "pH": -0.5667746912051241, "wall_seconds": 10.751389960001688}
{"shape": "K18", "r": 0.5125, "E": -7.572486446688415, "S": 1, "q": [0.0, 0.0], "pH": 0.7723177494930494, "wall_seconds": 13.979202813003212}
{"shape": "K18", "r": 0.36875, "E": -7.20331665380658, "S": 2, "q": [0.0, 0.0], "pH": 0.9606287884897623, "wall_seconds": 11.859652829996776}
{"shape": "K18", "r": 0.50625, "E": -7.543990918560221, "S": 1, "q": [0.0, 0.0], "pH": 0.765115824240051, "wall_seconds": 15.362911367003107}
{"shape": "K18", "r": 0.365625, "E": -7.195916133413717, "S": 2, "q": [0.0, 0.0], "pH": 0.9610971570421323, "wall_seconds": 11.746990396997717}
{"shape": "K18", "r": 0.503125, "E": -7.529797265586565, "S": 1, "q": [0.0, 0.0], "pH": 0.7609574122097946, "wall_seconds": 15.813166081999952}
{"shape": "K18", "r": 0.3640625, "E": -7.192959099054588, "S": 0, "q": [0.0, 0.5], "pH": -0.5663428336941768, "wall_seconds": 11.601988573995186}
{"shape": "K18", "r": 0.5015625, "E": -7.523016511417092, "S": 2, "q": [0.0, 0.0], "pH": 0.9403988244781118, "wall_seconds": 26.20984055300505}
{"shape": "K18", "r": 0.1125, "E": -7.027201374468421, "S": 0, "q": [0.0, 0.5], "pH": -0.5980070752695775, "wall_seconds": 15.582383187997038}
{"shape": "K18", "r": 0.125, "E": -7.030613286772205, "S": 0, "q": [0.0, 0.5], "pH": -0.5974980566553552, "wall_seconds": 13.954670771003293}
{"shape": "K18", "r": 0.1375, "E": -7.0344336036175115, "S": 0, "q": [0.0, 0.5], "pH": -0.5969199259094208, "wall_seconds": 12.121814949998225}
{"shape": "K18", "r": 0.1625, "E": -7.043347434884303, "S": 0, "q": [0.0, 0.5], "pH": -0.5955395373986192, "wall_seconds": 9.75154011700215}
{"shape": "K18", "r": 0.175, "E": -7.048466887966186, "S": 0, "q": [0.0, 0.5], "pH": -0.594727814710879, "wall_seconds": 11.210028950998094}
{"shape": "K18", "r": 0.1875, "E": -7.054046726126858, "S": 0, "q": [0.0, 0.5], "pH": -0.5938279411596776, "wall_seconds": 10.367686988996866}
{"shape": "K18", "r": 0.2125, "E": -7.066649162533911, "S": 0, "q": [0.0, 0.5], "pH": -0.5917388269161605, "wall_seconds": 10.28344358499453}
{"shape": "K18", "r": 0.2375, "E": -7.081288192855026, "S": 0, "q": [0.0, 0.5], "pH": -0.5892151491930314, "wall_seconds": 10.006209107996256}
{"shape": "K18", "r": 0.2625, "E": -7.098115300157267, "S": 0, "q": [0.0, 0.5], "pH": -0.58618582528316, "wall_seconds": 9.962289188006253}
{"shape": "K18", "r": 0.2875, "E": -7.1173026746335735, "S": 0, "q": [0.0, 0.5], "pH": -0.5825610982664161, "wall_seconds": 10.202097902998503}
{"shape": "K18", "r": 0.3125, "E": -7.13904637019006, "S": 0, "q": [0.0, 0.5], "pH": -0.578225905365091, "wall_seconds": 13.489356700003555}
{"shape": "K18", "r": 0.325, "E": -7.150945536587965, "S": 0, "q": [0.0, 0.5], "pH": -0.5757467158376933, "wall_seconds": 10.174849731003633}
{"shape": "K18", "r": 0.3375, "E": -7.163570042138846, "S": 0, "q": [0.0, 0.5], "pH": -0.5730302157917798, "wall_seconds": 10.578953813994303}
{"shape": "K18", "r": 0.3875, "E": -7.247845378797713, "S": 2, "q": [0.0, 0.0], "pH": 0.9578074789817753, "wall_seconds": 11.638310665999597}
{"shape": "K18", "r": 0.4125, "E": -7.307542438944073, "S": 2, "q": [0.0, 0.0], "pH": 0.9540209833724539, "wall_seconds": 11.031767888001923}
{"shape": "K18", "r": 0.425, "E": -7.337525824419729, "S": 2, "q": [0.0, 0.0], "pH": 0.9521193377482692, "wall_seconds": 8.308944402000634}
{"shape": "K18", "r": 0.4375, "E": -7.367596561949305, "S": 2, "q": [0.0, 0.0], "pH": 0.9502132246396049, "wall_seconds": 10.358604380002362}
{"shape": "K18", "r": 0.4625, "E": -7.427992997877978, "S": 2, "q": [0.0, 0.0], "pH": 0.9463904679013763, "wall_seconds": 9.879821529000765}
{"shape": "K18", "r": 0.475, "E": -7.4583153687495045, "S": 2, "q": [0.0, 0.0], "pH": 0.9444750643255434, "wall_seconds": 10.890357817996119}
{"shape": "K18", "r": 0.4875, "E": -7.4887184299591105, "S": 2, "q": [0.0, 0.0], "pH": 0.9425576842675583, "wall_seconds": 12.79876638299902}
{"shape": "K18", "r": 0.5375, "E": -7.687555158580644, "S": 1, "q": [0.0, 0.0], "pH": 0.7910132289494487, "wall_seconds": 9.79394324599707}
{"shape": "K18", "r": 0.5625, "E": -7.8037896135119365, "S": 1, "q": [0.0, 0.0], "pH": 0.8005577064815154, "wall_seconds": 8.011860313999932}
{"shape": "K18", "r": 0.575, "E": -7.862204109223952, "S": 1, "q": [0.0, 0.0], "pH": 0.8034495045582614, "wall_seconds": 8.949062530999072}
{"shape": "K18", "r": 0.5875, "E": -7.920773442926821, "S": 1, "q": [0.0, 0.0], "pH": 0.8055160835208319, "wall_seconds": 9.512815194997529}
{"shape": "K18", "r": 0.6125, "E": -8.038294057154364, "S": 1, "q": [0.0, 0.0], "pH": 0.8079137609256011, "wall_seconds": 10.728703287997632}
{"shape": "K18", "r": 0.625, "E": -8.097215761420435, "S": 1, "q": [0.0, 0.0], "pH": 0.808490883552399, "wall_seconds": 12.035259306998341}
{"shape": "K18", "r": 0.6375, "E": -8.156230573640157, "S": 1, "q": [0.0, 0.0], "pH": 0.8087654157200553, "wall_seconds": 18.507354352994298}
{"shape": "K18", "r": 0.6625, "E": -8.274508352041005, "S": 1, "q": [0.0, 0.0], "pH": 0.80862744035465, "wall_seconds": 11.162497516001167}
{"shape": "K18", "r": 0.6875, "E": -8.393078199784394, "S": 1, "q": [0.0, 0.0], "pH": 0.8078268112154274, "wall_seconds": 13.477299805999792}


```


### C18_scan.log

```text

{'number_of_index6_sublattices': 12, 'valid_K_compatible_shapes': 3, 'all_K_compatible_shapes_pointgroup_equivalent': True, 'control_T': [[1, 2], [0, 6]], 'geometry_metric': 'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.', 'caution': 'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'}
{"shape": "C18", "r": 0.1, "E": -6.944298333790638, "S": 0, "q": [0.0, 0.0], "pH": -0.5979653057711846, "wall_seconds": 9.773027766001178}
{"shape": "C18", "r": 0.15, "E": -6.961878773169014, "S": 0, "q": [0.0, 0.0], "pH": -0.59517770381404, "wall_seconds": 12.449051743002201}
{"shape": "C18", "r": 0.2, "E": -6.987448034996537, "S": 0, "q": [0.0, 0.0], "pH": -0.5909410188825451, "wall_seconds": 11.005403201001172}
{"shape": "C18", "r": 0.225, "E": -7.003426094983179, "S": 0, "q": [0.0, 0.0], "pH": -0.5881804218046687, "wall_seconds": 16.390853775999858}
{"shape": "C18", "r": 0.25, "E": -7.021655891667545, "S": 0, "q": [0.0, 0.0], "pH": -0.5849033022731984, "wall_seconds": 13.6924605910026}
{"shape": "C18", "r": 0.275, "E": -7.0422506291777625, "S": 0, "q": [0.0, 0.0], "pH": -0.5809939217805181, "wall_seconds": 8.44427527699736}
{"shape": "C18", "r": 0.3, "E": -7.065354738188504, "S": 0, "q": [0.0, 0.0], "pH": -0.5762544798419907, "wall_seconds": 6.893519569995988}
{"shape": "C18", "r": 0.35, "E": -7.138302381819909, "S": 2, "q": [0.0, 0.0], "pH": 0.9701577334583033, "wall_seconds": 9.993733410999994}
{"shape": "C18", "r": 0.4, "E": -7.258146619875335, "S": 2, "q": [0.0, 0.0], "pH": 0.9634810653827592, "wall_seconds": 8.2728408030016}
{"shape": "C18", "r": 0.45, "E": -7.379674924829612, "S": 2, "q": [0.0, 0.0], "pH": 0.956532682746168, "wall_seconds": 9.607045021999511}
{"shape": "C18", "r": 0.5, "E": -7.502775467884769, "S": 2, "q": [0.0, 0.0], "pH": 0.9493474704595691, "wall_seconds": 9.584738899997319}
{"shape": "C18", "r": 0.55, "E": -7.723052241296866, "S": 1, "q": [0.0, 0.0], "pH": 0.7868789874332696, "wall_seconds": 14.054238541000814}
{"shape": "C18", "r": 0.6, "E": -7.956924300913048, "S": 1, "q": [0.0, 0.0], "pH": 0.8031001274621625, "wall_seconds": 12.061304656002903}
{"shape": "C18", "r": 0.65, "E": -8.19304059916297, "S": 1, "q": [0.0, 0.0], "pH": 0.8063461258036289, "wall_seconds": 12.16067245999875}
{"shape": "C18", "r": 0.675, "E": -8.335941259932634, "S": 0, "q": [0.0, 0.0], "pH": 0.43365762226808613, "wall_seconds": 17.51349416800076}
{"shape": "C18", "r": 0.68, "E": -8.369237723424286, "S": 0, "q": [0.0, 0.0], "pH": 0.44250456387439896, "wall_seconds": 23.607270961998438}
{"shape": "C18", "r": 0.7, "E": -8.503566246669308, "S": 0, "q": [0.0, 0.0], "pH": 0.47211521642438187, "wall_seconds": 12.303600104998623}
{"shape": "C18", "r": 0.325, "E": -7.091167187883309, "S": 0, "q": [0.0, 0.0], "pH": -0.5703274161976458, "wall_seconds": 17.22855319899827}
{"shape": "C18", "r": 0.525, "E": -7.607644330654768, "S": 1, "q": [0.0, 0.0], "pH": 0.7653588043290396, "wall_seconds": 14.579771097996854}
{"shape": "C18", "r": 0.6625, "E": -8.253319204283123, "S": 0, "q": [0.0, 0.0], "pH": 0.40837075673773987, "wall_seconds": 12.14541293199727}
{"shape": "C18", "r": 0.3375, "E": -7.108620045676185, "S": 2, "q": [0.0, 0.0], "pH": 0.971778390708168, "wall_seconds": 11.39340555799572}
{"shape": "C18", "r": 0.5125, "E": -7.550630985090294, "S": 1, "q": [0.0, 0.0], "pH": 0.7468109107554037, "wall_seconds": 13.267521617999591}
{"shape": "C18", "r": 0.65625, "E": -8.222665698530562, "S": 1, "q": [0.0, 0.0], "pH": 0.8063362758098492, "wall_seconds": 11.74069615299959}
{"shape": "C18", "r": 0.33125, "E": -7.098075194617014, "S": 0, "q": [0.0, 0.0], "pH": -0.5685912694386152, "wall_seconds": 8.50228426099784}
{"shape": "C18", "r": 0.50625, "E": -7.522378758458252, "S": 1, "q": [0.0, 0.0], "pH": 0.7343017304509449, "wall_seconds": 11.124860322001041}
{"shape": "C18", "r": 0.659375, "E": -8.237486283606485, "S": 1, "q": [0.0, 0.0], "pH": 0.8063073656904872, "wall_seconds": 20.229337273995043}
{"shape": "C18", "r": 0.334375, "E": -7.101601169830126, "S": 0, "q": [0.0, 0.0], "pH": -0.5676758806936232, "wall_seconds": 15.975276888995722}
{"shape": "C18", "r": 0.503125, "E": -7.510519331896942, "S": 2, "q": [0.0, 0.0], "pH": 0.9488910012994459, "wall_seconds": 21.77297171799728}
{"shape": "C18", "r": 0.6609375, "E": -8.244898544595884, "S": 1, "q": [0.0, 0.0], "pH": 0.8062872535551672, "wall_seconds": 19.598476174003736}
{"shape": "C18", "r": 0.3359375, "E": -7.104917872526701, "S": 2, "q": [0.0, 0.0], "pH": 0.9719794866100372, "wall_seconds": 12.889647716998297}
{"shape": "C18", "r": 0.5046875, "E": -7.515348271111778, "S": 1, "q": [0.0, 0.0], "pH": 0.7307378288172554, "wall_seconds": 13.122981878994324}
{"shape": "C18", "r": 0.1125, "E": -6.947976527776314, "S": 0, "q": [0.0, 0.0], "pH": -0.5973914746197444, "wall_seconds": 14.436283791997994}
{"shape": "C18", "r": 0.125, "E": -6.952126789786114, "S": 0, "q": [0.0, 0.0], "pH": -0.5967376871801411, "wall_seconds": 15.033964987997024}
{"shape": "C18", "r": 0.1375, "E": -6.956757887508266, "S": 0, "q": [0.0, 0.0], "pH": -0.5960008705658901, "wall_seconds": 10.342618538001261}
{"shape": "C18", "r": 0.1625, "E": -6.967498618610788, "S": 0, "q": [0.0, 0.0], "pH": -0.5942645366905707, "wall_seconds": 11.09129938299884}
{"shape": "C18", "r": 0.175, "E": -6.973626862618976, "S": 0, "q": [0.0, 0.0], "pH": -0.5932572803771804, "wall_seconds": 8.618736784999783}
{"shape": "C18", "r": 0.1875, "E": -6.980273274001451, "S": 0, "q": [0.0, 0.0], "pH": -0.5921512604435424, "wall_seconds": 9.753419706998102}
{"shape": "C18", "r": 0.2125, "E": -6.995161851046011, "S": 0, "q": [0.0, 0.0], "pH": -0.5896200467955527, "wall_seconds": 13.272171224998601}
{"shape": "C18", "r": 0.2375, "E": -7.012252996493814, "S": 0, "q": [0.0, 0.0], "pH": -0.5866123132849647, "wall_seconds": 13.954368947001058}
{"shape": "C18", "r": 0.2625, "E": -7.031649553104247, "S": 0, "q": [0.0, 0.0], "pH": -0.5830374386653384, "wall_seconds": 3.8571613550011534}
{"shape": "C18", "r": 0.2875, "E": -7.053478232899767, "S": 0, "q": [0.0, 0.0], "pH": -0.5787452349031869, "wall_seconds": 4.499753598996904}
{"shape": "C18", "r": 0.3125, "E": -7.077906866928872, "S": 0, "q": [0.0, 0.0], "pH": -0.5734715287278088, "wall_seconds": 4.467489773996931}
{"shape": "C18", "r": 0.3625, "E": -7.168098475275254, "S": 2, "q": [0.0, 0.0], "pH": 0.9685167093033277, "wall_seconds": 4.663062468003773}
{"shape": "C18", "r": 0.375, "E": -7.198006001770059, "S": 2, "q": [0.0, 0.0], "pH": 0.9668563284907876, "wall_seconds": 4.380439396001748}
{"shape": "C18", "r": 0.3875, "E": -7.2280227504278205, "S": 2, "q": [0.0, 0.0], "pH": 0.9651775057674752, "wall_seconds": 3.9736959779984318}
{"shape": "C18", "r": 0.4125, "E": -7.288375614733887, "S": 2, "q": [0.0, 0.0], "pH": 0.9617677455790283, "wall_seconds": 3.712252976001764}
{"shape": "C18", "r": 0.425, "E": -7.318707842489365, "S": 2, "q": [0.0, 0.0], "pH": 0.9600382025936192, "wall_seconds": 4.63112741799705}
{"shape": "C18", "r": 0.4375, "E": -7.34914151074066, "S": 2, "q": [0.0, 0.0], "pH": 0.9582930141918756, "wall_seconds": 4.7142088990003685}
{"shape": "C18", "r": 0.4625, "E": -7.410306485860109, "S": 2, "q": [0.0, 0.0], "pH": 0.9547576378708273, "wall_seconds": 5.055237933003809}
{"shape": "C18", "r": 0.475, "E": -7.441034689119293, "S": 2, "q": [0.0, 0.0], "pH": 0.9529682386177586, "wall_seconds": 4.950719400003436}
{"shape": "C18", "r": 0.4875, "E": -7.471858122919285, "S": 2, "q": [0.0, 0.0], "pH": 0.9511647752319031, "wall_seconds": 5.036359489997267}
{"shape": "C18", "r": 0.5375, "E": -7.665163632602734, "S": 1, "q": [0.0, 0.0], "pH": 0.7780255079272791, "wall_seconds": 4.627496285000234}
{"shape": "C18", "r": 0.5625, "E": -7.781222618223524, "S": 1, "q": [0.0, 0.0], "pH": 0.7931619896211652, "wall_seconds": 6.225728306002566}
{"shape": "C18", "r": 0.575, "E": -7.839616383008513, "S": 1, "q": [0.0, 0.0], "pH": 0.7976515413025839, "wall_seconds": 4.424909667002794}
{"shape": "C18", "r": 0.5875, "E": -7.898193195203776, "S": 1, "q": [0.0, 0.0], "pH": 0.8008516423761075, "wall_seconds": 4.5662246759966365}
{"shape": "C18", "r": 0.6125, "E": -8.015788641101535, "S": 1, "q": [0.0, 0.0], "pH": 0.8046301268682446, "wall_seconds": 5.221180049004033}
{"shape": "C18", "r": 0.625, "E": -8.074770421994348, "S": 1, "q": [0.0, 0.0], "pH": 0.8056066807772656, "wall_seconds": 4.696355640000547}
{"shape": "C18", "r": 0.6375, "E": -8.133857550581048, "S": 1, "q": [0.0, 0.0], "pH": 0.8061492743904821, "wall_seconds": 4.759498667001026}
{"shape": "C18", "r": 0.6875, "E": -8.419411479344195, "S": 0, "q": [0.0, 0.0], "pH": 0.454610624688884, "wall_seconds": 5.433307260005677}


```


### finish18.log

```text

{'number_of_index6_sublattices': 12, 'valid_K_compatible_shapes': 3, 'all_K_compatible_shapes_pointgroup_equivalent': True, 'control_T': [[1, 2], [0, 6]], 'geometry_metric': 'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.', 'caution': 'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'}
{"shape": "K18", "r": 0.3646865486, "E": -7.193694908907828, "S": 2, "q": [0.0, 0.0], "pH": 0.9612376979882098, "wall_seconds": 4.198199275000661}
{"shape": "K18", "r": 0.501706589, "E": -7.523368452956724, "S": 2, "q": [0.0, 0.0], "pH": 0.9403766960439041, "wall_seconds": 5.240178157000628}
{"shape": "K18", "r": 1.0, "E": -10.459935635774283, "S": 0, "q": [0.0, 0.0], "pH": 0.7112177161643203, "wall_seconds": 4.493026425996504}
{"shape": "C18", "r": 0.3346873744, "E": -7.101957147023565, "S": 2, "q": [0.0, 0.0], "pH": 0.9721401363954603, "wall_seconds": 3.8616300799985765}
{"shape": "C18", "r": 0.5042126561, "E": -7.513215938866886, "S": 2, "q": [0.0, 0.0], "pH": 0.9487319269869249, "wall_seconds": 4.228693135999492}
{"shape": "C18", "r": 0.6619467984, "E": -8.249687192353582, "S": 1, "q": [0.0, 0.0], "pH": 0.806272316001412, "wall_seconds": 4.571706938004354}
{"shape": "C18", "r": 1.0, "E": -10.594618331987615, "S": 0, "q": [0.0, 0.0], "pH": 0.5859845150624994, "wall_seconds": 4.355918618006399}
PASS 18-spin endpoints/intersections, full-space controls, checkpoint reconstruction and r=0/r=1 sanity checks.


```


### audit27.log

```text

{'r': 0.2, 'file': 'j0.2_n12_r0.npz', 'cell_n': [1, 2], 'angular': 0, 'sha256': '8431b538cf6e432ca490dcaf4f4621454a357cf4aa6efc74034869f28c987e51', 'energy': -10.1859491284774, 'fresh_H_residual': 2.8171892745612456e-14, 'fresh_spin': {'S2': 0.7499999999999984, 'spin': 0.4999999999999992, 'nearest_half_integer_spin': 0.5, 'S2_eigen_residual': 5.799277077539944e-13}, 'H2_operator': 0.8565521192794763, 'R2_operator': 0.11103541754757891, 'H2_variance': 0.3244552426810764, 'R2_variance': 0.2275668157531657, 'H2_vs_ingested_G_error': 5.041522754822836e-13, 'R2_vs_ingested_G_error': 5.045311390894369e-12}
{'r': 0.35, 'file': 'p0.350000000_n00_r0_e0.npz', 'cell_n': [0, 0], 'angular': 0, 'sha256': '8357f49e3171a933805902b3fbaf6537c66a240916a48e7898d46e6f29ee0493', 'energy': -10.472017603235424, 'fresh_H_residual': 5.0007236876320586e-11, 'fresh_spin': {'S2': 15.75, 'spin': 3.5, 'nearest_half_integer_spin': 3.5, 'S2_eigen_residual': 2.582016572957432e-10}, 'H2_operator': 24.154209930468383, 'R2_operator': 2.006098314268092, 'H2_variance': 5.36598968979024, 'R2_variance': 0.10349926957098798, 'H2_vs_ingested_G_error': 2.4016344468691386e-12, 'R2_vs_ingested_G_error': 9.907630271754897e-13}
{'r': 0.5, 'file': 'j0.5_n00_r3.npz', 'cell_n': [0, 0], 'angular': 3, 'sha256': '4161855948d997f39112aab67de0d1349517f673c0f4af6ffa7ea4d7ff6dec27', 'energy': -11.083587020439092, 'fresh_H_residual': 2.988789654114859e-14, 'fresh_spin': {'S2': 8.750000000000002, 'spin': 2.5000000000000004, 'nearest_half_integer_spin': 2.5, 'S2_eigen_residual': 4.557377136942322e-13}, 'H2_operator': 22.83533605548076, 'R2_operator': 5.4093722577915955, 'H2_variance': 17.231446994407065, 'R2_variance': 2.2943713474130014, 'H2_vs_ingested_G_error': 1.55964130499342e-12, 'R2_vs_ingested_G_error': 9.308109838457312e-13}
{'r': 0.6, 'file': 'p0.600000000_n00_r0_e0.npz', 'cell_n': [0, 0], 'angular': 0, 'sha256': '0cc88450f6b7c88a421da96cdbfff890f86b4e1f04f73c92b8769f38673392bf', 'energy': -11.686276944519516, 'fresh_H_residual': 6.624621382213314e-12, 'fresh_spin': {'S2': 3.7499999999999982, 'spin': 1.4999999999999996, 'nearest_half_integer_spin': 1.5, 'S2_eigen_residual': 2.4840992260211862e-11}, 'H2_operator': 20.664009650567625, 'R2_operator': 9.588262313607583, 'H2_variance': 33.55242202586777, 'R2_variance': 12.209911837863615, 'H2_vs_ingested_G_error': 5.613287612504791e-13, 'R2_vs_ingested_G_error': 4.583000645652646e-13}
{'r': 0.7, 'file': 'p0.700000000_n00_r3_e0.npz', 'cell_n': [0, 0], 'angular': 3, 'sha256': '44857ffeee5bfa4bb3447831e5bee6b3ee409aa4be20f63336a0d275e41e2470', 'energy': -12.384408843832794, 'fresh_H_residual': 3.999233905988414e-12, 'fresh_spin': {'S2': 0.7500000000000002, 'spin': 0.5, 'nearest_half_integer_spin': 0.5, 'S2_eigen_residual': 1.0618211309628613e-11}, 'H2_operator': 18.23907874754623, 'R2_operator': 14.425800571638973, 'H2_variance': 46.45798954055692, 'R2_variance': 33.913742149688204, 'H2_vs_ingested_G_error': 7.105427357601002e-13, 'R2_vs_ingested_G_error': 6.004086117172847e-13}


```


### pilot36.log

```text

{'N': 36, 'r': 0.3, 'number_of_down_spins': 4, 'Sz_abs': 14, 'full_sector_dim': 58905, 'translation_Gamma_dim': 4947, 'basis_generation_seconds': 0.0009802590066101402, 'median_H_action_seconds': 0.0008149900022544898, 'H_action_samples_seconds': [0.00107136600126978, 0.0019513310035108589, 0.0007905609963927418, 0.0008149900022544898, 0.0007105070035322569], 'Hermiticity_error': 1.0007415106216802e-16, 'metadata_bytes_measured_arrays': 378624, 'complex_vector_bytes': 79152, 'peak_RSS_MiB': 48.203125, 'status': 'MEASURED_HIGH_SZ_PILOT; NOT_A_GROUND_STATE_CALCULATION'}
{'N': 36, 'r': 0.3, 'number_of_down_spins': 6, 'Sz_abs': 12, 'full_sector_dim': 1947792, 'translation_Gamma_dim': 162534, 'basis_generation_seconds': 0.05905545400310075, 'median_H_action_seconds': 0.048986677000357304, 'H_action_samples_seconds': [0.04832959300256334, 0.049606336004217155, 0.05064001400023699, 0.04750860399508383, 0.048986677000357304], 'Hermiticity_error': 5.612607933133757e-17, 'metadata_bytes_measured_arrays': 12502038, 'complex_vector_bytes': 2600544, 'peak_RSS_MiB': 72.3125, 'status': 'MEASURED_HIGH_SZ_PILOT; NOT_A_GROUND_STATE_CALCULATION'}


```


### tests.log

```text

{'status': 'PASS', 'independent_12spin_full_vs_all_translation_sectors': [{'r': 0.37, 'q': [0.0, 0.0], 'dimension': 246, 'full_vs_reduced_H_error': 2.906593245373189e-15, 'SU2_commutator_error': 9.590360962329457e-15, 'isometry_error': 1.2161883888976235e-15}, {'r': 0.37, 'q': [0.0, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.1793037611233078e-15, 'SU2_commutator_error': 6.618819534863394e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 0.37, 'q': [0.5, 0.0], 'dimension': 226, 'full_vs_reduced_H_error': 2.184837125596946e-15, 'SU2_commutator_error': 6.6566789665182284e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 0.37, 'q': [0.5, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.7010225238884088e-15, 'SU2_commutator_error': 7.691147285854097e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.0, 0.0], 'dimension': 246, 'full_vs_reduced_H_error': 3.5457679953657254e-15, 'SU2_commutator_error': 1.432492926783872e-14, 'isometry_error': 1.2161883888976235e-15}, {'r': 1.0, 'q': [0.0, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.878396735277038e-15, 'SU2_commutator_error': 1.040374687250249e-14, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.5, 0.0], 'dimension': 226, 'full_vs_reduced_H_error': 2.86795700357363e-15, 'SU2_commutator_error': 8.454111100346483e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.5, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 4.01981425250324e-15, 'SU2_commutator_error': 1.1755013129905261e-14, 'isometry_error': 7.021666937153402e-16}], 'r1_independent_primitive_NN_graph_checks': [{'N': 18, 'triangular_NN_edges': 54}, {'N': 18, 'triangular_NN_edges': 54}, {'N': 27, 'triangular_NN_edges': 81}, {'N': 36, 'triangular_NN_edges': 108}], 'full_spectrum_max_error_last_case': 4.618527782440651e-14}


```


### test36pilot.log

```text

PASS 1.2052890751028591e-14


```


### final_reproduction.log

```text

{'number_of_index6_sublattices': 12, 'valid_K_compatible_shapes': 3, 'all_K_compatible_shapes_pointgroup_equivalent': True, 'control_T': [[1, 2], [0, 6]], 'geometry_metric': 'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.', 'caution': 'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'}
PASS cached verified 27-spin sources, energy and spin-sum checks.
{'status': 'PASS', 'independent_12spin_full_vs_all_translation_sectors': [{'r': 0.37, 'q': [0.0, 0.0], 'dimension': 246, 'full_vs_reduced_H_error': 2.906593245373189e-15, 'SU2_commutator_error': 9.590360962329457e-15, 'isometry_error': 1.2161883888976235e-15}, {'r': 0.37, 'q': [0.0, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.1793037611233078e-15, 'SU2_commutator_error': 6.618819534863394e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 0.37, 'q': [0.5, 0.0], 'dimension': 226, 'full_vs_reduced_H_error': 2.184837125596946e-15, 'SU2_commutator_error': 6.6566789665182284e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 0.37, 'q': [0.5, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.7010225238884088e-15, 'SU2_commutator_error': 7.691147285854097e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.0, 0.0], 'dimension': 246, 'full_vs_reduced_H_error': 3.5457679953657254e-15, 'SU2_commutator_error': 1.432492926783872e-14, 'isometry_error': 1.2161883888976235e-15}, {'r': 1.0, 'q': [0.0, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 2.878396735277038e-15, 'SU2_commutator_error': 1.040374687250249e-14, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.5, 0.0], 'dimension': 226, 'full_vs_reduced_H_error': 2.86795700357363e-15, 'SU2_commutator_error': 8.454111100346483e-15, 'isometry_error': 7.021666937153402e-16}, {'r': 1.0, 'q': [0.5, 0.5], 'dimension': 226, 'full_vs_reduced_H_error': 4.01981425250324e-15, 'SU2_commutator_error': 1.1755013129905261e-14, 'isometry_error': 7.021666937153402e-16}], 'r1_independent_primitive_NN_graph_checks': [{'N': 18, 'triangular_NN_edges': 54}, {'N': 18, 'triangular_NN_edges': 54}, {'N': 27, 'triangular_NN_edges': 81}, {'N': 36, 'triangular_NN_edges': 108}], 'full_spectrum_max_error_last_case': 4.618527782440651e-14}
PASS 1.2052890751028591e-14
{'status': 'PASS', 'number_of_18spin_couplings_and_shapes': 122, 'number_of_18spin_momentum_minima': 732, 'number_of_18spin_returned_eigenpairs': 3660, 'max_ground_H_residual': 4.691651745880047e-14, 'max_all_eigenpair_H_residual': 2.2006568501437907e-12, 'max_ground_S2_residual': 4.907653602423959e-10, 'max_all_eigenpair_S2_residual': 1.1321672950292963e-09, 'max_18spin_total_spin_sum_rule_error': 1.857541898075965e-14, 'max_27spin_total_spin_sum_rule_error': 5.320854867818525e-12, 'max_18spin_energy_from_G_error': 1.4210854715202004e-14, 'note': 'Global ground-state coverage is at all sampled points. Numerical scans cannot rule out arbitrarily narrow re-entrant plateaus between unchanged endpoints.'}
PASS 122 fresh checkpoint checks and sector completeness.
PNG/SVG/PDF figures written.
Physics report written.
PASS reproducible outputs and fresh checkpoint verification.


```


### final_tables.log

```text

{'number_of_index6_sublattices': 12, 'valid_K_compatible_shapes': 3, 'all_K_compatible_shapes_pointgroup_equivalent': True, 'control_T': [[1, 2], [0, 6]], 'geometry_metric': 'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.', 'caution': 'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'}
PASS cached verified 27-spin sources, energy and spin-sum checks.
{'status': 'PASS', 'number_of_18spin_couplings_and_shapes': 122, 'number_of_18spin_momentum_minima': 732, 'number_of_18spin_returned_eigenpairs': 3660, 'max_ground_H_residual': 4.691651745880047e-14, 'max_all_eigenpair_H_residual': 2.2006568501437907e-12, 'max_ground_S2_residual': 4.907653602423959e-10, 'max_all_eigenpair_S2_residual': 1.1321672950292963e-09, 'max_18spin_total_spin_sum_rule_error': 1.857541898075965e-14, 'max_27spin_total_spin_sum_rule_error': 5.320854867818525e-12, 'max_18spin_energy_from_G_error': 1.4210854715202004e-14, 'note': 'Global ground-state coverage is at all sampled points. Numerical scans cannot rule out arbitrarily narrow re-entrant plateaus between unchanged endpoints.'}
PASS 122 fresh checkpoint checks and sector completeness.
PNG/SVG/PDF figures written.
Physics report written.
PASS reproducible outputs and fresh checkpoint verification.


```
