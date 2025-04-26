      function besy1 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension by1cs(14), bm1cs(21), bth1cs(24)
      external besj1, csevl, inits, r1mach
c
c series for by1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.87e-18
c                                         log weighted error  17.73
c                               significant figures required  17.83
c                                    decimal places required  18.30
c
      data by1 cs( 1) /    .0320804710 0611908629e0 /
      data by1 cs( 2) /   1.2627078974 33500450e0 /
      data by1 cs( 3) /    .0064999618 9992317500e0 /
      data by1 cs( 4) /   -.0893616452 8860504117e0 /
      data by1 cs( 5) /    .0132508812 2175709545e0 /
      data by1 cs( 6) /   -.0008979059 1196483523e0 /
      data by1 cs( 7) /    .0000364736 1487958306e0 /
      data by1 cs( 8) /   -.0000010013 7438166600e0 /
      data by1 cs( 9) /    .0000000199 4539657390e0 /
      data by1 cs(10) /   -.0000000003 0230656018e0 /
      data by1 cs(11) /    .0000000000 0360987815e0 /
      data by1 cs(12) /   -.0000000000 0003487488e0 /
      data by1 cs(13) /    .0000000000 0000027838e0 /
      data by1 cs(14) /   -.0000000000 0000000186e0 /
c
c series for bm1        on the interval  0.          to  6.25000d-02
c                                        with weighted error   5.61e-17
c                                         log weighted error  16.25
c                               significant figures required  14.97
c                                    decimal places required  16.91
c
      data bm1 cs( 1) /    .1047362510 931285e0 /
      data bm1 cs( 2) /    .0044244389 3702345e0 /
      data bm1 cs( 3) /   -.0000566163 9504035e0 /
      data bm1 cs( 4) /    .0000023134 9417339e0 /
      data bm1 cs( 5) /   -.0000001737 7182007e0 /
      data bm1 cs( 6) /    .0000000189 3209930e0 /
      data bm1 cs( 7) /   -.0000000026 5416023e0 /
      data bm1 cs( 8) /    .0000000004 4740209e0 /
      data bm1 cs( 9) /   -.0000000000 8691795e0 /
      data bm1 cs(10) /    .0000000000 1891492e0 /
      data bm1 cs(11) /   -.0000000000 0451884e0 /
      data bm1 cs(12) /    .0000000000 0116765e0 /
      data bm1 cs(13) /   -.0000000000 0032265e0 /
      data bm1 cs(14) /    .0000000000 0009450e0 /
      data bm1 cs(15) /   -.0000000000 0002913e0 /
      data bm1 cs(16) /    .0000000000 0000939e0 /
      data bm1 cs(17) /   -.0000000000 0000315e0 /
      data bm1 cs(18) /    .0000000000 0000109e0 /
      data bm1 cs(19) /   -.0000000000 0000039e0 /
      data bm1 cs(20) /    .0000000000 0000014e0 /
      data bm1 cs(21) /   -.0000000000 0000005e0 /
c
c series for bth1       on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.10e-17
c                                         log weighted error  16.39
c                               significant figures required  15.96
c                                    decimal places required  17.08
c
      data bth1cs( 1) /    .7406014102 6313850e0 /
      data bth1cs( 2) /   -.0045717556 59637690e0 /
      data bth1cs( 3) /    .0001198185 10964326e0 /
      data bth1cs( 4) /   -.0000069645 61891648e0 /
      data bth1cs( 5) /    .0000006554 95621447e0 /
      data bth1cs( 6) /   -.0000000840 66228945e0 /
      data bth1cs( 7) /    .0000000133 76886564e0 /
      data bth1cs( 8) /   -.0000000024 99565654e0 /
      data bth1cs( 9) /    .0000000005 29495100e0 /
      data bth1cs(10) /   -.0000000001 24135944e0 /
      data bth1cs(11) /    .0000000000 31656485e0 /
      data bth1cs(12) /   -.0000000000 08668640e0 /
      data bth1cs(13) /    .0000000000 02523758e0 /
      data bth1cs(14) /   -.0000000000 00775085e0 /
      data bth1cs(15) /    .0000000000 00249527e0 /
      data bth1cs(16) /   -.0000000000 00083773e0 /
      data bth1cs(17) /    .0000000000 00029205e0 /
      data bth1cs(18) /   -.0000000000 00010534e0 /
      data bth1cs(19) /    .0000000000 00003919e0 /
      data bth1cs(20) /   -.0000000000 00001500e0 /
      data bth1cs(21) /    .0000000000 00000589e0 /
      data bth1cs(22) /   -.0000000000 00000237e0 /
      data bth1cs(23) /    .0000000000 00000097e0 /
      data bth1cs(24) /   -.0000000000 00000040e0 /
c
c twodpi = 2.0/pi
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data nty1, ntm1, ntth1, xmin, xsml, xmax / 3*0, 3*0./
c
      if (nty1.ne.0) go to 10
      nty1 = inits (by1cs, 14, 0.1*r1mach(3))
      ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
      ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
      xmin = 1.571*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01)
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   if (x.le.0.) call seteru (29hbesy1   x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0) go to 20
c
      if (x.lt.xmin) call seteru (31hbesy1   x so small y1 overflows,
     1  31, 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besy1 = twodpi*alog(0.5*x)*besj1(x) +
     1  (0.5 + csevl (.125*y-1., by1cs, nty1))/x
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hbesy1   no precision because x is big, 37, 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(x)
      theta = x - 3.0*pi4 + csevl (z, bth1cs, ntth1) / x
      besy1 = ampl * sin (theta)
c
      return
      end


      function besy0 (x)
c august 1980 version.  w. fullerton, c3, los alamos scientific lab.
      dimension by0cs(13), bm0cs(21), bth0cs(24)
      external besj0, csevl, inits, r1mach
c
c series for by0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.20e-17
c                                         log weighted error  16.92
c                               significant figures required  16.15
c                                    decimal places required  17.48
c
      data by0 cs( 1) /   -.0112778393 92865573e0 /
      data by0 cs( 2) /   -.1283452375 6042035e0 /
      data by0 cs( 3) /   -.1043788479 9794249e0 /
      data by0 cs( 4) /    .0236627491 83969695e0 /
      data by0 cs( 5) /   -.0020903916 47700486e0 /
      data by0 cs( 6) /    .0001039754 53939057e0 /
      data by0 cs( 7) /   -.0000033697 47162423e0 /
      data by0 cs( 8) /    .0000000772 93842676e0 /
      data by0 cs( 9) /   -.0000000013 24976772e0 /
      data by0 cs(10) /    .0000000000 17648232e0 /
      data by0 cs(11) /   -.0000000000 00188105e0 /
      data by0 cs(12) /    .0000000000 00001641e0 /
      data by0 cs(13) /   -.0000000000 00000011e0 /
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
      data bm0 cs( 1) /    .0928496163 7381644e0 /
      data bm0 cs( 2) /   -.0014298770 7403484e0 /
      data bm0 cs( 3) /    .0000283057 9271257e0 /
      data bm0 cs( 4) /   -.0000014330 0611424e0 /
      data bm0 cs( 5) /    .0000001202 8628046e0 /
      data bm0 cs( 6) /   -.0000000139 7113013e0 /
      data bm0 cs( 7) /    .0000000020 4076188e0 /
      data bm0 cs( 8) /   -.0000000003 5399669e0 /
      data bm0 cs( 9) /    .0000000000 7024759e0 /
      data bm0 cs(10) /   -.0000000000 1554107e0 /
      data bm0 cs(11) /    .0000000000 0376226e0 /
      data bm0 cs(12) /   -.0000000000 0098282e0 /
      data bm0 cs(13) /    .0000000000 0027408e0 /
      data bm0 cs(14) /   -.0000000000 0008091e0 /
      data bm0 cs(15) /    .0000000000 0002511e0 /
      data bm0 cs(16) /   -.0000000000 0000814e0 /
      data bm0 cs(17) /    .0000000000 0000275e0 /
      data bm0 cs(18) /   -.0000000000 0000096e0 /
      data bm0 cs(19) /    .0000000000 0000034e0 /
      data bm0 cs(20) /   -.0000000000 0000012e0 /
      data bm0 cs(21) /    .0000000000 0000004e0 /
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
c
c twodpi = 2.0/pi
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data alnhaf / -0.6931471805 59945309e0 /
      data nty0, ntm0, ntth0, xsml, xmax / 3*0, 2*0./
c
      if (nty0.ne.0) go to 10
      nty0 = inits (by0cs, 13, 0.1*r1mach(3))
      ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
      ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   if (x.le.0.) call seteru (29hbesy0   x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besy0 = twodpi*(alnhaf+alog(x))*besj0(x) + .375 +
     1  csevl (.125*y-1., by0cs, nty0)
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hbesy0   no precision because x is big, 37, 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(x)
      theta = x - pi4 + csevl (z, bth0cs, ntth0) / x
      besy0 = ampl * sin (theta)
c
      return
      end

      function besj0 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bj0cs(13), bm0cs(21), bth0cs(24)
      external csevl, inits, r1mach
c
c series for bj0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   7.47e-18
c                                         log weighted error  17.13
c                               significant figures required  16.98
c                                    decimal places required  17.68
c
      data bj0 cs( 1) /    .1002541619 68939137e0 /
      data bj0 cs( 2) /   -.6652230077 64405132e0 /
      data bj0 cs( 3) /    .2489837034 98281314e0 /
      data bj0 cs( 4) /   -.0332527231 700357697e0 /
      data bj0 cs( 5) /    .0023114179 304694015e0 /
      data bj0 cs( 6) /   -.0000991127 741995080e0 /
      data bj0 cs( 7) /    .0000028916 708643998e0 /
      data bj0 cs( 8) /   -.0000000612 108586630e0 /
      data bj0 cs( 9) /    .0000000009 838650793e0 /
      data bj0 cs(10) /   -.0000000000 124235515e0 /
      data bj0 cs(11) /    .0000000000 001265433e0 /
      data bj0 cs(12) /   -.0000000000 000010619e0 /
      data bj0 cs(13) /    .0000000000 000000074e0 /
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
      data bm0 cs( 1) /    .0928496163 7381644e0 /
      data bm0 cs( 2) /   -.0014298770 7403484e0 /
      data bm0 cs( 3) /    .0000283057 9271257e0 /
      data bm0 cs( 4) /   -.0000014330 0611424e0 /
      data bm0 cs( 5) /    .0000001202 8628046e0 /
      data bm0 cs( 6) /   -.0000000139 7113013e0 /
      data bm0 cs( 7) /    .0000000020 4076188e0 /
      data bm0 cs( 8) /   -.0000000003 5399669e0 /
      data bm0 cs( 9) /    .0000000000 7024759e0 /
      data bm0 cs(10) /   -.0000000000 1554107e0 /
      data bm0 cs(11) /    .0000000000 0376226e0 /
      data bm0 cs(12) /   -.0000000000 0098282e0 /
      data bm0 cs(13) /    .0000000000 0027408e0 /
      data bm0 cs(14) /   -.0000000000 0008091e0 /
      data bm0 cs(15) /    .0000000000 0002511e0 /
      data bm0 cs(16) /   -.0000000000 0000814e0 /
      data bm0 cs(17) /    .0000000000 0000275e0 /
      data bm0 cs(18) /   -.0000000000 0000096e0 /
      data bm0 cs(19) /    .0000000000 0000034e0 /
      data bm0 cs(20) /   -.0000000000 0000012e0 /
      data bm0 cs(21) /    .0000000000 0000004e0 /
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
c
      data pi4 / 0.7853981633 9744831e0 /
      data ntj0, ntm0, ntth0, xsml, xmax / 3*0, 2*0./
c
      if (ntj0.ne.0) go to 10
      ntj0 = inits (bj0cs, 13, 0.1*r1mach(3))
      ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
      ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj0 = 1.0
      if (y.gt.xsml) besj0 = csevl (.125*y*y-1., bj0cs, ntj0)
      return
c
 20   if (y.gt.xmax) call seteru (
     1  42hbesj0   no precision because abs(x) is big, 42, 1, 2)
c
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(y)
      theta = y - pi4 + csevl (z, bth0cs, ntth0) / y
      besj0 = ampl * cos (theta)
c
      return
      end

      function besj1 (x)
c sept 1983 edition.  w. fullerton, c3, los alamos scientific lab.
      dimension bj1cs(12), bm1cs(21), bth1cs(24)
      external csevl, inits, r1mach
c
c series for bj1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   4.48e-17
c                                         log weighted error  16.35
c                               significant figures required  15.77
c                                    decimal places required  16.89
c
      data bj1 cs( 1) /   -.1172614151 3332787e0 /
      data bj1 cs( 2) /   -.2536152183 0790640e0 /
      data bj1 cs( 3) /    .0501270809 84469569e0 /
      data bj1 cs( 4) /   -.0046315148 09625081e0 /
      data bj1 cs( 5) /    .0002479962 29415914e0 /
      data bj1 cs( 6) /   -.0000086789 48686278e0 /
      data bj1 cs( 7) /    .0000002142 93917143e0 /
      data bj1 cs( 8) /   -.0000000039 36093079e0 /
      data bj1 cs( 9) /    .0000000000 55911823e0 /
      data bj1 cs(10) /   -.0000000000 00632761e0 /
      data bj1 cs(11) /    .0000000000 00005840e0 /
      data bj1 cs(12) /   -.0000000000 00000044e0 /
c
c series for bm1        on the interval  0.          to  6.25000d-02
c                                        with weighted error   5.61e-17
c                                         log weighted error  16.25
c                               significant figures required  14.97
c                                    decimal places required  16.91
c
      data bm1 cs( 1) /    .1047362510 931285e0 /
      data bm1 cs( 2) /    .0044244389 3702345e0 /
      data bm1 cs( 3) /   -.0000566163 9504035e0 /
      data bm1 cs( 4) /    .0000023134 9417339e0 /
      data bm1 cs( 5) /   -.0000001737 7182007e0 /
      data bm1 cs( 6) /    .0000000189 3209930e0 /
      data bm1 cs( 7) /   -.0000000026 5416023e0 /
      data bm1 cs( 8) /    .0000000004 4740209e0 /
      data bm1 cs( 9) /   -.0000000000 8691795e0 /
      data bm1 cs(10) /    .0000000000 1891492e0 /
      data bm1 cs(11) /   -.0000000000 0451884e0 /
      data bm1 cs(12) /    .0000000000 0116765e0 /
      data bm1 cs(13) /   -.0000000000 0032265e0 /
      data bm1 cs(14) /    .0000000000 0009450e0 /
      data bm1 cs(15) /   -.0000000000 0002913e0 /
      data bm1 cs(16) /    .0000000000 0000939e0 /
      data bm1 cs(17) /   -.0000000000 0000315e0 /
      data bm1 cs(18) /    .0000000000 0000109e0 /
      data bm1 cs(19) /   -.0000000000 0000039e0 /
      data bm1 cs(20) /    .0000000000 0000014e0 /
      data bm1 cs(21) /   -.0000000000 0000005e0 /
c
c series for bth1       on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.10e-17
c                                         log weighted error  16.39
c                               significant figures required  15.96
c                                    decimal places required  17.08
c
      data bth1cs( 1) /    .7406014102 6313850e0 /
      data bth1cs( 2) /   -.0045717556 59637690e0 /
      data bth1cs( 3) /    .0001198185 10964326e0 /
      data bth1cs( 4) /   -.0000069645 61891648e0 /
      data bth1cs( 5) /    .0000006554 95621447e0 /
      data bth1cs( 6) /   -.0000000840 66228945e0 /
      data bth1cs( 7) /    .0000000133 76886564e0 /
      data bth1cs( 8) /   -.0000000024 99565654e0 /
      data bth1cs( 9) /    .0000000005 29495100e0 /
      data bth1cs(10) /   -.0000000001 24135944e0 /
      data bth1cs(11) /    .0000000000 31656485e0 /
      data bth1cs(12) /   -.0000000000 08668640e0 /
      data bth1cs(13) /    .0000000000 02523758e0 /
      data bth1cs(14) /   -.0000000000 00775085e0 /
      data bth1cs(15) /    .0000000000 00249527e0 /
      data bth1cs(16) /   -.0000000000 00083773e0 /
      data bth1cs(17) /    .0000000000 00029205e0 /
      data bth1cs(18) /   -.0000000000 00010534e0 /
      data bth1cs(19) /    .0000000000 00003919e0 /
      data bth1cs(20) /   -.0000000000 00001500e0 /
      data bth1cs(21) /    .0000000000 00000589e0 /

      data bth1cs(22) /   -.0000000000 00000237e0 /
      data bth1cs(23) /    .0000000000 00000097e0 /
      data bth1cs(24) /   -.0000000000 00000040e0 /
c
c
      data pi4 / 0.7853981633 9744831e0 /
      data ntj1, ntm1, ntth1, xsml, xmin, xmax / 3*0, 3*0./
c
      if (ntj1.ne.0) go to 10
      ntj1 = inits (bj1cs, 12, 0.1*r1mach(3))
      ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
      ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (8.0*r1mach(3))
      xmin = 2.0*r1mach(1)
      xmax = 1.0/r1mach(4)
c
 10   y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj1 = 0.0
      if (y.eq.0.0) return
      if (y.le.xmin) call seteru (
     1  37hbesj1   abs(x) so small j1 underflows, 37, 1, 0)
      if (y.gt.xmin) besj1 = 0.5*x
      if (y.gt.xsml) besj1 = x * (.25 + csevl(.125*y*y-1., bj1cs, ntj1))
      return
c
 20   if (y.gt.xmax) call seteru (
     1  42hbesj1   no precision because abs(x) is big, 42, 2, 2)
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(y)
      theta = y - 3.0*pi4 + csevl (z, bth1cs, ntth1) / y
      besj1 = sign (ampl, x) * cos (theta)
c
      return
      end


      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C                 FOR FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
C                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR I1MACH.
C
      INTEGER IMACH(16), OUTPUT, SANITY, SMALL(2), I
C/6S
C/7S
      SAVE IMACH, SANITY
C/
      REAL RMACH
C
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   47 /
C      DATA IMACH(12) / -929 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   94 /
C      DATA IMACH(15) / -929 /
C      DATA IMACH(16) / 1069 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   47 /
C      DATA IMACH(12) / -929 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   94 /
C      DATA IMACH(15) / -929 /
C      DATA IMACH(16) / 1069 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA IMACH( 1) /            1 /
C      DATA IMACH( 2) /            1 /
C      DATA IMACH( 3) /            2 /
C      DATA IMACH( 4) /            1 /
C      DATA IMACH( 5) /           32 /
C      DATA IMACH( 6) /            4 /
C      DATA IMACH( 7) /            2 /
C      DATA IMACH( 8) /           31 /
C      DATA IMACH( 9) / :17777777777 /
C      DATA IMACH(10) /            2 /
C      DATA IMACH(11) /           23 /
C      DATA IMACH(12) /         -127 /
C      DATA IMACH(13) /         +127 /
C      DATA IMACH(14) /           47 /
C      DATA IMACH(15) /       -32895 /
C      DATA IMACH(16) /       +32637 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR VAX.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C  ***  ISSUE STOP 775 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532
     *           .AND. SMALL(2) .EQ. -448790528)
     *       .OR.     (SMALL(2) .EQ. 1117925532
     *           .AND. SMALL(1) .EQ. -448790528)) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,*)'Adjust autodoubled I1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine'
               WRITE(*,*)'and setting IMACH(I) = IMACH(I+3) for'
               WRITE(*,*)'I = 11, 12, and 13.'
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*               *** IBM ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SANITY = 987
            ELSE
               WRITE(*,*)'Adjust I1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
               STOP 777
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
         SANITY = 987
         END IF
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 30
C
      I1MACH = IMACH(I)
C/6S
C/7S
      IF (I .EQ. 6) I1MACH = 1
C/
      RETURN
C
 30   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
C
*     CALL FDUMP
C
      STOP
C
* /* C source for I1MACH -- remove the * in column 1 */
* /* Note that some values may need changing -- see the comments below. */
*#include <stdio.h>
*#include <float.h>
*#include <limits.h>
*#include <math.h>
*
*long i1mach_(long *i)
*{
*	switch(*i){
*	  case 1:  return 5;	/* standard input  unit -- may need changing */
*	  case 2:  return 6;	/* standard output unit -- may need changing */
*	  case 3:  return 7;	/* standard punch  unit -- may need changing */
*	  case 4:  return 0;	/* standard error  unit -- may need changing */
*	  case 5:  return 32;	/* bits per integer -- may need changing */
*	  case 6:  return 1;	/* Fortran 77 value: 1 character */
*	  			/*    per character storage unit */
*	  case 7:  return 2;	/* base for integers -- may need changing */
*	  case 8:  return 31;	/* digits of integer base -- may need changing */
*	  case 9:  return LONG_MAX;
*	  case 10: return FLT_RADIX;
*	  case 11: return FLT_MANT_DIG;
*	  case 12: return FLT_MIN_EXP;
*	  case 13: return FLT_MAX_EXP;
*	  case 14: return DBL_MANT_DIG;
*	  case 15: return DBL_MIN_EXP;
*	  case 16: return DBL_MAX_EXP;
*	  }
*
*	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
      data iunflo / 0 /
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end

      REAL FUNCTION R1MACH(I)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR R1MACH.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, I
C/6S
C/7S
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C      DATA SMALL(1) /     8388608 /
C      DATA LARGE(1) /  2139095039 /
C      DATA RIGHT(1) /   864026624 /
C      DATA DIVER(1) /   872415232 /
C      DATA LOG10(1) /  1050288283 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00564000000000000000B /
C      DATA RMACH(2) / 37767777777777777776B /
C      DATA RMACH(3) / 16414000000000000000B /
C      DATA RMACH(4) / 16424000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / O"00564000000000000000" /
C      DATA RMACH(2) / O"37767777777777777776" /
C      DATA RMACH(3) / O"16414000000000000000" /
C      DATA RMACH(4) / O"16424000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /, SC/987/
C
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532
     *          .AND. SMALL(2) .EQ. -448790528) THEN
*              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532
     *          .AND. SMALL(1) .EQ. -448790528) THEN
*              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,*)'Adjust autodoubled R1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*              *** IBM ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE
               WRITE(*,*)'Adjust R1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
               STOP 777
               END IF
            END IF
         SC = 987
         END IF
C
C  ***  ISSUE STOP 776 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (RMACH(4) .GE. 1.0) STOP 776
*C/6S
*C     IF (I .LT. 1  .OR.  I .GT. 5)
*C    1   CALL SETERR(24HR1MACH - I OUT OF BOUNDS,24,1,2)
*C/7S
*      IF (I .LT. 1  .OR.  I .GT. 5)
*     1   CALL SETERR('R1MACH - I OUT OF BOUNDS',24,1,2)
*C/
C
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)
      RETURN
C
* /* C source for R1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return FLT_MIN;
*	  case 2: return FLT_MAX;
*	  case 3: return FLT_EPSILON/FLT_RADIX;
*	  case 4: return FLT_EPSILON;
*	  case 5: return log10(FLT_RADIX);
*	  }
*
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END

      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end

      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end

      function inits (os, nos, eta)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c
      dimension os(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinits   number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinits   eta may be too small, 28,
     1  1, 2)
      inits = i
c
      return
      end

      function csevl (x, cs, n)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c             input arguments --
c x      value at which the series is to be evaluated.
c cs     array of n terms of a chebyshev series.  in eval-
c        uating cs, only half the first coef is summed.
c n      number of terms in array cs.
c
      dimension cs(1)
c
      if (n.lt.1) call seteru (28hcsevl   number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hcsevl   number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
     1  25hcsevl   x outside (-1,+1), 25, 1, 1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
      SUBROUTINE FDUMP
C  THIS IS A DUMMY ROUTINE TO BE SENT OUT ON
C  THE PORT SEDIT TAPE
C
      RETURN
      END
      subroutine eprint
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end

      subroutine e9rint(messg,nw,nerr,save)
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
      external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
 

