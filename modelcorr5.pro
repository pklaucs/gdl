pro modelcorr5, lambda0, deltalam, dt, trials, sxcorr, countrate, $
    FIXED = fixed, SEED = seed, $
    SNR_OUT = snr_out, RAND_OUT = rand_out, TEST=test

;git test
;git test
    ; run with modelcorr4, 600, 0.1, 0, 1000, xc, 1.0e9
    ; enter lambda0, deltalam in nm
    ; NOTE: for count rates, use mag 0 star --> 10300 photons/s/cm^2/nm
    ; (true for 540 nm, B-V=0)
    ; see B. Cameron Reed, JRAS Canada, Vol. 87, p. 123.

    ; In this version, each time step is 50ps, deadtime dt = 1540 (77ns)
    ; "FRAMES" ARE 409.6ns
    ; NOTE: routine plots xcorr shifted by 100 pix so that correlation peak
    ; can be clearly seen (not chopped at origin with left half wrapped
    ; around to pix 8192)

    ; commented out all dnf operations, because we are not using dnf for deadtime -paul (06/03/2020)
    ; note dnf stands for did not finish

    T = systime(1)
    tdim = 512 ; NOTE: use only 64, 512, or 1024 here for now.
    lambda0 = lambda0 * (1.0e-9)
    deltalam = deltalam * (1.0e-9)
    lamll = lambda0 - deltalam / 2.
    lamul = lambda0 + deltalam / 2.
    frequl = 3e8 / lamll ; 3e8 is speed of light
    freqll = 3e8 / lamul
    ;print, freqll, frequl
    beatpd = 1. / (frequl - freqll)
    beatpd = beatpd * 1.0e12
    if beatpd gt 50 then print, $ ; detector response time is 50 ps
        'Beat period > 50 ps! I am not prepared for that!!'
    if beatpd lt .1 then print, $
        'Beat period < .1 ps! I am not prepared for that!!'

    ;; countrate in Hz, now input above
    ;; if each time step is 50 ps, use the line below
    ;(factor has units of counts)
	   ;50e-12 is timescale of detector response
    factor = countrate*50.0e-12;;;
    ;; if each time step is 500 ps, use the line below
    ;factor = countrate*500.0e-12

    if not keyword_set(fixed) then begin
    ch1 = intarr(8192) ;; originally fltarr
    ch2 = intarr(8192) ; good # because power of 2 for fft
    ch3 = intarr(8192)
    s1 = dblarr(8192) ;; originally fltarr
    s2 = dblarr(8192)
    s2 = s1 + 1e-9
    endif

    l1 = [100, 100]
    xctmp_12 = dblarr(8192)
    xctmp_23 = dblarr(8192)
    xctmp_13 = dblarr(8192)

    mask = dblarr(tdim)

    if tdim eq 64 then width = round(tdim * 2. / beatpd)
    if tdim eq 512 then width = round(tdim / (beatpd * 10.0) * 2.) ; default
    if tdim eq 1024 then width = round(tdim / (beatpd * 20.0) * 2.)
    if width le 0 then width = 1
    if width ge tdim then begin
        width = tdim
        print, 'WARNING! tdim < width, so do not expect decent results!!'
    endif

    print, 'Beat Period = ', beatpd, ' ps,   mask width = ', width, ' pixels'
    mask(0:width - 1) = 1.0
    xcorr_12 = dblarr(8192)
    xcorr_23 = dblarr(8192)
    xcorr_13 = dblarr(8192)

    ;delay = abs(round(randomn(seed) * 819.2))
    ;delay = 800
    delay = 0
    print, delay

    ; not using dnf method for deadtime, so just commenting it out --paul
    ;dnf1 = 0
    ;dnf2 = 0
    ;dnf3 = 0
    countsdet_12 = 0
    countsdet_23 = 0
    countsdet_13 = 0

    ;; Set the seed to a fixed value. If no value is given, then set to 4357.
    if keyword_set(fixed) then begin
        if not keyword_set(seed) then seed = 4357
        ;; Initialize the random number generators.
        none = call_external('libmt19937ar.dylib', 'init_genrand_wrapper', $
                             seed)
    endif

    for j = 1, trials do begin

        if j / 100 eq j * 1.0 / 100. then print, 'Working on trial ', j, $ ; true when j=0, 100, 200
            ' of ', trials, ',  xctmp_12(100)= ', xctmp_12(100), ', xctmp_23(100)= ', xctmp_23(100),', xctmp_13(100)= ', xctmp_13(100)

        if keyword_set(fixed) then s1 = dblarr(8192)



        for j2 = 0, 8191 do begin

            ;; Use fixed numbers if fixed is set. Else, use random numbers.
            if keyword_set(fixed) then begin
                a1 = dblarr(tdim)
                none1 = call_external('libmt19937ar.dylib', $
                                      'genrand_real2_wrapper', a1, tdim)
                a1 *= (2 * !pi)
            endif $
            else a1 = randomu(seed, tdim) * 2 * !pi ; randomu(seed,tdim)--> nums between 0 and 1

            ;; Random Numbers Output
            if keyword_set(rand_out) then begin
                filename = '../data/fixed_randoms_idl_' + strtrim(j, 2) $
                           + '.dat'
                openw, lun, filename, /get_lun
                printf, lun, a1
                close, /all
            end

            ; a1 is 512 element array of random nums
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            if keyword_set(test) then begin
            ; 06/08/20 added by paul
            a2=a1
            ;inserts the first 51 elements of randomu
            a2(0:50)=randomu(seed,tdim)*2*!pi ;
            a2=complex(cos(a2),sin(a2))
            a2=fft(a2,-1)
            b2=a2*mask
            b2=fft(b2,1)
            b2=b2*conj(b2)
            endif
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


            a1 = complex(cos(a1), sin(a1)) ; uniform random numbers
            a1 = fft(a1, -1) ; -1 so forward
            b1 = a1*mask     ; get rid of high frequencies
            b1 = fft(b1, 1)  ; backward

            b1 = b1 * conj(b1) ; to get intensities
            ; the three lines below set integration interval relative to
            ; "speckle" lifetime
            if tdim eq 64 then s1(j2) = mean(b1(0:50 - 1))
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ; 06/08/20 added by Paul
            if tdim eq 512 then begin
              s1(j2) = mean(b1(0:500 - 1)) ; default tdim
              if keyword_set(test) then s2(j2) = mean(b2(0:500-1)) ; this is new
            endif
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            if tdim eq 1024 then s1(j2) = mean(b1(0:1000 - 1))
        endfor
        ; s1 is array of mean intensities at different time points
        s1 = s1 / mean(s1) * factor ; create light intensity array,
                                    ; scale s1 by mean and multiply by factor
                                    ; factor = countrate*50.0e-12

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; 06/08/20 added by Paul
        if keyword_set(test) then s2 = s2 / mean(s2) * factor
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;print, s1(0:49)
        ;s2=shift(s1,-delay)
        ;if tdim eq 64 then plot,findgen(tdim),b1
        ;if tdim eq 512 then plot,findgen(tdim)/10.,b1
        ;if tdim eq 1024 then plot,findgen(tdim)/20.,b1
        ;wait,1

        ;; Use fixed numbers if fixed is set. Otherwise, use random numbers.
        if keyword_set(fixed) then begin
            ch1 = intarr(8192)
            ch2 = intarr(8192)
            ch3 = intarr(8192)

            none2 = call_external('libmt19937ar.dylib', $
                                  'genrand_poisson_wrapper', $
                                  s1, ch1, ch2, 8192)
        endif $
        else begin
            for i = 0, 8191 do begin

              if keyword_set(test) then begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ; 06/08/20 added by Paul

                ch1(i) = randomu(seed, poisson = s1(i))

                ; if correlated signals are desired, use the following line
                ch2(i) = randomu(seed, poisson = s2(i))

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;3 Scopes;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;note using s1 here as did for ch1
                ;correlation between ch1 and ch3 should be highest
                ch3(i) = randomu(seed, poisson = s1(i))

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                ; if uncorrelated signals are desired, use the following line
                ; ch2(i) = randomu(seed,poisson=s2(i))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              endif else begin
                ch1(i) = randomu(seed, poisson = s1(i))

                ; if correlated signals are desired, use the following line
                ch2(i) = randomu(seed, poisson = s1(i))

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;3 Scopes;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ch3(i) = randomu(seed, poisson = s1(i))
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                ; if uncorrelated signals are desired, use the following line
                ; ch2(i) = randomu(seed,poisson=s2(i))

              endif
            endfor
        endelse

        ; if dnf1 gt 0 then ch1(0:dnf1 - 1) = 0
        ; if dnf2 gt 0 then ch2(0:dnf2 - 1) = 0
        ; if dnf3 gt 0 then ch3(0:dnf3 - 1) = 0

        ; dnf1 = 0
        ; dnf2 = 0
        ; dnf3 = 0

        ; no multiple counts
        for i = 0, 8191 do begin
            if ch1(i) gt 0 then ch1(i) = 1 ; no multiple counts in a bin
            if ch2(i) gt 0 then ch2(i) = 1 ; no multiple counts in a bin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;3 Scopes;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            if ch3(i) gt 0 then ch3(i) = 1 ; no multiple counts in a bin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        endfor

        ; Random Numbers Output
        if keyword_set(rand_out) then begin
            filename = '/Users/Paul/gdl/modelcorr5_output/fixed_channel1_idl_' + strtrim(j, 2) + '.dat'
            openw, lun, filename , /get_lun
            printf, lun, ch1
            close, /all

            filename = '/Users/Paul/gdl/modelcorr5_output/fixed_channel2_idl_' + strtrim(j, 2) + '.dat'
            openw, lun, filename , /get_lun
            printf, lun, ch2
            close, /all

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;3 Scopes;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            filename = '/Users/Paul/gdl/modelcorr5_output/fixed_channel3_idl_' + strtrim(j, 2) + '.dat'
            openw, lun, filename , /get_lun
            printf, lun, ch3
            close, /all
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        endif

        ; add in dead time, period of time when can't record photon arrivals
        ; after photon hits detector, some period of time till can detect again
        ; lowers snr
        if dt ge 1 then begin
            for i = 0, 8191 do begin
                ; channel 1
                ; essentially search for nonzero photon
                ; then set dt elements in array to 0 after the photon
                if ch1(i) gt 0 and i lt 8191 - dt then ch1(i + 1:i + dt) = 0
                if ch1(i) gt 0 and i ge 8191 - dt and i ne 8191 then begin
                    ch1(i + 1:8191) = 0
                    ; dnf1 = i + dt - 8191
                endif
                ; if ch1(i) gt 0 and i eq 8191 then dnf1 = dt

                ; channel 2
                if ch2(i) gt 0 and i lt 8191 - dt then ch2(i + 1:i + dt) = 0
                if ch2(i) gt 0 and i ge 8191 - dt and i ne 8191 then begin
                    ch2(i + 1:8191) = 0
                    ; dnf2 = i + dt - 8191
                endif
                ; if ch1(i) gt 0 and i eq 8191 then dnf1 = dt ; should this be dnf2?

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;3 Scopes;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ; channel 3
                 if ch3(i) gt 0 and i lt 8191 - dt then ch3(i + 1:i + dt) = 0
                 if ch3(i) gt 0 and i ge 8191 - dt and i ne 8191 then begin
                    ch3(i + 1:8191) = 0
                    ; dnf3 = i + dt - 8191
                 endif
                 if ch1(i) gt 0 and i eq 8191 then dnf3 = dt
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            endfor
        endif

        ;if j eq 1 then print, ch1(0:24)

        ;plot, ch1(0:2000)
        ;wait, 1
        d1 = fft(ch1, -1)
        d2 = fft(ch2, -1)
        d3 = fft(ch3, -1)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;correlations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;corr d1 & d2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        e_12 = d1 * conj(d2)
        ;if j eq 1 then print, real_part(e(4096:8191))
        e_12 = fft(e_12, 1) * 8192.
        ; stacking correlations
        xcorr_12 = xcorr_12 + double(e_12) ;; cross-correlation
        xctmp_12 = shift(xcorr_12, 100) ; shift it so that can see peak better
        ; shifted by 100 so should see peak at 100
        countsdet_12 = countsdet_12 + total(ch1)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;corr d2 & d3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        e_23 = d2 * conj(d3)
        ;if j eq 1 then print, real_part(e(4096:8191))
        e_23 = fft(e_23, 1) * 8192.
        ; stacking correlations
        xcorr_23 = xcorr_23 + double(e_23) ;; cross-correlation
        xctmp_23 = shift(xcorr_23, 100) ; shift it so that can see peak better
        ; shifted by 100 so should see peak at 100
        countsdet_23 = countsdet_23 + total(ch2)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;corr d1 & d3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        e_13 = d1 * conj(d3)
        ;if j eq 1 then print, real_part(e(4096:8191))
        e_13 = fft(e_13, 1) * 8192.
        ; stacking correlations
        xcorr_13 = xcorr_13 + double(e_13) ;; cross-correlation
        xctmp_13 = shift(xcorr_13, 100) ; shift it so that can see peak better
        ; shifted by 100 so should see peak at 100
        countsdet_13 = countsdet_13 + total(ch3)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        ;cgPlot, xctmp_12(0:500), Position=[0.10, 0.10, 0.30, 0.90]
        ;cgPlot, xctmp_23(0:500), Position=[0.40, 0.10, 0.60, 0.90], /noerase
        ;cgPlot, xctmp_13(0:500), Position=[0.70, 0.10, 0.90, 0.90], /noerase

        ;l2 = [xctmp(100), xctmp(100) * 5.0]
        ;oplot, l1, l2, linestyle = 1


        cgPlot, xctmp_12(0:500), Position=[0.10, 0.10, 0.90, 0.30]
        cgPlot, xctmp_23(0:500), Position=[0.10, 0.40, 0.90, 0.60], /noerase
        cgPlot, xctmp_13(0:500), Position=[0.10, 0.70, 0.90, 0.90], /noerase

        ;cgPlot, cgDemoData(17), Position=[0.40, 0.10, 0.60, 0.90], /noerase
        ;cgPlot, cgDemoData(17), Position=[0.70, 0.10, 0.90, 0.90], /noerase


        ;plot, xctmp(0:500)
        ;l2 = [xctmp(100), xctmp(100) * 5.0]
        ;oplot, l1, l2, linestyle = 1
        ;wait, .001
        ;countsdet = countsdet + total(ch1)
        ;s2 = s1 ; using s2 if uncorrelated signals

        ;; Output S/N, mean of signal around sd of area around peak
        ;if keyword_set(snr_out) then begin

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Scope 1 & 2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            z_12 = xctmp_12(0:200)
            z_12 = z_12-mean(z_12)
            sig_12 = dblarr(122)
            sig_12(0:60) = z_12(20:80) ; grabbing area around peak for sd
            sig_12(61:121) = z_12(120:180)
            snr_12 = z_12(100) / stddev(sig_12)

            ;; Open a new file during the first trial.
            if j eq 1 then begin
                openw, lun, '/Users/Paul/gdl/modelcorr5_output/snr_12_out_idl.txt', /get_lun
                printf, lun, j, snr_12
                Free_lun, lun
            endif
            ;; Append the rest of the S/N ratios to the new file.
            if j gt 1 then begin
                openu, lun, '/Users/Paul/gdl/modelcorr5_output/snr_12_out_idl.txt', /append
                printf, lun, j, snr_12
                Free_lun, lun
            endif
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Scope 2 & 3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        z_23 = xctmp_23(0:200)
        z_23 = z_23-mean(z_23)
        sig_23 = dblarr(122)
        sig_23(0:60) = z_23(20:80) ; grabbing area around peak for sd
        sig_23(61:121) = z_23(120:180)
        snr_23 = z_23(100) / stddev(sig_23)

        ;; Open a new file during the first trial.
        if j eq 1 then begin
            openw, lun, '/Users/Paul/gdl/modelcorr5_output/snr_23_out_idl.txt', /get_lun
            printf, lun, j, snr_23
            Free_lun, lun
        endif
        ;; Append the rest of the S/N ratios to the new file.
        if j gt 1 then begin
            openu, lun, '/Users/Paul/gdl/modelcorr5_output/snr_23_out_idl.txt', /append
            printf, lun, j, snr_23
            Free_lun, lun
        endif
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Scope 1 & 3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        z_13 = xctmp_13(0:200)
        z_13 = z_13-mean(z_13)
        sig_13 = dblarr(122)
        sig_13(0:60) = z_13(20:80) ; grabbing area around peak for sd
        sig_13(61:121) = z_13(120:180)
        snr_13 = z_13(100) / stddev(sig_13)

        ;; Open a new file during the first trial.
        if j eq 1 then begin
            openw, lun, '/Users/Paul/gdl/modelcorr5_output/snr_13_out_idl.txt', /get_lun
            printf, lun, j, snr_13
            Free_lun, lun
        endif
        ;; Append the rest of the S/N ratios to the new file.
        if j gt 1 then begin
            openu, lun, '/Users/Paul/gdl/modelcorr5_output/snr_13_out_idl.txt', /append
            printf, lun, j, snr_13
            Free_lun, lun
        endif



      ;  endif
    endfor

    ;countsdet = countsdet * 1.0 / trials

    ;print, 'count rate=', countrate, ',  detected rate=', countsdet / 409.6e-9

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Unecessary?;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; commented out because why are we calculting the shifted SNR?
    ; we can just use the non-shifted snr above

    ;sxcorr = shift(xcorr, 100)
    ;plot, sxcorr(0:500)
    ;l2 = [sxcorr(100), sxcorr(100) * 5.0]
    ;oplot, l1, l2, linestyle = 1

    ;z = sxcorr(0:200)
    ;z = z - mean(z)
    ;sig = dblarr(122)
    ;sig(0:60) = z(20:80)
    ;sig(61:121) = z(120:180)
    ;snr = z(100) / stddev(sig)
    ;print, 'Signal-to-Noise Ratio: ', snr
    ;print, systime(1) - T, ' seconds'
stop

end
