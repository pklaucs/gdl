pro modelcorr4_2, lambda0, deltalam, dt, frames, sxcorr, countrate, snr, rep

; 4/5/2020 paul
; added parameter snr for output

; enter lambda0, deltalam in nm
; NOTE: for count rates, use mag 0 star --> 10300 photons/s/cm^2/nm
; (true for 540 nm, B-V=0)
; see B. Cameron Reed, JRAS Canada, Vol. 87, p. 123.

; In this version, each time step is 50ps, deadtime dt = 1540 (77ns)
; "FRAMES" ARE 409.6ns
; NOTE: routine plots xcorr shifted by 100 pix so that correlation peak can
; be clearly seen (not chopped at origin with left half wrapped around to
; pix 8192)

tdim=512 ; NOTE: use only 64, 512, or 1024 here for now.
; enter wavelengths in nm
;lambda0 = 600.0 -- now entered above
;deltalam = 0.1 -- now entered above
lambda0 = lambda0*(1.0e-9) ; gives wavelength in m
deltalam = deltalam*(1.0e-9) ; line width of lambda?
lamll = lambda0-deltalam/2. ; lower limit
lamul = lambda0+deltalam/2. ; upper limit
frequl = 3e8/lamll ; upper limit
freqll = 3e8/lamul ; lower limit

; width of the filter defines the beat period?
;print, freqll, frequl
beatpd = 1./(frequl-freqll)
beatpd = beatpd*1.0e12
if beatpd gt 50 then print,'Beat period > 50 ps! I am not prepared for that!!'
if beatpd lt .1 then print,'Beat period < .1 ps! I am not prepared for that!!'

;countrate=1.0e8  ; in Hz, now input above
; if each time step is 50 ps, use the line below
factor = countrate*50.0e-12
; if each time step is 500 ps, use the line below
;factor = countrate*500.0e-12

l1=[100,100]
ch1=fltarr(8192)
ch2=fltarr(8192)
s1=fltarr(8192)
s2=fltarr(8192)
s2=s1+1e-9
xctmp=fltarr(8192)
mask = fltarr(tdim)
if tdim eq 64 then width=round(tdim*2./beatpd)
if tdim eq 512 then width=round(tdim/(beatpd*10.0)*2.)
if tdim eq 1024 then width=round(tdim/(beatpd*20.0)*2.)
if width le 0 then width = 1
if width ge tdim then begin
  width = tdim
  print,'WARNING! tdim < width, so do not expect decent results!!'
endif
print,'Beat Period = ',beatpd,' ps,   mask width = ',width,' pixels
mask(0:width-1)=1.0
xcorr=fltarr(8192)

;delay=abs(round(randomn(seed)*819.2))
;delay=800
delay=0
;print,delay

sig=fltarr(122)
sig1=fltarr(122)

dnf1=0
dnf2=0
countsdet=0
snrs=fltarr(frames+1)
time_out=fltarr(frames+1)

;z1_temp= MAKE_ARRAY(201, frames+1, /INTEGER, VALUE = 0)
;sig1_temp= MAKE_ARRAY(122, frames+1, /INTEGER, VALUE = 0)

for j=0,frames do begin
  tic
  if j/100 eq j*1.0/100. and j ne 0 then begin
    ;each new trial -- paul (every 100 frames)
    z=xctmp(0:200)
    ;print, 'xctmp(0:200): '
    ;print, z
    z=z-mean(z)
    sig(0:60)=z(20:80)
    sig(61:121)=z(120:180)
    snr=z(100)/stddev(sig)
;    print,'Signal-to-Noise Ratio: ',snr
    print,'Working on trial ',j,' of ',frames,   ',  S/N=',snr
  endif

  for j2=0,8191 do begin
    a1=randomu(seed,tdim)*2*!pi
    a1=complex(cos(a1),sin(a1))
    a1=fft(a1,-1)
    b1=a1*mask
    b1=fft(b1,1)
    b1=b1*conj(b1)
    ; the three lines below set integration interval relative to "speckle"
    ; lifetime
;    s1(j2) = mean(b1(0:50-1)) ; sets int. interval rel. to "speckle" lifetime
    if tdim eq 64 then s1(j2) = mean(b1(0:50-1))
    if tdim eq 512 then s1(j2) = mean(b1(0:500-1))
    if tdim eq 1024 then s1(j2) = mean(b1(0:1000-1))
  endfor

  s1=s1/mean(s1)*factor
;  s2=shift(s1,-delay)
;if tdim eq 64 then plot,findgen(tdim),b1
;if tdim eq 512 then plot,findgen(tdim)/10.,b1
;if tdim eq 1024 then plot,findgen(tdim)/20.,b1
;wait,1

  for i=0,8191 do begin
    ch1(i) = randomu(seed,poisson=s1(i))
    ; if correlated signals are desired, use the following line
    ch2(i) = randomu(seed,poisson=s1(i))
    ; if uncorrelated signals are desired, use the following line
;    ch2(i) = randomu(seed,poisson=s2(i))
  endfor

  if dnf1 gt 0 then ch1(0:dnf1-1)=0
  if dnf2 gt 0 then ch2(0:dnf2-1)=0
  dnf1=0
  dnf2=0

; no multiple counts
  for i=0,8191 do begin
    if ch1(i) gt 0 then ch1(i) = 1 ; no multiple counts in a bin
    if ch2(i) gt 0 then ch2(i) = 1 ; no multiple counts in a bin
  endfor

; add in dead time
  if dt ge 1 then begin
    for i=0,8191 do begin
      ; channel 1
      if ch1(i) gt 0 and i lt 8191-dt then ch1(i+1:i+dt)=0
      if ch1(i) gt 0 and i ge 8191-dt and i ne 8191 then begin
        ch1(i+1:8191)=0
        dnf1 = i+dt-8191
      endif
      if ch1(i) gt 0 and i eq 8191 then dnf1 = dt

      ; channel 2
      if ch2(i) gt 0 and i lt 8191-dt then ch2(i+1:i+dt)=0
      if ch2(i) gt 0 and i ge 8191-dt and i ne 8191 then begin
        ch2(i+1:8191)=0
        dnf2 = i+dt-8191
      endif
      if ch1(i) gt 0 and i eq 8191 then dnf1 = dt
    endfor
  endif

;plot,ch1(0:2000)
;wait,1

  d1=fft(ch1,-1)
  d2=fft(ch2,-1)
  e=d1*conj(d2)
  e=fft(e,1)*8192. ; cross correlation?

  xcorr = xcorr+float(e)
  xctmp=shift(xcorr,100)
  ;plot,xctmp(0:500) ;;;;;;;;;;;;;;;;;;;;;;;;;;;uncomment for plot
  l2=[xctmp(100),xctmp(100)*5.0]
  ;oplot,l1,l2,linestyle=1;;;;;;;;;;;;;;;;;;;;;;;;;;;uncomment for plot
  countsdet = countsdet + total(ch1)
  s2 = s1

  ;;;;;;;;;;;;;;;;;;I added this;;;;;;;;;;;;;;;;

  sxcorr = shift(xcorr,100)
  ;l2=[sxcorr(100),sxcorr(100)*5.0]
  z=sxcorr(0:200)
  ;print, 'sxcorr(0:200): '
  ;print, z
  z1=z-mean(z)
  sig1(0:60)=z1(20:80)
  sig1(61:121)=z1(120:180)
  ;print, 'sig1: ',sig1
  ;print, 'z1: ', z1
  ;z1_temp[*,j]=z1 ; added this to look at it after
  ;sig1_temp[*,j]=sig1; added this to look at it after
  snr=z1(100)/stddev(sig1)
  ;print, 'snr: ', snr
  ;print, 'j frame: '
  ;print, j
  ;print,'Signal-to-Noise Ratio: ',snr
  snrs[j]=snr
  ;print, 'frame #: ', j
  ;print, 'snr: ', snr
  ;print, snrs
  ;output snr and frame

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time_out[j]=toc()

endfor

; check if using modelcorr_rep
;if (strmatch(rep, 'rep')) then return

;;;;;;;;;;;;;;;;;;;;;;;;;MY FILE SAVE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
rep_out='null'
constant='null'
constant_value='null'
trials='null'
sd='null'
save_modelcorr_2, rep_out[0], constant[0],constant_value[0], trials[0], countrate[0], frames[0], $
  lambda0[0], deltalam[0], dt[0], snrs, sd[0], time_out
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;PART OF ORIGINAL CODE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; snr is already calced above

;countsdet = countsdet*1.0/frames
;print,'count rate=',countrate,',  detected rate=',countsdet/409.6e-9
;sxcorr = shift(xcorr,100)
;plot,sxcorr(0:500);;;;;;;;;;;;;;;;;;;;;;;;;;;uncomment for plot
;l2=[sxcorr(100),sxcorr(100)*5.0]
;oplot,l1,l2,linestyle=1;;;;;;;;;;;;;;;;;;;;;;;;;;;uncomment for plot
;z=sxcorr(0:200)
;print, 'sxcorr(0:200): '
;print, z
;z=z-mean(z)
;sig(0:60)=z(20:80)
;sig(61:121)=z(120:180)
;snr=z(100)/stddev(sig)
;print,'Signal-to-Noise Ratio: ',snr
;print,'Signal: ',z(100)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; free up MEMORY, idt it actually works when
;; looping over modelcorr- paul
;undefine, z, sig, sxcorr, ch1, ch2, s1, s2, xctmp, mask
;undefine, xcorr, z, a1, b1, d1, d2, e, l2
;undefine, dnf1, seed, sxcorr, countrate, countsdet, deltalam
;undefine, factor, frames, freqll, frequl, i, j, j2, l1, l2
;undefine, lambda0, lamll, lamul, lamul, mem, tdim, width

; calc the theoretical SNR
;n=   ;detector size (cm?)
; m=0
; F_0 = 10300 ;photons/s/cm^2/nm
; (2.512^-m)*10300*

end
