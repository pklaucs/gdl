pro save_modelcorr_2, rep, constant, constant_value, trials, countrate, frames, lam, dlam, dt, snrs, sd, time

; 05/25/2020 added frame_start as argument to function
  timestamp, ts; get date and time for saved file name
  ts=repstr(ts,':','_')

  if (strmatch(rep, 'rep')) then BEGIN
    append='_const_' + constant + '_' + strtrim(constant_value,2) + '_trials_'+ strtrim(trials[0],2)
  endif else BEGIN
    append='_maxframes_' +strtrim(frames,2) + '_countrate_' + strtrim(countrate,2)
  ENDELSE

  filename= ts[0] + append $
  + '_lam_' +strtrim(lam,2) + '_dlamfwhm_'+strtrim(dlam,2)+ $
  '_dt_'+strtrim(dt,2)

  filename=repstr(filename,'.','')
  filename=repstr(filename,'+','')
  filename=filename + '.csv'


  if (strmatch(rep, 'rep')) then BEGIN
    print, 'rep'
    column_headers=['trials','constant_var','frames','countrate',$
      'lambda','filter_fwhm_nm','deadtime', 'snr','sd', 'time']
    snrs_out=string(snrs)
    s=n_elements(frames)
    n_frames=transpose(string(frames))
    time_out=string(time)

    countrate_out=transpose(string(countrate))
    s=n_elements(countrate_out)
    lam_out=string(transpose(replicate(lam,s)))
    deltalamfwhm_out=string(transpose(replicate(dlam,s)))
    sd_out=string(sd)
    dt_out=string(transpose(replicate(dt[0],s)))
    trials_out=string(transpose(replicate(trials[0],s)))
    constant_out=string(transpose(replicate(constant[0],s)))

    data_out=[trials_out, constant_out, n_frames, countrate_out, lam_out, deltalamfwhm_out, dt_out, snrs_out, sd_out, time_out]

  endif else BEGIN
    print, 'non-rep'
    snrs_out=transpose(string(snrs))
    column_headers=['nframes','countrate',$
      'lambda','filter_fwhm_nm','deadtime','snr', 'time']
    n_frames=transpose(string(indgen(frames+1)))
    countrate_out=transpose(string(replicate(countrate,frames+1)))
    lam_out=transpose(string(replicate(lam,frames+1)))
    deltalamfwhm_out=transpose(string(replicate(dlam, frames+1)))
    time_out=transpose(string(time))

    dt_out=string(transpose(replicate(dt,frames+1)))

    ; no sd to output
    data_out=[n_frames, countrate_out, lam_out, deltalamfwhm_out,dt_out, snrs_out, time_out]

  endelse




  write_csv_data, data_out, column_headers, filename=filename



end
