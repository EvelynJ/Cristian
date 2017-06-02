;===============================================================

function S_N_calculator,x,y,cont_array,root,galaxy_ref,x_centre,y_centre,limit
  ;*** Calculate S/N for each spaxel. Using the cont_wavelength and cont_range,
  ;*** the Signal is calculated as the mean flux value within the defined wavelength
  ;*** range, while the nise is the standard deviation in that same range. This
  ;*** calculation assumes that the SD is dominated by noise and contains n
  ;*** significant spectral features
  
  
  ;prepare text files to save the result
  S_N_array=fltarr(x,y)
  close,01
  openw,01,root+galaxy_ref+'_S_N_array.txt'
  printf,01,'##########################################################################################'
  printf,01,'#          Xpix           Ypix               Signal              noise               S/N '
  printf,01,'##########################################################################################'
  
  close,11
  openw,11,root+galaxy_ref+'_pixel_array_full.txt'
  printf,11,'##########################################################################'
  printf,11,'#          Xpix           Ypix     '
  printf,11,'##########################################################################'
  
  close,21
  openw,21,root+galaxy_ref+'_S_N_array_full.txt'
  printf,21,'##########################################################################################'
  printf,21,'#          Xpix           Ypix               Signal              noise'
  printf,21,'##########################################################################################'
  
  off=4
  
  
  ;calculate signal and noise levels
  Signal=dblarr(x*y)
  noise=dblarr(x*y)
  S_N=dblarr(x*y)
  col_arr=dblarr(x*y)
  row_arr=dblarr(x*y)
  
  n=double(-1)
  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      n+=1
      Signal[n]=mean(cont_array[column,row,*])
      noise[n]=stddev(cont_array[column,row,*])
      if Signal[n] eq 0 then begin
        ;avoid Infinity values for pixels with value 0
        Signal[n]=0.0001
        noise[n]=1
      endif
      S_N[n]=Signal[n]/noise[n]
      S_N_array[column,row]=S_N[n]
      col_arr[n]=column
      row_arr[n]=row
    endfor
  endfor
 
  ;print out results to text file
  for m=0,n_elements(col_arr)-1,1 do begin
    if S_N[m] gt limit then printf,01,col_arr[m]-x_centre,row_arr[m]-y_centre,Signal[m],noise[m],Signal[m]/noise[m],format='(2f16.5,3f20.10)'
    if S_N[m] gt 0 then printf,11,col_arr[m]-x_centre,row_arr[m]-y_centre,format='(2f16.5)'
    printf,21,col_arr[m]-x_centre,row_arr[m]-y_centre,Signal[m],noise[m],format='(2f16.5,2f20.10)'
  endfor
  
  
  close,01
  close,11
  close,21
  
  return,S_N_array
end


;===============================================================
;===============================================================
pro bin2d_display_pixels, x, y, counts, pixelSize
  COMPILE_OPT IDL2, HIDDEN
  
  ; Plots colored pixels with the same pixels size, mapping out the bins
  
  PLOT, [MIN(x)-pixelSize,MAX(x)+pixelSize], [MIN(y)-pixelSize,MAX(y)+pixelSize], $
    /NODATA, /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec', /ISO, color=cgcolor('white'),$
    xrange=[MIN(x)-pixelSize,MAX(x)+pixelSize],yrange=[MIN(y)-pixelSize,MAX(y)+pixelSize]
  x1 = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
  y1 = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
  color = bytscl(counts)
  FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]
  
END


;##############################################################
;##############################################################
;##############################################################


pro measure_kinematics,input_file
;- first, download the IDL astron library, coyote codes, ppxf and voronoi binning
;- set paths:
;     e.g. export IDL_PATH=+/home/ejohnsto/IDLWorkspace82/astron/:$IDL_PATH
;- open IDL in the command line by typing:  idl
;- type:  measure_kinematics,'kinematics_input.txt'

read_input, input_file, setup



;*** Set up directory structure for all output files
root=setup.root
galaxy_ref=setup.galaxy_ref
file=setup.file
kinematics=setup.kinematics
stellib_dir=setup.stellib_dir

directory=root

;##############################################
;# First step is to log-rebin the datacube.
;# this step ensures that as you move through 
;# the datacube, the wavelength increases in 
;# logarithmic steps instead of linearly.
;# We need this logarithmic wavelength for when 
;# we measure the kinematics
;##############################################

fits_read,directory+file+'.fits',input,h

x=sxpar(h,'NAXIS1')
y=sxpar(h,'NAXIS2')
z=sxpar(h,'NAXIS3')

;calculate wavelength solution
wavelength0=sxpar(h,'CRVAL3')
step=sxpar(h,'CD3_3')
wavelength=fltarr(sxpar(h,'NAXIS3'))
for m=0,sxpar(h,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)



if setup.log_rebin_data eq 'y' then begin
  input_counts=input
  output=fltarr(x,y,z)
  lamRange2=[wavelength[0],wavelength[-1]]

  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      temp_spec=fltarr(z)
      temp_spec[*]=input_counts[column,row,*]
      ;log10_rebin rebins the datacube to log of base 10
      log10_rebin, lamRange2, temp_spec, new_spec, logLam2, VELSCALE=velScale
      output[column,row,*]=new_spec
    endfor
  endfor
  
  ;update header values and print out the log-binned datacube
  sxaddpar,h,'CRVAL3',logLam2[0]  
  sxaddpar,h,'CDELT3',logLam2[1]-logLam2[0]
  sxaddpar,h,'CD3_3',logLam2[1]-logLam2[0]
  
  fits_write,directory+file+'_LOG.fits',output,extname='FLUX'

  temp=mrdfits(directory+file+'.fits',0,h_temp)
  sxaddpar,h_temp,'CD1_1',sxpar(h,'CD1_1')
  sxaddpar,h_temp,'CD2_2',sxpar(h,'CD2_2')

  modfits,directory+file+'_LOG.fits',0,h_temp
  modfits,directory+file+'_LOG.fits',0,h,extname='FLUX'
endif


;simple stop message if the log-binned datacube doesn't exist
if setup.log_rebin_data eq 'n' and file_test(directory+file+'_LOG.fits') eq 0 then begin
  print,'You need to logarithmically bin the '
  print,'datacube before proceding. please set'
  print,'log_rebin_data to "y" in the input file'
  stop
endif
if file_test(directory+file+'_LOG.fits') eq 1 then file=file+'_LOG'



;##############################################
;# Now we need to bin the datacube. This step
;# provides information on which spaxels need
;# to be co-added to reach the desired S/N.
;##############################################
if setup.bin_data eq 'y' then begin
  print,'##############################'
  print,'#Voronoi binning the datacube#
  print,'##############################'
  
  
  ;========================================================
  ; *** Read in basic information about IFU data cube
  fits_read,root+file+'.fits',input_IFU,header_IFU
  x=sxpar(header_IFU,'NAXIS1')
  y=sxpar(header_IFU,'NAXIS2')
  z=sxpar(header_IFU,'NAXIS3')
  wavelength0=sxpar(header_IFU,'CRVAL3')
  step=sxpar(header_IFU,'CD3_3')
  
  wavelength=fltarr(sxpar(h,'NAXIS3'))
  for m=0,sxpar(h,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)
 
  ;========================================================
  ; *** identify the elements in the array corresponding to the region over
  ; *** which to measure the S/N of each spaxel
  
  cont_wavelength=setup.cont_wavelength
  cont_range=setup.cont_range
  x_centre=setup.x_centre
  y_centre=setup.y_centre
  targetSN=setup.targetSN
  
  
  cont_wavelength_log=alog10(cont_wavelength)
  cont_low_wavelength_log=alog10(cont_wavelength-cont_range/2)
  cont_high_wavelength_log=alog10(cont_wavelength+cont_range/2)
  
  sample=where(wavelength ge cont_low_wavelength_log and wavelength le cont_high_wavelength_log)
  
  result = FILE_TEST(root+kinematics,/DIRECTORY)
  if result eq 0 then spawn,'mkdir '+root+kinematics
  
  datacube=input_IFU
  
  ;========================================================
  ; *** Calculate S/N for each element in the array
  cont_array=datacube[*,*,sample]
  limit=setup.limit           ;only include pixels with S/N above this value in the binning
  
  
  S_N_array=S_N_calculator(x,y,cont_array,root,galaxy_ref,x_centre,y_centre,limit)
  
  ;========================================================
  ; *** Apply voronoi binning to the image
  readcol,root+galaxy_ref+'_S_N_array.txt',format='d,d,d,d',Xpix,Ypix,signal,noise,comment='#',/SILENT
  

  voronoi_2d_binning, xpix, ypix, signal, noise, targetSN, $
    binNum, xnde, ynde, xBar, yBar, sn, nPixels, scale, root,galaxy_ref, /QUIET;,/PLOT
    
  ; Save to a text file the initial coordinates of each pixel together
  ; with the corresponding bin number computed by this procedure.
  ; binNum uniquely specifies the bins and for this reason it is the only
  ; number required for any subsequent calculation on the bins.
  astrolib
  forprint, xpix, ypix, binNum, TEXTOUT=root+galaxy_ref+'_voronoi_2d_binning_output.txt', $
    COMMENT='          X"              Y"           BIN_NUM'
    
  ; Print out fits file with binned spectra
  n_bins=max(binNum)+1
  binned_spec=fltarr(z,n_bins)   ;2D array for binned spectra
  
  
  
  close,03
  openw,03,root+galaxy_ref+'_bin_centers.txt'
  printf,03,'#################################################################'
  printf,03,'#       Bin       Xcentre        Ycentre        S/N'
  printf,03,'#################################################################'
  
  ;collect together the spectra in each bin
  for bin=0,n_bins-1, 1 do begin
    printf,03,bin,xBar[bin],yBar[bin],sn[bin],format='(i9,3f15.3)'
    pixels=where(binNum eq bin)       ;identify elements of array to be used for each bin
    spec=fltarr(z)
    for n=0,n_elements(pixels)-1,1 do begin
      m=pixels[n]
      x_new=xpix[m]+x_centre
      y_new=ypix[m]+y_centre
      spec[*]+=datacube[x_new,y_new,*]    ;coadd relevant spectra
    endfor
    
    binned_spec[*,bin]=spec[*]
    
  endfor
  close,03
  
  
  ;update header file and print out binned spectra to a fits file
  sxaddpar,header_IFU,'CRVAL1',sxpar(header_IFU,'CRVAL3')
  sxaddpar,header_IFU,'CD1_1',sxpar(header_IFU,'CD3_3')
  sxaddpar,header_IFU,'CDELT1',sxpar(header_IFU,'CDELT3')
  sxaddpar,header_IFU,'CRVAL2',n_elements(binned_spec[*,0])
  sxaddpar,header_IFU,'CD2_2',1
  sxaddpar,header_IFU,'CDELT2',1
  sxdelpar,header_IFU,'CRVAL3'
  sxdelpar,header_IFU,'CDELT3'
  sxdelpar,header_IFU,'CD3_3'
  
  fits_write,root+kinematics+galaxy_ref+'_binned_spectra.fits',binned_spec,header_IFU,extname='FLUX'
  mkhdr,h0,binned_spec
  sxaddpar,h0,'CRVAL1',sxpar(header_IFU,'CRVAL1')
  sxaddpar,h0,'CD1_1',sxpar(header_IFU,'CD1_1')
  sxaddpar,h0,'CDELT1',sxpar(header_IFU,'CDELT1')
  modfits,root+kinematics+galaxy_ref+'_binned_spectra.fits',0,ho,exten_no=0
  
endif 


;##############################################
;# Now we have the binned spectra, we can 
;# start measuring the kinematics using ppxf
;##############################################


if setup.measure_kinematics eq 'y' then begin
  print,'############################'
  print,'# Measuring the kienmatics #
  print,'############################'
  
  Redshift=setup.Redshift
  PA=setup.PA
  
  fits_read,root+kinematics+galaxy_ref+'_binned_spectra.fits',binned_spec,header_IFU
  z=sxpar(header_IFU,'NAXIS1')
  n_bins=sxpar(header_IFU,'NAXIS2')
  c=3e5  ;speed of light

  wavelength0=sxpar(header_IFU,'CRVAL1')
  step=sxpar(header_IFU,'CD1_1')
  wavelength=fltarr(z)
  for n=0,z-1,1 do wavelength[n]=wavelength0+n*step
  
  
  a=(Redshift+1)*(Redshift+1)
  Velocity=c*((a-1)/(a+1))
  
  output=root+kinematics+galaxy_ref
  result = FILE_TEST(output+'_kinematics.txt')
  ;move kinematics results to a separate file in case you accidentally overwrite them
  if result eq 1 then begin
    spawn,'mv '+output+'_kinematics.txt '+output+'_old_kinematics.txt'
  endif
  
  
  ;read in first spectrum to get velocisty scale
  flux=binned_spec[*,0]
  
  
  ; Read the list of filenames from the Single Stellar Population library
  stellar_lib = file_search(stellib_dir+'*.fits',COUNT=nfiles_stellib)
  FWHM_tem = 2.51 ; MILES spectra have a resolution FWHM of 2.51A.
  
  
  ; Only use the wavelength range in common between galaxy and stellar library.
  ; Check user defined wavelength range, and if larger than that of stellar
  ; templates, resample for kinematics measurements
  if nfiles_stellib ne 0 then begin
    nfiles=nfiles_stellib
    stellib_header=headfits(stellar_lib[0])
    stellib_wave1=sxpar(stellib_header,'CRVAL1')
    stellib_step=sxpar(stellib_header,'CDELT1')
    stellib_length=sxpar(stellib_header,'NAXIS1')
    stellib_wave2=stellib_wave1+stellib_step*(stellib_length-1)
  endif else begin
    stellar_lib = file_search(stellib_dir+'*V',COUNT=nfiles_stellib_txt)
    nfiles=nfiles_stellib_txt
    readcol,stellar_lib[0],stellar_wave,stellar_flux,format='F,F'
    stellib_wave1=stellar_wave[0]
    stellib_step=stellar_wave[1]-stellar_wave[0]
    stellib_length=n_elements(stellar_wave)
    stellib_wave2=stellar_wave[stellib_length-1]
  endelse
  
  ;;;define sample region over which to measure the kinematics
  ;;;this range covers the optical absorption features
  sample = where(wavelength gt alog10(4830) and wavelength lt alog10(6800))
  
  galaxy = flux[sample]/median(flux[sample])  ; normalize spectrum to avoid numerical issues
  NAN=where(Finite(galaxy) EQ 0)
  galaxy[NAN]=0     ;catch NAN values
  wave = 10^wavelength[sample]
  noise = galaxy*0 + 1           ; Assume constant noise per pixel here
  
  
  ; Convert velocity step below to km/s
  velScale = alog(wave[1]/wave[0])*c
  FWHM_gal = setup.FWHM_gal
  
  
  ; Extract the wavelength range and logarithmically rebin one spectrum
  ; to the same velocity scale of the SAURON galaxy spectrum, to determine
  ; the size needed for the array which will contain the template spectra.
  ;
  
  if nfiles_stellib ne 0 then begin
    fits_read, stellar_lib[0], ssp, h2
    lamRange2 = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
  endif else begin
    readcol,stellar_lib[0],stellar_wave,ssp,format='F,F'
    lamRange2 = [stellar_wave[0],stellar_wave[stellib_length-1]]
  endelse
  log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
  stars_templates = dblarr(n_elements(sspNew),nfiles)
  
  
  ; Convolve the whole stellar library of spectral templates
  ; with the quadratic difference between the galaxy and the
  ; stellar instrumental resolution. Logarithmically rebin
  ; and store each template as a column in the array TEMPLATES.
  ;
  ; Quadratic sigma difference in pixels stellar_lib --> SDSS
  ; The formula below is rigorously valid if the shapes of the
  ; instrumental spectral profiles are well approximated by Gaussians.
  ;
  
  if abs(FWHM_tem) gt abs(FWHM_gal) then stop
  FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
  sigma = FWHM_dif/2.355/stellib_step;sxpar(h2,'CDELT1') ; Sigma difference in pixels
  
  
  ; IMPORTANT: To avoid spurious velocity offsets of the templates, the
  ; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below.
  ;
  lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
  for j=0,nfiles-1 do begin
    if nfiles_stellib ne 0 then fits_read, stellar_lib[j], ssp $
    else readcol,stellar_lib[j],ssp,format='X,F'
    
    
    ssp = convol(ssp,lsf) ; Degrade template to SDSS resolution
    log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
    stars_templates[*,j] = sspNew/median(sspNew) ; nrmalizes templates
    
  endfor
  
  
  
  ;run ppxf automatically on each binned spectrum
  for run=0,n_bins-1,1 do begin
    print,'*** now measuring the kinematics of binned spectrum number '+string(run,format='(i4.4)')+' out of '+string(n_bins-1,format='(i4.4)')
    delvarx,flux,galaxy,NAN,mult_factor,dv,goodpixels,start  ;clean up memory usage
    flux=binned_spec[*,run]
    
    galaxy = flux[sample]/median(flux[sample])  ; normalize spectrum to avoid numerical issues
    NAN=where(Finite(galaxy) EQ 0)
    galaxy[NAN]=0     ;catch NAN values
    wave = 10^wavelength[sample]
    noise = galaxy*0 + 1           ; Assume constant noise per pixel here
    mult_factor=median(flux[sample])
    
    ; Convert velocity step below to km/s
    velScale = alog(wave[1]/wave[0])*c
    FWHM_gal = setup.FWHM_gal
    

    
    
    
    ; The galaxy and the template spectra do not have the same starting wavelength.
    ; For this reason an extra velocity shift DV has to be applied to the template
    ; to fit the galaxy spectrum. We remove this artificial shift by using the
    ; keyword VSYST in the call to PPXF below, so that all velocities are
    ; measured with respect to DV. This assume the redshift is negligible.
    ; In the case of a high-redshift galaxy one should de-redshift its
    ; wavelength to the rest frame before using the line below (see above).
    ;
    
    
    dv = (logLam2[0] - alog(wave[0]))*c ; km/s
    
    ;identify regions to include in the fit by masking out 
    ;emission features and sky lines. Goodpixels1 and 2 will 
    ;contain those pixels to be included
    goodPixels1 = ppxf_determine_goodPixels(alog(wave),lamRange2,Velocity)  ;remove emission features from the fit
    goodPixels2=where(wave lt 5550 or wave gt 5600 )  ;remove the nasty sky line at ~5575Angstroms
    match, goodPixels1, goodPixels2, suba, subb    ;combine the two sets of good pixels
    goodpixels=goodPixels1[suba]
    
    ; Here the actual fit starts. The best fit is plotted on the screen.
    start = [Velocity, 2*velScale] ; (km/s), starting guess for [V,sigma]
    
    
    
    ;;;;Original verison of code, no gas measurements
    
    ppxf, stars_templates, galaxy, noise, velScale, start, sol,$
      GOODPIXELS=goodPixels, MOMENTS=2, DEGREE=4, $
      VSYST=dv, ERROR=error, BIAS=Bias, BESTFIT=bestfit
      
    output=root+kinematics+galaxy_ref
    close,50
    openw,50,output+'_kinematics.txt',/APPEND
    
    
    if run eq 0 then printf,50, '#', 'bin','V', 'sigma', 'h3', 'h4', 'h5', 'h6', FORMAT='(8A10)'
    printf,50, run,sol[0:5,0], FORMAT='(i5,f10.1,4f10.3,A10)'
    close,50
    
    set_plot,'ps'
    !p.multi=0
    device,file=output+'_kinematics_'+string(run,format='(i4.4)')+'.eps',/color,xoffset=0,yoffset=0,xsize=18,ysize=13
    mn = min(bestfit[goodPixels], MAX=mx)
    resid = mn + galaxy - bestfit
    
    cgplot, wave, galaxy, color='black', XTITLE='Observed Wavelength A', $
      YTITLE='Relative Flux', /XSTYLE, /YSTYLE, YRANGE=[-1, 1.2*max(galaxy)];, YRANGE=[-0.1, 2]
    cgoplot, wave, bestfit, color='orange'
    cgOPlot, wave[goodPixels], galaxy[goodPixels] - bestfit[goodPixels], PSYM=4,symsize=0.8, COLOR='limegreen'
    
    device,/close
    
    ;print out best fit spectrum if user wants to remove the emission features
    if run eq 0 then spawn,'mkdir '+root+kinematics+'best_fit_spectra/'
    mkhdr, h_spec, bestfit
    sxaddpar,h_spec,'CRVAL1',wavelength[0]
    sxaddpar,h_spec,'CDELT1',wavelength[1]-wavelength[0]
    sxaddpar,h_spec,'CD1_1',wavelength[1]-wavelength[0]
    fits_write,root+kinematics+'best_fit_spectra/stellar_bin'+string(run,format='(i4.4)')+'.fits',bestfit,h_spec
    
    vc=sol[0]/c
    print, 'Best-fitting redshift z:', sqrt((1+vc)/(1-vc))-1
    
    
  endfor
  
endif







;========================================================
; *** Use the kinematics measurements to plot the kinematics maps

if setup.plot_kinematics eq 'y' then begin
  print,'############################'
  print,'# Plotting the kinematics #
  print,'############################'
  
  
  ;read in kinematics
  readcol,root+kinematics+galaxy_ref+'_kinematics.txt',format='D,D,D,D,D,X,X',bin_n,vel_in,sigma_in,h3_in,h4_in,comment='#',/SILENT
  ;read in positions of each bin on the ccd image
  readcol,root+galaxy_ref+'_bin_centers.txt',format='D,D,D',bin_n_in,xbar,ybar,/silent
  
  ;read in pixel scale in x and y directions from file header
  fits_read,root+file+'.fits',input_IFU,header_IFU
  x_scale=abs(sxpar(header_IFU,'CD1_1')*3600)      ;*3600 to convert to arcsec
  y_scale=abs(sxpar(header_IFU,'CD2_2')*3600)
  
  
  ;calculate distance of each bin from the centre of the galaxy, taking
  ;into account the inclination
  n_bins=max(bin_n_in)+1
  radius=dblarr(n_bins)
  
  measure_circular_radius,bin_n_in,xbar*x_scale,ybar*x_scale,0,0,PA, radius
  
  
  ;Make a new array for the kinematics, and rearrange in order of min to max velocity.
  ; Then can remove outliers due to bad fits or foreground stars.
  central_bin=where(abs(radius) eq min(abs(radius)))
  vel_centre=vel_in[central_bin]
  vel_centre= vel_centre[0]       ;defines vel_centre as a float, and not an array
  
  kinematics_temp1=dblarr(8,n_bins)
  kinematics_temp1[0,*]=bin_n
  kinematics_temp1[1,*]=vel_in
  kinematics_temp1[2,*]=sigma_in
  kinematics_temp1[3,*]=h3_in
  kinematics_temp1[4,*]=h4_in
  kinematics_temp1[5,*]=radius
  
  sort_by=1   ;use 1 to sort by line-of-sight velocity
  kinematics_temp2=colsort(kinematics_temp1,sort_by)
  
  kinematics_temp3=dblarr(8,n_bins)
  new_central_bin=where(kinematics_temp2[1,*] eq vel_centre)   ;central bin after sorting
  
  ;identify any outlying points, and remove from fit
  ;outliers include those values where sigma=sigma_in or sigma_max (=1000)
  badbins=dblarr(100)
  badbins[*]=-1
  m=0
  l=0
  
  for n=0,new_central_bin[0],1 do begin
    if n gt new_central_bin[0]-10 then offset=-10 else offset=10
    if (kinematics_temp2[1,n+offset]-kinematics_temp2[1,n]) lt 200 then begin
      kinematics_temp3[*,m]=kinematics_temp2[*,n]
      m+=1
    endif
    if (kinematics_temp2[1,n+offset]-kinematics_temp2[1,n]) ge 200 then begin
      badbins[l]=kinematics_temp2[0,n]
      l+=1
    endif else if kinematics_temp2[2,n] eq 1000 AND badbins[l-1] ne kinematics_temp2[0,n] then begin
      badbins[l]=kinematics_temp2[0,n]
      l+=1
    endif else if kinematics_temp2[2,n] eq 80.0 AND badbins[l-1] ne kinematics_temp2[0,n] then begin
      badbins[l]=kinematics_temp2[0,n]
      l+=1
    endif
  endfor
  
  for n2=new_central_bin[0]+1,n_bins-1,1 do begin
    if (kinematics_temp2[1,n2]-kinematics_temp2[1,n2-10]) lt 200 then begin
      kinematics_temp3[*,m]=kinematics_temp2[*,n2]
      m+=1
    endif
    if (kinematics_temp2[1,n2]-kinematics_temp2[1,n2-10]) ge 200 then begin
      badbins[l]=kinematics_temp2[0,n2]
      l+=1
    endif else if kinematics_temp2[2,n2] eq 1000 AND badbins[l-1] ne kinematics_temp2[0,n2] then begin
      badbins[l]=kinematics_temp2[0,n2]
      l+=1
    endif else if kinematics_temp2[2,n2] eq 80.0 AND badbins[l-1] ne kinematics_temp2[0,n2] then begin
      badbins[l]=kinematics_temp2[0,n2]
      l+=1
    endif
    
  endfor
  non_zero=where(abs(kinematics_temp3[1,*]) gt 0,count)
  kinematics_sorted=kinematics_temp3[*,0:count-1]
  
  max_vel=max(kinematics_sorted[1,*])
  min_vel=min(kinematics_sorted[1,*])
  max_sigma=max(kinematics_sorted[2,*])+50
  
  
  if max_vel-vel_centre gt vel_centre-min_vel then y_range=max_vel-vel_centre+100
  if max_vel-vel_centre le vel_centre-min_vel then y_range=vel_centre-min_vel+100
  
  set_plot,'ps'
  !p.multi=0
  loadct,0
  device,file=root+kinematics+galaxy_ref+'_kinematics.eps',xoffset=0,yoffset=0,xsize=18,ysize=13
  !p.multi=[0,2,2]
  !p.symsize=0.7
  xmin=min((kinematics_sorted[5,*]))-5
  xmax=max((kinematics_sorted[5,*]))+5
  velmax=max((kinematics_sorted[1,*]))+50
  velmin=min((kinematics_sorted[1,*]))-50
  
  plot,kinematics_sorted[5,*],kinematics_sorted[1,*],psym=sym(1),yrange=[velmin,velmax],ytitle='V!ILOS!N (km/s)',$
    xtitle='distance (arcsec)',xrange=[xmin,xmax],xmargin=[10,3],ytickformat='(i5)',symsize=0.5,/ystyle,/xstyle
  ;  oplot,[-100,100],[velocity_NED,velocity_NED],linestyle=1
  xyouts,-24,(2*y_range)*0.9+vel_centre-y_range,galaxy_ref, charsize=0.7
  plot,kinematics_sorted[5,*],kinematics_sorted[2,*],psym=sym(1),yrange=[0,max_sigma],ytitle=greek('sigma')+' (km/s)',$
    xtitle='distance (arcsec)',xrange=[xmin,xmax],xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*(max_sigma),galaxy_ref, charsize=0.7
  
  plot,kinematics_sorted[5,*],kinematics_sorted[3,*],psym=sym(1),xrange=[xmin,xmax],yrange=[-0.3,0.3],ytitle='h!I3!N',$
    xtitle='distance (arcsec)',xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*0.6-0.3,galaxy_ref, charsize=0.7
  plot,kinematics_sorted[5,*],kinematics_sorted[4,*],psym=sym(1),xrange=[xmin,xmax],yrange=[-0.3,0.3],ytitle='h!I4!N',$
    xtitle='distance (arcsec)',xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*0.6-0.3,galaxy_ref, charsize=0.7
  
  device,/close
  
  
  ; *** Plot kinematics maps
  ; need to first create new kinematics arrays for each pixel
  readcol,root+galaxy_ref+'_voronoi_2d_binning_output.txt',format='D,D,D',x_in,y_in,bin_new,comment='#',/SILENT
  vel_new=dblarr(n_elements(x_in))
  sigma_new=dblarr(n_elements(x_in))
  h3_new=dblarr(n_elements(x_in))
  h4_new=dblarr(n_elements(x_in))
  x_new=dblarr(n_elements(x_in))
  y_new=dblarr(n_elements(x_in))
  m=double(0.)
  
  ;make array for whole image, and assign velocity values to each pixel
  for n=0,n_elements(x_in)-1,1 do begin
    new_element=where(kinematics_sorted[0,*] eq bin_new[n],count)
    new_element=new_element[0]
    
    if new_element ge 0 then begin
      vel_new[m]=kinematics_sorted[1,new_element]-vel_centre
      sigma_new[m]=kinematics_sorted[2,new_element]
      h3_new[m]=kinematics_sorted[3,new_element]
      h4_new[m]=kinematics_sorted[4,new_element]
      x_new[m]=x_in[n]*x_scale
      y_new[m]=y_in[n]*y_scale
      m+=1
    endif
  endfor
  
  count=double(where(sigma_new ne 0 and sigma_new ne 1000,countx))
  
  set_plot,'ps'
  cgloadct,33
  device,file=root+kinematics+galaxy_ref+'_kinematics_maps.eps',/color,xoffset=0,yoffset=0,xsize=18,ysize=13
  !p.multi=[0,2,1]
  !X.MARGIN=[10,0]
  pixelSize=x_scale;1
  
  
  
  bin2d_display_pixels, x_new[count], y_new[count], vel_new[count], pixelSize
  cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(vel_new),maxrange=max(vel_new),title='V!ILOS!N (km/s)',position=[0.40,0.61,0.42,0.945],format='(i5)'
  bin2d_display_pixels, x_new[count], y_new[count], sigma_new[count], pixelSize
  cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(sigma_new[count]),maxrange=max(sigma_new[count]),title=greek('sigma')+' (km/s)',position=[0.93,0.61,0.95,0.945],format='(i5)'
  !p.multi=0
  
  device,/close
  

  
  
endif





end