PRO read_input, setup_file, setup
;read in the input file for the IFU decomposition.
;
;Note- this code is based on the read_setup code distributed as part
;of the galapagos distribution (https://github.com/MegaMorph/galapagos, 
;Boris Haeussler, BorisHaeussler.astro@gmail.com)
;
;example usage
;   > read_input, setup_file, setup
;
;example:
  if file_test(setup_file) eq 0 then begin
     print, 'input file does not exist'
     stop
  ENDIF
  
  ON_IOERROR, bad_input
  
  len_num = 4                   ;length of the numbering scheme, e.g. 4 for 'A00)'
  
  
  ; 
  ; *** create structure to contail all the input information***
  ; create_struct works in pair, where the first is the name of the 
  ; variable and is followed by its value. the structure is created 
  ; using create_struct here, and then the input file is read in, 
  ; and the values put into the correct places in the structure using 
  ; the CASE keyword.
  ; 
  setup = create_struct('root', '', $
                        'galaxy_ref', '', $
                        'file', '', $
                        'kinematics', '', $
                        'stellib_dir','',$
                        'x_centre', 0., $
                        'y_centre', 0., $
                        'cont_wavelength', 0., $
                        'cont_range', 0., $
                        'targetSN', 0., $
                        'limit', 0., $
                        'Redshift', 0., $
                        'PA', 0., $
                        'FWHM_gal',0.,$
                        'start_wavelength', 0., $
                        'end_wavelength', 0., $
                        'no_bins', 0., $
                        'no_slices', 0., $
                        'log_rebin_data', '', $
                        'bin_data', '', $
                        'measure_kinematics', '', $
                        'plot_kinematics', '')


  
; check format for backwards compatibility
  block_bd = 0
  line = ''
  openr, 1, setup_file
  WHILE NOT eof(1) DO BEGIN
     readf, 1, line       
;get rid of leading and trailing blanks
     line = strtrim(line, 2)
;comment or empty line encountered?
     IF strmid(line, 0, 1) EQ '#' OR strlen(line) EQ 0 THEN CONTINUE       
;comment at end of line?
     pos = strpos(line, '#')
     IF pos EQ -1 THEN pos = strlen(line)
     content = strtrim(strmid(line, len_num, pos-len_num), 2)
     IF strupcase(strmid(line, 0, len_num)) eq 'G00)' then block_bd = 1
  ENDWHILE
  close, 1
  
  line = ''
  openr, 1, setup_file
  WHILE NOT eof(1) DO BEGIN
     readf, 1, line
     
;get rid of leading and trailing blanks
     line = strtrim(line, 2)
;    print, linestop
     
     
;comment or empty line encountered?
     IF strmid(line, 0, 1) EQ '#' OR strlen(line) EQ 0 THEN CONTINUE
     
;comment at end of line?
     pos = strpos(line, '#')
     IF pos EQ -1 THEN pos = strlen(line)
     
     content = strtrim(strmid(line, len_num, pos-len_num), 2)
     
     CASE strupcase(strmid(line, 0, len_num)) OF
;CHANGE=====have to trigger bad input / proper filenames
        
        'A00)': setup.root = content
        'A01)': setup.galaxy_ref = content
        'A02)': setup.file = content
        'A03)': setup.kinematics = content
        'A04)': setup.stellib_dir = content
        
        'C00)': setup.x_centre = float(content)
        'C01)': setup.y_centre = float(content)
        'C02)': setup.cont_wavelength = float(content)
        'C03)': setup.cont_range = float(content)
        'C04)': setup.targetSN = float(content)
        'C05)': setup.limit = float(content)
        'C06)': setup.Redshift = float(content)
        'C07)': setup.PA = float(content)
        'C08)': setup.FWHM_gal = float(content)

        'D02)': setup.start_wavelength = float(content)
        'D03)': setup.end_wavelength = float(content)
        'D04)': setup.no_bins = float(content)
        'D05)': setup.no_slices = float(content)
        
        'E00)': setup.log_rebin_data = content
        'E01)': setup.bin_data = content
        'E02)': setup.measure_kinematics = content
        'E03)': setup.plot_kinematics = content
        
     ENDCASE
  ENDWHILE
  close, 1
  
  
  return
  
bad_input:
  message, 'Invalid Entry in '+setup_file
END