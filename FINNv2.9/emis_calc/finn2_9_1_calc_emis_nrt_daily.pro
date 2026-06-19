;**************************************************************
; FINNv2.9.1 Emissions Calculation 
; Calculate emissions from pre-processed fire locations
; Reads fires produced from finn_py - pure python processor
; based on finn2_calc_emissions_v251.pro - uses v2.5 emissions factors
; called from daily processing cron job script
;**************************************************************

pro finn2_9_1_calc_emis_nrt_daily

  ; emissions date to process (YYYYJJJ), year of raster file
  args = COMMAND_LINE_ARGS()
  date_lab = args[0]
  year_rst = args[1]
            
  help,date_lab, year_rst

  today = bin_date(systime())
  todaystr = String(today[0:2],format='(i4,i2.2,i2.2)')  ;YYYYMMDD

  ; location of Fuel load, EFs, VOC speciation files
  path_inputs = '/glade/work/emmons/FINN_python/FINNv2.9nrt/emis_calc/finn_inputs'
  tag_fire = 'v2.9nrt'

  ; processed fires directory
  path_in = '/glade/derecho/scratch/emmons/finn2.9nrt_output'

  ; output emissions text files
  path_out = '/glade/derecho/scratch/emmons/finn2.9.1nrt_emis_txt'
  finnver = 'v2.9.1'

  ; processed fire files have dates labeled YYYYJJJ - year and julian day
  ; emissions date to process (YYYY-MM-DD)
  year = Fix(strmid(date_lab,0,4))
  doy = Fix(strmid(date_lab,4,3))
  help,year,doy
  jday = Julday(1,1,year)+doy-1
  caldat,jday,mm,dd,yy
  print,mm,dd,yy
  datestr = String(year,mm,dd,format='(i4,"-",i2.2,"-",i2.2)')

  ;parts = Strsplit(datestr,'-',/extract)
  ;year = Fix(parts[0])
  ;mm = Fix(parts[1])
  ;dd = Fix(parts[2])
  ;jday = Julday(mm,dd,year)
  ;doy = jday - Julday(1,1,year) + 1

  ;date_lab = String(year,doy,format='(i4,i3.3)')
  ;syrlc = String(year-1,format='(i4)')
  ;syrlc = args[1]
  year_emis = year

  ; date label for output files (YYYYMMDD)
  sdate_emis = String(year,mm,dd,format='(i4,i2.2,i2.2)')

  print,datestr,' ',date_lab,' ',sdate_emis

  simid = 'NRTmodvrs'

  file_in = path_in + '/out_'+tag_fire+'_'+date_lab+'_modlct_'+year_rst+'_modvcf_'+year_rst+'_regnum.csv'


  print,'--------- Starting processing of ',date_lab,' ',simid,' ---------'
  print,'Started on: ', SYSTIME()
  t0 = systime(1)               ;Procedure start time in seconds

  close, /all

  ; Open a log file
  logfile = path_out + '/logs/LOG_calcemis_FINN'+finnver+'_'+simid+'_'+sdate_emis+'_'+todaystr+'.txt'
  openw, ilun_log, logfile,/get_lun
  print, 'writing log file: ',logfile

  ; Log file for genveg assignment - a lot of output, comment out writing after checking
  log_genveg = path_out + '/logs/LOG_genveg_calcemis_FINN'+finnver+'_'+simid+'_'+sdate_emis+'_'+todaystr+'.txt'
  print,'log of genveg assignment: ',log_genveg
  openw, ilun_gv, log_genveg, /get_lun
  printf,ilun_gv, 'i, jday, fireid, polyid, lat, lon, lct_in, lct, tree_in, tree, herb_in, herb, bare_in, bare, flct, genveg'

  ;**************
  ; Standard input files (Fuel loads, EFs)
  file_fuelloads = path_inputs + '/Fuel_LOADS_NEW_022019.csv' 
  file_usfuel = path_inputs + '/LCTFuelLoad_fuel4_revisit20190521.csv'
  file_efs = path_inputs + '/EFs_byGenVeg_c20210601.csv'
  file_VOCsplit_M = path_inputs + '/NMOCfrac_byGenVeg_MOZ.csv'
  file_VOCsplit_S = path_inputs + '/NMOCfrac_byGenVeg_SAPRC.csv'
  file_VOCsplit_G = path_inputs + '/NMOCfrac_byGenVeg_GEOSCHEM.csv'

  ;  READ IN FUEL LOADING FILE: fuel loads for 5 veg types, for 13 regions
  ;  ALL FUEL INPUTS ARE IN g/m2 [-1 for missing values]
  fuel=read_csv(file_fuelloads, header=fueltypes)
  print,file_fuelloads,' Contains: ',fueltypes
  printf, ilun_log,file_fuelloads,' Contains: ',fueltypes
  ;print,'expected: region#, trop.for., temp.for., bor.for., woodySav, grassSav'
  if ((fueltypes[0] ne 'GlobalRegion') or (fueltypes[5] ne 'SavannaGrasslands')) then stop,'wrong format'
  ;   Set up fuel arrays
  regfuel = fuel.field1                 ;region index
  tffuel = fuel.field2  ;tropical forest fuels
  tefuel = fuel.field3  ;temperate forest fuels
  bffuel = fuel.field4  ;boreal forest fuels
  wsfuel = fuel.field5  ;woody savanna fuels
  grfuel = fuel.field6  ;grassland and savanna fuels

  ; READ in a secondary fuel loading file for use in US ONLY 
  ;   has fuel loads for tree and herb for each LCT 
  LCTfuel=read_csv(file_usfuel, header=usfueltypes)
  print,file_usfuel,' Contains: ',usfueltypes
  printf,ilun_log,file_usfuel,' Contains: ',usfueltypes
  ;print,'expected: Code,TREE,HERB'
  if (usfueltypes[2] ne 'HERB') then stop,'wrong format'
  lctfuelid = lctfuel.field1
  lcttree = lctfuel.field2
  lctherb = lctfuel.field3

  ; READ IN EMISSION FACTOR FILE [g compound emitted per kg dry biomass burned]
  ;   EFs for arbitrary number of species for GenVeg 1-6,9
  print,'Reading ',file_efs
  printf,ilun_log,'Reading ',file_efs
  ntype = 7  ;emission factors for genveg=1-6,9 
  openr, ilun,file_efs, /get_lun
  sdum = ' '
  readf,ilun,sdum  ;header/title 
  ;print,sdum
  readf,ilun,sdum  ;column names
  colnames_ef = strsplit(sdum,',',/extract)
  printf,ilun_log,file_efs,' contains: ',colnames_ef
  ;print,' expected: GenVegType,GenVegDescript, {species ...}'
  ncols = n_elements(colnames_ef)
  ef_genveg = intarr(ntype) ;genveg index for emission factors
  ;assumes first 2 cols are genveg, GenVegDesc, remaining columns are species EFs
  nspec = ncols - 2
  ef_species = colnames_ef[2:ncols-1]
  ;assumes 3rd row has molecular weights
  readf,ilun,sdum  ;MWs
  parts = strsplit(sdum,',',/extract)
  mws = Float(parts[2:ncols-1])
  ;read emission factors for all species
  emisfac = fltarr(ntype,nspec)
  for itype = 0,ntype-1 do begin
     readf,ilun,sdum
     cols = Strsplit(sdum,',',/extract)
     ef_genveg[itype] = fix(cols[0])
     emisfac[itype,*] = Float(cols[2:ncols-1])  ;emission factors in g/kg
  endfor
  free_lun,ilun

  print,'Have EFs for: ', ef_species
  print,' for GenVeg Types: ',ef_genveg
  for ispec = 0,nspec-1 do $
     printf,ilun_log,format='(a12,i4,7(f12.3))',ef_species[ispec], mws[ispec], emisfac[*,ispec]

  ;***************************************************************************************
  ; READ FIRE AND LAND COVER INPUT FILE (CREATED WITH PREPROCESSOR)
  ; **************************************************************************************
  ;  determine genveg
  ;  determine fuel loads and biomass burned
  ;  calculate emissions  
  ;--------------
  ; Read first line of fire file
  ; Set up arrays for saving emissions and other info
  ; -------------
  nfires = file_lines(file_in)-1L
  print,'# fires: ', nfires
  openr,ilun_in,file_in, /get_lun
  sdum=' '
  readf,ilun_in,sdum
  colnames_fires = strsplit(sdum,',',/extract)
  print,file_in,' contains: ',colnames_fires
  printf,ilun_log, file_in,' contains: ',colnames_fires

  ; v2.5.1 file:
  ;polyid,fireid,cen_lon,cen_lat,acq_date_utc,area_sqkm,alg_agg,v_lct,f_lct,v_tree,v_herb,v_bare,v_regnum,v_frp
  ; v2.9.1 file:
  ;polyid,fireid,cen_lon,cen_lat,acq_date_utc,area_sqkm,alg_agg,v_lct,f_lct,v_tree,v_herb,v_bare,v_regnum,v_frp
  ;0,0,52.94576,51.611920000000005,2026-05-30,1.0037730806999206,2,12,1.0,1.1666666666666667,55.5,43.333333333333336,7,33.48

  ;----------------
  ;  Read each line of input fire file 
  ;  Determine vegetation type, area, biomass burned
  ;  Save to arrays only valid fire points (correct date, valid vegetation, etc)
  ;----------------
  em_jday = intarr(nfires)
  em_date = lonarr(nfires)
  em_polyid = lonarr(nfires)
  em_fireid = lonarr(nfires)
  em_lat = fltarr(nfires)
  em_lon = fltarr(nfires)
  em_area = fltarr(nfires)
  em_bmass = fltarr(nfires)
  em_genveg = intarr(nfires)
  em_frp = fltarr(nfires)
  
  igood = 0L
  iskip_yr = 0L
  iskip_reg = 0L

  for i=0L,nfires-1 do begin
   readf,ilun_in,sdum
   parts = Strsplit(sdum,',',/extract,/preserve_null)
   if (n_elements(parts) lt 13) then begin
      printf,ilun_log,'input line wrong size: ',n_elements(parts),' ',sdum
      goto,skipfire
   endif

   polyid = long(parts[0])
   fireid = long(parts[1])
   lon = Float(parts[2])
   lat = Float(parts[3])
   area = Float(parts[5])
   lct = Fix(parts[7])
   flct = Float(parts[8])
   tree = Float(parts[9])
   herb = Float(parts[10])
   bare = Float(parts[11])
   globreg = Fix(parts[12])
   frp = Float(parts[13])

   dateparts = Strsplit(parts[4],'-',/extract)
   yy = Fix(dateparts[0])
   mm = Fix(dateparts[1])
   dd = Fix(dateparts[2])
   if (yy ne year_emis) then begin
     ;printf,ilun_log,'wrong year ',yy
     iskip_yr = iskip_yr+1
     goto,skipfire
   endif
   date = Long(Strjoin(dateparts))
   jday = Julday(mm,dd,yy) - Julday(1,1,yy) + 1.
   
   if (Size(globreg,/type) ne 2) then print,parts[12]
   if ((globreg lt 1) or (globreg gt 12)) then begin
      printf,ilun_log, format='("Fire ",i0," removed. global region:",i0," lon, lat: ",2f6.1)',i,globreg, lon,lat
      iskip_reg = iskip_reg+1
      goto, skipfire
   endif

   ; remove values of -9999 from VCF inputs
   if (tree lt 0.) then tree = 0.
   if (herb lt 0.) then herb = 0.
   if (bare lt 0.) then bare = 0.

   ; Calculate the total cover from the VCF product (CHECK TO MAKE SURE PERCENTAGES ADD TO 100%)
   totcov = tree+herb+bare

    ; Remove fires with no LCT assignment or in water bodies or snow/ice assigned by LCT
    ; LCT:
    ; 12, 14: cropland
    ; 13: urban
    ; 15: permanent snow or ice
    ; 16: barren
    ; 17: water
    ; 255: unclassified
    if ((lct ge 17) or (lct le 0) or (lct eq 15)) then begin 
      printf,ilun_log,format='("Fire ",i0," removed: lct = ",i0)',i,lct
      goto, skipfire
    endif
    if ((totcov ge 240.) or (totcov lt 1.)) then begin
       printf,ilun_log,format='("Fire ",i0," removed. totcov=",i0)',i,totcov
       goto, skipfire
    endif

    lct_in = lct
    tree_in = tree
    herb_in = herb
    bare_in = bare

    ; Scale VCF product to sum to 100. 
    if (totcov gt 101.) or (totcov lt 99.) then begin
      totcov_orig = totcov
      tree = tree_in*100./totcov
      herb = herb_in*100./totcov
      bare = bare_in*100./totcov
      totcov = bare + herb + tree
      printf,ilun_log,format='("Fire ",i0," had totcov adjusted from: ",i0," to: ",i0)',i,totcov_orig,totcov
    endif
    
    ; Fires with 100% bare cover reassign cover values based on LCT assignment
    if (bare ge 99.9) then begin
      if (lct le 5) then begin    ; Assign forest to the pixel
        tree = 60.
        herb = 40.
        bare = 0.
      endif
      if (lct ge 6 and lct le 8) or (lct eq 11) or (lct eq 14) then begin    ; Assign woody savanna to the pixel
        tree = 50.
        herb = 50.
        bare = 0.
      endif
      if (lct eq 9) or (lct eq 10) or (lct eq 12) or (lct eq 13) or (lct eq 16) then begin  ; Assign grassland to the pixel
        tree = 20.
        herb = 80.
        bare = 0.
     endif
     printf,ilun_log,format='("Fire ",i0," with LCT=",i0," had 100% bare adjusted to T/H/B: ",i0,"/",i0,"/",i0)',i,lct,tree,herb,bare
    endif

    ; ######################################################
    ; Assign Generic land cover to fire based on
    ;   global location and lct information
    ; ######################################################
    ;Generic land cover codes (genveg) are as follows:
    ;1 grassland
    ;2 shrub
    ;3 Tropical Forest
    ;4 Temperate Forest
    ;5 Boreal Forest
    ;6 Temperate Evergreen Forest
    ;7 Pasture
    ;8 Rice
    ;9 Crop (generic)
    ;10  Wheat
    ;11  Cotton
    ;12  Soy
    ;13  Corn
    ;14  Sorghum
    ;15  Sugar Cane
    genveg = -999

    case lct of
       1: begin   ; Evergreen Needleleaf Forest to Boreal or Temperate Evergreen
          if (lat gt 50.) then genveg = 5 else genveg = 6
       end 
       2: begin
           if (lat ge -23.5 and lat le 23.5) then begin
             genveg = 3 ; Tropical Forest
           endif else begin
             genveg = 4 ; Temperate Forest
           endelse
        end 
       3: begin  ; deciduous Needleleaf Forest to Boreal or Temperate forest
            if (lat gt 50.) then genveg = 5 else genveg = 4 
        end
       4: genveg = 4 ; Temperate Forest
       5: begin ; Mixed Forest, assign type by latitude
          if (lat gt 50.) then begin
             genveg = 5
          endif else if (lat ge -23.5 and lat le 23.5) then genveg = 3 else genveg = 4
        end 
       6: genveg = 2  ;Woody Savanna or Shrubs
       7: genveg = 2
       8: genveg = 2
       9: genveg = 1  ;Grasslands and Savanna
       10: genveg = 1
       11: genveg = 1
       12: genveg = 9  ;Croplands
       13: begin   ;Urban
           if (tree lt 40.) then begin
            genveg = 1        ; grasslands
            lct = 10       ; set to grassland
            endif else if (tree ge 40. and tree lt 60.) then begin
              genveg = 2  ; woody savannas
              lct = 8 ; set to woody savanna
            endif else if (tree ge 60.) then begin ; assign forest based on latitude
              if lat gt 50. then begin
                genveg = 5
                lct = 1  ; set to evergreen needleleaf forest
              endif else begin
                if lat ge -30. and lat le 30. then genveg = 3 else genveg = 4
                lct = 5 ; set to mixed forest
              endelse
            endif
         end
       14: genveg = 1
       16: genveg = 1
       else: genveg = -1
    endcase 

    if (genveg le 0) then begin
       printf,ilun_log,format='("Fire ",i0," does not have genveg set. lat,lon, LCT = ",i0)',i,lat,lon,lct
       goto,skipfire
    endif
    printf,ilun_gv, format='(i10,i5,2I10,2(f9.3),2i3,7(f6.1),i3)', $
        i, jday, fireid,polyid, lat, lon, lct_in, lct, tree_in, tree, herb_in, herb, bare_in, bare, flct, genveg 

    ; ####################################################
    ; Assign Fuel Loads based on Generic land cover
    ;   and global region location
    ;   units are in g dry mass/m2
    ; ####################################################
    ireg = globreg-1 

    ; Assign biomass density according to veg. type
    bmass1 = -1.
    case genveg of
       1: bmass1 = grfuel[ireg]
       2: bmass1 = wsfuel[ireg]
       3: bmass1 = tffuel[ireg]
       4: bmass1 = tefuel[ireg]
       5: bmass1 = bffuel[ireg]
       6: bmass1 = tefuel[ireg]
       9: bmass1 = 902.
    endcase
   ; Assign boreal forests in Southern Asia the biomass density of the temperate forest for the region
    if (genveg eq 5 and globreg eq 11) then bmass1 = tefuel[ireg]

    if (bmass1 lt 0.) then begin
      printf,ilun_log, 'BMASS1 < 0: Fire',i,' removed. genveg =',genveg,' globreg = ',globreg,' ireg = ',ireg
      goto, skipfire
    endif

    ; *****************************************************************************************
    ; Assign Burning Efficiencies based on Generic land cover (Hoezelmann et al. [2004] Table 5
    ; *****************************************************************************************
    ; ASSIGN CF VALUES (Combustion Factors)
    if (tree gt 60.) then begin      ;FOREST
      ; Values from Table 3 Ito and Penner [2004]
      CF1 = 0.30          ; Live Woody
      CF3 = 0.90          ; Leafy Biomass
      ;CF4 = 0.90          ; Herbaceous Biomass
      ;CF5 = 0.90          ; Litter Biomass
      ;CF6 = 0.30          ; Dead woody
    endif
    if (tree gt 40.) and (tree le 60.) then begin   ;WOODLAND
      CF3 = exp(-0.013*tree)       ; Apply to all herbaceous fuels
      CF1 = 0.30                   ; Apply to all coarse fuels in woodlands
      ; From Ito and Penner [2004]
    endif
    If (tree le 40.) then begin       ;GRASSLAND
      CF3 = 0.98     ;Range is between 0.44 and 0.98 - Assumed UPPER LIMIT!
    endif
    ; *******************************************************************************************
    ; Calculate the biomass burned of each classification (herbaceous, woody, and forest)
    ; These are in units of g dry matter/m2
    ; Bmass is the total burned biomass
    ; herbbm is the Herbaceous biomass burned
    ; coarsebm is the Woody biomass burned

    coarsebm = bmass1
    herbbm = grfuel[ireg]

    ; Determine if in North America and use updated fuel loading for North America (Global Region 1)
    if (globreg eq 1) then begin
      ; Assign coarse and herb biomass based on lct
      coarsebm = lcttree[lct]
      herbbm = lctherb[lct]
    endif

    ;  Grasslands
    if (tree le 40.) then begin
      Bmass = ((herb/100.)*herbbm*CF3)+((tree/100.)*herbbm*CF3)
      ; Assumed here that litter biomass = herbaceous biomass and that the percent tree
      ;   in a grassland cell contributes to fire fuels
      ; Assuming here that the duff and litter around trees burn
    endif
    ; Woodlands
    if (tree gt 40.) and (tree le 60.) then begin
      Bmass = ((herb/100.)*herbbm*CF3) + ((tree/100.)*(herbbm*CF3+coarsebm*CF1))
    endif
    ; Forests
    if (tree gt 60.) then begin
      Bmass = ((herb/100.)*herbbm*CF3) + ((tree/100.)*(herbbm*CF3+coarsebm*CF1))
    endif

    ; Convert units to be consistent; adjust area burned for vegetation and bare fraction
    bmass = bmass/1000. ; convert g-dm/m2 to kg-dm/m2
    areanow = area*flct*1.0e6    ; convert km2 to m2
    area_bare = areanow*(bare/100.)
    areanow = areanow - area_bare     ; remove bare area from being burned
    if (areanow lt 1.) then begin
       printf,ilun_log,'area = 0. area,flct,bare: ',areanow,flct,bare
       goto,skipfire
    endif

    em_jday[igood] = doy
    ;em_date[igood] = date
    em_polyid[igood] = polyid
    em_fireid[igood] = fireid
    em_lat[igood] = lat
    em_lon[igood] = lon
    em_area[igood] = areanow
    em_bmass[igood] = bmass
    em_genveg[igood] = genveg
    em_frp[igood] = frp

    igood = igood+1L

    skipfire:
endfor ; End loop over fires
free_lun,ilun_gv

print,'finished reading fire file'

ngood = igood
em_jday = em_jday[0:ngood-1]
em_date = em_date[0:ngood-1]
em_polyid = em_polyid[0:ngood-1]
em_fireid = em_fireid[0:ngood-1]
em_lat = em_lat[0:ngood-1]
em_lon = em_lon[0:ngood-1]
em_area = em_area[0:ngood-1]
em_bmass = em_bmass[0:ngood-1]
em_genveg = em_genveg[0:ngood-1]
em_frp = em_frp[0:ngood-1]

; Sort fires by day
;indsort = sort(em_jday)
;em_jday = em_jday[indsort]
;em_date = em_date[indsort]
;em_polyid = em_polyid[indsort]
;em_fireid = em_fireid[indsort]
;em_lat = em_lat[indsort]
;em_lon = em_lon[indsort]
;em_area = em_area[indsort]
;em_bmass = em_bmass[indsort]
;em_genveg = em_genveg[indsort]
;em_frp = em_frp[indsort]

;print,'finished sorting arrays of good points', ngood
print,' # fires with emissions: ',ngood
print,' % of total fires saved: ',float(ngood)/float(nfires)*100.

printf,ilun_log,'# fires skipped because wrong year: ',iskip_yr
printf,ilun_log,'# fires skipped because no region assigned: ',iskip_reg
printf,ilun_log,'# fires with emissions: ',ngood
printf,ilun_log,'% of total fires saved: ',float(ngood)/float(nfires)*100.
;
; Write base emissions with FRP
  outfile_txt = path_out + '/FINN'+finnver+'_'+simid + '_base_FRP_'+sdate_emis+'.txt'
  print,'Writing output to: ',outfile_txt
  openw,ilun_out,outfile_txt, /get_lun
  species_list = Strjoin(ef_species,',')

  printf,ilun_out, 'DAY,POLYID,FIREID,GENVEG,LATI,LONGI,AREA,BMASS,FRP,'+species_list
  format_out = '(I5,2(",",I10),",",I3,2(",",F9.3),2(",",E10.3),",",F9.3,'+string(nspec,format='(i0)')+'(",",E10.3))'

print,'Calculating base emissions... '
; CALCULATE EMISSIONS = area*BMASS*EF
; Units: EF[g-species/kg-dm]*[kg/g] * Area[m2] * Bmass[kg-dm/m2]
; Convert gas-phase species to [moles/fire/day] by scaling with molecular weight [kg/mole]
; Keep aerosols in [kg/fire/day] (MW=1 from EFs file)

emis_spec = fltarr(nspec)

for ifire = 0,ngood-1 do begin
    itype = where(ef_genveg eq em_genveg[ifire]) & itype = itype[0]
    if (itype lt 0) then begin
       printf,ilun_log,'no EF for this genveg: ',em_genveg[ifire]
       goto,skipfire1
    endif
    for ispec=0,nspec-1 do begin
       if (mws[ispec] ne 1.) then $
         emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire] /(mws[ispec]*1.e-3) $
       else emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire]
    endfor
    ; Write each fire to text file
    printf,ilun_out, format=format_out, em_jday[ifire],em_fireid[ifire],em_polyid[ifire], $
           em_genveg[ifire], em_lat[ifire],em_lon[ifire], em_area[ifire], em_bmass[ifire], $
           em_frp[ifire], emis_spec

    skipfire1:
endfor
free_lun,ilun_out


; Read factors to convert NMOC [kg/day] to each VOC [moles/day]
;  for MOZART, SAPRC and GEOS-Chem separately

print,'------ MOZART ------'
ind_NMOC = where(ef_species eq 'NMOC')
ind_NMOC = ind_NMOC[0]
have_file = File_test(file_VOCsplit_M)
if ((ind_NMOC  ge 0) and (have_file)) then begin
   print,'index of NMOC: ', ind_NMOC,' ',ef_species[ind_NMOC]
   print,'Reading ',file_VOCsplit_M
   openr, ilun_voc, file_VOCsplit_M, /get_lun
   sdum=''
   readf, ilun_voc, sdum  ;header line
   readf, ilun_voc, sdum  ;column labels
   ;GenVegIndex, EF_fuel_type, APIN,...
   colnames = strsplit(sdum,',', /extract)
   nvocs_split = n_elements(colnames)-2
   vocnames = colnames[2:nvocs_split+1]
   print,'MOZART VOC speciation for: '
   print,vocnames
   ntypes = 7
   voc_fraction = fltarr(nvocs_split,ntypes)
   genveg_voc = intarr(ntypes)
   for itype = 0,ntypes-1 do begin
    readf,ilun_voc, sdum
    cols = strsplit(sdum,',', /extract)
    genveg_voc[itype] = Fix(cols[0])
    for ivoc=0,nvocs_split-1 do voc_fraction[ivoc,itype] = float(cols[ivoc+2])
   endfor
   free_lun,ilun_voc
   ;for ivoc = 0,nvocs_split-1 do print,vocnames[ivoc],reform(voc_fraction[ivoc,*])
  endif else begin
     print,'Not calculating VOC speciation: NMOC not in Emission Factors file or cannot open ',file_VOCsplit_M
     nvocs_split = 0
  endelse

  ; SET UP OUTPUT TEXT FILE for base and MOZART species
  outfile_txt = path_out + '/FINN'+finnver+'_'+simid + '_MOZART_'+sdate_emis+'.txt'
  print,'Writing output to: ',outfile_txt
  openw,ilun_out,outfile_txt, /get_lun
  species_list = Strjoin(ef_species,',')
  if (nvocs_split gt 0) then vocs_list = ', '+Strjoin(vocnames,',') else vocs_list = ''

  printf,ilun_out, 'DAY,POLYID,FIREID,GENVEG,LATI,LONGI,AREA,BMASS,'+species_list+vocs_list
  format_out = '(I5,2(",",I10),",",I3,2(",",F9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3))'
  if (nvocs_split gt 0) then $
   format_out = '(I5,2(",",I10),",",I3,2(",",F9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3),'+string(nvocs_split,format='(i0)')+'(",",E10.3))'

print,'Calculating MOZART emissions... '

; CALCULATE EMISSIONS = area*BMASS*EF
; Units: EF[g-species/kg-dm]*[kg/g] * Area[m2] * Bmass[kg-dm/m2]
; Convert gas-phase species to [moles/fire/day] by scaling with molecular weight [kg/mole]
; Keep aerosols in [kg/fire/day] (MW=1 from EFs file)

emis_spec = fltarr(nspec)
if (nvocs_split gt 0) then emis_vocs = fltarr(nvocs_split)

for ifire = 0,ngood-1 do begin
    itype = where(ef_genveg eq em_genveg[ifire]) & itype = itype[0]
    if (itype lt 0) then begin
       printf,ilun_log,'no EF for this genveg: ',em_genveg[ifire]
       goto,skipfire2
    endif
    for ispec=0,nspec-1 do begin
       if (mws[ispec] ne 1.) then $
         emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire] /(mws[ispec]*1.e-3) $
       else emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire]
    endfor 
    ; Calculate emissions for VOCs (moles) as fraction of NMOC (kg-species)
    igen = where(genveg_voc eq em_genveg[ifire])
    if (igen[0] ge 0) then begin
     for ivoc = 0,nvocs_split-1 do begin
        emis_vocs[ivoc] = voc_fraction[ivoc,igen[0]] * emis_spec[ind_NMOC]
      endfor
    endif else emis_vocs[*] = 0.
    ; Write each fire to text file
    printf,ilun_out, format=format_out, em_jday[ifire],em_fireid[ifire],em_polyid[ifire], $
           em_genveg[ifire], em_lat[ifire],em_lon[ifire], em_area[ifire], em_bmass[ifire], $
           emis_spec, emis_vocs

    skipfire2:
endfor
free_lun,ilun_out

;---------------------------
; CALCULATE SAPRC speciation
; Read factors to convert NMOC [kg/day] to each VOC [moles/day]
;
; skip
;goto,skipsap
print,'------ SAPRC ------'
ind_NMOC = where(ef_species eq 'NMOC')
ind_NMOC = ind_NMOC[0]
have_file = File_test(file_VOCsplit_S)
if ((ind_NMOC  ge 0) and (have_file)) then begin
   print,'index of NMOC: ', ind_NMOC,' ',ef_species[ind_NMOC]
   print,'Reading ',file_VOCsplit_S
   openr, ilun_voc, file_VOCsplit_S, /get_lun
   sdum=''
   readf, ilun_voc, sdum  ;header line
   readf, ilun_voc, sdum  ;column labels
   ;GenVegIndex, EF_fuel_type, species...
   colnames = strsplit(sdum,',', /extract)
   nvocs_split = n_elements(colnames)-2
   vocnames = colnames[2:nvocs_split+1]
   print,'SAPRC VOC speciation for: '
   print,vocnames
   ntypes = 7
   voc_fraction = fltarr(nvocs_split,ntypes)
   genveg_voc = intarr(ntypes)
   for itype = 0,ntypes-1 do begin
    readf,ilun_voc, sdum
    cols = strsplit(sdum,',', /extract)
    genveg_voc[itype] = Fix(cols[0])
    for ivoc=0,nvocs_split-1 do voc_fraction[ivoc,itype] = float(cols[ivoc+2])
   endfor
   free_lun,ilun_voc
   ;for ivoc = 0,nvocs_split-1 do print,vocnames[ivoc],reform(voc_fraction[ivoc,*])
  endif else begin
     print,'Not calculating VOC speciation: NMOC not in Emission Factors file or cannot open ',file_VOCsplit_S
     nvocs_split = 0
  endelse
  print,'SAPRC #VOCs: ',nvocs_split
 
  ; SET UP OUTPUT TEXT FILE for base and SAPRC species
  outfile_txt = path_out + '/FINN'+finnver+'_'+simid + '_SAPRC_'+sdate_emis+'.txt'
  print,'Writing output to: ',outfile_txt
  openw,ilun_out,outfile_txt, /get_lun
  species_list = Strjoin(ef_species,',')
  if (nvocs_split gt 0) then vocs_list = ', '+Strjoin(vocnames,',') else vocs_list = ''

  printf,ilun_out, 'DAY,POLYID,FIREID,GENVEG,LATI,LONGI,AREA,BMASS,'+species_list+vocs_list
  format_out = '(I5,2(",",I10),",",I3,2(",",f9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3))'
  if (nvocs_split gt 0) then $
   format_out = '(I5,2(",",I10),",",I3,2(",",f9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3),'+string(nvocs_split,format='(i0)')+'(",",E10.3))'

print,'Calculating SAPRC emissions... '
emis_spec = fltarr(nspec)
if (nvocs_split gt 0) then emis_vocs = fltarr(nvocs_split)
for ifire = 0,ngood-1 do begin
    itype = where(ef_genveg eq em_genveg[ifire]) & itype = itype[0]
    if (itype lt 0) then begin
       printf,ilun_log,'no EF for this genveg: ',em_genveg[ifire]
       goto,skipfire3
    endif
    for ispec=0,nspec-1 do begin
       if (mws[ispec] ne 1.) then $
         emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire] /(mws[ispec]*1.e-3) $
       else emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire]
    endfor
    ; Calculate emissions for VOCs (moles) as fraction of NMOC (kg-species)
    igen = where(genveg_voc eq em_genveg[ifire])
    if (igen[0] ge 0) then begin
     for ivoc = 0,nvocs_split-1 do begin
        emis_vocs[ivoc] = voc_fraction[ivoc,igen[0]] * emis_spec[ind_NMOC]
      endfor
    endif else emis_vocs[*] = 0.
    ; Write each fire to text file
    printf,ilun_out, format=format_out, em_jday[ifire],em_fireid[ifire],em_polyid[ifire], $
           em_genveg[ifire], em_lat[ifire],em_lon[ifire], em_area[ifire], em_bmass[ifire], $
           emis_spec, emis_vocs
    skipfire3:
endfor
free_lun,ilun_out

skipsap:

;---------------------------
; CALCULATE GEOSCHEM speciation
; Read factors to convert NMOC [kg/day] to each VOC [moles/day]
;
; skip because not sure anyone is using these
goto, skipgc 

print,'------ GEOSCHEM ------'
ind_NMOC = where(ef_species eq 'NMOC')
ind_NMOC = ind_NMOC[0]
have_file = File_test(file_VOCsplit_G)
if ((ind_NMOC  ge 0) and (have_file)) then begin
   print,'index of NMOC: ', ind_NMOC,' ',ef_species[ind_NMOC]
   print,'Reading ',file_VOCsplit_G
   openr, ilun_voc, file_VOCsplit_G, /get_lun
   sdum=''
   readf, ilun_voc, sdum  ;header line
   readf, ilun_voc, sdum  ;column labels
   ;GenVegIndex, EF_fuel_type, species...
   colnames = strsplit(sdum,',', /extract)
   nvocs_split = n_elements(colnames)-2
   vocnames = colnames[2:nvocs_split+1]
   print,'GEOS-Chem VOC speciation for: '
   print,vocnames
   ntypes = 7
   voc_fraction = fltarr(nvocs_split,ntypes)
   genveg_voc = intarr(ntypes)
   for itype = 0,ntypes-1 do begin
    readf,ilun_voc, sdum
    cols = strsplit(sdum,',', /extract)
    genveg_voc[itype] = Fix(cols[0])
    for ivoc=0,nvocs_split-1 do voc_fraction[ivoc,itype] = float(cols[ivoc+2])
   endfor
   free_lun,ilun_voc
   ;for ivoc = 0,nvocs_split-1 do print,vocnames[ivoc],reform(voc_fraction[ivoc,*])
  endif else begin
     print,'Not calculating VOC speciation: NMOC not in Emission Factors file or cannot open ',file_VOCsplit_G
     nvocs_split = 0
  endelse
  print,'GC #VOCS:',nvocs_split

  ; SET UP OUTPUT TEXT FILE for base and GC species
  outfile_txt = path_out + '/FINN'+finnver+'_'+simid + '_GEOSCHEM_'+sdate_emis+'.txt'
  print,'Writing output to: ',outfile_txt
  openw,ilun_out,outfile_txt, /get_lun
  species_list = Strjoin(ef_species,',')
  if (nvocs_split gt 0) then vocs_list = ', '+Strjoin(vocnames,',') else vocs_list = ''

  printf,ilun_out, 'DAY,POLYID,FIREID,GENVEG,LATI,LONGI,AREA,BMASS,'+species_list+vocs_list
  format_out = '(I5,2(",",I10),",",I3,2(",",f9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3))'
  if (nvocs_split gt 0) then $
   format_out = '(I5,2(",",I10),",",I3,2(",",f9.3),2(",",E10.3),'+string(nspec,format='(i0)')+'(",",E10.3),'+string(nvocs_split,format='(i0)')+'(",",E10.3))'

print,'Calculating GEOSCHEM emissions... '
emis_spec = fltarr(nspec)
if (nvocs_split gt 0) then emis_vocs = fltarr(nvocs_split)
for ifire = 0,ngood-1 do begin
    itype = where(ef_genveg eq em_genveg[ifire]) & itype = itype[0]
    if (itype lt 0) then begin
       printf,ilun_log,'no EF for this genveg: ',em_genveg[ifire]
       goto,skipfire4
    endif
    for ispec=0,nspec-1 do begin
       if (mws[ispec] ne 1.) then $
         emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire] /(mws[ispec]*1.e-3) $
       else emis_spec[ispec] = emisfac[itype,ispec]*1.e-3 * em_area[ifire] * em_bmass[ifire]
    endfor
    ; Calculate emissions for VOCs (moles) as fraction of NMOC (kg-species)
    igen = where(genveg_voc eq em_genveg[ifire])
    if (igen[0] ge 0) then begin
     for ivoc = 0,nvocs_split-1 do begin
        emis_vocs[ivoc] = voc_fraction[ivoc,igen[0]] * emis_spec[ind_NMOC]
      endfor
    endif else emis_vocs[*] = 0.
    ; Write each fire to text file
    printf,ilun_out, format=format_out, em_jday[ifire],em_fireid[ifire],em_polyid[ifire], $
           em_genveg[ifire], em_lat[ifire],em_lon[ifire], em_area[ifire], em_bmass[ifire], $
           emis_spec, emis_vocs
    skipfire4:
endfor
free_lun,ilun_out
skipgc: 

t1 = systime(1)-t0
print,'Running time: '+ strtrim(string(fix(t1)/60,t1 mod 60,format='(i3,1h:,i2.2)'),2)
print,'Completed at: ', SYSTIME()
printf,ilun_log,'Running time: '+ strtrim(string(fix(t1)/60,t1 mod 60,format='(i3,1h:,i2.2)'),2)
printf,ilun_log,'Completed at: ', SYSTIME()

free_lun,ilun_log
close,/all                      ;close all text files
 
end

