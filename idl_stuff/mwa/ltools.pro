; ltools.pro
; A small collection of tools useful for dealing with "L-file" data
; By CLMW 05-May-2010

;********************************************;

; reorder_coarse
; This function will reorder the coarse channels in a 768 fine channel
; dataset in order to correct for the Rx coarse PFB output order.
; It accepts an array of dimensionality [768,*,*] and returns an array
; of the same dimensionality

function reorder_coarse,data,center


  minchan=center-12 ; figure out the minimum channel
  nbank1=0
  nbank2=0
  correct_chan=intarr(24)

  ; Figure out the  order of the channels based on the 2 PFB banks
  for i=minchan,minchan+23 do if i le 128 then nbank1++ else nbank2++
  for i=0,nbank1-1 do correct_chan[i]=i
  for i=0,nbank2-1 do correct_chan[i+nbank1]=23-i

  outdata=data*0 ;make the output array
  
  for i=0,23 do begin
     outdata[32*correct_chan[i]:32*correct_chan[i]+31,*,*]=data[32*i:32*i+31,*,*]
  end

  return,outdata
end


;********************************************;

; read_file_list
; This is a dumb wrapper around read_ascii to make an array of
; strings based on the lines of an ascii input file

pro read_file_list,filename,list

myTemplate={version:1,datastart:0,delimiter:',',missingvalue:0,fieldcount:1,fieldtypes:7,fieldnames:"field1",fieldlocations:0,fieldgroups:0,commentsymbol:'#'}
text=read_ascii(filename,template=myTemplate)
list=text.field1

end

;********************************************;

; ccindex
; This function determines which cross-correlations in an L-file array
; include a certain input number (input).  "arr" is the array which is
; to be subscripted (it just uses it to get size information).  It
; will return an array of the array subscripts to use to get these
; baselines.  Optionally, "labels" will list which other inputs are
; included in each baseline.
; An example would be: index=ccindex(7,arr,labels=lab).  To access
; only the baselines with input 7 then you could say arr[*,index,*]


function ccindex,input,arr,labels=labels

  s=size(arr,/dim) 
  ncorrs=long(s[1]) ; get the second dimension of the array, which corresponds to the number of correlations
  ninps=1/2.0*(1.0+sqrt(1.0+double(ncorrs)*8.0))
  if ninps ne fix(ninps) then begin & print,"Bad array size" & return,-1 & end

  ninps=fix(ninps)
  ;print,ninps
  if input ge ninps then begin & print,"Input out of range" & return,-1 & end
   
  bls=intarr(ninps-1)
  labels=intarr(ninps-1)


  index=0L
  t1=0
  t2=0

  ; couldn't think of a better way easily than to just loop through...
  for i=0,ncorrs-1 do begin
     t2++
     if t2 ge ninps then begin
        t1++
        t2=t1+1
     end
     if t1 lt ninps and t2 lt ninps then begin
        if t1 eq input or t2 eq input then begin
           bls[index]=i
           if t1 eq input then labels[index]=t2
           if t2 eq input then labels[index]=t1
           index++
        end
     end

 
  end

  return,bls
end

;********************************************;

; waterfall
; This program will make waterfall plots of arrays.  It is
; meant to be used with arrays of L-files.  min or max will set the
; range of the plot.  If "cchan" is given (the center coarse channel
; of the data), then the channel axis be converted to frequency.

pro waterfall,arr,min=mina,max=maxa,title=tit,cchan=cchan


  p=reform(arr)
  s=size(p,/dim)
  if n_elements(s) ne 2 then begin
     print,"Not a 2D array"
     return
  end

  if keyword_set(cchan) and s[0] ne 768 then begin
     print,"Currently will only convert to frequency for a 768-channel array"
     return
  end



  if keyword_set(cchan) then begin
     minfreq=(cchan-12.5)*1.28
     maxfreq=minfreq+30.72
     plot,[0],[0],/nodata,xrange=[minfreq,maxfreq],yrange=[0,s[1]],xtit="Freq (MHz)",ytit="Time",xsty=1,ysty=1,title=tit
  end else begin
  plot,[0],[0],/nodata,xrange=[0,s[0]],yrange=[0,s[1]],xtit="Channel",ytit="Time",xsty=1,ysty=1,title=tit
end

  c=!p.clip ; get the plot clipping range

  xsize=c[2]-c[0]-1
  ysize=c[3]-c[1]-1


  f=congrid(p,xsize,ysize,cubic=-0.5) ; resize the array

;  if keyword_set(mina) then begin
;     m=where(f lt mina)
;     if m[0] ne -1 then f[m]=mina
;  end
;
;  if keyword_set(maxa) then begin
;     m=where(f gt maxa)
;     if m[0] ne -1 then f[m]=maxa
;  end

  device,decomposed=0
  if keyword_set(maxa) then mx=maxa else mx=max(f)
  if keyword_set(mina) then mn=mina else mn=min(f)
  fs=(f-mn)*255/(mx-mn)
  tv,fs,c[0]+1,c[1]+1

  ;replot the axes
  if keyword_set(cchan) then begin
     minfreq=(cchan-12.5)*1.28
     maxfreq=minfreq+30.72
     plot,[0],[0],/nodata,xrange=[169.6,200.32],yrange=[0,s[1]],xtit="Freq (MHz)",ytit="Time",/noerase,xsty=1,ysty=1,title=tit
  end else begin
  plot,[0],[0],/nodata,xrange=[0,s[0]],yrange=[0,s[1]],xtit="Channel",ytit="Time",/noerase,xsty=1,ysty=1,title=tit
end

end

;********************************************;

; dtvcombine
; This will read in arrays of data for analyzing the DTV interference
; and combine them into a massive array for plotting/analyzing.

pro dtvcombine,infile,auto,cross=cross,ant=ant,onebl=onebl,bls=bls

  read_File_list,infile,list

  if keyword_set(onebl) and not keyword_set(ant) then begin
     print,"Need to specify an antenna for the baseline"
     return
  end

  if keyword_set(ant) then begin
     dummy=fltarr(768,2016)+!values.f_nan
     ccsel=ccindex(ant,dummy,lab=bls)
  
     if keyword_set(onebl) then begin

        w=where(bls eq onebl)
        if w eq -1 then begin
           print,"Invalid baseline"
           return
        end
        ccsel=ccsel[w]
     end
  end

  scantimes=intarr(n_elements(list))

  for i=0,n_elements(list)-1 do begin

     if list[i] eq "gap" then begin
        scantimes[i]=30
     end else begin
        openr,lun,'/r3/X13/mwadas1/xraid/raw_data/'+list[i]+'/'+list[i]+'_das1.lacspc',/get
        info=fstat(lun)
        fsize=info.size
        free_lun,lun
        scantimes[i]=fsize/(192L*64L*4L)
        if scantimes[i] eq 0 then begin & scantimes[i]=30 & print,"No data for file "+list[i] & list[i]="gap" & end
     end
  end


  totscans=total(scantimes)

  if keyword_set(ant) then ninps=1 else ninps=64
  auto=fltarr(768,ninps,totscans)
  
  if arg_present(cross) then begin
     if keyword_set(onebl) then ncross=1 else if keyword_set(ant) then ncross=63 else ncross=2016
     cross=complexarr(768,ncross,totscans)
  end


  index=0L
  for i=0L,n_elements(list)-1 do begin

     if list[i] eq "gap" then begin
        data=fltarr(768,64,scantimes[i])+!values.f_nan
        if keyword_set(ant) then data=data[*,ant,*]
     end else begin

        print,"Reading file "+string(i+1,format='(I3)')+" of "+string(n_elements(list),format='(I3)')
        ntime=scantimes[i]

        da1=fltarr(192,64,ntime,/nozero)
        da2=fltarr(192,64,ntime,/nozero)
        da3=fltarr(192,64,ntime,/nozero)
        da4=fltarr(192,64,ntime,/nozero)

        openr,lun,'/r3/X13/mwadas1/xraid/raw_data/'+list[i]+'/'+list[i]+'_das1.lacspc',/get
        readu,lun,da1
        free_lun,lun
        ;d1=read_binary('/r3/X13/mwadas1/xraid/raw_data/'+list[i]+'/'+list[i]+'_das1.lacspc',data_type=4,data_dims=[192,64,ntime])
     
        openr,lun,'/r3/X13/mwadas2/xraid/raw_data/'+list[i]+'/'+list[i]+'_das2.lacspc',/get
        readu,lun,da2
        free_lun,lun
        ;d2=read_binary('/r3/X13/mwadas2/xraid/raw_data/'+list[i]+'/'+list[i]+'_das2.lacspc',data_type=4,data_dims=[192,64,ntime])

        openr,lun,'/r3/X13/mwadas3/xraid/raw_data/'+list[i]+'/'+list[i]+'_das3.lacspc',/get
        readu,lun,da3
        free_lun,lun
        ;d3=read_binary('/r3/X13/mwadas3/xraid/raw_data/'+list[i]+'/'+list[i]+'_das3.lacspc',data_type=4,data_dims=[192,64,ntime])
     
        openr,lun,'/r3/X13/mwadas4/xraid/raw_data/'+list[i]+'/'+list[i]+'_das4.lacspc',/get
        readu,lun,da4
        free_lun,lun
        ;d4=read_binary('/r3/X13/mwadas4/xraid/raw_data/'+list[i]+'/'+list[i]+'_das4.lacspc',data_type=4,data_dims=[192,64,ntime])

        if keyword_set(ant) then begin
           data=fltarr(768,1,ntime)
           data[0:191,0,*]=da1[*,ant,*]
           data[192:383,0,*]=da2[*,ant,*]
           data[384:575,0,*]=da3[*,ant,*]
           data[576:767,0,*]=da4[*,ant,*]
        end else begin
           data=fltarr(768,64,ntime)
           data[0:191,*,*]=da1
           data[192:383,*,*]=da2
           data[384:575,*,*]=da3
           data[576:767,*,*]=da4
        end
        da1=0
        da2=0
        da3=0
        da4=0
        data=reorder_coarse(data,145)
     end

     ;if i eq 0 then auto=data else auto=[[[auto]],[[data]]]
     auto[*,*,index:index+scantimes[i]-1]=data
     data=0

     if arg_present(cross) then begin
        if list[i] eq "gap" then begin
           data=complexarr(768,2016,30)+complex(!values.f_nan,+!values.f_nan)
           if keyword_set(ant) then data=data[*,ccsel,*]
        end else begin

           dc1=complexarr(192,2016,ntime,/nozero)
           dc2=complexarr(192,2016,ntime,/nozero)
           dc3=complexarr(192,2016,ntime,/nozero)
           dc4=complexarr(192,2016,ntime,/nozero)

           openr,lun,'/r3/X13/mwadas1/xraid/raw_data/'+list[i]+'/'+list[i]+'_das1.lccspc',/get
           readu,lun,dc1
           free_lun,lun
           ;d1=read_binary('/r3/X13/mwadas1/xraid/raw_data/'+list[i]+'/'+list[i]+'_das1.lccspc',data_type=6,data_dims=[192,2016,ntime])
           
           openr,lun,'/r3/X13/mwadas2/xraid/raw_data/'+list[i]+'/'+list[i]+'_das2.lccspc',/get
           readu,lun,dc2
           free_lun,lun
           ;d2=read_binary('/r3/X13/mwadas2/xraid/raw_data/'+list[i]+'/'+list[i]+'_das2.lccspc',data_type=6,data_dims=[192,2016,ntime])
           
           openr,lun,'/r3/X13/mwadas3/xraid/raw_data/'+list[i]+'/'+list[i]+'_das3.lccspc',/get
           readu,lun,dc3
           free_lun,lun
           ;d3=read_binary('/r3/X13/mwadas3/xraid/raw_data/'+list[i]+'/'+list[i]+'_das3.lccspc',data_type=6,data_dims=[192,2016,ntime])
           
           openr,lun,'/r3/X13/mwadas4/xraid/raw_data/'+list[i]+'/'+list[i]+'_das4.lccspc',/get
           readu,lun,dc4
           free_lun,lun
           ;d4=read_binary('/r3/X13/mwadas4/xraid/raw_data/'+list[i]+'/'+list[i]+'_das4.lccspc',data_type=6,data_dims=[192,2016,ntime])

           data=complexarr(768,ncross,ntime)
           if keyword_set(ant) then begin
              
              data[0:191,*,*]=dc1[*,ccsel,*]
              data[192:383,*,*]=dc2[*,ccsel,*]
              data[384:575,*,*]=dc3[*,ccsel,*]
              data[576:767,*,*]=dc4[*,ccsel,*]
           end else begin
              data[0:191,*,*]=dc1
              data[192:383,*,*]=dc2
              data[384:575,*,*]=dc3
              data[576:767,*,*]=dc4
           end
           ;data=[dc1,dc2,dc3,dc4]
           ;if keyword_set(ant) then data=data[*,ccsel,*]
           data=reorder_coarse(data,145)
           dc1=0
           dc2=0
           dc3=0
           dc4=0
        end

        ;if i eq 0 then cross=data else cross=[[[cross]],[[data]]]
        cross[*,*,index:(index+scantimes[i]-1)]=data
        data=0
     end
     index+=scantimes[i]
  end

end

;********************************************;

; read_auto
; This function will read in a .LACSPC file and
; return it as an array


function read_auto,infile,nchans=nchans,ninps=ninps

  if not keyword_set(nchans) then nchans=192
  if not keyword_set(ninps) then ninps=64

  openr,lun,infile,/get
  info=fstat(lun)
  fsize=info.size
  nscans=fsize/(long(nchans)*ninps*4L)
  free_lun,lun
  r=read_binary(infile,data_type=4,data_dims=[nchans,ninps,nscans])
  return,r
end

;********************************************;

; read_cross
; This function will read in a .LCCSPC file and
; return it as an array

function read_cross,infile,nchans=nchans,ninps=ninps

  if not keyword_set(nchans) then nchans=192
  if not keyword_set(ninps) then ninps=64

  openr,lun,infile,/get
  info=fstat(lun)
  fsize=info.size
  nscans=fsize/(long(nchans)*ninps*(ninps-1)*4L)
  free_lun,lun
  r=read_binary(infile,data_type=6,data_dims=[nchans,ninps*(ninps-1)/2,nscans])
  return,r
end
