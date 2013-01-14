
FUNCTION color_image,image,r,g,b,missing=missing
IF N_Elements(missing) NE 1 THEN missing=0
dimension=(size(Image,/dimension))[0]
elements=(size(Image,/dimension))[1]
image_c=fltarr(3,dimension,elements)

image_c[0,*,*]=Interpolate(r,Image,missing=missing)
image_c[1,*,*]=Interpolate(g,Image,missing=missing)
image_c[2,*,*]=Interpolate(b,Image,missing=missing)

RETURN,image_c
END


PRO ImageFast,Image,filename=filename,data_directory=data_directory,no_colorbar=no_colorbar,right_colorbar=right_colorbar,left_colorbar=left_colorbar,$
    tick_spacing=tick_spacing,units=units,logarithmic_color=logarithmic_color,reverse_color=reverse_color,hist_equalize=hist_equalize,low=low,high=high,$
    background=background,significant_figures=significant_figures,zero_black=zero_black,resize=resize,color_table=color_table,$
    project_name=project_name,full_range=full_range,transparent=transparent

;Keywords not yet supported:tick_spacing=tick_spacing,units=units,logarithmic_color=logarithmic_color,reverse_color=reverse_color,hist_equalize=hist_equalize
;writes to a png image by default.
;set keyword transparent to make zero values in the final png image transparent

DefsysV,'!black',exists=astr ;check if astronomy system variables have been defined
IF astr EQ 0 THEN astrolib
IF N_Elements(background) EQ 0 THEN background=!black
IF Keyword_Set(full_range) THEN BEGIN high=max(image) & low=min(image) & ENDIF
Image_use=Image
dimension=(size(Image_use,/dimension))[0]
elements=(size(Image_use,/dimension))[1]
IF Keyword_Set(resize) THEN BEGIN
    Image_use=Rebin(Image_use,dimension*resize,elements*resize)
    dimension=(size(Image_use,/dimension))[0]
    elements=(size(Image_use,/dimension))[1]
ENDIF
scale=2.^Ceil(Alog(64./(dimension<elements))/Alog(2))
IF scale GE 2 THEN BEGIN
    dimension*=scale
    elements*=scale
    Image_use=Rebin(Image_use,dimension,elements)
ENDIF ELSE Image_use=Image_use
xsize=dimension
ysize=elements
colorbar_scale=(2.^Ceil(Alog((dimension<elements)/256.)/Alog(2)))>1
colorbar_width=90.*colorbar_scale
IF Keyword_Set(significant_figures) THEN IF significant_figures GT 2 THEN colorbar_width+=(significant_figures-2)*colorbar_scale
xstart=0
ystart=0
;xoffset = 0.5 & yoffset = 0.5
CASE 1 OF
    Keyword_Set(no_colorbar):
    Keyword_Set(left_colorbar):BEGIN
        xsize+=colorbar_width
;       ysize=ysize>256
        xstart+=colorbar_width
    END
    Keyword_Set(right_colorbar): BEGIN
        xsize+=colorbar_width
;       ysize=ysize>256
    END
    ELSE:BEGIN
        ysize+= colorbar_width
;       xsize=xsize>256
        ystart+=colorbar_width
    ENDELSE
ENDCASE

IF N_Elements(filename) EQ 0 THEN filename_use='ImageFast.png' ELSE filename_use=filename
UPNAME=StrUpCase(filename_use)
ppng=strpos(UPNAME,'.PNG')
IF ppng EQ -1 THEN filename_use=filename_use+String('.png')
RootDirectory=rootdir(project_name)

IF Keyword_Set(data_directory) THEN filename_full=filepath(filename_use,Root_dir=RootDirectory,subdir=data_directory) $
    ELSE filename_full=filepath(filename_use,Root_dir=RootDirectory)
;openps, filename_full, xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset

;use XYOUTS to write text
;use CONTOUR for images?
;set keyword /FILL or use /CELL_FILL followed by a second call using /OVERPLOT for the namesake contour plot
;also: LEVELS MAX_VALUE MIN_VALUE C_COLORS
;The MIN_CURVE_SURF function can be used to smoothly interpolate both regularly and irregularly sampled surfaces before contouring
;use TV for basic image output keyword TRUE={1|2|3} interleaves the color over the dimension specified-> set to one as a standard convention
;manipulate !P.Multi for multiple plots/images
;TVBOX - Draw a box of specified size on the image display
;TVCIRCLE - Draw a circle of specified radius on the image display
;TVELLIPSE - Draw an ellipse of specified axes on the image display

nsigma=3.
IF Keyword_Set(logarithmic_color) THEN Image_use[where(Image_use)]=ALog10(Abs(Image_use[where(Image_use)]))
image_hist=Simple_histogram(Image_use,center=meanval,sigma=sigma)
IF N_Elements(low) EQ 0 THEN low=(meanval-nsigma*sigma)>Min(Image_use)
IF N_Elements(high) EQ 0 THEN high=(meanval+nsigma*sigma)<Max(Image_use)
zero_i=where(Image_use EQ 0,n_zero)
low_i=where(Image_use LT low,n_low)
high_i=where(Image_use GT high,n_high)
image_use=Image_use-low
image_use*=255./(high-low)
IF n_low GT 0 THEN image_use[low_i]=0
IF n_high GT 0 THEN image_use[high_i]=255.
r=255.*Sin(1.5*!Pi*findgen(256.)/255.+!Pi/4.)>0
r[0]=0 & r[255]=255.
g=255.*Sin(1.5*!Pi*findgen(256.)/255.-!Pi/4.)>0
g[255]=255.
b=255.*Sin(1.5*!Pi*findgen(256.)/255.+5*!Pi/4)>0
b[255]=255.
endramp=(findgen(64)+1)*4.-1
r[192:*]=r[192:*]>endramp
g[192:*]=g[192:*]>endramp
b[192:*]=b[192:*]>endramp
missing=0
IF Keyword_Set(zero_black) THEN IF n_zero GT 0 THEN image_use[zero_i]=0 

IF Keyword_Set(color_table) THEN BEGIN 
    loadct,color_table,rgb=rgb 
    r=Float(rgb[*,0])
    g=Float(rgb[*,1])
    b=Float(rgb[*,2])
ENDIF
image_c=color_image(image_use,r,g,b,missing=missing)
;image_c is a (3,dimension,elements) array, with the first dimension encoding color
set_plot,'win'
;device, filename=filename_full, /color, bits=8, /times, /isolatin1, xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset, /inches,landscape=0

free=1
;!P.background=background
window,free=free,xsize=xsize,ysize=ysize,retain=2
IF background NE !black THEN BEGIN
    bg=lonarr(3,xsize,ysize)
    bg_color=intarr(3)
    bg[0,*,*]=(bg_color[0]=background mod 256)
    bg[1,*,*]=(bg_color[1]=Floor((background-bg_color[0])/256) mod 256)
    bg[2,*,*]=(bg_color[2]=Floor(((background-bg_color[0])/256 -bg_color[1])/256))
    tv,bg,0,0,true=1
ENDIF

tv,image_c,xstart,ystart,xsize=xsize*2,ysize=ysize*2,true=1
IF not Keyword_Set(no_colorbar) THEN BEGIN
    IF N_Elements(significant_figures) EQ 0 THEN significant_figures=2.
    n_labels=5<(Floor(Alog(dimension<elements)/Alog(2)-3))    
    labels0=findgen(n_labels)*(high-low)/(n_labels-1)+low
    IF Keyword_Set(logarithmic_color) THEN labels0=10.^labels0
    range_test=Round(Alog10(high-low))
    label_precision=10.^(range_test-significant_figures)
    labels1=Round(labels0/label_precision)*label_precision
    neg_test=intarr(n_labels) & FOR i=0,n_labels-1 DO IF labels1[i] LT 0 THEN neg_test[i]=1
    labels=Strarr(n_labels)
    IF (range_test-significant_figures GE 2) OR (range_test LE -2) THEN BEGIN
        FOR i=0,n_labels-1 DO BEGIN
            IF labels1[i] EQ 0 THEN BEGIN 
;                format_code=String(format='("(I",I1,")")',2+significant_figures) 
                labels[i]=Strmid(StrTrim(String(format='(F)',0),2),0,2+significant_figures)
                CONTINUE
            ENDIF
            int_length=8+significant_figures+neg_test[i]
            IF int_length LT 10 THEN format_code=String(format='("(E",I1,".",I1,")")',int_length,significant_figures) $
                ELSE format_code=String(format='("(E",I2,".",I1,")")',int_length,significant_figures)
            labels[i]=String(format=format_code,labels1[i])
        ENDFOR
    ENDIF ELSE BEGIN
        FOR i=0,n_labels-1 DO BEGIN
            IF labels1[i] EQ 0 THEN BEGIN 
;                format_code=String(format='("(I",I1,")")',2+significant_figures) 
                labels[i]=Strmid(StrTrim(String(format='(F)',0),2),0,2+significant_figures)
                CONTINUE
            ENDIF
            int_length=Floor(Alog10(Abs(labels1[i])))
            IF int_length LT 0 THEN int_length=2+Abs(int_length)+significant_figures ELSE int_length+=1
            int_length=int_length>(2+significant_figures)
            int_length+=neg_test[i]
;            IF int_length LT 10 THEN format_code=String(format='("(F",I1,")")',int_length) $
;                ELSE format_code=String(format='("(F",I2,")")',int_length)
            labels[i]=Strmid(StrTrim(String(format='(F)',labels1[i]),2),0,int_length)
        ENDFOR
    
    ENDELSE
    labels=Strtrim(labels,2)
    
;    range_test=Max(Abs(Alog10(Abs(labels[where(labels)]))))
;    IF range_test GE 5 THEN labels=StrTrim(String(format='(G10.2e2)' ,labels),2) ELSE labels=StrTrim(String(format='(G10.7)' ,labels),2)
;    
    
    c_h=16*colorbar_scale
    CASE 1 OF
       Keyword_Set(left_colorbar): BEGIN
         c_w=(256*colorbar_scale)<(ysize-16*colorbar_scale)
           ystart_c=ystart
           IF c_w LT ysize THEN ystart_c+=(ysize-c_w)/2.
           xstart_c=colorbar_width-c_h

           label_x=intarr(n_labels)
           label_y=c_w*indgen(n_labels)/(n_labels-1)+ystart_c-4*colorbar_scale
           orientation=0
       END
       Keyword_Set(right_colorbar):BEGIN
         c_w=(256*colorbar_scale)<(ysize-16*colorbar_scale)
           ystart_c=ystart
           IF c_w LT ysize THEN ystart_c+=(ysize-c_w)/2.
           xstart_c=dimension
           label_x=intarr(n_labels)+dimension+c_h
           label_y=c_w*indgen(n_labels)/(n_labels-1)+ystart_c-4*colorbar_scale

           orientation=0
       END
       ELSE: BEGIN
         c_w=(256*colorbar_scale)<(xsize-16*colorbar_scale)
           xstart_c=xstart
           IF c_w LT xsize THEN xstart_c+=(xsize-c_w)/2.
           ystart_c=colorbar_width-c_h

           label_x=c_w*indgen(n_labels)/(n_labels-1)+xstart_c-4*colorbar_scale
           label_y=ystart_c+intarr(n_labels)
           orientation=-90
       ENDELSE
    ENDCASE

    colorbar=indgen(c_w)#Replicate(1,c_h)*255./(c_w-1.)
    colorbar[0,*]=0
    colorbar[*,0]=0
    colorbar[c_w-1,*]=0
    colorbar[*,c_h-1]=0
    IF Keyword_Set(right_colorbar) OR Keyword_Set(left_colorbar) THEN colorbar=rotate(colorbar,1)
    colorbar_c=color_image(colorbar,r,g,b,missing=missing)
    tv,colorbar_c,xstart_c,ystart_c,true=1
    label_color=Abs((background-!white) mod 256.^3.)
    XYouts,label_x,label_y,labels,orientation=orientation,/device,color=label_color,charsize=colorbar_scale,charthick=colorbar_scale
ENDIF
result=tvrd(true=1)

;IF Keyword_Set(transparent) THEN BEGIN
;    IF n_zero GT 0 THEN BEGIN
;        zero_x=zero_i mod dimension
;        zero_y=Floor(zero_i/dimension)
;        zero_x2=zero_x+xstart
;        zero_y2=zero_y+ystart
;        zero_i2=zero_x2+zero_y2*xsize
;        write_png,filename_full,result,transparent=zero_i2
;    ENDIF
;ENDIF ELSE 
write_png,filename_full,result
wdelete

END


