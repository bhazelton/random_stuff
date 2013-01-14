FUNCTION clean_source_condense,source_array,radius=radius,simple=simple
IF N_Elements(radius) EQ 0 THEN radius=0.5
simple=1 ;only this form is written yet. Should be decent, though

sx=Reform(source_array[0,*])
sy=Reform(source_array[1,*])
sf=Reform(source_array[3,*])
sf2=Reform(source_array[2,*])
sf_ratio=sf[0]/sf2[0]
ns=N_Elements(sx)

group_id=fltarr(ns)-1
IF Keyword_Set(simple) THEN BEGIN
    g_id=0
    FOR si=0,ns-1 DO BEGIN 
        IF group_id[si] GE 0 THEN CONTINUE ;skip sources already grouped
        si_use=where(group_id EQ -1,n_use)        
        dx=sx[si]-sx[si_use]
        dy=sy[si]-sy[si_use]
        dr=sqrt(dx^2.+dy^2.)
        group_i=where(dr LE radius,n_group) ;guaranteed at least one
        group_id[si_use[group_i]]=g_id
        g_id+=1
    ENDFOR
ENDIF ELSE BEGIN
    ;use some form of he astro idl library star grouping
ENDELSE

hgroup=histogram(group_id,binsize=1,min=0,reverse_ind=gri)
ng=max(group_id)
source_array_condensed=fltarr(4,ng)
FOR gi=0,ng-1 DO BEGIN
    si_g=gri[gri[gi]:gri[gi+1]-1]; guaranteed at least one source per group
    sx_weight=Total(sx[si_g]*sf[si_g])/Total(sf[si_g])
    sy_weight=Total(sy[si_g]*sf[si_g])/Total(sf[si_g])
    sf_tot=Total(sf2[si_g])
    
;    sf_weight=Total(sf[si_g])
    source_array_condensed[*,gi]=[sx_weight,sy_weight,sf_tot,sf_tot*sf_ratio]
ENDFOR

order=Reverse(sort(Reform(source_array_condensed[3,*])))
source_array_condensed=source_array_condensed[*,order]
RETURN,source_array_condensed
END