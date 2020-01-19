PRO find_gr_matchvol

;domain='brisbane'
;domain='sydney'
domain='darwin'

; Set the radars to compare
if domain eq 'brisbane' then begin
   ;rad1='IDR66'
   ;rad1='IDR50'
   rad1='IDR28'
   ;rad2='IDR50'
   ;rad2='IDR28'
   rad2='IDR08'
endif
if domain eq 'sydney' then begin
   ;rad1='IDR71'
   rad1='IDR03'
   ;rad2='IDR03'
   rad2='IDR04'
endif
if domain eq 'darwin' then begin
   rad1='IDR59'
   rad2='IDR63'
endif

; Set the radar scans
vcp1=154
vcp2=154

; Set parameters
rmax=115.    ; maximum range (km)
zmin=2000.   ; minimum height (m)
zmax=12000.  ; maximum height (m)
dmax=500.    ; maximum distance between bins (m)
vfmax=0.1    ; maximum fractional difference in bin volumes
dr=250.      ; range gate spacing (m)
naz=360      ; number of azimuths (m)
daz=1.       ; azimuth spacing (deg)

; Set all the radar properties
if domain eq 'brisbane' then begin
   nrad=4
   rads='IDR'+['50','66','28','08']
   lonr=[152.5390,153.2400,152.9510,152.5770]
   latr=[-27.6080,-27.7178,-29.6220,-25.9574]
   zr=[372.0,174.0,40.0,375.0]
   rname=['Marburg','Brisbane','Grafton','Gympie']
   bwth=[1.9,1.0,1.9,2.0]
endif
if domain eq 'sydney' then begin
   nrad=4
   rads='IDR'+['03','04','40','71']
   rname=['Wollongong','Newcastle','Canberra','Sydney']
   lonr=[150.8752,152.0250,149.5122,151.2094]
   latr=[-34.2625,-32.7300,-35.6614,-33.7008]
   zr=[449.,84.,1384.,222.]
   bwth=[2.0,1.9,1.9,1.0]
endif
if domain eq 'darwin' then begin
   nrad=2
   rads='IDR'+['59','63']
   rname=['CPOL','Berrimah']
   lonr=[131.044,130.925]
   latr=[-12.249,-12.457]
   zr=[50.,51.]
   bwth=[1.0,1.0]
endif

; Set the VCP properties
nvcp=7
vcps=[140,150,151,152,153,154,155]
tilts=STRARR(nvcp)
tilts[0]='0.5,0.9,1.3,1.8,2.4,3.1,4.2,5.6,7.4,10.0,13.3,17.9,23.9,32.0'
tilts[1]='0.5,1.2,1.9,2.7,3.5,4.7,6.0,7.5,9.2,11.0,13.0,16.0,20.0,25.0,32.0'
tilts[2]='0.5,0.8,1.1,1.4,1.9,2.5,3.3,4.4,5.8,7.7,10.3,13.6,18.1,24.1,32.0'
tilts[3]='0.5,0.9,1.1,1.4,1.9,2.5,3.3,4.4,5.8,7.7,10.3,13.6,18.1,24.1,32.0'
tilts[4]='0.5,0.7,0.9,1.2,1.6,2.2,3.0,4.0,5.4,7.2,9.8,13.1,17.7,23.8,32.0'
tilts[5]='0.5,0.9,1.3,1.8,2.4,3.1,4.2,5.6,7.4,10.0,13.3,17.9,23.9,32.0,43.1'
tilts[6]='0.500000,0.890625,1.29688,1.79688,2.39062,3.09375,4.18750,5.59375,7.39062,10.0000,13.2969,17.8906,23.8906,32.0000,43.0938'

; Note the indicies of the two radars
ir1=(WHERE(rads eq rad1))[0]
ir2=(WHERE(rads eq rad2))[0]

; Extract parameters for the two radars
lonr1=lonr[ir1]
latr1=latr[ir1]
h1=zr[ir1]
b1=bwth[ir1]
rmax1=rmax*1000
lonr2=lonr[ir2]
latr2=latr[ir2]
h2=zr[ir2]
b2=bwth[ir2]
rmax2=rmax*1000

; Note the vcp indices
iv1=(WHERE(vcps eq vcp1))[0]
iv2=(WHERE(vcps eq vcp2))[0]

; Set the range values
nr1=ROUND(rmax1/dr)
nr2=ROUND(rmax2/dr)
r1=dr/2.+FINDGEN(nr1)*dr
r2=dr/2.+FINDGEN(nr2)*dr

; Set the azimuth angle values
az=FINDGEN(naz)*daz

; Note the elevation angles and number of tilts for each radar
el1=FLOAT(STRSPLIT(tilts[iv1],',',/extract))
nel1=N_ELEMENTS(el1)
el2=FLOAT(STRSPLIT(tilts[iv2],',',/extract))
nel2=N_ELEMENTS(el2)

; Adjust rmax2 to account for maximum volume of radar 1
; (at elevation angle of 0deg)
vmax1=(rmax1^2)*dr*!dtor*(b1+daz)*!dtor*b1
vmax2=vmax1*(1+vfmax)
rmax2=MIN([rmax2,SQRT(vmax2/(dr*!dtor*(b1+daz)*!dtor*b1))])

PRINT,rmax1,rmax2

; Set the initial domain centre as the midpoint between the two radars
; in lon-lat space
lon0=0.5*(lonr1+lonr2)
lat0=0.5*(latr1+latr2)

; Iterate to get the true domain centre
done=0
x0p=0
y0p=0
while done eq 0 do begin

   ; Set the projection
   smap=MAP_PROJ_INIT('Transverse Mercator',ellipsoid=24, $
                      center_longitude=lon0,center_latitude=lat0)

   ; Convert radar coordinates to x,y
   res=MAP_PROJ_FORWARD(lonr1,latr1,map_structure=smap)
   xr1=res[0]
   yr1=res[1]
   res=MAP_PROJ_FORWARD(lonr2,latr2,map_structure=smap)
   xr2=res[0]
   yr2=res[1]

   ; Find the domain centre
   x0=(xr1+xr2)/2.
   y0=(yr1+yr2)/2.

   ; Convert to lat-lon
   res=MAP_PROJ_INVERSE(x0,y0,map_structure=smap)
   lon0=res[0]
   lat0=res[1]

   if MIN(ABS([x0-x0p,y0-y0p])) lt 1. then done=1

   xp0=x0
   yp0=y0

endwhile

; Note the effective Earth radius
a=smap.a
e2=smap.e2
b=a*SQRT(1-e2)
invf=1/(1-b/a)
ac=SQRT((a^4*COS(!dtor*lat0)^2+b^4*SIN(!dtor*lat0)^2)/ $
        (a^2*COS(!dtor*lat0)^2+b^2*SIN(!dtor*lat0)^2))
ae=(4/3.)*ac

; Note the coordinates of the two radars
res=MAP_PROJ_FORWARD(lonr1,latr1,map_structure=smap)
xr1=res[0]
yr1=res[1]
res=MAP_PROJ_FORWARD(lonr2,latr2,map_structure=smap)
xr2=res[0]
yr2=res[1]

print,SQRT((xr1-xr2)^2+(yr1-yr2)^2)/2.

; Create 3D coordinate arrays for each radar
r11=REBIN(r1,nr1,naz,nel1,/sample)
az11=REBIN(REFORM(az,1,naz,1),nr1,naz,nel1,/sample)
el11=REBIN(REFORM(el1,1,1,nel1),nr1,naz,nel1,/sample)
r22=REBIN(r2,nr2,naz,nel2,/sample)
az22=REBIN(REFORM(az,1,naz,1),nr2,naz,nel2,/sample)
el22=REBIN(REFORM(el2,1,1,nel2),nr2,naz,nel2,/sample)

; Compute the approximate volume of each bin for each radar
vol11=dr*r11*!dtor*(b1+daz*COS(!dtor*el11))*r11*!dtor*b1
vol22=dr*r22*!dtor*(b2+daz*COS(!dtor*el22))*r22*!dtor*b2

; Compute the x,y,z coordinates of radar 1
z11=SQRT(r11^2+(ae+h1)^2+2*r11*(ae+h1)*SIN(!dtor*el11))-ae
s11=ae*ASIN(r11*COS(!dtor*el11)/(ae+z11))
x11=xr1+s11*COS(!dtor*(90-az11))
y11=yr1+s11*SIN(!dtor*(90-az11))

; Compute the x,y,z coordinates of radar 2
z22=SQRT(r22^2+(ae+h2)^2+2*r22*(ae+h2)*SIN(!dtor*el22))-ae
s22=ae*ASIN(r22*COS(!dtor*el22)/(ae+z22))
x22=xr2+s22*COS(!dtor*(90-az22))
y22=yr2+s22*SIN(!dtor*(90-az22))

; Compute the coordinates of radar 1 wrt radar 2
x12=x11-xr2
y12=y11-yr2
z12=z11
s12=SQRT(x12^2+y12^2)
el12=!radeg*ATAN((COS(s12/ae)-(ae+h2)/(ae+z12))/SIN(s12/ae))
r12=(ae+z12)*SIN(s12/ae)/COS(!dtor*el12)
az12=90-!radeg*ATAN(y12,x12)
az12=(az12+360) mod 360  ; make all values +ve

; Compute the coordinates of radar 2 wrt radar 1
x21=x22-xr1
y21=y22-yr1
z21=z22
s21=SQRT(x21^2+y21^2)
el21=!radeg*ATAN((COS(s21/ae)-(ae+h1)/(ae+z21))/SIN(s21/ae))
r21=(ae+z21)*SIN(s21/ae)/COS(!dtor*el21)
az21=90-!radeg*ATAN(y21,x21)
az21=(az21+360) mod 360  ; make all values +ve

; Compute the volume of radar 2 pulses at radar 1 pixels and the
; volume of radar 1 pulses at radar 2 pixels
vol12=dr*(r12^2)*!dtor*(b2+daz*COS(!dtor*el12))*!dtor*b2
vol21=dr*(r21^2)*!dtor*(b1+daz*COS(!dtor*el21))*!dtor*b1

; Compute the volume fractional differences
vf12=ABS(vol12-vol11)/(0.5*(vol12+vol11))
vf21=ABS(vol22-vol21)/(0.5*(vol22+vol21))

; Find all radar 1 pixels which fall within the coverage of radar 2
; and which have a comparable volume
i1=WHERE(el12 ge MIN(el2)-b2/2. and el12 le MAX(el2)+b2/2. and $
         r11 le rmax1+dr/2. and $
         r12 le rmax2+dr/2. and $
         z12 ge zmin and z12 le zmax and $
         vf12 le vfmax,n1)

; Find all radar 2 pixels which fall within the coverage of radar 1
; and which have a comparable volume
i2=WHERE(el21 ge MIN(el1)-b1/2. and el21 le MAX(el1)+b1/2. and $
         r22 le rmax2+dr/2. and $
         r21 le rmax1+dr/2. and $
         z21 ge zmin and z21 le zmax and $
         vf21 le vfmax,n2)

PRINT,n1,n2

; Note the min and max ranges
PRINT,MIN(r11[i1]),MAX(r11[i1]),MIN(r22[i2]),MAX(r22[i2])

; Open a file to store the data
;dir='/media/robwarr/Elements/data/radar/brisbane/lookup/volmatch'
if domain ne 'darwin' then begin
   dir='/media/robwarr/MyBook/data/radar/'+domain+'/lookup/volmatch'
endif else begin
   dir='/media/robwarr/Elements3/data/radar/'+domain+'/lookup/volmatch'
endelse
;dir=dir+'_new'
file=dir+'/'+rad1+'-'+SCROP(vcp1)+'_vs_'+rad2+'-'+SCROP(vcp2)+ $
     '_rmax'+SCROP(rmax,ndecs=0)+'_dr'+SCROP(dr,ndecs=0)+ $
     '_dmax'+SCROP(dmax,ndecs=0)+'.txt'
OPENW,lun,file,/get_lun

count=0L
minv=1000.
maxv=0.
mind1=1000.
maxd1=0.
mind2=1000.
maxd2=0.

sec1=SYSTIME(/seconds)

; Loop over the radar 1 pixels
for i=0,n1-1 do begin

   ;if i mod 1000 eq 0 then print,i

   j1=i1[i]

   ; Determine the distance between this pixel and the radar 2 pixels
   d=SQRT((x22[i2]-x11[j1])^2+(y22[i2]-y11[j1])^2+(z22[i2]-z11[j1])^2)

   ; Determine the differences in volume between this pixel and the 
   ; radar 2 pixels as a fraction of the mean pixel volume
   vf=ABS(vol22[i2]-vol11[j1])/(0.5*(vol22[i2]+vol11[j1]))

   ; Check for overlapping and volume-matched pixels
   imatch=WHERE(d le dmax and vf le vfmax,nmatch)
   if nmatch gt 0 then begin
      if nmatch gt 1 then begin
         tmp=d[imatch]
         iclose=(WHERE(tmp eq MIN(tmp)))[0]
         imatch=imatch[iclose]
      endif else imatch=imatch[0]
      j2=i2[imatch]
      ;print,''
      print,count+1,i+1,n1
      ;print,d[imatch],vf[imatch]
      ;print,r11[j1]/1000.,az11[j1],el11[j1],z11[j1],vol1[j1]/1e9
      ;print,r22[j2]/1000.,az22[j2],el22[j2],z22[j2],vol2[j2]/1e9
      ;print,z11[j1],z22[j2],vol11[j1]/1e9,vol22[j2]/1e9
      count=count+1
      PRINTF,lun,j1,j2
      minv=MIN([minv,vol11[j1]/1e9,vol22[j2]/1e9])
      maxv=MAX([maxv,vol11[j1]/1e9,vol22[j2]/1e9])
      mind1=MIN([mind1,s11[j1]/1000.])
      mind2=MIN([mind2,s22[j2]/1000.])
      maxd1=MAX([maxd1,s11[j1]/1000.])
      maxd2=MAX([maxd2,s22[j2]/1000.])
   endif

endfor

print,count

print,minv,maxv,ROUND(mind1),ROUND(maxd1),ROUND(mind2),ROUND(maxd2)

FREE_LUN,lun

sec2=SYSTIME(/seconds)

print,(sec2-sec1)/60.

END
