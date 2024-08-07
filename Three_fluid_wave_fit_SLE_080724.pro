pro Three_fluid_wave_fit_SLE_080724,SH1,SH0,coll,gplanet,rho_obs_ratio,windowme,titleme,$
latitude,rotrate,nugas0,kzgas0val,gas1gas0phasediff,returned_lx,returned_lz,decfeb,$
returnval,SH2,coll2,rho_obs_ratio2,gas2gas0phasediff,gravityacoustic
;SLE 08/07/24
;This code fits to the identified wave perturbations (density
;pertrubation amplituderatio and phase offset)
;for 3 gas species - CO2, N2 and O.
;It returns the best fit gravity wave parameters and
;an estimation of the goodness of fit. 
;The code is derived from Cui, J., Y. Lian, and I. C. F. Müller-Wodarg
;(2013), Compositional effects in Titan’s thermospheric gravity waves,
;Geophys. Res. Lett., 40, 43–47, doi:10.1029/2012GL054621
;and adapted to Mars by England, S. L., G. Liu, E. Yiğit,
;P. R. Mahaffy, M. Elrod, M. Benna, H. Nakagawa, N. Terada, and
;B. Jakosky (2017), MAVEN NGIMS observations of atmospheric gravity
;waves in the Martian thermosphere, J. Geophys. Res. Space Physics,
;122, doi:10.1002/ 2016JA023475.
  
;Inputs: Scale height for 3 species (SH0, 1, 2) in m
;CO2-N2 and CO2-O collision rates, in s-1 (coll and coll2)
;acceleration due to gravity at the relevant altitude for the
;observation (gplanet) in m/s2
;The ratio of the density amplitudes in CO2-N2 and CO2-O
;(rho_obs_ratio and rho_obs_ratio2)
;Inputs for the headings for plots generated (windowme, titleme and decfeb)
;The latitude of the observation in degrees (latitude)
;The rotation rate of the planet (rotate) in radians/s
;The viscosity coefficient for CO2 (nugas0)
;Phase difference between CO2 and N2, CO2 and O(gas1gas0phasediff and
;gas2gas0phasediff) in degrees
;gravityacoustic a toggle that should be set to 'gravity'

;Outputs: 
;best fit horizontal, vertical wavelength in m (returned_lx,returned_lz)
;residual % of fit vs model (returnval)
;plots of the model parameters - amplitude ratios and phase
;differences vs horizontal and
;vertical wavelengths


set_plot,'x'
!p.multi=[0,2,2,0]

;Start with the observed values...
SHgas0=sh0
SHgas1=sh1
SHgas2=sh2
collgas1=coll
collgas2=coll2

;finally, sounds speeds squared..
sound2_gas0=1.3*gplanet*SHgas0
sound2_gas1=1.4*gplanet*SHgas1
sound2_gas2=1.*gplanet*SHgas2;atomic

NBgas0=(gplanet*(1.3-1.)^0.5)/(sound2_gas0^0.5) ;per Hines, Cui 2013. 

gas1_gas0_obs_ratio=rho_obs_ratio
gas2_gas0_obs_ratio=rho_obs_ratio2
;from latitude, calculate coriolis_f
coriolis_f=2.*rotrate*sin(latitude/180.*!pi)

;Run CO2-N2 first...

shortk=2.*!pi/(2000.*1000.)
longk=2.*!pi/(10.*2000.)

kres=shortk/5.
nk=fix(longk/kres)+1
lx_array=fltarr(nk,nk)
lz_array=fltarr(nk,nk)
ltot_array=fltarr(nk,nk)
ktot_array=fltarr(nk,nk)
kz_array=fltarr(nk,nk)
kx_array=fltarr(nk,nk)
omega_array=fltarr(nk,nk)

;full-solution approach....
real_ratio_array=fltarr(nk,nk)
imag_ratio_array=fltarr(nk,nk)
mag_ratio_Array=fltarr(nk,nk)
im_over_re_ratio_array=fltarr(nk,nk)
for ik=0,nk-1 do begin
  kx=ik*kres+shortk
  lx=2.*!pi/kx
  for jk=0,nk-1 do begin
    ;find k, lambda etc
    kz=jk*kres+shortk
    lz=2.*!pi/kz
    lx_array(ik,jk)=lx
    lz_array(ik,jk)=lz
    ktot=(kx^2+kz^2)^0.5
    ktot_Array(ik,jk)=ktot
    ltot=1./((1./(lx^2)+1./(lz^2))^0.5)
    ltot_array(ik,jk)=ltot

;If gravity waves, use this
if gravityacoustic eq 'gravity' then begin
    ;find omega - from the shorter disp rel
;    omega2=(nbgas0^2*kx^2+coriolis_F^2*(kz^2+1./(4.*SHgas0^2)))/(kx^2+kz^2+1./(4.*SHgas0^2))  ;per Cui 2014
    omega2=(nbgas0^2*kx^2)/(kx^2+kz^2+1./(4.*SHgas0^2))  ;per Cui 2013
;    omega2=(nbgas0^2*ktot^2)/(kx^2+kz^2+1./(4.*SHgas0^2))  ;per Cui 2013
    omega=omega2^0.5
    omega_array(ik,jk)=omega
endif


    kz=-1.*kz  
    kx=-1.*kx

  kz_Array(ik,jk)=kz
  kx_array(ik,jk)=kx

;all of these should be stored as complex numbers
;trying gas1
;a
rval=0.
ival=omega
justa=dcomplex(rval,ival)
;b
rval=0.
ival=-1.*kx
justb=dcomplex(rval,ival)
;c
rval=-1./(2.*SHgas0)
ival=-1.*kz
justc=dcomplex(rval,ival)
;d
rval=0.
ival=-1.*kx*sound2_gas0/1.3 ;Gamma is 1.3 for CO2
justd=dcomplex(rval,ival)
;e
rval=(-1./(2.*SHgas0))*sound2_gas0/1.3
ival=(-1.*kz)*sound2_gas0/1.3
juste=dcomplex(rval,ival)
;g
rval=gplanet
ival=0.
justg=dcomplex(rval,ival)
;h
rval=1.3/SHgas0-1./SHgas0
ival=0.
justh=dcomplex(rval,ival)
;i
rval=0.
ival=-1.*1.3*omega
justi=dcomplex(rval,ival)
;k
rval=0.
ival=-1.*kx
justk=dcomplex(rval,ival)
;l
rval=-1./(2.*SHgas1)
ival=-1.*kz
justl=dcomplex(rval,ival)
;n
rval=collgas1
ival=omega
justn=dcomplex(rval,ival)
;m
rval=0.
ival=-1.*kx*sound2_gas1/1.4
justm=dcomplex(rval,ival)
;X
rval=-1l*collgas1
ival=0.
justX=dcomplex(rval,ival)
;q
rval=(-1./(2.*SHgas1))*sound2_gas1/1.4
ival=-1.*kz*sound2_gas1/1.4
justq=dcomplex(rval,ival)
;Y
rval=-1l*collgas1
ival=0.
justY=dcomplex(rval,ival)
;r
rval=1.4/SHgas1-1./SHgas1
ival=0.
justr=dcomplex(rval,ival)
;s
rval=0.
ival=-1.*1.4*omega
justs=dcomplex(rval,ival)


;write this out equation by equation, using values for everything - just like I have in excel...
;define each equation as...
;r0, p0, u0, w0, r1, p1, u1, w1
;testing an alternative method
eq1=[justa,0.,justb,justc,0.,0.,0.,0.]
eq2=[0.,justd,justa,0.,0.,0.,0.,0.]
eq3=[justg,juste,0.,justa,0.,0.,0.,0.]
eq4=[justi,justa,0.,justh,0.,0.,0.,0.]
eq5=[0.,0.,0.,0.,justa,0,justk,justl]
eq6=[0.,0.,justx,0.,0.,justm,justn,0.]
eq7=[0.,0.,0.,justy,justg,justq,0.,justn]
eq8=[0.,0.,0.,0.,justs,justa,0.,justr]

;now, follow my solution...
;take eq5 and multiply it by 6(6)/5(6)
eq9=eq5*eq6(6)/eq5(6)
;eq6-eq9
eq10=eq6-eq9
;take eq8 and *7(5)/8(5)
eq11=eq8*eq7(5)/eq8(5)
;eq7-eq11
eq12=eq7-eq11
;take eq8 and multiply by  10(5)/8(5)
eq13=eq8*eq10(5)/eq8(5)
;eq10-eq13
eq14=eq10-eq13
;take eq14 and 12(7)/14(7)
eq15=eq14*eq12(7)/eq14(7)
;eq12-eq15
eq16=eq12-eq15
;take eq1 and * 16(2)/1(2)
eq17=eq1*eq16(2)/eq1(2)
;eq16-17
eq18=eq16-eq17
;take eq3 and * 4(1)/3(1)
eq19=eq3*eq4(1)/eq3(1)
;eq4-eq19
eq20=eq4-eq19
;eq20*18(3)/20(3)
eq21=eq20*eq18(3)/eq20(3)
;21-18
eq22=eq21-eq18
;now, the only significant values left are in rho0 and rho1.
complex_rho0=-1.*eq22(0)
complex_rho1=eq22(4)
;stop

complex_rho_ratio=complex_rho0/complex_rho1
;save the real & imaginary components
real_ratio_array(ik,jk)=1l/real_part(complex_rho_ratio)
imag_ratio_array(ik,jk)=imaginary(complex_rho_ratio)
mag_ratio_array(ik,jk)=abs(complex_rho_Ratio)
im_over_re_ratio_array(ik,jk)=imaginary(complex_rho_ratio)/real_part(complex_rho_ratio)
  endfor
endfor

;There's now an issue if real_ratio_array = 0., where the value has gone to infinite...
if min(real_ratio_array) le 0. then real_ratio_array(where(real_ratio_array le 0.))=10.


;now do O-CO2
;set up arrays to save into - to not overwrite O
real_ratio_array2=fltarr(nk,nk)
imag_ratio_array2=fltarr(nk,nk)
mag_ratio_Array2=fltarr(nk,nk)
im_over_re_ratio_array2=fltarr(nk,nk)
;now copy all of O over N2
SHgas1=SHgas2
collgas1=collgas2
sound2_gas1=sound2_gas2
;NBgas1=NBgas2
gas1_gas0_obs_ratio=gas2_gas0_obs_ratio
;now run existing code
for ik=0,nk-1 do begin
  kx=ik*kres+shortk
  lx=2.*!pi/kx
  for jk=0,nk-1 do begin
    ;find k, lambda etc
    kz=jk*kres+shortk
    lz=2.*!pi/kz
    lx_array(ik,jk)=lx
    lz_array(ik,jk)=lz
    ktot=(kx^2+kz^2)^0.5
    ktot_Array(ik,jk)=ktot
    ltot=1./((1./(lx^2)+1./(lz^2))^0.5)
    ltot_array(ik,jk)=ltot

    ;If gravity waves, use this
    if gravityacoustic eq 'gravity' then begin
      omega2=(nbgas0^2*kx^2)/(kx^2+kz^2+1./(4.*SHgas0^2))  ;per Cui 2013
      omega=omega2^0.5
      omega_array(ik,jk)=omega
    endif


    kz=-1.*kz  
    kx=-1.*kx

    kz_Array(ik,jk)=kz
    kx_array(ik,jk)=kx

    ;now.... try calculating all of the individual coefficients that are used....
    ;all of these should be stored as complex numbers
    ;trying gas1
    ;a
    rval=0.
    ival=omega
    justa=dcomplex(rval,ival)
    ;b
    rval=0.
    ival=-1.*kx
    justb=dcomplex(rval,ival)
    ;c
    rval=-1./(2.*SHgas0)
    ival=-1.*kz
    justc=dcomplex(rval,ival)
    ;d
    rval=0.
    ival=-1.*kx*sound2_gas0/1.3 ;Gamma is 1.3 for CO2
    justd=dcomplex(rval,ival)
    ;e
    rval=(-1./(2.*SHgas0))*sound2_gas0/1.3
    ival=(-1.*kz)*sound2_gas0/1.3
    juste=dcomplex(rval,ival)
    ;g
    rval=gplanet
    ival=0.
    justg=dcomplex(rval,ival)
    ;h
    rval=1.3/SHgas0-1./SHgas0
    ival=0.
    justh=dcomplex(rval,ival)
    ;i
    rval=0.
    ival=-1.*1.3*omega
    justi=dcomplex(rval,ival)
    ;k
    rval=0.
    ival=-1.*kx
    justk=dcomplex(rval,ival)
    ;l
    rval=-1./(2.*SHgas1)
    ival=-1.*kz
    justl=dcomplex(rval,ival)
    ;n
    rval=collgas1
    ival=omega
    justn=dcomplex(rval,ival)
    ;m
    rval=0.
    ival=-1.*kx*sound2_gas1/1.4
    justm=dcomplex(rval,ival)
    ;X
    rval=-1l*collgas1
    ival=0.
    justX=dcomplex(rval,ival)
    ;q
    rval=(-1./(2.*SHgas1))*sound2_gas1/1.4
    ival=-1.*kz*sound2_gas1/1.4
    justq=dcomplex(rval,ival)
    ;Y
    rval=-1l*collgas1
    ival=0.
    justY=dcomplex(rval,ival)
    ;r
    rval=1.4/SHgas1-1./SHgas1
    ival=0.
    justr=dcomplex(rval,ival)
    ;s
    rval=0.
    ival=-1.*1.4*omega
    justs=dcomplex(rval,ival)

 
    ;write this out equation by equation, using values for everything - 
    ;define each equation as...
    ;r0, p0, u0, w0, r1, p1, u1, w1
    eq1=[justa,0.,justb,justc,0.,0.,0.,0.]
    eq2=[0.,justd,justa,0.,0.,0.,0.,0.]
    eq3=[justg,juste,0.,justa,0.,0.,0.,0.]
    eq4=[justi,justa,0.,justh,0.,0.,0.,0.]
    eq5=[0.,0.,0.,0.,justa,0,justk,justl]
    eq6=[0.,0.,justx,0.,0.,justm,justn,0.]
    eq7=[0.,0.,0.,justy,justg,justq,0.,justn]
    eq8=[0.,0.,0.,0.,justs,justa,0.,justr]


    ;take eq5 and multiply it by 6(6)/5(6)
    eq9=eq5*eq6(6)/eq5(6)
    ;eq6-eq9
    eq10=eq6-eq9
    ;take eq8 and *7(5)/8(5)
    eq11=eq8*eq7(5)/eq8(5)
    ;eq7-eq11
    eq12=eq7-eq11
    ;take eq8 and multiply by  10(5)/8(5)
    eq13=eq8*eq10(5)/eq8(5)
    ;eq10-eq13
    eq14=eq10-eq13
    ;take eq14 and 12(7)/14(7)
    eq15=eq14*eq12(7)/eq14(7)
    ;eq12-eq15
    eq16=eq12-eq15
    ;take eq1 and * 16(2)/1(2)
    eq17=eq1*eq16(2)/eq1(2)
    ;eq16-17
    eq18=eq16-eq17
    ;take eq3 and * 4(1)/3(1)
    eq19=eq3*eq4(1)/eq3(1)
    ;eq4-eq19
    eq20=eq4-eq19
    ;eq20*18(3)/20(3)
    eq21=eq20*eq18(3)/eq20(3)
    ;21-18
    eq22=eq21-eq18
    ;now, the only significant values left are in rho0 and rho1.
    complex_rho0=-1.*eq22(0)
    complex_rho1=eq22(4)

    
    complex_rho_ratio=complex_rho0/complex_rho1
    ;save the real & imaginary components
    real_ratio_array2(ik,jk)=1l/real_part(complex_rho_ratio);SLE 07/15/16
    imag_ratio_array2(ik,jk)=imaginary(complex_rho_ratio)
    mag_ratio_array2(ik,jk)=abs(complex_rho_Ratio)
    im_over_re_ratio_array2(ik,jk)=imaginary(complex_rho_ratio)/real_part(complex_rho_ratio)

  endfor
endfor

;as above
if min(real_ratio_array2) le 0. then real_ratio_array2(where(real_ratio_array2 le 0.))=10.

;end O-CO2
;reset N2 values...
SHgas1=sh1
collgas1=coll
sound2_gas1=1.4*gplanet*SHgas1
NBgas1=(gplanet*(1.4-1.)^0.5)/(sound2_gas0^0.5)
gas1_gas0_obs_ratio=rho_obs_ratio


haxis=fltarr(nk)
vaxis=fltarr(nk)
for ik=0,nk-1 do begin
  haxis(ik)=(ik*kres+shortk);/2.2
  vaxis(ik)=(ik*kres+shortk)
endfor


delta_Array=fltarr(nk,nk)

phase_diff_angles=abs(imag_ratio_array/!pi*180.)
phase_diff_angles2=abs(imag_ratio_array2/!pi*180.)

delta_phase_Array=fltarr(nk,nk)

;Update Delta array to include O as well
delta_3_Array=fltarr(nk,nk)
for ik=0,nk-1 do begin
  for jk=0,nk-1 do begin
    sigma_val=real_ratio_array(ik,jk)

    angle=phase_diff_angles(ik,jk)/180.*!pi
    gas1gas0phasediff_rad=gas1gas0phasediff/180.*!pi
;O
sigma_val2=real_ratio_array2(ik,jk)

angle2=phase_diff_angles2(ik,jk)/180.*!pi
gas2gas0phasediff_rad=gas2gas0phasediff/180.*!pi
;
;per Cui et al., these should be *1/4 as we're summing 4 values.
    delta_3_array(ik,jk)=(1./4.)*(((sigma_val-gas1_gas0_obs_ratio)/gas1_gas0_obs_ratio)^2.+$
      ((angle-gas1gas0phasediff_rad)/gas1gas0phasediff_rad)^2.+$
      ((sigma_val2-gas2_gas0_obs_ratio)/gas2_gas0_obs_ratio)^2.+$
      ((angle2-gas2gas0phasediff_rad)/gas2gas0phasediff_rad)^2.$
      )^0.5
  endfor
endfor



returnval=min(delta_3_array)*100.
mindelta=where(delta_3_array eq min(delta_3_array))
returned_lx=-1.*(2.*!pi/kx_array(mindelta))
returned_lz=-1.*(2.*!pi/kz_array(mindelta))





!p.multi=[0,2,3,0]
set_plot,'ps'
loadct,39
xtickme=[2.*!pi/2000000.,2.*!pi/1000000.,2.*!pi/500000.,2.*!pi/200000.,2.*!pi/100000.,2.*!pi/50000.,2.*!pi/20000.]
xtickmenames=['2000','1000','500','200','100','50','20']
device,filename='Plots/'+strcompress(decfeb,/remove_all)+'_Cui_mag_ratio_figure_'+gravityacoustic+'.ps',/color

contour,real_Ratio_array,haxis,vaxis,title='Real part of ratio CO2/N2',/fill,levels=findgen(200)/199.*5.,$
  xticks=5,xtickv=xtickme,xtickn=xtickmenames,/xlog,yticks=6,ytickv=xtickme,ytickn=xtickmenames,/ylog,xtitle='lx, km',ytitle='lz, km'
contour,real_Ratio_array2,haxis,vaxis,title='Real part of ratio CO2/O',/fill,levels=findgen(200)/199.*5.,$
    xticks=5,xtickv=xtickme,xtickn=xtickmenames,/xlog,yticks=6,ytickv=xtickme,ytickn=xtickmenames,/ylog,xtitle='lx, km',ytitle='lz, km'

contour,phase_diff_angles,haxis,vaxis,title='CO2-N2 Phase Difference, degrees'+titleme,/fill,$
  xticks=5,xtickv=xtickme,xtickn=xtickmenames,/xlog,yticks=6,ytickv=xtickme,ytickn=xtickmenames,/ylog,xtitle='lx, km',ytitle='lz, km',levels=findgen(91)/90.*40.
contour,phase_diff_angles2,haxis,vaxis,title='CO2-O Phase Difference, degrees'+titleme,/fill,$
    xticks=5,xtickv=xtickme,xtickn=xtickmenames,/xlog,yticks=6,ytickv=xtickme,ytickn=xtickmenames,/ylog,xtitle='lx, km',ytitle='lz, km',levels=findgen(91)/90.*40.

contour,delta_3_array,haxis,vaxis,title='Delta parameter, with O'+titleme,/fill,$
  xticks=5,xtickv=xtickme,xtickn=xtickmenames,/xlog,yticks=6,ytickv=xtickme,ytickn=xtickmenames,/ylog,xtitle='lx, km',ytitle='lz, km',levels=findgen(256)/255.*6.

plot,findgen(10),/nodata,xs=4,ys=4
device,/close

set_plot,'x'

!p.position=0
!p.multi=0
;stop
return
stop
end
