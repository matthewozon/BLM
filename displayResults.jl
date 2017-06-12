    using PyPlot

function convertToFloat64(str::Array{ByteString,1},len::Int64)
    #compute sizes
    numTimeStep=round(Int64,size(str,1))
    numLevel=round(Int64,length(str[1])/len)
    #allocate a 2D array for the data
    strDouble=Array(Cdouble,numLevel,numTimeStep)

    #convert each element
    for i in 1:numTimeStep
        for j in 1:numLevel
            idxS=1+len*(j-1)
            idxE=len*j
            strDouble[j,i]=parse(Cdouble,str[i][idxS:idxE])
        end
    end

    #return
    strDouble
end


#load the files
f=open("output/theta.dat","r")
thetaSTR=readlines(f)
close(f)

f=open("output/ua.dat","r")
uaSTR=readlines(f)
close(f)

f=open("output/va.dat","r")
vaSTR=readlines(f)
close(f)

f=open("output/time.dat","r")
tSTR=readlines(f)
close(f)

f=open("output/h.dat","r")
hSTR=readlines(f)
close(f)

f=open("output/kt.dat","r")
ktSTR=readlines(f)
close(f)

f=open("output/kmt.dat","r")
kmtSTR=readlines(f)
close(f)

f=open("output/kht.dat","r")
khtSTR=readlines(f)
close(f)

f=open("output/Ri.dat","r")
riSTR=readlines(f)
close(f)

f=open("output/emi_iso.dat","r")
isoSTR=readlines(f)
close(f)

f=open("output/emi_alp.dat","r")
alpSTR=readlines(f)
close(f)

f=open("output/OH.dat","r")
ohSTR=readlines(f)
close(f)

f=open("output/HO2.dat","r")
ho2STR=readlines(f)
close(f)

f=open("output/H2SO4.dat","r")
h2so4STR=readlines(f)
close(f)

f=open("output/isoprene.dat","r")
isopreneSTR=readlines(f)
close(f)

f=open("output/alpha.dat","r")
alphapSTR=readlines(f)
close(f)

f=open("output/ELVOC.dat","r")
elvocSTR=readlines(f)
close(f)

#convert strings to arrays of float64
theta=convertToFloat64(thetaSTR,25)
ua=convertToFloat64(uaSTR,25)
va=convertToFloat64(vaSTR,25)
t=convertToFloat64(tSTR,8)
kt=convertToFloat64(ktSTR,25)
kmt=convertToFloat64(kmtSTR,25)
kht=convertToFloat64(khtSTR,25)
ri=convertToFloat64(riSTR,25)
iso=convertToFloat64(isoSTR,25)
alp=convertToFloat64(alpSTR,25)
oh=convertToFloat64(ohSTR,25)
ho2=convertToFloat64(ho2STR,25)
h2so4=convertToFloat64(h2so4STR,25)
isoprene=convertToFloat64(isopreneSTR,25)
alphap=convertToFloat64(alphapSTR,25)
elvoc=convertToFloat64(elvocSTR,25)
h= [0.0,   10,   20,   30,   40, 50,   60,   70,   80,   90, 100,  120,  140,  160,  180,  200,  230,  260,  300,  350,  400,  450,  500,  550,  600,  650,  700,  800,  900, 1000,  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]


#plotting variables
locs = [0.0, 10.0, 20.0, 30.0, 40.0, 49.0]
labels = (string(h[50]),string(h[41]),string(h[31]),string(h[21]),string(h[11]),string(h[1]))

if false
#now we can plot
figure(1)
#imshow(flipdim(theta,1)-273.15)
contourf(repmat(t,50,1),repmat(h,1,121),theta-273.15,1000)
colorbar()
title("potential temperature (C)")
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)

figure(2)
#imshow(flipdim(ua,1))
contourf(repmat(t,50,1),repmat(h,1,121),ua,1000)
colorbar()
title("wind speed component u (m s^{-1})")
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)

figure(3)
#imshow(flipdim(va,1))
contourf(repmat(t,50,1),repmat(h,1,121),va,1000)
colorbar()
title("wind speed component v (m s^{-1})")
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)



######## K-values
figure(4)
#imshow(flipdim(kt,1))
contourf(repmat(t,49,1),repmat(h[1:end-1],1,121),kt,1000)
colorbar()
s = @sprintf "K-value (m^{2} s^{-1}) (%3.2f,%5.2f)" minimum(kt) maximum(kt)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)

figure(5)
#imshow(flipdim(kmt,1))
contourf(repmat(t,49,1),repmat(h[1:end-1],1,121),kmt,1000)
colorbar()
s = @sprintf "K_m-value (m^{2} s^{-1}) (%3.2f,%5.2f)" minimum(kmt) maximum(kmt)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)

figure(6)
#imshow(flipdim(kht,1))
contourf(repmat(t,49,1),repmat(h[1:end-1],1,121),kht,1000)
colorbar()
s = @sprintf "K_h-value (m^{2} s^{-1}) (%3.2f,%5.2f)" minimum(kht) maximum(kht)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)
end

##### Richardson number
figure(7)
#imshow(log(1.0+log(1.0-minimum(ri)+flipdim(ri,1))))
#imshow(flipdim(ri,1),vmin=-1.0,vmax=300.0)
contourf(repmat(t,49,1),repmat(h[1:end-1],1,121),ri,1000,vmin=-1.0,vmax=300.0)
minri=minimum(ri)
maxri=maximum(ri)
#dRi=(maxri-minri)/9.0
#ticksColor=(minri,minri+dRi,minri+2.0dRi,minri+3.0dRi,minri+4.0dRi,minri+5.0dRi,minri+6.0dRi,minri+7.0dRi,minri+8.0dRi,minri+9.0dRi)#, format="%7.2f"
colorbar()#(norm=LogNorm(vmin=minri, vmax=maxri))#(extend="max")#(ticks=ticksColor,format="%4.2f")
s = @sprintf "Richardson number (%3.2f,%7.2f)" minimum(ri) maximum(ri)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)



######### wind speed
### component U
#figure(8)
#plot(ua[:,1],h, linewidth=2,marker="o")
#plot(ua[:,25],h, linewidth=2,marker="o")
#plot(ua[:,49],h, linewidth=2,marker="o")
#plot(ua[:,73],h, linewidth=2,marker="o")
#plot(ua[:,97],h, linewidth=2,marker="o")
#plot(ua[:,121],h, linewidth=2,marker="o")
#title("wind speed component U")
#legend(["0d","1d","2d","3d","4d","5d"])
#xlabel("velocity (m s^{-1})")
#ylabel("Height (m)")


figure(9)
plot(ua[:,1],h, linewidth=2,marker="o")
plot(ua[:,5],h, linewidth=2,marker="o")
plot(ua[:,9],h, linewidth=2,marker="o")
plot(ua[:,13],h, linewidth=2,marker="o")
plot(ua[:,17],h, linewidth=2,marker="o")
plot(ua[:,21],h, linewidth=2,marker="o")
title("wind speed component U: day 1")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity (m s^{-1})")
ylabel("Height (m)")



figure(10)
plot(ua[:,1+100],h, linewidth=2,marker="o")
plot(ua[:,5+100],h, linewidth=2,marker="o")
plot(ua[:,9+100],h, linewidth=2,marker="o")
plot(ua[:,13+100],h, linewidth=2,marker="o")
plot(ua[:,17+100],h, linewidth=2,marker="o")
plot(ua[:,21+100],h, linewidth=2,marker="o")
title("wind speed component U: day 5")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity (m s^{-1})")
ylabel("Height (m)")





######### wind speed
### component V
#figure(11)
#plot(va[:,1],h, linewidth=2,marker="o")
#plot(va[:,25],h, linewidth=2,marker="o")
#plot(va[:,49],h, linewidth=2,marker="o")
#plot(va[:,73],h, linewidth=2,marker="o")
#plot(va[:,97],h, linewidth=2,marker="o")
#plot(va[:,121],h, linewidth=2,marker="o")
#title("wind speed component V")
#legend(["0d","1d","2d","3d","4d","5d"])
#xlabel("velocity (m s^{-1})")
#ylabel("Height (m)")


figure(12)
plot(va[:,1],h, linewidth=2,marker="o")
plot(va[:,5],h, linewidth=2,marker="o")
plot(va[:,9],h, linewidth=2,marker="o")
plot(va[:,13],h, linewidth=2,marker="o")
plot(va[:,17],h, linewidth=2,marker="o")
plot(va[:,21],h, linewidth=2,marker="o")
title("wind speed component V: day 1")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity (m s^{-1})")
ylabel("Height (m)")


figure(13)
plot(va[:,1+100],h, linewidth=2,marker="o")
plot(va[:,5+100],h, linewidth=2,marker="o")
plot(va[:,9+100],h, linewidth=2,marker="o")
plot(va[:,13+100],h, linewidth=2,marker="o")
plot(va[:,17+100],h, linewidth=2,marker="o")
plot(va[:,21+100],h, linewidth=2,marker="o")
title("wind speed component V: day 5")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity (m s^{-1})")
ylabel("Height (m)")




##### potential temperature
#figure(14)
#plot(theta[:,1]-273.15,h, linewidth=2,marker="o")
#plot(theta[:,25]-273.15,h, linewidth=2,marker="o")
#plot(theta[:,49]-273.15,h, linewidth=2,marker="o")
#plot(theta[:,73]-273.15,h, linewidth=2,marker="o")
#plot(theta[:,97]-273.15,h, linewidth=2,marker="o")
#plot(theta[:,121]-273.15,h, linewidth=2,marker="o")
#title("once a day (potential temperature theta)")
#legend(["0d","1d","2d","3d","4d","5d"])
#xlabel("potential temperature theta (C)")
#ylabel("Height (m)")

figure(15)
plot(theta[:,1]-273.15,h, linewidth=2,marker="o")
plot(theta[:,5]-273.15,h, linewidth=2,marker="o")
plot(theta[:,9]-273.15,h, linewidth=2,marker="o")
plot(theta[:,13]-273.15,h, linewidth=2,marker="o")
plot(theta[:,17]-273.15,h, linewidth=2,marker="o")
plot(theta[:,21]-273.15,h, linewidth=2,marker="o")
title("first day (potential temperature theta)")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("potential temperature theta (C)")
ylabel("Height (m)")

figure(16)
plot(theta[:,1+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,5+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,9+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,13+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,17+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,21+100]-273.15,h, linewidth=2,marker="o")
title("fifth day (potential temperature theta)")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("potential temperature theta (C)")
ylabel("Height (m)")




##### temperature at the ground  boundary
figure(17)
plot(24.*reshape(t,121),reshape(theta[1,:],121)-273.15, linewidth=2)
title("boundary temperature")
ylabel("potential temperature theta (C)")
xlabel("time (h)")



figure(18)
plot(squeeze(t,1),squeeze(iso,1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(alp,1), linewidth=2,marker="o")
title("emission rates")
legend(["isoprene","alpha"])
xlabel("time")
ylabel("emission rate")



figure(19)
contourf(repmat(t,50,1),repmat(h,1,121),oh,1000)
colorbar()
s = @sprintf "OH (cm^{-3}) (%3.2f,%10.2f)" minimum(oh) maximum(oh)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)


figure(20)
contourf(repmat(t,50,1),repmat(h,1,121),ho2,1000)
colorbar()
s = @sprintf "HO2 (cm^{-3}) (%3.2f,%10.2f)" minimum(ho2) maximum(ho2)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)


figure(21)
contourf(repmat(t,50,1),repmat(h,1,121),h2so4,1000)
colorbar()
s = @sprintf "H2SO4 (cm^{-3}) (%3.2f,%10.2f)" minimum(h2so4) maximum(h2so4)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)


figure(22)
contourf(repmat(t,50,1),repmat(h,1,121),isoprene,1000)
colorbar()
s = @sprintf "isoprene (cm^{-3}) (%3.2f,%10.2f)" minimum(isoprene) maximum(isoprene)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)


figure(23)
contourf(repmat(t,50,1),repmat(h,1,121),alphap,1000)
colorbar()
s = @sprintf "alpha-pinene (cm^{-3}) (%3.2f,%10.2f)" minimum(alphap) maximum(alphap)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)


figure(24)
contourf(repmat(t,50,1),repmat(h,1,121),elvoc,1000)
colorbar()
s = @sprintf "ELVOC (cm^{-3}) (%3.2f,%10.2f)" minimum(elvoc) maximum(elvoc)
title(s)
xlabel("time (h)")
ylabel("height (m)")
#yticks(locs,labels)



figure(25)
s = @sprintf "OH (cm^{-3}) (%3.2f,%10.2f)" minimum(oh) maximum(oh)
title(s)
plot(squeeze(t,1),squeeze(oh[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(oh[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("oh [molec cm^{-3}]")
legend(["10 m","50 m"])


figure(26)
s = @sprintf "HO2 (cm^{-3}) (%3.2f,%10.2f)" minimum(ho2) maximum(ho2)
title(s)
plot(squeeze(t,1),squeeze(ho2[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(ho2[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("oh [molec cm^{-3}]")
legend(["10 m","50 m"])


figure(27)
s = @sprintf "H2SO4 (cm^{-3}) (%3.2f,%10.2f)" minimum(h2so4) maximum(h2so4)
title(s)
plot(squeeze(t,1),squeeze(h2so4[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(h2so4[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("oh [molec cm^{-3}]")
legend(["10 m","50 m"])


figure(28)
s = @sprintf "isoprene (cm^{-3}) (%3.2f,%10.2f)" minimum(isoprene) maximum(isoprene)
title(s)
plot(squeeze(t,1),squeeze(isoprene[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(isoprene[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("isoprene [molec cm^{-3}]")
legend(["10 m","50 m"])

figure(29)
s = @sprintf "alphap (cm^{-3}) (%3.2f,%10.2f)" minimum(alphap) maximum(alphap)
title(s)
plot(squeeze(t,1),squeeze(alphap[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(alphap[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("isoprene [molec cm^{-3}]")
legend(["10 m","50 m"])


figure(30)
s = @sprintf "elvoc (cm^{-3}) (%3.2f,%10.2f)" minimum(elvoc) maximum(elvoc)
title(s)
plot(squeeze(t,1),squeeze(elvoc[2,:],1), linewidth=2,marker="o")
plot(squeeze(t,1),squeeze(elvoc[6,:],1), linewidth=2,marker="o")
xlabel("time")
ylabel("isoprene [molec cm^{-3}]")
legend(["10 m","50 m"])


#isoprene=convertToFloat64(isopreneSTR,25)
#alphap=convertToFloat64(alphapSTR,25)
#elvoc=convertToFloat64(elvocSTR,25)
