    using PyPlot

# define some useful functions

# conversion of bulk string to double
function convertToFloat64(str::Array{ByteString,1},len::Int64)
    # str: string data
    # len: length of each datum in character
    
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

# load data from file and return a 2D array of double
function loadData(fileName::ASCIIString,len::Int64)
    # fileName: file name containing the organized data
    # len:      length of each datum in character

    # load the data in a string array
    f=open(fileName,"r") # should check if the file is open (isopen)
    dataSTR=readlines(f)           # should check if there actually are data (isempty)
    close(f)

    # convert the string data to Float64 array
    convertToFloat64(dataSTR,len)
end

# display data on the actual grid (taking into account the pixel coordinates)
function imshowData(figNum::Int64,_t::Array{Cdouble,1},_h::Array{Cdouble,1},Z::Array{Cdouble,2};_norm=:NoNorm,_vmin=0.0,_vmax=1.0,_edgecolors="face",_shading="None")
    # get the focus on the figure or create a figure
    figure(figNum)

    # get the current axis
    ax = axes()

    # display image if possible
    if (length(_t),length(_h))==size(Z) || (length(_h),length(_t))==size(Z)
        pcolormesh(repmat(_t,1,length(_h))',repmat(_h,1,length(_t)),Z,edgecolors=_edgecolors,shading=_shading,norm=matplotlib[:colors][:Normalize](vmin=_vmin,vmax=_vmax)) #shading="gouraud"
    else
        @sprintf "Cannot display the image because of a size mismatch: (%i,%i)!=%i,%i" size(Z,1) size(Z,2) length(_t) length(_h)
    end

    #return axis handler
    ax
end



###############################################
##### script for display purposes: meteo ######
###############################################



#load the data
theta = loadData("output/theta.dat",25)
ua = loadData("output/ua.dat",25)
va = loadData("output/va.dat",25)
t  = loadData("output/time.dat",8)
t=squeeze(t,1)
kt = loadData("output/kt.dat",25)
kmt = loadData("output/kmt.dat",25)
kht = loadData("output/kht.dat",25)
ri = loadData("output/Ri.dat",25)
h= [0.0,   10,   20,   30,   40, 50,   60,   70,   80,   90, 100,  120,  140,  160,  180,  200,  230,  260,  300,  350,  400,  450,  500,  550,  600,  650,  700,  800,  900, 1000,  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]


#now we can plot
figure(1)
imshowData(1,t,h,theta-273.15,_norm=:Normalize,_vmin=minimum(theta)-273.15,_vmax=maximum(theta)-273.15,_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("temperature [C]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
title("potential temperature \$\\theta\$ [C]")
xlabel("time [h]")
ylabel("height [m]")

figure(2)
imshowData(2,t,h,ua,_norm=:Normalize,_vmin=minimum(ua),_vmax=maximum(ua),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("velocity [m s\$^{-1}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
title("wind speed component \$u\$ [m s\$^{-1}\$]")
xlabel("time [d]")
ylabel("height [m]")

figure(3)
imshowData(3,t,h,va,_norm=:Normalize,_vmin=minimum(va),_vmax=maximum(va),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("velocity [m s\$^{-1}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
title("wind speed component \$v\$ [m s\$^{-1}\$]")
xlabel("time [d]")
ylabel("height [m]")



######## K-values
figure(4)
imshowData(4,t,0.5*(h[1:end-1]+h[2:end]),kt,_norm=:Normalize,_vmin=minimum(kt),_vmax=maximum(kt),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("K [m\$^{2}\$ s\$^{-1}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "K-value (%1.2e,%1.2e) m\$^{2}\$ s\$^{-1}\$" minimum(kt) maximum(kt)
title(s)
xlabel("time [d]")
ylabel("height [m]")

figure(5)
imshowData(5,t,0.5*(h[1:end-1]+h[2:end]),kmt,_norm=:Normalize,_vmin=minimum(kmt),_vmax=maximum(kmt),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("K\$_m\$ [m\$^{2}\$ s\$^{-1}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "K\$_m\$-value (%1.2e,%1.2e) m\$^{2}\$ s\$^{-1}\$" minimum(kmt) maximum(kmt)
title(s)
xlabel("time [d]")
ylabel("height [m]")

figure(6)
imshowData(6,t,0.5*(h[1:end-1]+h[2:end]),kht,_norm=:Normalize,_vmin=minimum(kht),_vmax=maximum(kht),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("K\$_h\$ [m\$^{2}\$ s\$^{-1}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "K\$_h\$-value (%1.2e,%1.2e) m\$^{2}\$ s\$^{-1}\$" minimum(kht) maximum(kht)
title(s)
xlabel("time [d]")
ylabel("height [m]")

##### Richardson number
figure(7)
minri=minimum(ri)
maxri=maximum(ri)
imshowData(7,t,0.5*(h[1:end-1]+h[2:end]),ri,_norm=:Normalize,_vmin=-1.0,_vmax=300.0,_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("Ri [\$\\tilde{}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "Richardson number (%1.2e,%1.2e)" minimum(ri) maximum(ri)
title(s)
xlabel("time [d]")
ylabel("height [m]")


######### wind speed
### component U
figure(8)
plot(ua[:,1],h, linewidth=2,marker="o")
plot(ua[:,5],h, linewidth=2,marker="o")
plot(ua[:,9],h, linewidth=2,marker="o")
plot(ua[:,13],h, linewidth=2,marker="o")
plot(ua[:,17],h, linewidth=2,marker="o")
plot(ua[:,21],h, linewidth=2,marker="o")
title("wind speed component U: day 1")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity [m s\$^{-1}\$]")
ylabel("Height [m]")



figure(9)
plot(ua[:,1+100],h, linewidth=2,marker="o")
plot(ua[:,5+100],h, linewidth=2,marker="o")
plot(ua[:,9+100],h, linewidth=2,marker="o")
plot(ua[:,13+100],h, linewidth=2,marker="o")
plot(ua[:,17+100],h, linewidth=2,marker="o")
plot(ua[:,21+100],h, linewidth=2,marker="o")
title("wind speed component \$u\$: day 5")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity [m s\$^{-1}\$]")
ylabel("Height [m]")





######### wind speed
### component V
figure(10)
plot(va[:,1],h, linewidth=2,marker="o")
plot(va[:,5],h, linewidth=2,marker="o")
plot(va[:,9],h, linewidth=2,marker="o")
plot(va[:,13],h, linewidth=2,marker="o")
plot(va[:,17],h, linewidth=2,marker="o")
plot(va[:,21],h, linewidth=2,marker="o")
title("wind speed component V: day 1")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity [m s\$^{-1}\$]")
ylabel("Height [m]")


figure(11)
plot(va[:,1+100],h, linewidth=2,marker="o")
plot(va[:,5+100],h, linewidth=2,marker="o")
plot(va[:,9+100],h, linewidth=2,marker="o")
plot(va[:,13+100],h, linewidth=2,marker="o")
plot(va[:,17+100],h, linewidth=2,marker="o")
plot(va[:,21+100],h, linewidth=2,marker="o")
title("wind speed component \$v\$: day 5")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("velocity [m s\$^{-1}\$]")
ylabel("Height [m]")




##### potential temperature
figure(12)
plot(theta[:,1]-273.15,h, linewidth=2,marker="o")
plot(theta[:,5]-273.15,h, linewidth=2,marker="o")
plot(theta[:,9]-273.15,h, linewidth=2,marker="o")
plot(theta[:,13]-273.15,h, linewidth=2,marker="o")
plot(theta[:,17]-273.15,h, linewidth=2,marker="o")
plot(theta[:,21]-273.15,h, linewidth=2,marker="o")
title("first day (potential temperature \$\\theta\$)")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("potential temperature \$\\theta\$ [C]")
ylabel("Height [m]")

figure(13)
plot(theta[:,1+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,5+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,9+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,13+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,17+100]-273.15,h, linewidth=2,marker="o")
plot(theta[:,21+100]-273.15,h, linewidth=2,marker="o")
title("fifth day (potential temperature \$\\theta\$)")
legend(["0h","4h","8h","12h","16h","20h"])
xlabel("potential temperature \$\\theta\$ [C]")
ylabel("Height [m]")




##### temperature at the ground  boundary
figure(14)
plot(24.*reshape(t,121),reshape(theta[1,:],121)-273.15, linewidth=2)
title("boundary temperature")
ylabel("potential temperature \$\\theta\$ [C]")
xlabel("time [h]")
