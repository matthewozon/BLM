    using PyPlot

function convertToFloat64(str::Array{ByteString,1},len::Int64)
    #compute sizes
    numTimeStep=round(Int64,size(str,1))
    #println(summary(length(str[1])/len))
    numLevel=round(Int64,length(str[1])/len)
    #allocate a 2D array for the data
    strDouble=Array(Cdouble,numLevel,numTimeStep)

    #println(summary(numTimeStep))
    #println(summary(numLevel))
    #println(summary(strDouble))
    #println(numTimeStep," ",numLevel)
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

#convert strings to arrays of float64
theta=convertToFloat64(thetaSTR,25)
ua=convertToFloat64(uaSTR,25)
va=convertToFloat64(vaSTR,25)
t=convertToFloat64(tSTR,8)
h= [0.0,   10,   20,   30,   40, 50,   60,   70,   80,   90, 100,  120,  140,  160,  180,  200,  230,  260,  300,  350,  400,  450,  500,  550,  600,  650,  700,  800,  900, 1000,  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]


#now we can plot
figure(1)
imshow(theta)
title("potential temperature")


figure(2)
imshow(ua)
title("wind speed component u")


figure(3)
imshow(va)
title("wind speed component v")


figure(4)
plot(ua[:,1],h)
plot(ua[:,25],h)
plot(ua[:,49],h)
plot(ua[:,73],h)
plot(ua[:,97],h)
plot(ua[:,121],h)
title("once a day (wind speed component U)")

figure(5)
plot(ua[:,1],h)
plot(ua[:,5],h)
plot(ua[:,9],h)
plot(ua[:,13],h)
plot(ua[:,17],h)
plot(ua[:,21],h)
title("first day (wind speed component U)")


figure(6)
plot(ua[:,1+100],h)
plot(ua[:,5+100],h)
plot(ua[:,9+100],h)
plot(ua[:,13+100],h)
plot(ua[:,17+100],h)
plot(ua[:,21+100],h)
title("fifth day (wind speed component U)")



figure(7)
plot(va[:,1],h)
plot(va[:,25],h)
plot(va[:,49],h)
plot(va[:,73],h)
plot(va[:,97],h)
plot(va[:,121],h)
title("once a day (wind speed component V)")

figure(8)
plot(va[:,1],h)
plot(va[:,5],h)
plot(va[:,9],h)
plot(va[:,13],h)
plot(va[:,17],h)
plot(va[:,21],h)
title("first day (wind speed component V)")


figure(9)
plot(va[:,1+100],h)
plot(va[:,5+100],h)
plot(va[:,9+100],h)
plot(va[:,13+100],h)
plot(va[:,17+100],h)
plot(va[:,21+100],h)
title("fifth day (wind speed component V)")



figure(10)
plot(theta[:,1],h)
plot(theta[:,25],h)
plot(theta[:,49],h)
plot(theta[:,73],h)
plot(theta[:,97],h)
plot(theta[:,121],h)
title("once a day (potential temperature theta)")

figure(11)
plot(theta[:,1],h)
plot(theta[:,5],h)
plot(theta[:,9],h)
plot(theta[:,13],h)
plot(theta[:,17],h)
plot(theta[:,21],h)
title("first day (potential temperature theta)")


figure(12)
plot(theta[:,1+100],h)
plot(theta[:,5+100],h)
plot(theta[:,9+100],h)
plot(theta[:,13+100],h)
plot(theta[:,17+100],h)
plot(theta[:,21+100],h)
title("fifth day (potential temperature theta)")
