function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end


@testset "parEstimation.jl" begin
    light_scans = [1, 2, 3, 4, 5, 6]
    heavy_scans = [1, 2, 3, 4, 5, 6]
    light_transitions = UnorderedDictionary(
        ["y3+1","y4+1","y5+1"],
        [
            Float32[1.1273232f7, 8.631513f6, 5.4651805f6, 3.0412648f6, 1.6092901f6, 1.024025f6],
            Float32[4.8260845f6, 3.7873248f6, 2.3309442f6, 1.3229029f6, 704660.44, 436613.62],
            Float32[5.183467f6, 3.9857915f6, 2.6672412f6, 1.4236808f6, 772209.56, 519827.88]
        ]
    )
    heavy_transitions = UnorderedDictionary{String, Vector{Float32}}()
    for (key, value) in pairs(light_transitions)
        insert!(heavy_transitions, key, 100*value)
    end

    par, goodness_of_fit = fitPAR(light_scans, heavy_scans, light_transitions, heavy_transitions)
    @test Tol(par, 100.0)
    @test abs(goodness_of_fit-0.0) < 0.001

    X = Float32[]
    for (key, value) in pairs(light_transitions)
        append!(X, value)
    end
    X = reshape(X, (length(X), 1))
    y = Float32[]
    for (key, value) in pairs(heavy_transitions)
        append!(y, value)
    end

    #Robust model should be equivalent to OLS model
    @test Tol((X\y)[1], par)

    #Add interference to the light transitions and try again
    light_transitions["y3+1"] .+= maximum(light_transitions["y3+1"])    
    par, goodness_of_fit = fitPAR(light_scans, heavy_scans, light_transitions, heavy_transitions)
    
    X = Float32[]
    for (key, value) in pairs(light_transitions)
        append!(X, value)
    end
    X = reshape(X, (length(X), 1))
    y = Float32[]
    for (key, value) in pairs(heavy_transitions)
        append!(y, value)
    end

    #In the presence of interference, the robust model should no longer be equivalent to OLS model
    @test Tol((X\y)[1], par) != true
    
    @test size(quant) == (5, 5)
    @test Set(quant.sequence) == Set(["GLPNVQR"
    "LQAVTDDHIR"
    "VAVWGNK"
    "VAIDAGYR"
    #"VEYLDDR"
    "VGVNGFGR"
    #"GYILQAK"
    ])
end
#end
#=
X = zeros(6*3, 1)
for (i, pair) in enumerate(pairs(light_transitions))
    X[((i-1)*6 + 1):i*6, 1] = pair[2]
end
X[:,2] .= 1.0
y = Float32[]
for (i, pair) in enumerate(pairs(heavy_transitions))
    for value in pair[2]
        if i == 2
            if length(y)>9
                push!(y, 1000000000 + (value))
            else
                push!(y, value)
            end
        else
            push!(y, (value) + Float32(rand([-1, 1])[1]*rand(1)[1]*100000000))
        end
    end
end
model = rlm(X/mean(X), y/mean(X), MMEstimator{TukeyLoss}(), initial_scale=:mad)
#model = rlm(X[:,1]/mean(X[:,1]), y/mean(X[:,1]), MMEstimator{TukeyLoss}(), initial_scale=:mad)
#model = rlm(X[:,1]/mean(X[:,1]), y/mean(X[:,1]), MMEstimator{TukeyLoss}(), initial_scale=:mad)
plot(X[X[:,1].!=0.0,1], y[X[:,1].!=0.0], seriestype = :scatter)

model = rlm(X[X[:,1].!=0.0,:], y[X[:,1].!=0.0], MMEstimator{TukeyLoss}(), initial_scale=:mad)

model = rlm(X, y, MMEstimator{TukeyLoss}(), initial_scale=:mad)

plot(X[:,1], y, seriestype = :scatter)
plot!([0.0, maximum(X)], [0.0, coef(model)[1]*maximum(X)])

X = [342.54638671875; 0.0; 368.5027160644531; 0.0; 0.0; 507.86566162109375; 375.87689208984375; 0.0; 0.0; 400.209716796875; 0.0; 0.0; 339.0276794433594; 0.0; 320.4367980957031; 0.0; 362.7376403808594; 380.0841369628906; 0.0; 419.9341735839844; 562.9187622070312; 406.54083251953125; 0.0; 0.0; 291.9933776855469; 0.0;;]
y = [52852.3984375, 232431.03125, 222513.171875, 375204.53125, 358891.84375, 325395.21875, 291480.21875, 359446.84375, 361143.71875, 303827.21875, 317379.5625, 266394.40625, 281501.34375, 37803.0390625, 190922.625, 155942.828125, 316945.34375, 296424.46875, 262111.359375, 221025.9375, 293807.75, 340294.5, 220664.671875, 256177.0625, 182155.015625, 202120.890625]


X = [0.0; 2492.72412109375; 6826.6474609375; 8723.4248046875; 8205.0224609375; 3938.034423828125; 1363.9224853515625; 572.5691528320312; 2752.986328125; 7992.8544921875; 10192.841796875; 7889.11083984375; 4371.8017578125; 2252.360595703125; 0.0; 1310.2606201171875; 4996.7021484375; 5562.02685546875; 4165.08544921875; 3060.96923828125; 1642.3353271484375; 0.0; 1231.5587158203125; 2408.059814453125; 4880.72265625; 3892.404541015625; 2066.1025390625; 1085.94482421875;;]
y = [44749.25, 343252.71875, 1.131552e6, 1.582338e6, 1.416228e6, 953178.25, 470382.0, 46541.828125, 458083.875, 1.378756875e6, 1.791954875e6, 1.652650375e6, 1.052758e6, 513532.25, 0.0, 286253.09375, 860109.1875, 1.119175125e6, 1.019547375e6, 675746.625, 323832.5, 15461.53125, 184265.875, 574117.5, 774431.4375, 740790.75, 460860.5625, 218290.078125]

X = [1164.6656494140625; 1306.4122314453125; 840.1860961914062; 579.1719970703125; 641.4588623046875; 0.0; 656.5302734375; 0.0; 0.0; 2382.042724609375; 2023.518310546875; 1246.87255859375; 1000.4410400390625; 1437.4708251953125; 697.4054565429688; 1390.4549560546875; 929.9027709960938; 485.04931640625; 135289.375; 128612.375; 91589.2890625; 67215.9765625; 48425.70703125; 32135.6484375; 21423.4609375; 12018.896484375; 7678.91015625;;]
y = [57729.890625, 65582.0703125, 64046.47265625, 91722.453125, 87955.3046875, 99730.6953125, 122570.796875, 112786.265625, 88998.125, 59671.453125, 63266.63671875, 84103.0703125, 79927.40625, 112120.9140625, 132137.109375, 117590.78125, 98385.4921875, 85419.0859375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
a= X[(X[:,1].!=0.0) .& (y.!=0.0),:]
b=y[(X[:,1].!=0.0) .& (y.!=0.0)]
function getModel(X::Matrix{T}, y::Vector{T}, I::BitVector) where T <: AbstractFloat
    rlm(@view(X[I,:])./mean(y), @view(y[I])./mean(y), MMEstimator{TukeyLoss}(), initial_scale=:mad)
end

function getModel(X::Matrix{T}, y::Vector{T}, I::BitVector) where T <: AbstractFloat
    rlm(X[I,:]./mean(y),y[I]./mean(y), TauEstimator{TukeyLoss}(), initial_scale=:mad)
end

X = [6309.74658203125; 8303.84765625; 9002.759765625; 10999.056640625; 14785.7138671875; 16989.86328125; 14809.75390625; 13639.9296875; 9815.6142578125; 6352.00537109375; 4893.71435546875; 735.6868286132812; 1012.072265625; 1440.4442138671875; 2637.63671875; 3835.33544921875; 4302.349609375; 2882.954345703125; 2374.57958984375; 1236.435546875; 1032.617431640625; 629.5667114257812; 835.488525390625; 1980.774169921875; 1508.6524658203125; 4485.21337890625; 5193.28759765625; 7666.80322265625; 6466.359375; 4543.419921875; 2555.9970703125; 1549.7501220703125; 1086.2591552734375; 8203.455078125; 14674.314453125; 23357.873046875; 40105.94140625; 56674.41796875; 67405.2421875; 59418.9140625; 46138.6328125; 31431.806640625; 17982.40625; 9409.5576171875;;]
y = [258337.765625, 414211.40625, 603681.3125, 826679.375, 962697.75, 952656.75, 716837.5, 511010.875, 324589.6875, 170046.828125, 108235.5703125, 107293.984375, 187567.84375, 276593.6875, 351463.75, 406741.59375, 385348.3125, 282297.875, 182287.3125, 122238.65625, 47484.8125, 30582.515625, 182713.578125, 299544.8125, 467758.9375, 617538.6875, 738394.875, 716880.875, 519617.125, 373223.875, 223759.953125, 119972.7421875, 77893.4453125, 1.76037625e6, 2.8511565e6, 4.511994e6, 6.1785615e6, 7.277979e6, 7.0853065e6, 5.321215e6, 3.5316765e6, 2.09896975e6, 1.242870125e6, 766873.0]

X = [0.0; 0.0; 0.0; 334.046630859375; 828.4470825195312; 1024.75830078125; 1771.9224853515625; 2123.4033203125; 2626.49951171875; 3022.205322265625; 2921.033203125; 1922.810546875; 2831.904541015625; 2247.67236328125; 1660.0445556640625; 1060.215087890625; 1094.4190673828125; 767.7047729492188; 0.0; 0.0; 1070.746826171875; 1347.78564453125; 1995.390869140625; 3087.507080078125; 3632.519287109375; 5708.24755859375; 5969.76318359375; 5800.08349609375; 6947.0166015625; 6381.3173828125; 4996.50830078125; 5773.0771484375; 3796.20703125; 3707.3037109375; 1537.9156494140625; 2055.7587890625; 592.7316284179688; 1160.4761962890625; 1504.987548828125; 2809.416259765625; 6089.0625; 6555.42431640625; 10353.548828125; 12548.5361328125; 15966.859375; 15739.80078125; 14726.83203125; 14730.931640625; 15703.1591796875; 11881.478515625; 11566.8486328125; 9399.46875; 7345.91552734375; 5046.32080078125; 0.0; 0.0; 469.2434387207031; 1382.074951171875; 1472.045166015625; 2698.220703125; 4312.33984375; 3984.9765625; 5009.8154296875; 5868.24462890625; 5343.43212890625; 5279.10009765625; 5242.61474609375; 3815.791259765625; 3521.82568359375; 2133.60791015625; 3119.4384765625; 1560.6461181640625;;]
y = [14888.1142578125, 50975.8515625, 100736.3359375, 141623.796875, 250636.3125, 388673.8125, 474415.0625, 596670.75, 612769.8125, 639153.75, 666881.6875, 633754.375, 623379.4375, 497927.625, 463788.15625, 408077.8125, 354678.53125, 272698.625, 31860.34765625, 132546.234375, 241878.359375, 408140.75, 613362.1875, 897029.375, 1.139099625e6, 1.3976125e6, 1.634739125e6, 1.703326625e6, 1.705717625e6, 1.615322e6, 1.462416375e6, 1.2691205e6, 1.170907375e6, 1.059360375e6, 789649.1875, 630949.5625, 96186.9296875, 340087.25, 622924.0625, 960776.875, 1.567353875e6, 2.23006875e6, 2.86151425e6, 3.64290375e6, 4.115078e6, 4.3186265e6, 4.27756e6, 4.03195525e6, 3.88872925e6, 3.44737125e6, 2.9884135e6, 2.631653e6, 2.033596875e6, 1.665974375e6, 45646.1875, 104915.4140625, 209863.5625, 356309.4375, 611866.0, 820286.9375, 993218.125, 1.28606925e6, 1.433952875e6, 1.468773375e6, 1.488338125e6, 1.470022e6, 1.37718075e6, 1.199327125e6, 1.0408444375e6, 961842.5625, 730695.9375, 566071.6875]


X = [0.0; 0.0; 593.3040161132812; 680.9676513671875; 533.0647583007812; 742.5764770507812; 835.7925415039062; 1277.613525390625; 465.0248718261719; 977.941162109375; 453.80657958984375; 524.1534423828125; 400.8550109863281; 1025.7659912109375; 3729.7734375; 7977.93408203125; 33058.66796875; 58342.8359375; 96970.9765625; 125303.8515625; 105264.4453125; 68260.9375; 40450.44921875; 27465.583984375; 12540.0400390625; 10328.353515625;;]
y = [170845.453125, 373183.40625, 489773.0625, 640513.0, 896556.25, 894648.0625, 1.0691845e6, 1.056674e6, 1.2729035e6, 1.0171798125e6, 937198.375, 593592.5625, 694190.0625, 0.0, 19976.91796875, 0.0, 0.0, 32533.732421875, 0.0, 42338.94140625, 47968.13671875, 54009.5546875, 40316.9609375, 40895.29296875, 0.0, 33167.5078125]


X = [5647.8369140625; 6628.3076171875; 507.0470886230469; 8355.84765625; 13309.650390625; 14522.13671875; 13654.9013671875; 16488.599609375; 14314.76953125; 11612.3046875; 11920.208984375; 8753.75390625; 6444.62646484375; 5001.57177734375; 3521.316650390625; 2028.278564453125; 570.1524658203125; 1184.244873046875; 2045.870849609375; 2765.747802734375; 3292.171142578125; 4070.550048828125; 4552.0849609375; 4060.694091796875; 3315.963134765625; 2447.675048828125; 2127.295654296875; 1659.0201416015625; 1650.6153564453125; 1300.4232177734375; 747.2733764648438; 0.0; 1013.890380859375; 2284.882568359375; 3022.373291015625; 4095.805419921875; 6542.93017578125; 6977.5380859375; 7079.29296875; 7680.25732421875; 7975.08154296875; 6442.70263671875; 5676.60107421875; 3636.54833984375; 3423.8330078125; 2126.964599609375; 1613.52392578125; 774.7852783203125; 12436.3994140625; 23137.39453125; 36733.5390625; 56046.68359375; 72690.8671875; 90197.34375; 94232.328125; 91408.609375; 91542.34375; 77891.6484375; 65229.36328125; 48593.80859375; 34813.515625; 23575.701171875; 14340.1806640625; 9955.935546875;;]
y = [188529.84375, 283520.4375, 383212.6875, 485976.5, 573218.4375, 626788.75, 659732.0625, 594011.25, 536762.875, 536762.875, 466817.59375, 342485.84375, 257796.359375, 158859.578125, 76235.3359375, 76235.3359375, 54419.17578125, 91854.9140625, 133258.5, 160239.09375, 202187.578125, 201098.421875, 208030.0625, 203628.34375, 184173.78125, 184173.78125, 136529.640625, 118192.3671875, 88580.0625, 66881.6640625, 0.0, 0.0, 91172.4375, 178197.59375, 257583.625, 333822.0, 392629.25, 408956.53125, 432181.75, 430384.5, 321428.0, 321428.0, 295049.40625, 242776.71875, 197271.125, 120689.015625, 51242.890625, 51242.890625, 1.193580375e6, 2.05123275e6, 2.96587e6, 3.97473275e6, 4.633026e6, 5.0743135e6, 5.3432595e6, 5.04211e6, 4.2099095e6, 4.2099095e6, 3.54649325e6, 2.74803625e6, 2.06750025e6, 1.4481585e6, 549959.375, 549959.375]


I = (X[:,1].!=0.0) .& (y.!=0.0)
rlm(X[I,:]./mean(y),y[I]./mean(y), loss, initial_scale=:mad, maxiter = 200)
model = getModel(X, y, I)
plot(X[I,:], y[I], seriestype = :scatter)
plot!([0.0, maximum(X)], [0.0, coef(model)[1]*maximum(X)])

plot!([0.0, maximum(X)], [0.0, 0.6911979570792074*maximum(X)])


model = rlm(@view(X[(X[:,1].!=0.0) .& (y.!=0.0),:]), @view(y[(X[:,1].!=0.0) .& (y.!=0.0)]), MMEstimator{TukeyLoss}(), initial_scale=:mad)
model = rlm(X[X[:,1].!=0.0,:], y[X[:,1].!=0.0], MMEstimator{TukeyLoss}(), initial_scale=:mad)



plot(a[:,1], seriestype = :scatter)
plot!([0.0, maximum(a)], [0.0, coef(model)[1]*maximum(a)])


function foo(x)
    if x == 2
        @warn "you have been warned"
    end
    x
end
@test_log foo(2)

df = DataFrame(A=X[I,1], B=y[I])

using GLM
ols = lm(@formula(A  ~ B), df)


heavy = [1.0, 2.0, 3.0, 4.0]
light = [1.4, 2.6, 3.3, 3.8, 4.1]
getScanPairs(heavy, light)
end

light_scans = [1, 2, 3, 4, 5, 6]
heavy_scans = [1, 2, 3, 4, 5, 6]
light_transitions = UnorderedDictionary(
    ["y3+1","y4+1","y5+1"],
    [
        Float32[1.024025f6
        1.6092901f6
        3.0412648f6
        5.4651805f6
        8.631513f6
        1.1273232f7
        1.1273232f7
        8.631513f6
        5.4651805f6
        3.0412648f6
        1.6092901f6
        1.024025f6],
        Float32[ 436613.62
        704660.44
             1.3229029f6
             2.3309442f6
             3.7873248f6
             4.8260845f6
             4.8260845f6
             3.7873248f6
             2.3309442f6
             1.3229029f6
        704660.44
        436613.62],
        Float32[ 519827.88
        772209.56
             1.4236808f6
             2.6672412f6
             3.9857915f6
             5.183467f6
             5.183467f6
             3.9857915f6
             2.6672412f6
             1.4236808f6
        772209.56
        519827.88]
    ]
)
heavy_transitions = UnorderedDictionary{String, Vector{Float32}}()
for (key, value) in pairs(light_transitions)
    insert!(heavy_transitions, key, 100*value)
end
function plotPairedFragmentIonChromatogram(light_transitions::UnorderedDictionary{String, Vector{Float32}}, 
    heavy_transitions::UnorderedDictionary{String, Vector{Float32}}, light_rts::Vector{Float32}, heavy_rts::Vector{Float32}, prot_name::String,
    peptide_sequence::String, 
    par::Real, goodness_of_fit::Real, out_path::String)
    prot_name = replace(prot_name, r"/" => "|")
    p = plot(title = prot_name*"-"*peptide_sequence, fontfamily="helvetica", legend=:outertopright)
    if (length(heavy_rts) == 0) | (length(light_rts) == 0)
        return 
    end
    max_light = 0.0
    max_heavy = 0.0

    for (color, t) in enumerate(keys(heavy_transitions))
        if isassigned(light_transitions, t)
            if maximum(light_transitions[t])>max_light
                max_light = maximum(light_transitions[t])
            end
        end
        if maximum(heavy_transitions[t])>max_heavy
            max_heavy = maximum(heavy_transitions[t])
        end
    end

    for (color, t) in enumerate(keys(heavy_transitions))
        if isassigned(light_transitions, t)
            I = (light_transitions[t].>=1e-12)
            #I = (light_transitions[t].>=1e-12)
            plot!(p, light_rts[I], -1*(max_heavy/max_light)*light_transitions[t][I], color = color, legend = false, label = nothing)
            plot!(p, light_rts[I], -1*(max_heavy/max_light)*light_transitions[t][I], seriestype=:scatter, legend = false, color = color, label = nothing)
        end
        I = (heavy_transitions[t].>=1e-12)
        #I = (heavy_transitions[t].>=1e-12)
        plot!(p, heavy_rts[I], heavy_transitions[t][I], color = color, legend = :outertopright, label = t)
        plot!(p, heavy_rts[I], heavy_transitions[t][I], seriestype=:scatter, color = color, label = nothing)
    end

    ticks = [1*max_heavy, 0, -1*(max_heavy/max_light)*max_light]
    ticklabels = [ @sprintf "%.2E" x for x in [max_heavy, 0, max_light]]
    #println("par ", par)
    #par = string(par)
    par = @sprintf "%.2E" par
    goodness_of_fit= @sprintf "%.2E" 100.0*goodness_of_fit
    annotate!(minimum(heavy_rts), max_heavy, fontfamily = "helvetica", text("PAR $par \nstderr/mean % $goodness_of_fit", :black, :top, :left, 8))
    yticks!((ticks, 
             ticklabels), 
             axis = 1)
    #savefig(joinpath(out_path, prot_name*"-"*peptide_sequence*".pdf"))
    #plot!(legend=:outertopright)
end

light_transitions = UnorderedDictionary(
    ["y3+1","y4+1","y5+1"],
    [
        Float32[1.024025f6
        1.6092901f6
        3.0412648f6
        5.4651805f6
        8.631513f6
        1.1273232f7
        1.1273232f7
        8.631513f6
        5.4651805f6
        3.0412648f6
        1.6092901f6
        1.024025f6],
        Float32[ 436613.62
        704660.44
             1.3229029f6
             2.3309442f6
             3.7873248f6
             4.8260845f6
             4.8260845f6
             3.7873248f6
             2.3309442f6
             1.3229029f6
        704660.44
        436613.62],
        Float32[ 519827.88
        772209.56
             1.4236808f6
             2.6672412f6
             3.9857915f6
             5.183467f6
             5.183467f6
             3.9857915f6
             2.6672412f6
             1.4236808f6
        772209.56
        519827.88]
    ]
)
heavy_transitions = UnorderedDictionary{String, Vector{Float32}}()
heavy_transitions = UnorderedDictionary{String, Vector{Float32}}()
for (key, value) in pairs(light_transitions)
    insert!(heavy_transitions, key, 100*value)
end
times = [Float32(x) for x in 1:12]
plotPairedFragmentIonChromatogram(light_transitions, heavy_transitions, times, times, "Examaple","", 0.0, 0.0, "test")
savefig("/Users/n.t.wamsley/Desktop/chromatogram_no_interference.pdf")
par, goodness_of_fit = fitPAR(light_scans, heavy_scans, light_transitions, heavy_transitions)
X = Float32[]
for (key, value) in pairs(light_transitions)
    append!(X, value)
end
X = reshape(X, (length(X), 1))
y = Float32[]
for (key, value) in pairs(heavy_transitions)
    append!(y, value)
end
p = plot()
for (color, key) in enumerate(keys(light_transitions))
    plot!(p, light_transitions[key], heavy_transitions[key], seriestype=:scatter,color = color)
end
plot!(p, [0.0, maximum(X)], [0.0, (X\y)[1]*maximum(X)], label = "OLS")
plot!([0.0, maximum(X)], [0.0, par*maximum(X)], label = "Robust")
p
savefig("/Users/n.t.wamsley/Desktop/ModelFit_no_interference.pdf")

light_transitions["y4+1"] .*= 2
plotPairedFragmentIonChromatogram(light_transitions, heavy_transitions, times, times, "test","test", 0.0, 0.0, "test")
savefig("/Users/n.t.wamsley/Desktop/chromatogram_interference_a.pdf")
par, goodness_of_fit = fitPAR(light_scans, heavy_scans, light_transitions, heavy_transitions)
X = Float32[]
for (key, value) in pairs(light_transitions)
    append!(X, value)
end
X = reshape(X, (length(X), 1))
y = Float32[]
for (key, value) in pairs(heavy_transitions)
    append!(y, value)
end
p = plot()
for (color, key) in enumerate(keys(light_transitions))
    plot!(p, light_transitions[key], heavy_transitions[key], seriestype=:scatter,color = color)
end
plot!(p, [0.0, maximum(X)], [0.0, (X\y)[1]*maximum(X)], label = "OLS")
plot!([0.0, maximum(X)], [0.0, par*maximum(X)], label = "Robust")
p
savefig("/Users/n.t.wamsley/Desktop/ModelFit_interference_a.pdf")
light_transitions["y4+1"] .= 0.5e7
plotPairedFragmentIonChromatogram(light_transitions, heavy_transitions, times, times, "test","test", 0.0, 0.0, "test")
savefig("/Users/n.t.wamsley/Desktop/chromatogram_interference_b.pdf")
par, goodness_of_fit = fitPAR(light_scans, heavy_scans, light_transitions, heavy_transitions)
X = Float32[]
for (key, value) in pairs(light_transitions)
    append!(X, value)
end
X = reshape(X, (length(X), 1))
y = Float32[]
for (key, value) in pairs(heavy_transitions)
    append!(y, value)
end
p = plot()
for (color, key) in enumerate(keys(light_transitions))
    plot!(p, light_transitions[key], heavy_transitions[key], seriestype=:scatter,color = color)
end
plot!(p, [0.0, maximum(X)], [0.0, (X\y)[1]*maximum(X)], label = "OLS")
plot!([0.0, maximum(X)], [0.0, par*maximum(X)], label = "Robust")
p
savefig("/Users/n.t.wamsley/Desktop/ModelFit_interference_b.pdf")
=#