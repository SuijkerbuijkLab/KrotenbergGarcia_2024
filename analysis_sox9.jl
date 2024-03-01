using Polynomials
using StatsKit
using CSV
using Statistics
using LinearAlgebra
using Plots
using DelimitedFiles
using DataFrames
using XLSX
using StatsPlots
path_plate0="C:/Users/Ledes002/Desktop/secondSet/2_SOX9_JV_Exp025_AKGExp141_142_156/Overview positions_SOX9.xlsx"
data_exp025=DataFrame(XLSX.readtable(path_plate0, "JV_Exp025"));
select!(data_exp025, Not("Useful?-Crop info"))
data_exp141=DataFrame(XLSX.readtable(path_plate0, "AKG_Exp141"));
select!(data_exp141, Not("Staining"))
select!(data_exp141, Not("Notes for cropping"))
select!(data_exp141, Not("Discard"))
data_exp156=DataFrame(XLSX.readtable(path_plate0, "AKG_Exp156"));
select!(data_exp156, Not("Staining"))
select!(data_exp156, Not("Notes for cropping"))
select!(data_exp156, Not("Discard"))
data_exp142=DataFrame(XLSX.readtable(path_plate0, "AKG_Exp142_diff"));
select!(data_exp142, Not("Staining"))
select!(data_exp142, Not("Notes for cropping"))
select!(data_exp142, Not("Discard"))
#index of summaries for each image 
indexX=1; indexY=2; indexZ=3; indexC=4; indexA=5; indexB=6; indexAng=7; 
indexVx=8; indexVx=9; indexVz=10; indexVol=11; indexvolEll=12; indexRaw=13; 
indexAverageOscar=14; indexXB=15; indexYB=16; indexWB=17; indexHB=18; 
indexXF=19; indexYF=20; indexFMax=21; indexFmin=22; indexFang=23; indexXM=24; 
indexYM=25; indexSkew=26; indexKurt=27; indexMean=28; indexSTD=29; indexMedian=30; 
indexCh1_vol=31; indexCh1_Int=32; indexCh2_vol=33; indexCh2_Int=34; indexCh3_vol=35; indexCh3_Int=36;
path="C:\\Users\\Ledes002\\Desktop\\secondSet\\2_SOX9_JV_Exp025_AKGExp141_142_156\\proc"
pathSummariesImage="C:\\Users\\Ledes002\\Desktop\\secondSet\\2_SOX9_JV_Exp025_AKGExp141_142_156\\proc\\OSCAR_output_borders\\OSCAR_summaries"
# pathSummariesImage=path*"\\OSCAR_output_borders"; 
pathSummaryWhole=path*"\\InternMeasurements_TLM.txt";
pathO=pathSummariesImage*"\\nucleiOut-test"; 

varSearch="Summary_Exp";
dataFolder = CSV.read(pathSummaryWhole, DataFrame); # output path.
#
nodenames = ["name","Exp", "scene", "org", "ExpGroup", "GroupID", "ID", "Phenotype", "Phenotype2", "slice", "xc", "yc", 
"Intensity_Marker", "Intensity_DAPI", "sox9ONdapi", "normInt", "volume"];
df_nuclei=DataFrame([[] for _ = nodenames] , nodenames)
alpha=0.25; alpha2=0.25;
newMax=65535; newMin=0;
for i in 1:size(dataFolder,1)
    # name, minR,2minZ and maxZ (aqui medZ de IJ es minZ)
    minR=dataFolder[i,2]
    medZProc=dataFolder[i,3]
    minZProc=dataFolder[i,4]
    minZ=dataFolder[i,5]
    medZ=dataFolder[i,6]
    maxZ=dataFolder[i,7]
    TotalNuclei=dataFolder[i,8]
    file=dataFolder[i,1]; #file=file[1:end-4]; #cambiar.tif por .txt (por el pipe IJ > julia)    
    #
    nameFile=file; 
    println("Processing the file ", nameFile);
    indExp=collect(findall.("Exp", file)[1]);      
    indT2=collect(findall.("_-Scene-", file)[1]); 
    indT3=collect(findall.("_-Org_-", file)[1]);
    scene=file[indT2[end]+1:end]; scene=parse(Int64,scene);
    org=file[indT3[end]+1:indT2[1]-1]; 
    exp=file[indExp[end]+1:indT3[1]-1]; 
    println("Processing the file ", file);  file=file*".txt"; file="Summary_"*file; 
    #
    if exp=="025"
        indT=findall(x -> x==scene, data_exp025[:,:Scene])[1]
        grp_sample=data_exp025[indT,:Group]
    elseif exp=="141"
        indT=findall(x -> x==scene, data_exp141[:,:Scene])[1]
        grp_sample=data_exp141[indT,:Group]
    elseif exp=="156"
        indT=findall(x -> x==scene, data_exp156[:,:Scene])[1]
        grp_sample=data_exp156[indT,:Group]
    elseif exp=="142"
        indT=findall(x -> x==scene, data_exp142[:,:Scene])[1]
        grp_sample=data_exp142[indT,:Group]
    else
        println("at index: ", i, "is ", exp)
    end
    if lowercase(grp_sample)=="crc" || lowercase(grp_sample)=="cancer"
        grp_sample="C"
    elseif lowercase(grp_sample)=="mix"
        grp_sample="M"
    elseif lowercase(grp_sample)=="wt"
        grp_sample="W"
    end
    # opening the summary of each image to get numb of objc for eah Channel
    dataH=CSV.read(pathSummariesImage*"\\$file", DataFrame);  
    numbWT=0;
    intenV_WT=Vector{}(); intenV_Cancer=Vector{}(); 
    #
    sceneVt=Vector{}(); orgVt=Vector{}(); groupVt=Vector{}(); labelPropioV=Vector{}();
    IDsV=Vector{}(); filesV=Vector{}(); expVt=Vector{}();
    phenV=Vector{}(); sliceV=Vector{}(); phenV2=Vector{}();
    xcV=Vector{}(); ycV=Vector{}(); volumeV=Vector{}();
    intenV=Vector{}(); intenDapiV=Vector{}();
    #
    volsV=dataH[:, indexVol]; avgIntenMarker=dataH[:, indexCh3_Int]./dataH[:, indexVol]
    for j in 1:size(dataH, 1)
        VolTot=dataH[j, indexVol]; #dataH[j, indexCh1]dataH[j, indexCh2] + dataH[j, indexCh3];
        VolCh1=dataH[j, indexCh1_vol];
        VolCh2=dataH[j, indexCh2_vol];
        avgCh3=avgIntenMarker[j]
        # avgCh3=dataH[j, indexCh3_Int]/dataH[j, indexVol];
        avgChDAPI=dataH[j, indexAverageOscar];
        #
        push!(groupVt, grp_sample)
        push!(sceneVt,scene)
        push!(orgVt,org)
        push!(filesV, nameFile)
        push!(IDsV, ("nuc-"*string(j)))
        push!(expVt, exp)
        push!(sliceV, ceil(dataH[j, indexZ]))
        push!(xcV, ceil(dataH[j, indexX]))
        push!(ycV, ceil(dataH[j, indexY]))
        push!(volumeV, VolTot)
        #index of summaries for each image 
        if isnan(dataH[j, indexCh3_Int])
            push!(intenV, 0)
        else
            push!(intenV, avgCh3)
        end 
        if isnan(avgChDAPI)
            push!(intenDapiV, 0)
        else
            push!(intenDapiV, avgChDAPI)
        end  
        phT2=nothing; phT=nothing
        if (grp_sample)=="W"
            numbWT=numbWT+1
            # phT="WT_pure"
            phT="WT"
            phT2="PureWT"
        elseif (grp_sample)=="C" 
            # phT="Cancer_pure"
            phT="Cancer"
            phT2="PureCancer" 
        else (grp_sample)=="M"            
            if VolCh2>alpha*VolTot && VolCh2>alpha2*VolCh1
                numbWT=numbWT+1;  
                phT="WT" 
            else
                phT="Cancer"
            end
            phT2=phT
        end
        push!(phenV, phT)  
        push!(phenV2, phT2)     
        #
        push!(labelPropioV, "E-"*exp*"_Gr-"*phT2)
    end
    #
    col_labels = nodenames;
    rel=intenV./intenDapiV; 
    oldMax=maximum(intenV)
    oldMin=minimum(intenV)
    varT2=newMax.*(intenV.-oldMin)./(oldMax-oldMin)
    println(all(phenV.==phenV2))
    dataV=hcat(filesV, expVt, sceneVt, orgVt, groupVt,labelPropioV, IDsV, phenV,phenV2,sliceV,xcV, ycV, intenV, intenDapiV, rel, varT2,volumeV)
    dfmain=DataFrame(Tables.table(dataV, header = col_labels))
    dfmain=filter(:volume => n -> n > 200, dfmain); 
    dfmain=filter(:Intensity_Marker => n -> n > 0, dfmain); 
    df_nuclei=vcat(df_nuclei, dfmain)   
    # CSV.write(pathO*"\\"*nameFile*".csv", dfmain, delim = '\t')
    #
end
#
grpdInfo=groupby(df_nuclei, :Exp); namesK=keys(grpdInfo);
normInt2=Vector{}(); dfh=DataFrame(); labelV=Vector{}();
namesTV=Vector{}();
limits_Group1_wt=Vector{}(); limits_Group2_wt=Vector{}(); limits_Group3_wt=Vector{}();
limits_Group1_c=Vector{}(); limits_Group2_c=Vector{}(); limits_Group3_c=Vector{}();
for i in 1:size(keys(grpdInfo),1)
    dataT=DataFrame(grpdInfo[i])
    varT=dataT[:, :Intensity_Marker]
    maxL=varT./maximum(varT);
    k=namesK[i]; TAG="Exp-"*k.Exp #*"Ph-"*k.Phenotype
    stringsL = [TAG for i in 1:size(dataT,1)]
    labelV=vcat(labelV, stringsL)
    normInt2=vcat(normInt2, maxL)
    dfh=vcat(dfh,dataT); 
    #
    push!(namesTV, k.Exp)
    dataT_wt=filter(:Phenotype => n -> n == "WT", dataT);
    limswt=percentile(dataT_wt[:,:Intensity_Marker], [25 50 75])
    push!(limits_Group1_wt, limswt[1])
    push!(limits_Group2_wt, limswt[2])
    push!(limits_Group3_wt, limswt[3])
    dataT_c=filter(:Phenotype => n -> n == "Cancer", dataT); 
    limsc=percentile(dataT_c[:,:Intensity_Marker], [25 50 75])
    push!(limits_Group1_c, limsc[1])
    push!(limits_Group2_c, limsc[2])
    push!(limits_Group3_c, limsc[3])
    #
end
dfh[!, :relInt]=normInt2;
# dfh[!, :IDs]=labelV;
#
dfh_g=groupby(dfh, :Phenotype)
df_temp=DataFrame(dfh_g[("WT",)])
limsWT=percentile(df_temp[:,:relInt], [25 50 75])
df_temp=DataFrame(dfh_g[("Cancer",)])
limsCancer=percentile(df_temp[:,:relInt], [25 50 75])
#
dfh_g=groupby(dfh, :name)
namesK=keys(dfh_g);
namesV=Vector{}(); sceneV=Vector{}(); orgV=Vector{}(); expV=Vector{}(); 
totalNumb=Vector{}(); percWT=Vector{}();  groupV=Vector{}(); groupV2=Vector{}();  
IDpropioV=Vector{}();
avgMarker_wt=Vector{}(); avgMarker_cancer=Vector{}();
avgMarker_wt2=Vector{}(); avgMarker_cancer2=Vector{}();
fraction1_wt=Vector{}(); fraction2_wt=Vector{}(); fraction3_wt=Vector{}();
fraction1_c=Vector{}(); fraction2_c=Vector{}(); fraction3_c=Vector{}();
fraction1_wt_2=Vector{}(); fraction2_wt_2=Vector{}(); fraction3_wt_2=Vector{}();
fraction1_c_2=Vector{}(); fraction2_c_2=Vector{}(); fraction3_c_2=Vector{}();
for i in 1:size(namesK,1)
    k=namesK[i]; TAG=k.name #*"Ph-"*k.Phenotype
    #calcular relInt para WT y para Cancer
    dataT=DataFrame(dfh_g[i])
    push!(namesV, TAG);
    push!(sceneV,dataT[1,:scene])
    push!(orgV,dataT[1,:org])
    push!(expV,dataT[1,:Exp])
    push!(groupV, dataT[1,:Phenotype])
    push!(groupV2, dataT[1,:Phenotype2])

    push!(IDpropioV,"Exp_"*dataT[1,:Exp]*"Gr_"*dataT[1,:Phenotype2]);

    push!(totalNumb,size(dataT,1))
    dataT_wt=filter(:Phenotype => n -> n == "WT", dataT); 
    push!(percWT, size(dataT_wt,1)/size(dataT,1));
    vtWT=0; vtWT2=0; p1_wt=0; p2_wt=0; p3_wt=0; 
    p1_wt_2=0; p2_wt_2=0; p3_wt_2=0; 
    dataT_cancer=filter(:Phenotype => n -> n == "Cancer", dataT); 
    vtC=0; vtC2=0; p1_c=0; p2_c=0; p3_c=0;    
    p1_c_2=0; p2_c_2=0; p3_c_2=0;   
    indexSearch_limitsExp=findall(x -> x==dataT[1,:Exp], namesTV)[1]
    if size(dataT_wt,1)>0
        vtWT=mean(dataT_wt[:,:Intensity_Marker])
        tempTT2=dataT_wt[:,:Intensity_Marker]; vtWT=mean(tempTT2)
        tempTT=dataT_wt[:,:relInt]
        vtWT2=mean(tempTT)
        p1_wt=length(filter((x) -> 0<=x<limsWT[1], tempTT))/length(tempTT)
        p2_wt=length(filter((x) -> limsWT[1]<=x<=limsWT[3], tempTT))/length(tempTT)
        p3_wt=length(filter((x) -> limsWT[3]<x<=maximum(tempTT), tempTT))/length(tempTT)

        p1_wt_2=length(filter((x) -> 0<=x<limits_Group1_wt[indexSearch_limitsExp], tempTT2))/length(tempTT2)
        p2_wt_2=length(filter((x) -> limits_Group1_wt[indexSearch_limitsExp]<=x<=limits_Group3_wt[indexSearch_limitsExp], tempTT2))/length(tempTT2)
        p3_wt_2=length(filter((x) -> limits_Group3_wt[indexSearch_limitsExp]<x<=maximum(dataT_wt[:,:Intensity_Marker]), tempTT2))/length(tempTT2)
    end
    if size(dataT_cancer,1)>0
        tempTT2=dataT_cancer[:,:Intensity_Marker]; vtC=mean(tempTT2)
        tempTT=dataT_cancer[:,:relInt]
        vtC2=mean(tempTT)
        p1_c=length(filter((x) -> 0<=x<limsCancer[1], tempTT))/length(tempTT)
        p2_c=length(filter((x) -> limsCancer[1]<=x<=limsCancer[3], tempTT))/length(tempTT)
        p3_c=length(filter((x) -> limsCancer[3]<x<=maximum(tempTT), tempTT))/length(tempTT)

        p1_c_2=length(filter((x) -> 0<=x<limits_Group1_c[indexSearch_limitsExp], tempTT))/length(tempTT2)
        p2_c_2=length(filter((x) -> limits_Group1_c[indexSearch_limitsExp]<=x<=limits_Group3_c[indexSearch_limitsExp], tempTT))/length(tempTT2)
        p3_c_2=length(filter((x) -> limits_Group3_c[indexSearch_limitsExp]<x<=maximum(dataT_cancer[:,:Intensity_Marker]), tempTT))/length(tempTT2)
    end
    push!(avgMarker_wt, vtWT); push!(avgMarker_wt2, vtWT2)
    push!(fraction1_wt, p1_wt); push!(fraction2_wt, p2_wt); push!(fraction3_wt, p3_wt)
    push!(fraction1_wt_2, p1_wt_2); push!(fraction2_wt_2, p2_wt_2); push!(fraction3_wt_2, p3_wt_2)
    push!(avgMarker_cancer, vtC); push!(avgMarker_cancer2, vtC2)
    push!(fraction1_c, p1_c); push!(fraction2_c, p2_c); push!(fraction3_c, p3_c)
    push!(fraction1_c_2, p1_c_2); push!(fraction2_c_2, p2_c_2); push!(fraction3_c_2, p3_c_2)
    #
end
col_labels = ["name", "Exp", "scene", "org", "ExpGroup","ExpGroup2", "IDpropio", "numbTotal", 
"wt/total","avgMarker_wt", "avgMarker_cancer", "avgRelInt_wt", "avgRelInt_cancer", 
"WT_p25_RelInt", "WT_p25p75_RelInt", "WT_p75_RelInt","C_p25_RelInt", "C_p25p75_RelInt", "C_p75_RelInt",
"WT_p25_byExp", "WT_p25p75_byExp", "WT_p75_byExp","C_p25_byExp", "C_p25p75_byExp", "C_p75_byExp"];
col_labels = Symbol.(col_labels);
dataV=hcat(namesV,expV,sceneV,orgV,groupV,groupV2, IDpropioV, 
totalNumb,percWT, avgMarker_wt, avgMarker_cancer, avgMarker_wt2, avgMarker_cancer2, 
fraction1_wt, fraction2_wt, fraction3_wt, fraction1_c, fraction2_c, fraction3_c,
fraction1_wt_2, fraction2_wt_2, fraction3_wt_2, fraction1_c_2, fraction2_c_2, fraction3_c_2)
dfmain=DataFrame(Tables.table(dataV, header = col_labels))


# real intensities of organoids (average) for wt and cancer populations for each Experiment
# This shows the variability in experimental conditions. 
gr(size=(2048,1024))
@df dfmain violin(string.(:IDpropio), :avgMarker_wt, linewidth=0,  legend=false)
@df dfmain dotplot!(string.(:IDpropio), :avgMarker_wt, 
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot)
@df dfmain boxplot!(string.(:IDpropio), :avgMarker_wt, fillalpha=0.75, linewidth=2, legend=false)
# plot!(xtickfont = font(12, "Comic Sans"))
# plot!(ytickfont = font(12, "Comic Sans"))

plot!(formatter=:plain)
gr(size=(2048,1024))
@df dfmain violin(string.(:IDpropio), :avgMarker_cancer, linewidth=0,  legend=false)
@df dfmain dotplot!(string.(:IDpropio), :avgMarker_cancer, 
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot)
@df dfmain boxplot!(string.(:IDpropio), :avgMarker_cancer, fillalpha=0.75, linewidth=2, legend=false)
# plot!(xtickfont = font(12, "Comic Sans"))
# plot!(ytickfont = font(12, "Comic Sans"))
plot!(formatter=:plain)

#this plot is about nuclei. To get the relative intensity, I will show two plots, 
#one with the intensitites for each experiment of all nuclei and another one with relint

# nodenames = ["name","Exp", "scene", "org", "ExpGroup", "ID", "Phenotype", "Phenotype2", "slice", "xc", "yc", 
# "Intensity_Marker", "Intensity_DAPI", "sox9ONdapi", "normInt", "volume"];

gr(size=(5000,3000))
@df dfh violin(string.(:GroupID), :Intensity_Marker, linewidth=0,  legend=false)
@df dfh dotplot!(string.(:GroupID), :Intensity_Marker, 
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot)
@df dfh boxplot!(string.(:GroupID), :Intensity_Marker, fillalpha=0.75, linewidth=2, legend=false)
# plot!(xtickfont = font(16, "Comic Sans"))
# plot!(ytickfont = font(16, "Comic Sans"))
plot!(formatter=:plain)

gr(size=(2048,1024))
# @df dfh violin(string.(:Phenotype2), :relInt, linewidth=0,  legend=false)
@df dfh dotplot!(string.(:Phenotype2), :relInt, 
markershape = :hexagon,
markersize = 5,
markeralpha = 0.25,
markercolor = :orange,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot)
@df dfh boxplot!(string.(:Phenotype2), :relInt, fillalpha=0.75, linewidth=2, legend=false)
# plot!(xtickfont = font(16, "Comic Sans"))
# plot!(ytickfont = font(16, "Comic Sans"))
plot!(formatter=:plain)


gr(leg = false, bg = :lightgrey)
fig = plot(layout=(2,2), legend= false)
# @df dfmain violin!(string.(:ExpGroup), :avgRelInt_wt, linewidth=0,  legend=false)
@df dfmain dotplot!(string.(:ExpGroup), :avgRelInt_wt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=1)
@df dfmain boxplot!(string.(:ExpGroup), :avgRelInt_wt, fillalpha=0.75, linewidth=2, legend=false,
subplot=1)

# @df dfmain violin(string.(:ExpGroup), :WT_p25_RelInt, linewidth=0,  legend=false)
@df dfmain dotplot!(string.(:ExpGroup), :WT_p25_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=2)
@df dfmain boxplot!(string.(:ExpGroup), :WT_p25_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=2)


@df dfmain dotplot!(string.(:ExpGroup), :WT_p25p75_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=3)
@df dfmain boxplot!(string.(:ExpGroup), :WT_p25p75_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=3)


@df dfmain dotplot!(string.(:ExpGroup), :WT_p75_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=4)
@df dfmain boxplot!(string.(:ExpGroup), :WT_p75_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=4)
plot!(formatter=:plain)


gr(leg = false, bg = :lightgrey)
fig = plot(layout=(2,2), legend= false)
@df dfmain dotplot!(string.(:ExpGroup), :avgRelInt_cancer,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=1)
@df dfmain boxplot!(string.(:ExpGroup), :avgRelInt_cancer, fillalpha=0.75, linewidth=2, legend=false,
subplot=1)

@df dfmain dotplot!(string.(:ExpGroup), :C_p25_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=2)
@df dfmain boxplot!(string.(:ExpGroup), :C_p25_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=2)


@df dfmain dotplot!(string.(:ExpGroup), :C_p25p75_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=3)
@df dfmain boxplot!(string.(:ExpGroup), :C_p25p75_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=3)


@df dfmain dotplot!(string.(:ExpGroup), :C_p75_RelInt,
markershape = :hexagon,
markersize = 5,
markeralpha = 0.75,
markercolor = :black,
markerstrokewidth = 3,
markerstrokealpha = 0.2,
markerstrokecolor = :black,
markerstrokestyle = :dot,
subplot=4)
@df dfmain boxplot!(string.(:ExpGroup), :C_p75_RelInt, fillalpha=0.75, linewidth=2, legend=false,
subplot=4)
plot!(formatter=:plain)