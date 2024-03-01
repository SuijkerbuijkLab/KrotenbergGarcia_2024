using Polynomials
using StatsKit
using CSV
using Statistics
using LinearAlgebra
using Plots
using DelimitedFiles
using DataFrames
using XLSX
path_plate0="C:\\Users\\Ledes002\\Desktop\\secondSet\\3_Differentiations_d1VSd3_Exp142_149_150_155_161\\Overview positions_Diff_d1_d3.xlsx"
data_exp142_d3=DataFrame(XLSX.readtable(path_plate0, "Exp142_d3"));
replace!(data_exp142_d3.Scene, missing => -1);
select!(data_exp142_d3, Not("Notes for cropping"))
select!(data_exp142_d3, Not("Discard"))
data_exp149_d1=DataFrame(XLSX.readtable(path_plate0, "Exp149_d1"));
replace!(data_exp149_d1.Scene, missing => -1);
select!(data_exp149_d1, Not("Notes for cropping"))
select!(data_exp149_d1, Not("Discard"))
data_exp149_d3=DataFrame(XLSX.readtable(path_plate0, "Exp149_d3"));
replace!(data_exp149_d3.Scene, missing => -1);
select!(data_exp149_d3, Not("Notes for cropping"))
select!(data_exp149_d3, Not("Discard"))
data_exp150_d1=DataFrame(XLSX.readtable(path_plate0, "Exp150_d1"));
replace!(data_exp150_d1.Scene, missing => -1);
select!(data_exp150_d1, Not("Notes for cropping"))
select!(data_exp150_d1, Not("Discard"))
data_exp150_d3=DataFrame(XLSX.readtable(path_plate0, "Exp150_d3"));
replace!(data_exp150_d3.Scene, missing => -1);
select!(data_exp150_d3, Not("Notes for cropping"))
select!(data_exp150_d3, Not("Discard"))
data_exp155_d1=DataFrame(XLSX.readtable(path_plate0, "Exp155_d1"));
replace!(data_exp155_d1.Scene, missing => -1);
select!(data_exp155_d1, Not("Notes for cropping"))
select!(data_exp155_d1, Not("Discard"))
data_exp155_d3=DataFrame(XLSX.readtable(path_plate0, "Exp155_d3_with green channel"));
replace!(data_exp155_d3.Scene, missing => -1);
select!(data_exp155_d3, Not("Notes for cropping"))
select!(data_exp155_d3, Not("Discard"))
data_exp161_d1=DataFrame(XLSX.readtable(path_plate0, "Exp161_d1"));
replace!(data_exp161_d1.Scene, missing => -1);
select!(data_exp161_d1, Not("Notes for cropping"))
select!(data_exp161_d1, Not("Discard"))	
data_exp161_d3=DataFrame(XLSX.readtable(path_plate0, "Exp161_d3"));
replace!(data_exp161_d3.Scene, missing => -1);
select!(data_exp161_d3, Not("Notes for cropping"))
select!(data_exp161_d3, Not("Discard"))	
#index of summaries for each image 
indexX=1; indexY=2; indexZ=3; indexC=4; indexA=5; indexB=6; indexAng=7; 
indexVx=8; indexVx=9; indexVz=10; indexVol=11; indexvolEll=12; indexRaw=13; 
indexAverageOscar=14; indexXB=15; indexYB=16; indexWB=17; indexHB=18; 
indexXF=19; indexYF=20; indexFMax=21; indexFmin=22; indexFang=23; indexXM=24; 
indexYM=25; indexSkew=26; indexKurt=27; indexMean=28; indexSTD=29; indexMedian=30; 
indexCh1_vol=31; indexCh1_Int=32; indexCh2_vol=33; indexCh2_Int=34; indexCh3_vol=35; indexCh3_Int=36;
path="C:\\Users\\Ledes002\\Desktop\\secondSet\\3_Differentiations_d1VSd3_Exp142_149_150_155_161\\proc"
pathSummariesImage="C:\\Users\\Ledes002\\Desktop\\secondSet\\3_Differentiations_d1VSd3_Exp142_149_150_155_161\\proc\\OSCAR_output_normal"
# pathSummariesImage=path*"\\OSCAR_output_borders"; 
pathSummaryWhole=path*"\\InternMeasurements_TLM_final.txt";
pathO=path*"\\Summaries_Files"; 
varSearch="Summary_Exp";
dataFolder = CSV.read(pathSummaryWhole, DataFrame); # output path.
#
namesV=Vector{}(); sceneV=Vector{}(); orgV=Vector{}(); expV=Vector{}(); 
totalNumb=Vector{}(); WT=Vector{}();  percWT=Vector{}();  
groupV=Vector{}(); dayV=Vector{}(); conditionV=Vector{}();
#
alpha=0.25; alpha2=0.25;
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
    nameFile=file; push!(namesV, nameFile);
    indExp=collect(findall.("Exp", file)[1]);    
    indTime=collect(findall.("_d", file)[1]);    
    indT2=collect(findall.("_-Scene-", file)[1]); 
    indT3=collect(findall.("_-Org_-", file)[1]);
    scene=file[indT2[end]+1:end]; scene=parse(Int64,scene); push!(sceneV,scene)
    org=file[indT3[end]+1:indT2[1]-1]; push!(orgV,org)
    exp=file[indExp[end]+1:indTime[1]-1]; push!(expV,exp)
    time=file[indTime[end]+1:indT3[1]-1]; push!(dayV, time)

    println("Processing the file ", i);  file=file*".txt"; file="Summary_"*file; 
    #
    if "Exp"*exp*"_d"*time=="Exp142_d3"
        indT=findall(x -> x==scene, data_exp142_d3[:,:Scene])[1]
        grp_sample=data_exp142_d3[indT,:Group]
        condition_sample=data_exp142_d3[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp149_d1"
        indT=findall(x -> x==scene, data_exp149_d1[:,:Scene])[1]
        grp_sample=data_exp149_d1[indT,:Group]
        condition_sample=data_exp149_d1[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp149_d3"
        indT=findall(x -> x==scene, data_exp149_d3[:,:Scene])[1]
        grp_sample=data_exp149_d3[indT,:Group]
        condition_sample=data_exp149_d3[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp150_d1"
        indT=findall(x -> x==scene, data_exp150_d1[:,:Scene])[1]
        grp_sample=data_exp150_d1[indT,:Group]
        condition_sample=data_exp150_d1[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp150_d3"
        indT=findall(x -> x==scene, data_exp150_d3[:,:Scene])[1]
        grp_sample=data_exp150_d3[indT,:Group]
        condition_sample=data_exp150_d3[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp155_d1"
        indT=findall(x -> x==scene, data_exp155_d1[:,:Scene])[1]
        grp_sample=data_exp155_d1[indT,:Group]
        condition_sample=data_exp155_d1[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp155_d3"
        indT=findall(x -> x==scene, data_exp155_d3[:,:Scene])[1]
        grp_sample=data_exp155_d3[indT,:Group]
        condition_sample=data_exp155_d3[indT,:Condition]
    elseif "Exp"*exp*"_d"*time=="Exp161_d1"
        indT=findall(x -> x==scene, data_exp161_d1[:,:Scene])[1]
        grp_sample=data_exp161_d1[indT,:Group]
        condition_sample=data_exp161_d1[indT,:Differentiation]
    elseif "Exp"*exp*"_d"*time=="Exp161_d3"
        indT=findall(x -> x==scene, data_exp161_d3[:,:Scene])[1]
        grp_sample=data_exp161_d3[indT,:Group]
        condition_sample=data_exp161_d3[indT,:Differentiation]
    end
    push!(groupV, lowercase(grp_sample))
    push!(conditionV,lowercase(condition_sample)) 
    # opening the summary of each image to get numb of objc for eah Channel
    dataH=CSV.read(pathSummariesImage*"\\$file", DataFrame);  
    numbWT=0;
    intenV_WT=Vector{}(); intenV_Cancer=Vector{}(); 
    for j in 1:size(dataH, 1)
        VolTot=dataH[j, indexVol]; #dataH[j, indexCh1]dataH[j, indexCh2] + dataH[j, indexCh3];
        VolCh1=dataH[j, indexCh1_vol];
        VolCh2=dataH[j, indexCh2_vol];
        avgChDAPI=dataH[j, indexAverageOscar];
        #
        #index of summaries for each image 
        if VolCh2>alpha*VolTot && VolCh2>alpha2*VolCh1
            #
            numbWT=numbWT+1;
        end
        #
    end
    push!(totalNumb, size(dataH, 1));
    push!(WT, numbWT);
    push!(percWT, numbWT/size(dataH, 1))
end
col_labels = ["name", "Exp", "Day", "scene", "org", "ExpGroup", "condition", "numbTotal", 
"wt","percWT"];
col_labels = Symbol.(col_labels);
dataV=hcat(namesV,expV,dayV,sceneV,orgV, groupV,conditionV, totalNumb, WT, percWT)
dfmain=DataFrame(Tables.table(dataV, header = col_labels))
XLSX.openxlsx(pathO*"\\Data_normal.xlsx", mode="w") do xf
    nodenames = string.(names(dfmain))
    XLSX.addsheet!(xf, "Data")
    sheet=xf["Data"]
    for r in 1:length(nodenames)
        sheet[XLSX.CellRef(1, r)] = nodenames[r]
        for ii in 1:size(dfmain,1)
            coordX=ii+1;
            sheet[XLSX.CellRef(coordX, r)]=dfmain[ii,r]; 
        end
    end
        #
end
# CSV.write(pathO*"\\Images_GrouppedInfo.txt", dfmain, delim = '\t')

XLSX.openxlsx(pathO*"\\Data_GrouppedByExpCondition_normal.xlsx", mode="w") do xf
    grpdInfo=groupby(dfmain, :ExpGroup)
    namesK=keys(grpdInfo);
    for i in 1:size(namesK,1)
        #
        k=namesK[i]; TAG=k.ExpGroup
        #
        dataT=DataFrame(grpdInfo[i])
        nodenames = string.(names(dataT))
        XLSX.addsheet!(xf, TAG)
        sheet=xf[TAG]
        for r in 1:length(nodenames)
            sheet[XLSX.CellRef(1, r)] = nodenames[r]
            for ii in 1:size(dataT,1)
                coordX=ii+1;
                sheet[XLSX.CellRef(coordX, r)]=dataT[ii,r]; 
            end
        end
        #
    end
end