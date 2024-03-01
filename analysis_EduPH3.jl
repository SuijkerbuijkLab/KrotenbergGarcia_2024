using Polynomials
using StatsKit
using CSV
using Statistics
using LinearAlgebra
using Plots
using DelimitedFiles
using DataFrames
using XLSX
#index of summaries for each image. Por el processing, channel 1 es cancer y channel 2 es WT
indexX=1; indexY=2; indexZ=3; indexC=4; indexA=5; indexB=6; indexAng=7; 
indexVx=8; indexVx=9; indexVz=10; indexVol=11; indexvolEll=12; indexRaw=13; 
indexAverageOscar=14; indexXB=15; indexYB=16; indexWB=17; indexHB=18; 
indexXF=19; indexYF=20; indexFMax=21; indexFmin=22; indexFang=23; indexXM=24; 
indexYM=25; indexSkew=26; indexKurt=27; indexMean=28; indexSTD=29; indexMedian=30; 
index_volPH3=31; index_IntPH3=32; index_volPH3_high=33; index_IntPH3_high=34; 
index_volEDU=35; index_IntEDU=36; index_volCherry=37; index_IntCherry=38;
#
path_plate0="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\Overview positions Exp152_Mario.xlsx"
data_plate0=DataFrame(XLSX.readtable(path_plate0, "Sheet1"));
data_plate0=rename(data_plate0, "scene nr"=>:scene)
select!(data_plate0, Not("notes for cropping"))
path_plate1="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\Overview positions_Exp153_Plate1.xlsx"
data_plate1=DataFrame(XLSX.readtable(path_plate1, "Sheet1"));
data_plate1=rename(data_plate1, "scene nr"=>:scene)
select!(data_plate1, Not("notes for cropping"))
path_plate2="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\Overview positions_Exp153_Plate2.xlsx"
data_plate2=DataFrame(XLSX.readtable(path_plate2, "Sheet1"));
data_plate2=rename(data_plate2, "scene nr"=>:scene)
select!(data_plate2, Not("notes for cropping"))
#
path="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\proc\\OSCAR_output_channels\\OSCAR_summaries_2"
pathO="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\proc\\OSCAR_output_channels\\OSCAR_summaries_2\\INFO"
if !ispath(pathO)
    mkdir(pathO);
end
pathXZ="C:\\Users\\Ledes002\\Desktop\\secondSet\\EdU_pulse_Exp152_Exp153\\proc\\out\\InternMeasurements_TLM_second.txt"
dataFolder = CSV.read(pathXZ, DataFrame);
dataFolder=sort(dataFolder, :name)
#
alpha_mcherry=0.1; 
alpha_edu=0.1; 
alpha_ph3=0.1;
namesV=Vector{}(); 
sceneV=Vector{}();
orgV=Vector{}(); 
expV=Vector{}();
groupV=Vector{}();  
timeEDU=Vector{}();
totalNumb=Vector{}(); 
wtV=Vector{}();  
cancerV=Vector{}(); 
ph3_totalV=Vector{}(); 
edu_totalV=Vector{}(); 
ph3_wtV=Vector{}(); 
ph3_cancerV=Vector{}();
edu_wtV=Vector{}(); 
edu_cancerV=Vector{}();        
ph3_eduV=Vector{}(); 
ph3_edu_wtV=Vector{}(); 
ph3_edu_cancerV=Vector{}(); 
for main_index in 1:size(dataFolder,1)
    file=dataFolder[main_index,:name]; nameFile="Summary_"*file*".txt"
    push!(namesV, file); 
    println("Processing the file2 ", nameFile);
    indExp=collect(findall.("Exp", file)[1]);      
    indT2=collect(findall.("_-Scene-", file)[1]); 
    indT3=collect(findall.("_-Org_-", file)[1]);
    scene=file[indT2[end]+1:end]; scene=parse(Int64,scene); push!(sceneV,scene)
    org=file[indT3[end]+1:indT2[1]-1]; push!(orgV,org)
    exp=file[indExp[end]+1:indT3[1]-1]; push!(expV,exp)
    #
    if exp=="153Plate0"
        indT=findall(x -> x==scene, data_plate0[:,:scene])[1]
        grp_sample=data_plate0[indT,:group]
        time_sample=data_plate0[indT,:EdU]
    elseif exp=="153Plate1"
        indT=findall(x -> x==scene, data_plate1[:,:scene])[1]
        grp_sample=data_plate1[indT,:group]
        time_sample=data_plate1[indT,:EdU]
    elseif exp=="153Plate2"
        indT=findall(x -> x==scene, data_plate2[:,:scene])[1]
        grp_sample=data_plate2[indT,:group]
        time_sample=data_plate2[indT,:EdU]
    else
        println("at index: ", rowI, "is ", exp_sample)
    end
    push!(groupV, grp_sample)
    push!(timeEDU, time_sample) 
    #
    dataH=CSV.read(path*"\\$nameFile", DataFrame); 
    TotalNuclei=size(dataH, 1);
    wt=0; cancer=0; ph3_total=0; edu_total=0; 
    ph3_wt=0;
    ph3_cancer=0;
    edu_wt=0;
    edu_cancer=0;        
    ph3_edu=0;
    ph3_edu_wt=0;
    ph3_edu_cancer=0;
    vec_wt=Vector{}(); vec_ph3High=Vector{}();
    vec_ph3=Vector{}(); vec_edu=Vector{}();
    for j in 1:Int.(TotalNuclei)
        VolTot=dataH[j, indexVol]; #dataH[j, indexCh1]dataH[j, indexCh2] + dataH[j, indexCh3];
        #edu
        VolEdU=dataH[j, index_volEDU];
        avgEdu=dataH[j, index_IntEDU]/dataH[j, indexVol];
        #ph3
        Volph3=dataH[j, index_volPH3];
        avgph3=dataH[j, index_IntPH3]/dataH[j, indexVol];
        #ph3_high
        Volph3high=dataH[j, index_volPH3_high];
        avgph3high=dataH[j, index_volPH3_high]/dataH[j, indexVol];
        #mcherry==WT
        Volwt=dataH[j, index_volCherry];
        avgwt=dataH[j, index_IntCherry]/dataH[j, indexVol];
        #
        vwt=0;
        if Volwt>alpha_mcherry*VolTot 
            vwt=1;
        end
        push!(vec_wt, vwt)
        vph3=0;
        if Volph3>alpha_ph3*VolTot 
            vph3=1;
        end
        push!(vec_ph3, vph3)
        vedu=0;
        if VolEdU>alpha_edu*VolTot 
            vedu=1;
        end
        push!(vec_edu, vedu)
        vph3H=0
        if Volph3high>alpha_ph3*VolTot                 
            vph3H=1
        end
        push!(vec_ph3High, vph3H)
        # no bias
        if grp_sample=="WT"
            vwt=1
        elseif grp_sample=="CRC"
            vwt=0
        end
        if vwt+vph3==2
            wt=wt+1;
            ph3_total=ph3_total+1;
            ph3_wt=ph3_wt+1;
            if vwt+vph3+vedu==3
                edu_total=edu_total+1;
                edu_wt=edu_wt+1;
                ph3_edu=ph3_edu+1;
                ph3_edu_wt=ph3_edu_wt+1;
            end
        elseif vph3==1
            cancer=cancer+1;
            ph3_total=ph3_total+1;
            ph3_cancer=ph3_cancer+1;
            if vedu==1
                edu_total=edu_total+1;
                edu_cancer=edu_cancer+1;
                ph3_edu=ph3_edu+1;
                ph3_edu_cancer=ph3_edu_cancer+1;
            end
        elseif vwt==1   
            wt=wt+1;
            if vedu==1
                edu_wt=edu_wt+1;
                edu_total=edu_total+1;
            end
        else
            cancer=cancer+1;
            if vedu==1
                edu_cancer=edu_cancer+1;
                edu_total=edu_total+1;
            end
        end
        #
    end
    dataH[!, "WT_tag"]=vec_wt
    dataH[!, "PH3_tag"]=vec_ph3
    dataH[!, "PH3high_tag"]=vec_ph3High
    dataH[!, "EDU_tag"]=vec_edu
    #
    fileOUT=dirOut*"\\Categories_"*nameFile;
    # py"genrate_result_image"(summary,path,tamaño imagen x y z)
    CSV.write("$fileOUT",  dataH, delim = '\t', use_mmap=false)  
    push!(totalNumb, TotalNuclei)
    push!(wtV, wt)
    push!(cancerV, cancer)
    push!(ph3_totalV, ph3_total)
    push!(edu_totalV, edu_total)
    push!(ph3_wtV, ph3_wt)
    push!(ph3_cancerV, ph3_cancer)
    push!(edu_wtV, edu_wt)
    push!(edu_cancerV, edu_cancer)
    push!(ph3_eduV, ph3_edu)
    push!(ph3_edu_wtV, ph3_edu_wt)
    push!(ph3_edu_cancerV, ph3_edu_cancer)     
end
dataFolder[!,:num3Dcells]=totalNumb;
CSV.write(pathO*"\\OSCAR_img_conditions.txt", dataFolder, delim = '\t', use_mmap=false)
#
col_labels = ["name", "Exp", "scene", "org", "ExpGroup", "EdUTime", "numbTotal", 
"wt", "cancer","PH3", "EdU",
"PH3+_WT", "PH3+_CANCER","EdU+_WT", "EdU+_CANCER",
"EdU+_PH3+", "EdU+_PH3+_WT", "EdU+_PH3+_CANCER"];
col_labels = Symbol.(col_labels);

dataV=hcat(namesV,expV,sceneV,orgV,groupV,timeEDU,totalNumb,
wtV, cancerV, ph3_totalV, edu_totalV, 
ph3_wtV, ph3_cancerV, edu_wtV, edu_cancerV,
ph3_eduV, ph3_edu_wtV, ph3_edu_cancerV);

dfmain=DataFrame(Tables.table(dataV, header = col_labels))
CSV.write(pathO*"\\data_noBias.txt", dfmain, delim = '\t', use_mmap=false)
dfmain_sorted=sort(dfmain,[:ExpGroup, :EdUTime, :scene])
