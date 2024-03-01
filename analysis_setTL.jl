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
indexCh1_vol=31; indexCh1_Int=32; indexCh2_vol=33; indexCh2_Int=34;
#
path="C:\\Users\\Ledes002\\Desktop\\time-lapse\\ultimo\\kernelXY_fixedXZ\\OSCAR_output_KernelByTL"
cd(path); files = readdir(path)
pathO="C:\\Users\\Ledes002\\Desktop\\time-lapse\\ultimo\\kernelXY_fixedXZ\\OSCAR_output_KernelByTL\\o"


varSearch="Summary_Exp";
files_Sel=filter(x -> occursin(varSearch, x), files);
alpha=0; alpha2=0.25;
IDv=Vector{}(); namesV=Vector{}(); frameV=Vector{}(); sceneV=Vector{}(); orgV=Vector{}(); expV=Vector{}(); 
numbCancer=Vector{}(); numbWT=Vector{}(); totalNumb=Vector{}();
for main_index in 1:length(files_Sel)
        #
        file = files_Sel[main_index];   
        indSumm=collect(findall.(varSearch, file)[1]); 
        nameFile="Exp"*file[indSumm[end]+1:end]; push!(namesV, nameFile);
        println("Processing the file2 ", nameFile);
        indExp=collect(findall.("Exp", file)[1]);      
        indT=collect(findall.("Frame-", file)[1]);      
        indT2=collect(findall.("_Scene-", file)[1]); 
        indT3=collect(findall.("_Org-", file)[1]);
        indT_end=collect(findall.(".", file)[1]);
        frame=file[indT[end]+1:indT2[1]-1]; frame=parse(Int64,frame); push!(frameV,frame)
        scene=file[indT2[end]+1:indT3[1]-1]; scene=parse(Int64,scene); push!(sceneV,scene)
        org=file[indT3[end]+1:indT_end[1]-1]; org=parse(Int64,org); push!(orgV,org)
        exp=file[indExp[end]+1:indT[1]-2]; exp=parse(Int64,exp); push!(expV,exp)
        ID="E-"*string(exp)*"_S-"*string(scene)*"_O-"*string(org); push!(IDv, ID)
        #
        dataH=CSV.read(path*"\\$file", DataFrame); 
        TotalNuclei=size(dataH, 1); 
        ch1n=0; ch2n=0; 
        for j in 1:Int.(TotalNuclei)
            VolTot=dataH[j, indexVol]; #dataH[j, indexCh1]dataH[j, indexCh2] + dataH[j, indexCh3];

            Ael=dataH[j, indexA]/2; 
            Bel=dataH[j, indexB]/2;
            Cel=dataH[j, indexC];
            volTot=(4/3)*pi*Ael*Bel*Cel;
            VolCh1=dataH[j, indexCh1_vol];
            VolCh2=dataH[j, indexCh2_vol];
            avgCh1=dataH[j, indexCh1_Int]/dataH[j, indexVol];
            avgCh2=dataH[j, indexCh2_Int]/dataH[j, indexVol];
            #
            if VolCh2>alpha*VolTot && VolCh2>alpha2*VolCh1 #for this dataset WT channel was saved as chnl 2
                #
                ch2n=ch2n+1;
            else
                ch1n=ch1n+1;
            end
        end
        push!(numbCancer, ch1n); push!(numbWT, ch2n); push!(totalNumb,TotalNuclei);
        #
end
col_labels = ["ID","name", "Exp", "frame", "scene", "org", "numbTotal", "numbCancer", "numbWT"];
col_labels = Symbol.(col_labels);
dataV=hcat(IDv,namesV,expV,frameV,sceneV,orgV,totalNumb,numbCancer,numbWT)
objects3D_df = DataFrame(Tables.table(dataV, header = col_labels))
dataF_sorted=sort(objects3D_df,[:Exp,:scene,:org,:frame])
# CSV.write(path*"\\Summary_TL_kernelsbyTL.tsv", dataF_sorted, delim = '\t')
##
names_to_exclude=["E-133_S-123_O-2","E-133_S-122_O-1", "E-133_S-22_O-1", "E-133_S-12_O-1",
"E-133_S-73_O-1", "E-133_S-68_O-1", "E-133_S-85_O-1", "E-137_S-41_O-1", 
"E-133_S-101_O-1","E-133_S-103_O-1","E-133_S-107_O-1","E-133_S-27_O-1","E-133_S-32_O-1","E-133_S-36_O-1"];
data_filtered=deepcopy(objects3D_df); cc=0;
for ind2 in 1:size(dataF_sorted,1)
    ID=dataF_sorted[ind2,:ID]
    a=collect(findall(x->x==ID,names_to_exclude))
    if length(a)>0
        cc=cc+1
        deleteat!(data_filtered,ind2)
        println("venga",ID)
    end
end
println("Han sido eliminados ", cc+1, " archivos, porque no valen!")
# CSV.write(path*"\\Summary_TL_kernelsbyTL_filtered.tsv", data_filtered, delim = '\t')
dataF_3=groupby(data_filtered, [:Exp,:scene,:org])
##
struct item
    name
    experiment
    scene
    organoid
    measurements
end
# wt 133>  15 wt, 20 MIX, 9 CANCER
scenes_Exp133_wt=[1,5,6,7,8,9,10,11,12,13,15,16,18,19,21,22,23,24,25,118,120,121,122,123,124,125,126,127,128,129]
scenes_Exp133_mix=[26,27,30,31,32,35,36,39,40,42,46,47,49,52,54,55,59,61,62,92,93,94,95,101,102,103,107,110,114,116]
scenes_Exp133_cancer=[64,65,66,67,68,73,75,77,78,80,81,82,84,85,86,90,91]
# 137> 8WT, 3 mix, 8 cancer
scenes_Exp137_wt=[66,67,69,70,71,73,74,75,76,78]
scenes_Exp137_mix=[47,48,52]
scenes_Exp137_cancer=[39,41,89,90,91,92]
# 139> 14 wt, 11 mix, 13 cancer
scenes_Exp139_wt=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]
scenes_Exp139_mix=[15,17,18,19,22,25,26,28,30,31,33]
scenes_Exp139_cancer=[40,42,43,44,45,46,47,49,50,51,52,53,54]
# 140>12 WT, 16 mix, 12 cancer
scenes_Exp140_wt=[1,2,3,4,5,6,7,8,9,10,11,14]
scenes_Exp140_mix=[17,18,25,26,27,29,36,39,43,44,45,46,47,48,50,51]
scenes_Exp140_cancer=[54,55,56,57,58,59,60,61,62,63,64,65]
length_wt=length(scenes_Exp133_wt)+length(scenes_Exp137_wt)+length(scenes_Exp139_wt)+length(scenes_Exp140_wt)
length_mix=length(scenes_Exp133_mix)+length(scenes_Exp137_mix)+length(scenes_Exp139_mix)+length(scenes_Exp140_mix)
length_cancer=length(scenes_Exp133_cancer)+length(scenes_Exp137_cancer)+length(scenes_Exp139_cancer)+length(scenes_Exp140_cancer)

WT_pure=Vector{}(); 
cancer_pure=Vector{}(); 
mix=Vector{}(); 
samplesV=Vector{}(); colLabels=Vector{}(); 
WT_df=DataFrame()
cancer_df=DataFrame()
mix_df=DataFrame()
colLabels_wt=Vector{}(); 
colLabels_cancer=Vector{}(); 
colLabels_mix=Vector{}(); 

maxt=15; #max num of frames exp139
WTpure_bytime_total=ones(maxt,length(dataF_3))*-1
WTpure_bytime_cancer=ones(maxt,length(dataF_3))*-1
WTpure_bytime_wt=ones(maxt,length(dataF_3))*-1
WTpure_bytime_percWT=ones(maxt,length(dataF_3))*-1
WTpure_bytime_percCancer=ones(maxt,length(dataF_3))*-1
ccWT=0; 
mix_bytime_total=ones(maxt,length(dataF_3))*-1
mix_bytime_cancer=ones(maxt,length(dataF_3))*-1
mix_bytime_wt=ones(maxt,length(dataF_3))*-1
mix_bytime_percWT=ones(maxt,length(dataF_3))*-1
mix_bytime_percCancer=ones(maxt,length(dataF_3))*-1
ccMIX=0;
cancerPure_bytime_total=ones(maxt,length(dataF_3))*-1
cancerPure_bytime_cancer=ones(maxt,length(dataF_3))*-1
cancerPure_bytime_wt=ones(maxt,length(dataF_3))*-1
cancerPure_bytime_percWT=ones(maxt,length(dataF_3))*-1
cancerPure_bytime_percCancer=ones(maxt,length(dataF_3))*-1
ccCANCER=0;
for main_index in 1:length(dataF_3)
    dataT=dataF_3[main_index];
    nameFile=dataT[1,:name];
    exp=dataT[1,:Exp]
    scn=dataT[1,:scene]
    orgn=dataT[1,:org]
    ID=dataT[1,:ID]
    dataT=sort(select(dataT, [:ID, :frame, :numbTotal, :numbCancer, :numbWT]),:frame)
    push!(samplesV,item(nameFile, exp, scn, orgn, dataT))
    ##
    println(nameFile)
    push!(colLabels, dataT[1,:ID])
    ##
    phenotype=NaN
    if exp==133
        a=findall(x -> x==scn, scenes_Exp133_wt)
        c=findall(x -> x==scn, scenes_Exp133_mix)
        b=findall(x -> x==scn, scenes_Exp133_cancer)
        if length(a)>0
            phenotype="WT"
            push!(WT_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(b)>0
            phenotype="cancer"
            push!(cancer_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(c)>0
            phenotype="mix"
            push!(mix,item(nameFile, exp, scn, orgn, dataT))
        end
    elseif exp==137
        a=findall(x -> x==scn, scenes_Exp137_wt)
        c=findall(x -> x==scn, scenes_Exp137_mix)
        b=findall(x -> x==scn, scenes_Exp137_cancer)
        if length(a)>0
            phenotype="WT"
            push!(WT_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(b)>0
            phenotype="cancer"
            push!(cancer_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(c)>0
            phenotype="mix"
            push!(mix,item(nameFile, exp, scn, orgn, dataT))
        end
    elseif exp==139
        a=findall(x -> x==scn, scenes_Exp139_wt)
        c=findall(x -> x==scn, scenes_Exp139_mix)
        b=findall(x -> x==scn, scenes_Exp139_cancer)
        if length(a)>0
            phenotype="WT"
            push!(WT_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(b)>0
            phenotype="cancer"
            push!(cancer_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(c)>0
            phenotype="mix"
            push!(mix,item(nameFile, exp, scn, orgn, dataT))
        end
    elseif exp==140
        a=findall(x -> x==scn, scenes_Exp140_wt)
        c=findall(x -> x==scn, scenes_Exp140_mix)
        b=findall(x -> x==scn, scenes_Exp140_cancer)
        if length(a)>0
            phenotype="WT"
            push!(WT_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(b)>0
            phenotype="cancer"
            push!(cancer_pure,item(nameFile, exp, scn, orgn, dataT))
        elseif length(c)>0
            phenotype="mix"
            push!(mix,item(nameFile, exp, scn, orgn, dataT))
        end
    end
    if phenotype=="WT"
        ccWT=ccWT+1; 
        WT_df=vcat(WT_df, dataT);
        for indT in 1:size(dataT,1)
            WTpure_bytime_total[indT,ccWT]=dataT[indT,:numbTotal]
            WTpure_bytime_cancer[indT,ccWT]=dataT[indT,:numbCancer]
            WTpure_bytime_wt[indT,ccWT]=dataT[indT,:numbWT]   
            #
            WTpure_bytime_percWT[indT,ccWT]=dataT[indT,:numbWT]/dataT[indT,:numbTotal]
            WTpure_bytime_percCancer[indT,ccWT]=dataT[indT,:numbCancer]/dataT[indT,:numbTotal]
            #
        end
        push!(colLabels_wt, dataT[1,:ID])
    elseif phenotype=="cancer"
        ccCANCER=ccCANCER+1;      
        cancer_df=vcat(cancer_df, dataT)
        for indT in 1:size(dataT,1)
            cancerPure_bytime_total[indT,ccCANCER]=dataT[indT,:numbTotal]
            cancerPure_bytime_cancer[indT,ccCANCER]=dataT[indT,:numbCancer]
            cancerPure_bytime_wt[indT,ccCANCER]=dataT[indT,:numbWT]      
            #
            cancerPure_bytime_percWT[indT,ccCANCER]=dataT[indT,:numbWT]/dataT[indT,:numbTotal]
            cancerPure_bytime_percCancer[indT,ccCANCER]=dataT[indT,:numbCancer]/dataT[indT,:numbTotal]
            #  
        end
        push!(colLabels_cancer, dataT[1,:ID])
    elseif phenotype=="mix"
        ccMIX=ccMIX+1;
        mix_df=vcat(mix_df, dataT)
        for indT in 1:size(dataT,1)
            mix_bytime_total[indT,ccMIX]=dataT[indT,:numbTotal]
            mix_bytime_cancer[indT,ccMIX]=dataT[indT,:numbCancer]
            mix_bytime_wt[indT,ccMIX]=dataT[indT,:numbWT] 
            #
            mix_bytime_percWT[indT,ccMIX]=dataT[indT,:numbWT]/dataT[indT,:numbTotal]
            mix_bytime_percCancer[indT,ccMIX]=dataT[indT,:numbCancer]/dataT[indT,:numbTotal]
            #  
        end
        push!(colLabels_mix, dataT[1,:ID])
    else
        println("estoque es?", scn, "exp",exp); 
    end
    ##
end
WTpure_bytime_total=WTpure_bytime_total[:,1:ccWT]
WTpure_bytime_cancer=WTpure_bytime_cancer[:,1:ccWT]
WTpure_bytime_wt=WTpure_bytime_wt[:,1:ccWT]
WTpure_bytime_percWT=WTpure_bytime_percWT[:,1:ccWT]
WTpure_bytime_percCancer=WTpure_bytime_percCancer[:,1:ccWT]

mix_bytime_total=mix_bytime_total[:,1:ccMIX]
mix_bytime_cancer=mix_bytime_cancer[:,1:ccMIX]
mix_bytime_wt=mix_bytime_wt[:,1:ccMIX]
mix_bytime_percWT=mix_bytime_percWT[:,1:ccMIX]
mix_bytime_percCancer=mix_bytime_percCancer[:,1:ccMIX]

cancerPure_bytime_total=cancerPure_bytime_total[:,1:ccCANCER]
cancerPure_bytime_cancer=cancerPure_bytime_cancer[:,1:ccCANCER]
cancerPure_bytime_wt=cancerPure_bytime_wt[:,1:ccCANCER]
cancerPure_bytime_percWT=cancerPure_bytime_percWT[:,1:ccCANCER]
cancerPure_bytime_percCancer=cancerPure_bytime_percCancer[:,1:ccCANCER]

totalM= DataFrame(Tables.table(WTpure_bytime_total, header = colLabels_wt))
cancerM= DataFrame(Tables.table(WTpure_bytime_cancer, header = colLabels_cancer))
wtM= DataFrame(Tables.table(WTpure_bytime_wt, header = colLabels_mix))
CSV.write(pathO*"\\total_wtpure.txt", totalM, delim = '\t')
CSV.write(pathO*"\\cancer_wtpure.txt", cancerM, delim = '\t')
CSV.write(pathO*"\\wt_wtpure.txt", wtM, delim = '\t')

totalM= DataFrame(Tables.table(cancerPure_bytime_total, header = colLabels_wt))
cancerM= DataFrame(Tables.table(cancerPure_bytime_cancer, header = colLabels_cancer))
wtM= DataFrame(Tables.table(cancerPure_bytime_wt, header = colLabels_mix))
CSV.write(pathO*"\\total_cancerpure.txt", totalM, delim = '\t')
CSV.write(pathO*"\\cancer_cancerpure.txt", cancerM, delim = '\t')
CSV.write(pathO*"\\wt_cancerpure.txt", wtM, delim = '\t')

totalM= DataFrame(Tables.table(mix_bytime_total, header = colLabels_wt))
cancerM= DataFrame(Tables.table(mix_bytime_cancer, header = colLabels_cancer))
wtM= DataFrame(Tables.table(mix_bytime_wt, header = colLabels_mix))
CSV.write(pathO*"\\total_mix.txt", totalM, delim = '\t')
CSV.write(pathO*"\\cancer_mix.txt", cancerM, delim = '\t')
CSV.write(pathO*"\\wt_mix.txt", wtM, delim = '\t')







