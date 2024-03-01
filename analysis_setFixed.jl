# Directory paths (assuming these have been removed as per your request)
# pathSummariesImage = "C:\\Users\\mario\\OneDrive\\Escritorio\\images\\sub\\out\\xy_files_channels_base2D\\channels_new";
# pathSummariesImage = "C:\\Users\\mario\\OneDrive\\Escritorio\\images\\sub\\out\\OSCAR_output_nuclei"; 
pathSummariesImage = "E:\\wholeorganoid\\preprocessed\\images\\out\\xy_files_base_xzMean\\channels_V2nuclei"
pathSummaryWhole = "E:\\wholeorganoid\\preprocessed\\images\\out\\InternMeasurements_OSCAR.txt";

# List files in directory
files = readdir(pathSummariesImage);

# Indices for summary data of each image 
indexX = 1; indexY = 2; indexZ = 3; indexC = 4; indexA = 5; indexB = 6; indexAng = 7; 
indexVx = 8; indexVx = 9; indexVz = 10; indexVol = 11; indexvolEll = 12; indexRaw = 13; 
indexAverageOscar = 14; indexXB = 15; indexYB = 16; indexWB = 17; indexHB = 18; 
indexXF = 19; indexYF = 20; indexFMax = 21; indexFmin = 22; indexFang = 23; indexXM = 24; 
indexYM = 25; indexSkew = 26; indexKurt = 27; indexMean = 28; indexSTD = 29; indexMedian = 30; 
indexCh1_vol = 31; indexCh1_Int = 32; indexCh2_vol = 33; indexCh2_Int = 34; indexType2D = 35;

# Initialize vectors
numbCh1 = Vector{}(); numbCh2 = Vector{}(); totalNumb = Vector{}();
names = Vector{}(); organoid = Vector{}(); set = Vector{}(); day = Vector{}(); scene = Vector{};

# Filtering parameter
alpha = 0.5;

for i in 1:length(files)
    file = files[i];
    if file[end-3:end] == ".txt" && file[1:8] == "Summary_"
        println("Processing the file ", file);
        
        # Extracting metadata from file name
        organoidT = parse(Int64, file[18]);
        setT = file[30];
        dayT = parse(Int64, file[35]);
        sceneT = parse(Int64, file[end-5:end-4]);
        
        # Reading summary data
        dataH = CSV.read(pathSummariesImage*"\\$file", DataFrame); 
        
        TotalNuclei = size(dataH, 1);
        ch1n = 0; ch2n = 0;
        
        for j in 1:Int.(TotalNuclei)
            VolTot = dataH[j, indexVol];
            VolCh1 = dataH[j, indexCh1_vol];
            VolCh2 = dataH[j, indexCh2_vol];
            avgCh1 = dataH[j, indexCh1_Int] / dataH[j, indexVol];
            avgCh2 = dataH[j, indexCh2_Int] / dataH[j, indexVol];
            
            if VolCh1 > alpha * VolTot #for this dataset,  WT channel == number 1
                ch2n += 1;
            else
                ch1n += 1;
            end
        end
        
        push!(numbCh1, ch1n);
        push!(numbCh2, ch2n);
        push!(totalNumb, TotalNuclei);
        push!(names, file[12:end-4]);
        push!(organoid, organoidT);
        push!(set, setT);
        push!(day, dayT);
        push!(scene, sceneT);
    end
end

col_labels = ["name", "Organoid", "set", "day", "scene", "numbTotal", "numbCancer", "numbWT"];
col_labels = Symbol.(col_labels);
dataV = hcat(names, organoid, set, day, scene, totalNumb, numbCh1, numbCh2);
objects3D_df = DataFrame(Tables.table(dataV, header = col_labels));

CSV.write(pathSummariesImage*"\\Channels_AlphaFixToNEWW_$alpha.txt", objects3D_df, delim = '\t');
