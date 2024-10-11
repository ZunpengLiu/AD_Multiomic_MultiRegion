#!/bin/sh

# Define input and output directories
int=~/Input
out=~/Output

# Define the path to ChromHMM.jar
chromhmm_jar="/net/bmc-lab5/data/kellis/users/zunpeng/04_Softwares/ChromHMM/n/ChromHMM/ChromHMM.jar"

# Define the list of brain tissues and corresponding sample IDs, defined from the Epimap data
brain_samples=(
    "BSS00369:BRN.FTL.CTX"
    "BSS00371:BRN.FTL.CTX"
    "BSS01125:BRN.HPC"
    "BSS01126:BRN.HPC"
    "BSS01124:BRN.HPC"
    "BSS01271:BRN.MID.FTL"
    "BSS01272:BRN.MID.FTL"
    "BSS01273:BRN.MID.FTL.GYR"
    "BSS00077:BRN.AG"
    "BSS00078:BRN.AG"
    "BSS00089:BRN.AST"
    "BSS00090:BRN.AST.CRBLLM"
    "BSS00091:BRN.AST.HPC"
)

# Loop through each brain sample
for sample in "${brain_samples[@]}"; do
    # Split the sample into ID and name
    IFS=":" read -r sample_id sample_name <<< "$sample"
    
    # Define the segments file for this sample
    segments="~/04_Softwares/ChromHMM/ChromHMM/epiMap/${sample_id}_18_CALLS_segments.bed"
    
    # Check if the segments file exists
    if [ -f "$segments" ]; then
        # Run the ChromHMM OverlapEnrichment command
        echo "Running ChromHMM OverlapEnrichment for $sample_id ($sample_name) using $segments"
        java -mx1000000M -jar $chromhmm_jar OverlapEnrichment -labels -b 1 $segments $int $out/${sample_id}_BRAIN_$sample_name
    else
        echo "Segments file $segments not found, skipping $sample_id ($sample_name)"
    fi
done

echo "All ChromHMM OverlapEnrichment jobs completed!"
