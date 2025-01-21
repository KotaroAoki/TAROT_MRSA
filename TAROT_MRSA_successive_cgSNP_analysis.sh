#!/bin/bash
#TAROT_MRSA_successive_cgSNP_analysis

# Create output directory and get its path
output_dir=`date +%Y%m%d%H%M`
mkdir $output_dir
cd  $output_dir
output_dir_fullpath=`pwd` 
cd ..

### Set variables ###
start_dir_fullpath=`pwd` 
start_dir_name=`basename $start_dir_fullpath`
reference_dir="/tmp/Saureus_reference_strain_genome" 
bam_dir="/tmp/Saureus_bam"
TOOLSPATH="/tmp/program"
################

# Current directory
cfullname=`pwd`
cname=`basename $cfullname`

#.fastq decompression
unpigz -f *.gz

# Concatenate fastq files
cat  {"fastq_runid_"*,*"barcode"*,*"RBK114"*}.fastq > cat_reads.fastq
# Nanoplot reads QC (before trim)
NanoPlot --fastq cat_reads.fastq --loglength -o cat_reads_nanoplot_results
# Nanofilt filtering
NanoFilt -q 10 --headcrop 75 -l 1000 < cat_reads.fastq > cat_reads_trimmed_q10_minlen1k.fastq
# Nanoplot reads QC (after trim)
NanoPlot --fastq cat_reads_trimmed_q10_minlen1k.fastq --loglength -o cat_reads_trimmed_q10_minlen1k_nanoplot_results

# Downsampling (auto-calculated to about 100x depth for a 5Mb genome)
# Step 1: Calculate the average read length
total_read_nucleotide=`grep "Total bases" cat_reads_trimmed_q10_minlen1k_nanoplot_results/NanoStats.txt | awk -F":" '{print $2}' | perl -pe 's/ |\,//g'`
# Step 2: Calculate the number of reads needed to achieve 100x coverage
# Formula: (genome size * coverage) / (2 * average read length)
genome_size=5000000
coverage=100
target_nucleotide=$((genome_size * coverage))
sampling_rate=$(echo "scale=5; $target_nucleotide / $total_read_nucleotide" | bc)
# Step 3: Downsample the input FASTQ files
if (( $(echo "$sampling_rate >= 0 && $sampling_rate < 1" | bc -l) )); then
    # Execute if sampling_rate is between 0 and 1
    seqkit sample -p $sampling_rate cat_reads_trimmed_q10_minlen1k.fastq > cat_reads_trimmed_q10_minlen1k_depth100.fastq
elif (( $(echo "$sampling_rate >= 1" | bc -l) )); then
    # Execute if sampling_rate is 1 or greater
    cat cat_reads_trimmed_q10_minlen1k.fastq > cat_reads_trimmed_q10_minlen1k_depth100.fastq
fi

# Self-error correction
dechat -i cat_reads_trimmed_q10_minlen1k_depth100.fastq -o cat_reads_trimmed_q10_minlen1k_depth100_dechat -t 20

# Convert fasta to fastq using wgsim
wgsim -h -1 300 -2 300 -N1000000 -e0 -r0 -R0 -X0 -d0 cat_reads_trimmed_q10_minlen1k_depth100_dechat.ec.fa cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq

# Search for the best reference (specifying all genome files ".fasta" in $reference_dir)
# MLST only S. aureus just in case
rm -f cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_stringMLST.txt
stringMLST.py --predict -P /tmp/program/stringMLST/datasets/Staphylococcus_aureus/Staphylococcus_aureus -1 cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq -2 cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq -o cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_stringMLST.txt

MLST=`awk 'NR>1 {gsub(/\*/, "", $1); print $9}' cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_stringMLST.txt`

db_STs=`ls $reference_dir | cut -d "." -f 1 | sort | uniq`

found=0
for db_ST in $db_STs; do
    if [[ "$db_ST" == ST"$MLST" ]]; then
        found=1
        break
    fi
done

if [[ $found -eq 1 ]]; then
    bestref=`echo ST$MLST`
    match_ref="mlst"
else
    # Search for the best reference (specifying all genome files ".fasta" in $reference_dir)
    reference_genome=`find $reference_dir -maxdepth 1 -type f -name "*.fasta" ! -name ".*" | tr '\n' ',' | sed 's/,$//'`
    bbsplit.sh -Xmx20g ref=$reference_genome in1=cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq in2=cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq outu1=clean1.fq.gz outu2=clean2.fq.gz refstats=refstats.out basename=out_%.fq.gz

    # Extract the line with the maximum %unambiguousReads
    bestref=`awk 'NR==1 {max=$2; max_line=$0} $2>max {max=$2; max_line=$0} END{print max_line}' refstats.out | awk '{print $1}'`

    match_ref="bbmap"
fi

# Mapping
cat cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq > cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R12.fastq
bwa index $reference_dir/$bestref".fasta"
bwa bwasw -t 16 $reference_dir/$bestref".fasta" cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R12.fastq > out12.sam
samtools view -@ 16 -S -b out12.sam > out12.bam
samtools sort -@ 16 out12.bam -o $cname"_"$bestref"_out12.sorted.bam"
samtools index $cname"_"$bestref"_out12.sorted.bam"

# Branch based on the value of match_ref
if [ "$match_ref" = "mlst" ]; then
    # Store BAM copies by reference
    mkdir -p "$bam_dir" # Create directory if $bam_dir does not exist
    mkdir -p "$bam_dir/$bestref" # Create directory if $bestref directory does not exist
    cp $cname"_"$bestref"_out12.sorted.bam" $bam_dir/$bestref
    cp $cname"_"$bestref"_out12.sorted.bam.bai" $bam_dir/$bestref

    # Create bam_list
    find $bam_dir/$bestref -name "*sorted.bam" -type f | sort > $output_dir_fullpath/bam_list.txt

    # Create name_list
    find $bam_dir/$bestref -name "*sorted.bam" -type f |sort | sed 's|.*/||;s/\.sorted\.bam//' | sed 's/^/>/' | sed 's/_out12//' | awk -F "_" '{print $1}' > $output_dir_fullpath/name_list.txt

    # Delete unnecessary files
    rm *.sam out12.bam cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq clean1.fq.gz clean2.fq.gz cat_reads.fastq cat_reads_trimmed_q10_minlen1k.fastq cat_reads_trimmed_q10_minlen1k_depth100.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R12.fastq out_*.fq.gz 

    # Compress .fastq files
    pigz -f *.fastq

    #### Do not execute SNP extraction if the number of strains is less than 4 ####
    # Count the number of lines in the file
    linecount=$(wc -l < "$output_dir_fullpath/name_list.txt")

    if [ "$linecount" -ge 4 ]; then

        # SNP extraction
        # Move to output directory
        cd $output_dir_fullpath
        # mpileup and consensus calling
        samtools mpileup -B -f $reference_dir/$bestref".fasta" -b bam_list.txt | java -jar $TOOLSPATH/VarScan.v2.3.9.jar mpileup2cns > mpileup2cns.txt
        # Create core genome multifasta
        grep -v "N:" mpileup2cns.txt | grep -v "/" | grep "Pass" | perl -pe 's/\t/ /g' | cut -d " " -f 11-300 | perl -pe 's/:/ /g' | awk '{print $1,$7,$13,$19,$25,$31,$37,$43,$49,$55,$61,$67,$73,$79,$85,$91,$97,$103,$109,$115,$121,$127,$133,$139,$145,$151,$157,$163,$169,$175,$181,$187,$193,$199,$205,$211,$217,$223,$229,$235,$241,$247,$253,$259,$265,$271,$277,$283,$289,$295,$301,$307,$313,$319,$325,$331,$337,$343,$349,$355,$361,$367,$373,$379,$385,$391,$397,$403,$409,$415,$421,$427,$433,$439,$445,$451,$457,$463,$469,$475,$481,$487,$493,$499,$505,$511,$517,$523,$529,$535,$541,$547,$553,$559,$565,$571,$577,$583,$589,$595,$601,$607,$613,$619,$625,$631,$637,$643,$649,$655,$661,$667,$673,$679,$685,$691,$697,$703,$709,$715,$721,$727,$733,$739,$745,$751,$757,$763,$769,$775,$781,$787,$793,$799,$805,$811,$817,$823,$829,$835,$841,$847,$853,$859,$865,$871,$877,$883,$889}' | grep -v -e Q -e W -e R -e Y -e U -e I -e O -e P -e S -e D -e F -e H -e J -e K -e L -e Z -e X -e V -e B -e N -e M -e q -e w -e r -e y -e u -e i -e o -e p -e s -e d -e f -e h -e j -e k -e l -e z -e x -e v -e b -e n -e m | perl -pe 's/ /\t/g' | perl $TOOLSPATH/transpose1_v2.pl | perl -pe 's/\t//g' > mpileup2cns_tran.txt
        paste name_list.txt mpileup2cns_tran.txt | perl -pe 's/\t/\n/g' > mpileup2cns_tran.fasta
        # File transformation (creating original)
        grep -v "N:" mpileup2cns.txt | grep -v "/" | grep "Pass" | perl -pe 's/\t/ /g' | cut -d " " -f 11-300 | perl -pe 's/:/ /g' | awk '{print $1,$7,$13,$19,$25,$31,$37,$43,$49,$55,$61,$67,$73,$79,$85,$91,$97,$103,$109,$115,$121,$127,$133,$139,$145,$151,$157,$163,$169,$175,$181,$187,$193,$199,$205,$211,$217,$223,$229,$235,$241,$247,$253,$259,$265,$271,$277,$283,$289,$295,$301,$307,$313,$319,$325,$331,$337,$343,$349,$355,$361,$367,$373,$379,$385,$391,$397,$403,$409,$415,$421,$427,$433,$439,$445,$451,$457,$463,$469,$475,$481,$487,$493,$499,$505,$511,$517,$523,$529,$535,$541,$547,$553,$559,$565,$571,$577,$583,$589,$595,$601,$607,$613,$619,$625,$631,$637,$643,$649,$655,$661,$667,$673,$679,$685,$691,$697,$703,$709,$715,$721,$727,$733,$739,$745,$751,$757,$763,$769,$775,$781,$787,$793,$799,$805,$811,$817,$823,$829,$835,$841,$847,$853,$859,$865,$871,$877,$883,$889}' | grep -v -e Q -e W -e R -e Y -e U -e I -e O -e P -e S -e D -e F -e H -e J -e K -e L -e Z -e X -e V -e B -e N -e M -e q -e w -e r -e y -e u -e i -e o -e p -e s -e d -e f -e h -e j -e k -e l -e z -e x -e v -e b -e n -e m | perl -pe 's/ /\t/g' > mpileup2cns_original.txt
        # File transformation (creating common, comparison expressions change depending on number of strains)
        compare_strains=`grep ">" name_list.txt | wc -l`   # Pick up column descriptions based on the number of strains to compare
        conpare_strains_object=`cat $TOOLSPATH/core_genome_SNPs_strain_number_database.txt | awk -F "\t" -v number="${compare_strains}" '$1==number {print $3}'`
        eval ${conpare_strains_object} # Inserts comparison text corresponding to the number of strains into awk command. Supports up to 300 strains.
        # Compare original and common to leave only SNP positions
        diff mpileup2cns_original.txt mpileup2cns_common.txt | grep "<" | cut -d" " -f "2-100" > mpileup2cns_nondupli.txt
        # Generate multifasta with only SNPs
        perl $TOOLSPATH/transpose1_v2.pl mpileup2cns_nondupli.txt | perl -pe 's/\t//g' > mpileup2cns_nondupli_tran.txt
        paste name_list.txt mpileup2cns_nondupli_tran.txt | perl -pe 's/\t/\n/g' > mpileup2cns_nondupli_tran.fasta
        # Convert fasta to phylip format
        python $TOOLSPATH/convert_fasta_to_phylip.py -i mpileup2cns_nondupli_tran.fasta -o mpileup2cns_nondupli_tran.phy
        # PhyML (create starting tree for ClonalFrameML)
        phyml -i mpileup2cns_nondupli_tran.phy -m GTR -b 100
        # RAxML (tree creation for reference before ClonalFrameML)
        raxmlHPC-PTHREADS -f a -s mpileup2cns_nondupli_tran.fasta -n mpileup2cns_nondupli_tran_RAxML -m GTRGAMMA -T 4 -x 12345 -p12345 -# 100
        # Delete unnecessary files and compress
        rm mpileup2cns_tran.txt mpileup2cns_original.txt mpileup2cns_common.txt mpileup2cns_nondupli.txt mpileup2cns_nondupli_tran.txt
        pigz mpileup2cns.txt
        # Detection of homologous recombination regions (ClonalFrameML)
        # Calculate kappa
        AG=`grep "A <-> G" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        CT=`grep "C <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        AC=`grep "A <-> C" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        AT=`grep "A <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        CG=`grep "C <-> G" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        GT=`grep "G <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        kappa=`echo "scale=5; (($AG + $CT) / 2) /  (($AC + $AT + $CG + $GT) / 4)" | bc`
           
        # Execute ClonalFrameML
        ClonalFrameML mpileup2cns_nondupli_tran.phy_phyml_tree.txt mpileup2cns_tran.fasta mpileup2cns_tran -kappa $kappa -emsim 100
        # Exclude homologous recombination regions
        python $TOOLSPATH/remove_rec_from_ClonalFrameML_output_v1.1.py -a mpileup2cns_tran.fasta -o mpileup2cns_tran.importation_status.txt
        # Format sequences with only SNPs
        grep -v ">" mpileup2cns_tran.recRemoved.fasta | perl -pe 's/\n/_/g' | fold -1 | perl -pe 's/\n/\t/g' | perl -pe 's/\t_\t/\n/g' | perl $TOOLSPATH/transpose1_v2.pl  | sed '$d' > mpileup2cns_tran.recRemoved.tran
            
        # Create common file (File transformation for common creation, comparison expressions depend on number of strains)
        compare_strains_aftre_clonalframeml=`grep ">" name_list.txt | wc -l`   # Pick up column descriptions based on the number of strains to compare
        conpare_strains_object_clonalframeml=`cat $TOOLSPATH/core_genome_SNPs_strain_number_database_after_clonalframeml.txt | awk -F "\t" -v number="${compare_strains}" '$1==number {print $3}'`
        eval ${conpare_strains_object_clonalframeml} 
        diff mpileup2cns_tran.recRemoved.tran mpileup2cns_tran.recRemoved.tran.common | grep "<" | cut -d " " -f 2-100 > mpileup2cns_tran.recRemoved.tran.common.snp
        perl $TOOLSPATH/transpose1_v2.pl mpileup2cns_tran.recRemoved.tran.common.snp | perl -pe 's/\t//g' > mpileup2cns_tran.recRemoved.tran.common.snp.prefas
        paste name_list.txt mpileup2cns_tran.recRemoved.tran.common.snp.prefas | perl -pe 's/\t/\n/g' > mpileup2cns_tran.recRemoved.tran.common.snp.fasta
        rm mpileup2cns_tran.recRemoved.tran mpileup2cns_tran.recRemoved.tran.common mpileup2cns_tran.recRemoved.tran.common.snp mpileup2cns_tran.recRemoved.tran.common.snp.prefas
        # Molecular phylogenetic analysis (RAxML)
        rm -f *RAxML*
        raxmlHPC-PTHREADS -f a -s mpileup2cns_tran.recRemoved.tran.common.snp.fasta -n mpileup2cns_tran.recRemoved.tran.common.snp_RAxML -m GTRGAMMA -T 4 -x 12345 -p12345 -# 1000
        # Count SNPs
        snp-dists mpileup2cns_tran.recRemoved.tran.common.snp.fasta > mpileup2cns_tran.recRemoved.tran.common.snp_distance.tsv
    fi

else
    # Matching by bbsplit
    # Store BAM copies by reference
    mkdir -p "$bam_dir" # Create directory if $bam_dir does not exist
    mkdir -p "$bam_dir/${bestref}_bbsplit" # Create directory if ${bestref}_bbsplit does not exist
    cp $1_$bestref"_out12.sorted.bam" $bam_dir/${bestref}_bbsplit
    cp $1_$bestref"_out12.sorted.bam.bai" $bam_dir/${bestref}_bbsplit

    # Create bam_list
    find $bam_dir/${bestref}_bbsplit -name "*sorted.bam" -type f | sort > $output_dir_fullpath/bam_list.txt

    # Create name_list
    find $bam_dir/${bestref}_bbsplit -name "*sorted.bam" -type f |sort | sed 's|.*/||;s/\.sorted\.bam//' | sed 's/^/>/' | sed 's/_out12//' | awk -F "_" '{print $1}' > $output_dir_fullpath/name_list.txt

    # Delete unnecessary files
    rm *.sam out12.bam cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R1.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R2.fastq clean1.fq.gz clean2.fq.gz cat_reads.fastq cat_reads_trimmed_q10_minlen1k.fastq cat_reads_trimmed_q10_minlen1k_depth100.fastq cat_reads_trimmed_q10_minlen1k_depth100_dechat_wgsim_R12.fastq out_*.fq.gz 

    # Compress .fastq files
    pigz -f *.fastq

    #### Do not execute SNP extraction if the number of strains is less than 4 ####
    # Count the number of lines in the file
    linecount=$(wc -l < "$output_dir_fullpath/name_list.txt")

    if [ "$linecount" -ge 4 ]; then

        # SNP extraction
        # Move to output directory
        cd $output_dir_fullpath
        # mpileup and consensus calling
        samtools mpileup -B -f $reference_dir/$bestref"_bbsplit.fasta" -b bam_list.txt | java -jar $TOOLSPATH/VarScan.v2.3.9.jar mpileup2cns > mpileup2cns.txt
        # Create core genome multifasta
        grep -v "N:" mpileup2cns.txt | grep -v "/" | grep "Pass" | perl -pe 's/\t/ /g' | cut -d " " -f 11-300 | perl -pe 's/:/ /g' | awk '{print $1,$7,$13,$19,$25,$31,$37,$43,$49,$55,$61,$67,$73,$79,$85,$91,$97,$103,$109,$115,$121,$127,$133,$139,$145,$151,$157,$163,$169,$175,$181,$187,$193,$199,$205,$211,$217,$223,$229,$235,$241,$247,$253,$259,$265,$271,$277,$283,$289,$295,$301,$307,$313,$319,$325,$331,$337,$343,$349,$355,$361,$367,$373,$379,$385,$391,$397,$403,$409,$415,$421,$427,$433,$439,$445,$451,$457,$463,$469,$475,$481,$487,$493,$499,$505,$511,$517,$523,$529,$535,$541,$547,$553,$559,$565,$571,$577,$583,$589,$595,$601,$607,$613,$619,$625,$631,$637,$643,$649,$655,$661,$667,$673,$679,$685,$691,$697,$703,$709,$715,$721,$727,$733,$739,$745,$751,$757,$763,$769,$775,$781,$787,$793,$799,$805,$811,$817,$823,$829,$835,$841,$847,$853,$859,$865,$871,$877,$883,$889}' | grep -v -e Q -e W -e R -e Y -e U -e I -e O -e P -e S -e D -e F -e H -e J -e K -e L -e Z -e X -e V -e B -e N -e M -e q -e w -e r -e y -e u -e i -e o -e p -e s -e d -e f -e h -e j -e k -e l -e z -e x -e v -e b -e n -e m | perl -pe 's/ /\t/g' | perl $TOOLSPATH/transpose1_v2.pl | perl -pe 's/\t//g' > mpileup2cns_tran.txt
        paste name_list.txt mpileup2cns_tran.txt | perl -pe 's/\t/\n/g' > mpileup2cns_tran.fasta
        # File transformation (creating original)
        grep -v "N:" mpileup2cns.txt | grep -v "/" | grep "Pass" | perl -pe 's/\t/ /g' | cut -d " " -f 11-300 | perl -pe 's/:/ /g' | awk '{print $1,$7,$13,$19,$25,$31,$37,$43,$49,$55,$61,$67,$73,$79,$85,$91,$97,$103,$109,$115,$121,$127,$133,$139,$145,$151,$157,$163,$169,$175,$181,$187,$193,$199,$205,$211,$217,$223,$229,$235,$241,$247,$253,$259,$265,$271,$277,$283,$289,$295,$301,$307,$313,$319,$325,$331,$337,$343,$349,$355,$361,$367,$373,$379,$385,$391,$397,$403,$409,$415,$421,$427,$433,$439,$445,$451,$457,$463,$469,$475,$481,$487,$493,$499,$505,$511,$517,$523,$529,$535,$541,$547,$553,$559,$565,$571,$577,$583,$589,$595,$601,$607,$613,$619,$625,$631,$637,$643,$649,$655,$661,$667,$673,$679,$685,$691,$697,$703,$709,$715,$721,$727,$733,$739,$745,$751,$757,$763,$769,$775,$781,$787,$793,$799,$805,$811,$817,$823,$829,$835,$841,$847,$853,$859,$865,$871,$877,$883,$889}' | grep -v -e Q -e W -e R -e Y -e U -e I -e O -e P -e S -e D -e F -e H -e J -e K -e L -e Z -e X -e V -e B -e N -e M -e q -e w -e r -e y -e u -e i -e o -e p -e s -e d -e f -e h -e j -e k -e l -e z -e x -e v -e b -e n -e m | perl -pe 's/ /\t/g' > mpileup2cns_original.txt
        # File transformation (creating common, comparison expressions depend on number of strains)
        compare_strains=`grep ">" name_list.txt | wc -l`   # Pick up column descriptions based on the number of strains to compare
        conpare_strains_object=`cat $TOOLSPATH/core_genome_SNPs_strain_number_database.txt | awk -F "\t" -v number="${compare_strains}" '$1==number {print $3}'`
        eval ${conpare_strains_object}
        # Compare original and common to leave only SNP positions
        diff mpileup2cns_original.txt mpileup2cns_common.txt | grep "<" | cut -d" " -f "2-100" > mpileup2cns_nondupli.txt
        # Generate multifasta with only SNPs
        perl $TOOLSPATH/transpose1_v2.pl mpileup2cns_nondupli.txt | perl -pe 's/\t//g' > mpileup2cns_nondupli_tran.txt
        paste name_list.txt mpileup2cns_nondupli_tran.txt | perl -pe 's/\t/\n/g' > mpileup2cns_nondupli_tran.fasta
        # Convert fasta to phylip format
        python $TOOLSPATH/convert_fasta_to_phylip.py -i mpileup2cns_nondupli_tran.fasta -o mpileup2cns_nondupli_tran.phy
        # PhyML (create starting tree for ClonalFrameML)
        phyml -i mpileup2cns_nondupli_tran.phy -m GTR -b 100
        # RAxML (tree creation for reference before ClonalFrameML)
        raxmlHPC-PTHREADS -f a -s mpileup2cns_nondupli_tran.fasta -n mpileup2cns_nondupli_tran_RAxML -m GTRGAMMA -T 4 -x 12345 -p12345 -# 100
        # Delete unnecessary files and compress
        rm mpileup2cns_tran.txt mpileup2cns_original.txt mpileup2cns_common.txt mpileup2cns_nondupli.txt mpileup2cns_nondupli_tran.txt
        pigz mpileup2cns.txt
        # Detection of homologous recombination regions (ClonalFrameML)
        # Calculate kappa
        AG=`grep "A <-> G" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        CT=`grep "C <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        AC=`grep "A <-> C" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        AT=`grep "A <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        CG=`grep "C <-> G" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        GT=`grep "G <-> T" mpileup2cns_nondupli_tran.phy_phyml_stats.txt | awk '{print $4}'`
        kappa=`echo "scale=5; (($AG + $CT) / 2) /  (($AC + $AT + $CG + $GT) / 4)" | bc`
           
        # Execute ClonalFrameML
        ClonalFrameML mpileup2cns_nondupli_tran.phy_phyml_tree.txt mpileup2cns_tran.fasta mpileup2cns_tran -kappa $kappa -emsim 100
        # Exclude homologous recombination regions
        python $TOOLSPATH/remove_rec_from_ClonalFrameML_output_v1.1.py -a mpileup2cns_tran.fasta -o mpileup2cns_tran.importation_status.txt
        # Format sequences with only SNPs
        grep -v ">" mpileup2cns_tran.recRemoved.fasta | perl -pe 's/\n/_/g' | fold -1 | perl -pe 's/\n/\t/g' | perl -pe 's/\t_\t/\n/g' | perl $TOOLSPATH/transpose1_v2.pl  | sed '$d' > mpileup2cns_tran.recRemoved.tran
            
        # Create common file (File transformation for common creation, comparison expressions depend on number of strains)
        compare_strains_aftre_clonalframeml=`grep ">" name_list.txt | wc -l`   # Pick up column descriptions based on the number of strains to compare
        conpare_strains_object_clonalframeml=`cat $TOOLSPATH/core_genome_SNPs_strain_number_database_after_clonalframeml.txt | awk -F "\t" -v number="${compare_strains}" '$1==number {print $3}'`
        eval ${conpare_strains_object_clonalframeml}
        diff mpileup2cns_tran.recRemoved.tran mpileup2cns_tran.recRemoved.tran.common | grep "<" | cut -d " " -f 2-100 > mpileup2cns_tran.recRemoved.tran.common.snp
        perl $TOOLSPATH/transpose1_v2.pl mpileup2cns_tran.recRemoved.tran.common.snp | perl -pe 's/\t//g' > mpileup2cns_tran.recRemoved.tran.common.snp.prefas
        paste name_list.txt mpileup2cns_tran.recRemoved.tran.common.snp.prefas | perl -pe 's/\t/\n/g' > mpileup2cns_tran.recRemoved.tran.common.snp.fasta
        rm mpileup2cns_tran.recRemoved.tran mpileup2cns_tran.recRemoved.tran.common mpileup2cns_tran.recRemoved.tran.common.snp mpileup2cns_tran.recRemoved.tran.common.snp.prefas
        # Molecular phylogenetic analysis (RAxML)
        rm -f *RAxML*
        raxmlHPC-PTHREADS -f a -s mpileup2cns_tran.recRemoved.tran.common.snp.fasta -n mpileup2cns_tran.recRemoved.tran.common.snp_RAxML -m GTRGAMMA -T 4 -x 12345 -p12345 -# 1000
        # Count SNPs
        snp-dists mpileup2cns_tran.recRemoved.tran.common.snp.fasta > mpileup2cns_tran.recRemoved.tran.common.snp_distance.tsv
    fi

fi