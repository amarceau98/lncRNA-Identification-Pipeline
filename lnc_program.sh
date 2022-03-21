#!/bin/bash
#make directory to keep data in 

start=$(pwd)
echo Give your project a name:
read project_path
mkdir $project_path
cd $project_path

mkdir lncRNA_Identification
mkdir softwares

cd lncRNA_Identification

mkdir hisat2_mapping1
mkdir hisat2_mapping2

cd ..

#Get info 
echo Please format answer as path name to folder containing fasta files for each chromosome, found at hgdownload.soe.ucsc.edu/downloads.html
echo Where are chromosome files located?
read -e index

echo Please format answer as path name
echo Where is your data located?
read -e fastq_dir

echo How many samples are present?
read -e file_number

echo Please format answer as path name
echo Where is the reference genome, .gtf? Recommended: NCBI refGene
read -e reference

echo Please format answer as path name
echo Where is the a second reference genome, .gtf? Recommended: Ensembl 
read -e reference_2

echo Please format answer as path name
echo Where is the reference genome, .fa?
read -e reference_fa

echo Please format answer as the full path name, file included
echo Where is the RefSeq RNA database, .faa.gz? Can be found at https://ftp.ncbi.nih.gov/refseq
read -e hmmer_RefSeq

echo Please format answer as the full path name, file included
echo Where is the RefSeq protein database, .fna.gz? Can be found at https://ftp.ncbi.nih.gov/refseq
read -e blast_RefSeq


#Build index
gunzip $index/*.gz

cd softwares
wget --content-disposition https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
unzip -q *.*
chmod 777 hisat2-2.2.1/*

cd ..
cp $index/*.*  $start/$project_path
list=$(echo $start/$project_path/*.fa | tr ' ' ',')

$start/$project_path/softwares/hisat2-2.2.1/hisat2-build -f $list index

mkdir $start/$project_path/lncRNA_Identification/Index
mv index.* $start/$project_path/lncRNA_Identification/Index
rm *.fa

#Align twice, first to generate junction file, second to implement for highest possible accuracy

cd $fastq_dir

for i in {1..$file_number}
do 
	echo Are sequencing files paired, yes-answer 1, no-answer 2?
	read -e paired_split

	echo What format are the input files, fastq-answer 1, qseq-answer 2, or fasta files-answer 3?
	read -e input_format


	if [ $paired_split == 1 ] && [ $input_format == 1 ]
		then
		echo Input first file of pair:
		read -e paired_names_1
		echo Input second file of pair:
		read -e paired_names_2
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
		
	else 
		echo Loading...
	fi

	if [ $paired_split == 1 ] && [ $input_format == 2 ]
	then
		echo Input first file of pair:
		read -e paired_names_1
		echo Input second file of pair:
		read -e paired_names_2
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
	else 
		echo Loading...
	fi

	if [ $paired_split == 1 ] && [ $input_format == 3 ]
	then
		echo Input first file of pair:
		read -e paired_names_1
		echo Input second file of pair:
		read -e paired_names_2
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
	else 
		echo Loading...
	fi

	if [ $paired_split == 2 ] && [ $input_format == 1 ]
	then
		echo Input files name:
		read -e unpaired_names
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
	
	else 
		echo Loading...
	fi

	if [ $paired_split == 2 ] && [ $input_format == 2 ]
	then
		echo Input files name:
		read -e unpaired_names
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
	
	else 
		echo Loading...
	fi

	if [ $paired_split == 2 ] && [ $input_format == 3 ]
	then
		echo Input files name:
		read -e unpaired_names
		echo What wold you like to name the output file?
		read -e output_name
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
	else 
		echo Unable to align files.
	fi
done

#Convert .sam to .bam

#wget -P $project_path/softwares https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
#wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
#tar xvjf $project_path/softwares/samtools-1.14.tar.bz2 
#cd samtools-1.14
#make
#alias samtools_location = "$project_path/softwares/samtools-1.14"

#wget -P $project_path/softwares https://sourceforge.net/projects/bowtie-bio/files/latest/download
#gunzip -k $project_path/softwares/bowtie-1.3.1-src.zip
#alias bowtie_location = "$project_path/softwares/bowtie-1.3.1-src"

#echo How many samples are present?
#read -e file_number_2

#cd lncRNA_Identification/hisat2_mapping2

for numbers in $file_number_2
do 
	echo What is the file name, be sure to use the same name as used for mapping?
	read -e file_name
	
	samtools view -S -b $file_name.sam > $file_name.bam
	samtools sort $file_name.bam > $file_name.sorted.bam
done

#Stringtie
echo Please format answer as path name
echo Where is the Stringtie software located?
read -e Stringtie_location

cd lncRNA_Identification
mkdir Stringtie
cd hisat2_mapping2

echo How many samples are present?
read -e file_number_2

for numbers in $file_number_2
do 
	echo What is the file name, be sure to use the same name as used for mapping?
	read -e file_name
	
	
	$stringtie_location/stringtie lncRNA_Identification/hisat2_mapping2/$file_name.sam -o lncRNA_Identification/Stringtie/$file_name.gtf -G $reference -A $file_name.tab
done

#Determine consensus sequence
cd $tophat_maping_dir

echo How many conditions are present?
read -e Number_Of_Conditions

for numbers in $Number_Of_Conditions
do 
	echo Condition name?
	read -e Name
	
	echo List samples in this condition, separated by a space, including .gtf.
	read -e $sample_list

	
	$stringtie_location/stringtie --merge $sample_list -o lncRNA_Identification/Stringtie/$Name -G $reference 

done


#Isolate intergenic--CuffCompare
wget -P $project_path/softwares http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.tar.gz
gunzip -k $project_path/softwares/cufflinks-2.2.1.tar.gz
tar -tvf $project_path/softwares/cufflinks-2.2.1.tar
alias Cuffcompare_location = "$project_path/softwares/cufflinks-2.2.1"

cd lncRNA_Identification/Stringtie
mkdir lncRNA_Identification/CuffCompare

echo Where has this program been downloaded to?
read -e $software_location
chmod 777 $software_location/get_intergenicLoci_list.pl
chmod 777 $software_location/get_intergenicGTF.pl
chmod 777 $software_location/process_fa.pl
chmod 777 $software_location/summary_gtf.pl
chmod 777 $software_location/sixFrameTranslate.pl
	

#Get softwares
cd  $start/$project_path/softwares
wget --content-disposition https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-src.zip
wget --content-disposition http://eddylab.org/software/hmmer/hmmer.tar.gz
unzip -q \*.zip
gunzip -q *.gz
tar -xf *.tar
blast=$($start/$project_path/softwares/ncbi*)
hmmer=$($start/$project_path/softwares/hmmer*)


cd $start/$project_path/lncRNA_Identification
mkdir BLAST
mkdir Pfam

wget -P $project_path/softwares http://biopython.org/DIST/biopython-1.79.zip
gunzip -k $project_path/softwares/biopython-1.79.zip
alias biopython_location = "$project_path/softwares/biopython-1.79"

wget -P $project_path/softwares http://cpc2.gao-lab.org/data/CPC2-beta.tar.gz
gunzip -k $project_path/softwares/CPC2-beta.tar.gz
tar -tvf $project_path/softwares/CPC2-beta.tar
alias CPC_location = "$project_path/softwares/CPC2-beta"

mkdir $project_path/lncRNA_Identification/CPC 
mkdir $project_path/lncRNA_Identification/Pfam
	
wget -P $project_path/softwares http://eddylab.org/software/hmmer/hmmer.tar.gz
gunzip -k $project_path/softwares/hmmer.tar.gz
tar -tvf $project_path/softwares/hmmer.tar
alias hmmer_location = "$project_path/softwares/hmmer"



#$blast/makeblastdb -in $blast_RefSeq -dbtype nucl -o $start/$project_path/lncRNA_Identification/BLAST
#$hmmer/hmmbuild -o $project_path/lncRNA_Identification/Pfam -n hmmbuild --amino $hmmer_RefSeq



for numbers in $Number_Of_Conditions
do
	echo Condition name, be sure it matches to those above?
	read -e Name_2
	
	$Cuffcompare_location/cuffcompare -s $reference_fa -r $reference -o lncRNA_Identification/Stringtie/$Name_2.ref.gtf lncRNA_Identification/Stringtie/$Name.gtf
	
	$Cuffcompare_location/cuffcompare -s $reference_fa -r $reference_2 -o lncRNA_Identification/Stringtie/$Name_2.ens.gtf lncRNA_Identification/Stringtie/$Name.gtf

	
	$software_location/get_intergenicLoci_list.pl
	$software_location/get_intergenicGTF.pl
	
	#Make reports

	head -n1 $Name_2.summary.intergenic_loci.txt > $Name_2.summary.intergenic_loci.sorted.txt
	tail -n+2 $Name_2.summary.intergenic_loci.sorted.txt | sort -k3,3V -k5,5n >> $Name_2.summary.intergenic_loci.sorted.txt


	#Remove annotated genes and loci
	#May be stringtie??
	$CuffCompare_location/gffread -e -w $lncRNA_Identification/CuffCompare/$Name_2.intergenic_loci.step1.fa -g $lncRNA_Identification/CuffCompare/$Name_2.intergenic_loci.gtf

	$software_location/process_fa.pl
	$software_location/summary_gtf.pl
	
	#Remove those with coding potential
	cd $project_path/lncRNA_Identification/CPC 

	python $CPC_location/CPC2.py -i $Name_2.intergenic_loci_2.fa -o $Name_2.CPC.txt
	#Need to add reverse compliment

	awk -F, '$7 == "noncoding\r"' $Name_2.CPC.txt > $Name_2.CPC_removed.txt

	#Remove size
	awk -F, '$2 >= "200\r"' $Name_2.CPC_removed.txt > $Name_2.size_filtered.txt

	#Remove BLAST hits
	cd $project_path/lncRNA_Identification/BLAST
	$blast_location/blastn –db $project_path/lncRNA_Identification/BLAST –query (((((make fasta from txt file)  –out $project_path/lncRNA_Identification/BLAST/$Name_2.BLAST.out -outfmt 6

	#Remove those with hits
	awk -F, '$11 > "10\r"' $Name_2.BLAST.out > $Name_2.BLAST.removed.out
	
	#Turn back to .fa
	cd $project_path/lncRNA_Identification/Pfam
	
	#Convert to protein sequence
	$software_location/sixFrameTranslate.pl
	
	#Remove Pfam hits
	$hmmer_location/hmmsearch -tformat $hmmer_RefSeq lncRNA_Identification/Pfam/$Name_2.BLAST.removed.ff -o lncRNA_Identification/Pfam/$Name_2.Pfam_table.txt
	#Change this: awk -F, '$11 > "10\r"' $Name_2.Pfam_table.txt > $Name_2.Pfam_removed.txt
	
	#Final list output:
	
done
