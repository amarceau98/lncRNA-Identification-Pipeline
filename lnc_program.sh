#!/bin/bash
#make directory to keep data in 

start=$(pwd)
echo Give your project a name:
read project_path
#mkdir $project_path
#cd $project_path

#mkdir lncRNA_Identification
#mkdir softwares

#cd $start/$project_path/lncRNA_Identification
#mkdir hisat2_mapping1
#mkdir hisat2_mapping2
#mkdir CPC 
#mkdir BLAST
#cd $start/$project_path/lncRNA_Identification/BLAST
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
#cd ..
#mkdir Pfam


#cd ..

#Get info 
#echo Please format answer as path name to folder containing fasta files for each chromosome, found at hgdownload.soe.ucsc.edu/downloads.html
#echo Where are chromosome files located?
#read -e index

#echo Please format answer as path name
#echo Where is your data located?
#read -e fastq_dir

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

echo Please format answer as path name
echo Where is the fasta file for un-scaffolded transcripts,in a .fa format?
read -e unscaffolded

#Build index
#gunzip $index/*.gz

#cd softwares
#wget -q --content-disposition https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
#unzip -q *.*
#chmod 777 hisat2-2.2.1/*

#cd ..
#cp $index/*.*  $start/$project_path
#list=$(echo $start/$project_path/*.fa | tr ' ' ',')

#$start/$project_path/softwares/hisat2-2.2.1/hisat2-build -f $list index

#mkdir $start/$project_path/lncRNA_Identification/Index
#mv index.* $start/$project_path/lncRNA_Identification/Index
#rm *.fa

#Align twice, first to generate junction file, second to implement for highest possible accuracy

#cd $fastq_dir

#for i in {1..$file_number}
#do 
#	echo Are sequencing files paired, yes-answer 1, no-answer 2?
#	read -e paired_split
#
#	echo What format are the input files, fastq-answer 1, qseq-answer 2, or fasta files-answer 3?
#	read -e input_format
#
#
#	if [ $paired_split == 1 ] && [ $input_format == 1 ]
#		then
#		echo Input first file of pair:
#		read -e paired_names_1
#		echo Input second file of pair:
#		read -e paired_names_2
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#		
#	else 
#		echo Loading...
#	fi
#
#	if [ $paired_split == 1 ] && [ $input_format == 2 ]
#	then
#		echo Input first file of pair:
#		read -e paired_names_1
#		echo Input second file of pair:
#		read -e paired_names_2
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#	else 
#		echo Loading...
#	fi
#
#	if [ $paired_split == 1 ] && [ $input_format == 3 ]
#	then
#		echo Input first file of pair:
#		read -e paired_names_1
#		echo Input second file of pair:
#		read -e paired_names_2
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -1 $fastq_dir/$paired_names_1 -2 $fastq_dir/$paired_names_2 -S $output_name.hisat2-1.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#	else 
#		echo Loading...
#	fi
#
#	if [ $paired_split == 2 ] && [ $input_format == 1 ]
#	then
#		echo Input files name:
#		read -e unpaired_names
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -q -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1		
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#	
#	else 
#		echo Loading...
#	fi
#
#	if [ $paired_split == 2 ] && [ $input_format == 2 ]
#	then
#		echo Input files name:
#		read -e unpaired_names
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -qseq -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#	
#	else 
#		echo Loading...
#	fi
#
#	if [ $paired_split == 2 ] && [ $input_format == 3 ]
#	then
#		echo Input files name:
#		read -e unpaired_names
#		echo What wold you like to name the output file?
#		read -e output_name
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-1.sam --novel-splicesite-outfile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-1.summary.txt
#		$start/$project_path/softwares/hisat2-2.2.1/hisat2 -f -x $start/$project_path/lncRNA_Identification/Index/index -U $fastq_dir/$unpaired_names -S $output_name.hisat2-2.sam --novel-splicesite-infile $output_name.hisat2-1.junctions --summary-file $output_name.hisat2-2.summary.txt
#		mv $output_name.hisat2-1* $start/$project_path/lncRNA_Identification/hisat2_mapping1
#		mv $output_name.hisat2-2* $start/$project_path/lncRNA_Identification/hisat2_mapping2
#	else 
#		echo Unable to align files.
#	fi
#done

#Convert .sam to .bam
#cd $start/$project_path/softwares
#wget -q https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
#bzip2 -d samtools*
#tar -xf $start/$project_path/softwares/samtools*
#rm $start/$project_path/softwares/*.tar
#cd samtools-1.15
#./configure --disable-lzma
#make

#cd ..
#wget -q --content-disposition https://sourceforge.net/projects/bowtie-bio/files/latest/download
#unzip bowtie2*.zip

#cd $start/$project_path/softwares
#cd $start/$project_path/softwares/samtools*
#chmod 777 *
#cd hts*
#chmod 777 *
#cd .. 

#cd $start/$project_path/softwares
 

#echo Running Samtools

#for files in $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.sam
#do
#$start/$project_path/softwares/samtools*/samtools view -S -b $files > $files.bam
#$start/$project_path/softwares/samtools*/samtools sort $files.bam > $files.sorted.bam
#echo Done $files
#done

#echo Running Stringtie


#Stringtie
#cd $start/$project_path/softwares
#wget -q --content-disposition http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz
#gunzip stringtie*
#tar -xf stringtie*
#rm stringtie*.tar
#cd stringtie*
#chmod 777 *

#cd $start/$project_path/lncRNA_Identification
#mkdir Stringtie
#cd Stringtie

#for files in $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.sorted.bam
#do
#	$start/$project_path/softwares/stringtie*/stringtie $files -o $files.ref1.gtf -G $reference -A $files.ref1.tab
#	mv -n $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.gtf $start/$project_path/lncRNA_Identification/Stringtie/
#	mv -n $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.tab $start/$project_path/lncRNA_Identification/Stringtie/
#done

#for files in $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.sorted.bam
#do
#	$start/$project_path/softwares/stringtie*/stringtie $files -o $files.ref2.gtf -G $reference_2 -A $files.ref2.tab
#	mv -n $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.gtf $start/$project_path/lncRNA_Identification/Stringtie/
#	mv -n $start/$project_path/lncRNA_Identification/hisat2_mapping2/*.tab $start/$project_path/lncRNA_Identification/Stringtie/
#done


#Isolate intergenic--CuffCompare
#cd $start/$project_path/softwares
#wget -q --content-disposition http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
#gunzip cufflinks*.gz
#tar -xf cufflinks*.tar
#rm cufflinks*.tar
#cd cufflinks*
#chmod 777 *


#cd $start/$project_path/softwares
#wget -q --content-disposition http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz
#gunzip gffread*.gz
#tar -xf gffread*.tar
#rm gffread*.tar

#cd $start/$project_path/softwares
#wget -q --content-disposition https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
#mv $start/$project_path/softwares/bedtools.static.binary bedtools
#chmod a+x $start/$project_path/softwares/bedtools



#cd $start/$project_path/lncRNA_Identification
#mkdir CuffCompare

#Pull dependencies from GitHub
#cd $start/$project_path/softwares
#wget -q --content-disposition https://github.com/amarceau98/lncRNA-Identification-Pipeline/archive/refs/heads/main.zip
#unzip lncRNA*
#rm lncRNA*.zip
#cd lncRNA*
#chmod 777 *.pl	



#cd $start/$project_path/softwares
#wget -q --content-disposition http://biopython.org/DIST/biopython-1.79.tar.gz
#wget -q --content-disposition https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
#gunzip CPC*.gz
#gunzip biopython*  
#tar -xf CPC*.tar
#tar -xf biopython*
#rm CPC*.tar
#rm biopython*.tar
#cd $start/$project_path/softwares/CPC*
#cd libs/libsvm/
#gzip -dc libsvm-*.tar.gz | tar xf -
#cd libsvm-*/
#make clean && make

#cd $start/$project_path/softwares
#wget -q --content-disposition https://github.com/lh3/seqtk/archive/refs/heads/master.zip
#unzip seqtk-master.zip
#rm seqtk-master.zip
#cd $start/$project_path/softwares/seqtk-master
#make


#cd $start/$project_path/softwares
#wget -q --content-disposition https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
#gunzip ncbi*  
#tar -xf ncbi*.tar
#rm ncbi*.tar

#make database 
#cd $start/$project_path/lncRNA_Identification/BLAST/
#wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
#gunzip nt.gz
#$start/$project_path/softwares/ncbi*/bin/makeblastdb -in nt -dbtype nucl -parse_seqids 

#cd

#Get hmmer
#cd $start/$project_path/softwares
#wget -q --content-disposition http://eddylab.org/software/hmmer/hmmer.tar.gz
#gunzip hmmer*
#tar -xf hmmer*
#rm hmmer*.tar
#cd hmmer*
#./configure
#make

#cd $start/$project_path/lncRNA_Identification/Pfam
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
#gunzip Pfam-A.hmm.gz
#$start/$project_path/softwares/hmmer*/src/hmmpress $start/$project_path/lncRNA_Identification/Pfam/Pfam-A.hmm


#sed -i 's/*hr/chr/g' $reference
#sed -i 's/*hr/chr/g' $reference_2
#sed -i 's/*hr/chr/g' $reference_fa

#for for i in {1..$file_number}
#do
#echo Sample name?
#read output_name
#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref1.gtf 
#do
#grep 'chrUn*' $files > $start/$project_path/lncRNA_Identification/Stringtie/$output_name.ref1.Unknown.gtf
#grep -v 'chrUn*' $files > $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref1.Known.gtf
#tail -n +2 $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref1.Known.gtf > $start/$project_path/lncRNA_Identification/Stringtie/$output_name.ref1.Known.gtf
#rm $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref1.Known.gtf
#done

#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref2.gtf 
#do
#grep 'chrUn*' $files > $start/$project_path/lncRNA_Identification/Stringtie/$output_name.ref2.Unknown.gtf
#grep -v 'chrUn*' $files > $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref2.Known.gtf
#tail -n +2 $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref2.Known.gtf > $start/$project_path/lncRNA_Identification/Stringtie/$output_name.ref2.Known.gtf
#rm $start/$project_path/lncRNA_Identification/Stringtie/1.$output_name.ref2.Known.gtf
#done

#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref1.Known.gtf
#do
# 	$start/$project_path/softwares/cufflink*/cuffcompare $files -s $reference_fa -r $reference -o $files.known.cc.ref1
#done

#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref1.Unknown.gtf
#do
#	$start/$project_path/softwares/cufflink*/cuffcompare $files -s $unscaffolded -r $reference -o $files.unknown.cc.ref1
#done

#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref2.Known.gtf
#do
#	$start/$project_path/softwares/cufflink*/cuffcompare $files -s $reference_fa -r $reference_2 -o $files.known.cc.ref2
#done

#for files in $start/$project_path/lncRNA_Identification/Stringtie/*.ref2.Unknown.gtf
#do
# 	$start/$project_path/softwares/cufflink*/cuffcompare $files -s $unscaffolded -r $reference_2 -o $files.unknown.cc.ref2
#done

#cd $start/$project_path/lncRNA_Identification/Stringtie
#mkdir Keep
#mv *.sorted.* Keep/
#mv *Unknown.gtf Keep/
#mv *Known.gtf Keep/
#mv *.* $start/$project_path/lncRNA_Identification/CuffCompare
#cd Keep/
#mv *.* $start/$project_path/lncRNA_Identification/Stringtie/
#cd ..
#rm -r Keep

#cat $start/$project_path/lncRNA_Identification/CuffCompare/*.known.cc.ref1.combined.gtf $start/$project_path/lncRNA_Identification/CuffCompare/*.unknown.cc.ref1.combined.gtf > $start/$project_path/lncRNA_Identification/CuffCompare/$output_name.ref1.total.gtf
#cat $start/$project_path/lncRNA_Identification/CuffCompare/*.known.cc.ref2.combined.gtf $start/$project_path/lncRNA_Identification/CuffCompare/*.unknown.cc.ref2.combined.gtf > $start/$project_path/lncRNA_Identification/CuffCompare/$output_name.ref2.total.gtf
#done 

#Determine consensus sequence

echo How many conditions are present?
read Number_Of_Conditions

chmod 777 $start/$project_path/lncRNA_Identification/Stringtie/*
mkdir $start/$project_path/lncRNA_Identification/Filtering

for i in {1..$Number_Of_Conditions}
do 
	echo Condition name?
	read Name
	
	echo How many samples within this condition? 
	read samples_in_condition
for i in {1..$samples_in_condition}
do
	echo Sample name?
	read Sample

#	ref1="$start/$project_path/lncRNA_Identification/CuffCompare/$Sample.ref1.total.gtf"
#	ref1_1=`echo $ref1`	
#	ref2="$start/$project_path/lncRNA_Identification/CuffCompare/$Sample.ref2.total.gtf"
#	ref2_1=`echo $ref2`
#	output1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.list_intergenic.txt"
#	output1_1=`echo $output1`
	
#	cd $start/$project_path/softwares/lnc*/
#	cp get_intergenicLoci_list.pl get_intergenicLoci_list_1.pl	

#	sed -i "s,loci_replace_ref1,${ref1_1},g" $start/$project_path/softwares/lnc*/get_intergenicLoci_list_1.pl
#	sed -i "s,loci_replace_ref2,${ref2_1},g" $start/$project_path/softwares/lnc*/get_intergenicLoci_list_1.pl
#	sed -i "s,list_replace,${output1_1},g" $start/$project_path/softwares/lnc*/get_intergenicLoci_list_1.pl
	
#	perl $start/$project_path/softwares/lnc*/get_intergenicLoci_list_1.pl

#	rm $start/$project_path/softwares/lnc*/get_intergenicLoci_list_1.pl

#	loci1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.list_intergenic.txt"
#	loci1_1=`echo $loci1`	
#	gtf1="$start/$project_path/lncRNA_Identification/CuffCompare/$Sample.ref1.total.gtf"
#	gtf1_1=`echo $gtf1`
#	out1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.gtf"
#	out1_1=`echo $out1`

#	loci1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.list_intergenic.txt"
#	loci1_1=`echo $loci1`	
#	gtf2="$start/$project_path/lncRNA_Identification/CuffCompare/$Sample.ref2.total.gtf"
#	gtf2_1=`echo $gtf2`
#	out2="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.gtf"
#	out2_1=`echo $out2`

#	cd $start/$project_path/softwares/lnc*/
#	cp get_intergenicGTF.pl get_intergenicGTF_1.pl	
#	cp get_intergenicGTF.pl get_intergenicGTF_2.pl	

#	sed -i "s,loci_replace,${loci1_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_1.pl	
#	sed -i "s,gtf_replace,${gtf1_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_1.pl
#	sed -i "s,out_replace,${out1_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_1.pl	
#	sed -i "s,loci_replace,${loci1_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_2.pl	
#	sed -i "s,gtf_replace,${gtf2_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_2.pl	
#	sed -i "s,out_replace,${out2_1},g" $start/$project_path/softwares/lnc*/get_intergenicGTF_2.pl


#	perl $start/$project_path/softwares/lnc*/get_intergenicGTF_1.pl
#	perl $start/$project_path/softwares/lnc*/get_intergenicGTF_2.pl

#	rm $start/$project_path/softwares/lnc*/get_intergenicGTF_1.pl 
#	rm $start/$project_path/softwares/lnc*/get_intergenicGTF_2.pl 


#	grep 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.gtf > $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.ChrUn.gtf	
#	grep -v 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.gtf > $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref1.intergenic.ChrKnown.gtf
#	tail -n +2 $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref1.intergenic.ChrKnown.gtf > $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.ChrKnown.gtf
#	rm $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref1.intergenic.ChrKnown.gtf

#	grep 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.gtf > $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.ChrUn.gtf	
#	grep -v 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.gtf > $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref2.intergenic.ChrKnown.gtf
#	tail -n +2 $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref2.intergenic.ChrKnown.gtf > $start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.ChrKnown.gtf
#	rm $start/$project_path/lncRNA_Identification/Filtering/1.$Sample.ref2.intergenic.ChrKnown.gtf

	
#	gtf1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.ChrUn.gtf"
#	gtf1_1=`echo $gtf1`	
#	out1="$start/$project_path/lncRNA_Identification/Filtering/$Sample.summary.ref1.ChrUn.txt"
#	out1_1=`echo $out1`

	
#	gtf2="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref1.intergenic.ChrKnown.gtf"
#	gtf2_1=`echo $gtf2`	
#	out2="$start/$project_path/lncRNA_Identification/Filtering/$Sample.summary.ref1.ChrKnown.txt"
#	out2_1=`echo $out2`

#	gtf3="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.ChrUn.gtf"
#	gtf3_1=`echo $gtf3`	
#	out3="$start/$project_path/lncRNA_Identification/Filtering/$Sample.summary.ref2.ChrUn.txt"
#	out3_1=`echo $out3`

#	gtf4="$start/$project_path/lncRNA_Identification/Filtering/$Sample.ref2.intergenic.ChrKnown.gtf"
#	gtf4_1=`echo $gtf4`	
#	out4="$start/$project_path/lncRNA_Identification/Filtering/$Sample.summary.ref2.ChrKnown.txt"
#	out4_1=`echo $out4`

#	perl $start/$project_path/softwares/lnc*/summary_gtf.pl -i $gtf1_1 -o $out1_1
#	perl $start/$project_path/softwares/lnc*/summary_gtf.pl -i $gtf2_1 -o $out2_1
#	perl $start/$project_path/softwares/lnc*/summary_gtf.pl -i $gtf3_1 -o $out3_1
#	perl $start/$project_path/softwares/lnc*/summary_gtf.pl -i $gtf4_1 -o $out4_1
done

#	sort -u $start/$project_path/lncRNA_Identification/Filtering/*.ref*.txt >> $start/$project_path/lncRNA_Identification/Filtering/$Name.summary.txt
#	sed -e '1i\TCONS\tXLOC\tchr\tstrand\tstart\tend\tnum_exons\tlength\tstarts\tends' $start/$project_path/lncRNA_Identification/Filtering/$Name.summary.txt >> $start/$project_path/lncRNA_Identification/Filtering/$Name.overall.summary.txt
#	rm $start/$project_path/lncRNA_Identification/Filtering/$Name.summary.txt
#	tail -n+2 $start/$project_path/lncRNA_Identification/Filtering/$Name.overall.summary.txt | sort -k3,3V -k5,5n >> $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic_loci.sorted.txt

#	awk ' (NR==1) || ($8 > 199 ) ' $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic_loci.sorted.txt > $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic_loci.size.txt	
#	awk '{ print $3 "\t" $5 "\t" $6 "\t" $2}' $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic_loci.size.txt > $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic.bed

	
#	grep 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic.bed > $start/$project_path/lncRNA_Identification/Filtering/$Name.UnknownChr.bed
#	grep -v 'chrUn*' $start/$project_path/lncRNA_Identification/Filtering/$Name.intergenic.bed > $start/$project_path/lncRNA_Identification/Filtering/1.$Name.KnownChr.bed
#	tail -n +2 $start/$project_path/lncRNA_Identification/Filtering/1.$Name.KnownChr.bed > $start/$project_path/lncRNA_Identification/Filtering/$Name.KnownChr.bed
#	rm $start/$project_path/lncRNA_Identification/Filtering/1.$Name.KnownChr.bed

#	$start/$project_path/softwares/bedtools getfasta -fi $reference_fa -bed $start/$project_path/lncRNA_Identification/Filtering/$Name.KnownChr.bed -fo $start/$project_path/lncRNA_Identification/Filtering/$Name.KnownChr.fa
#	$start/$project_path/softwares/bedtools getfasta -fi $unscaffolded -bed $start/$project_path/lncRNA_Identification/Filtering/$Name.UnknownChr.bed -fo $start/$project_path/lncRNA_Identification/Filtering/$Name.UnknownChr.fa
#	cat $start/$project_path/lncRNA_Identification/Filtering/$Name.KnownChr.fa $start/$project_path/lncRNA_Identification/Filtering/$Name.UnknownChr.fa > $start/$project_path/lncRNA_Identification/CPC/$Name.Intergenic.fa
#	$start/$project_path/softwares/seqtk-master/seqtk seq -r $start/$project_path/lncRNA_Identification/CPC/$Name.Intergenic.fa > $start/$project_path/lncRNA_Identification/CPC/$Name.Intergenic.RC.fa

#	$start/$project_path/softwares/CPC*/bin/CPC2.py -i $start/$project_path/lncRNA_Identification/CPC/$Name.Intergenic.fa -o $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Out
#	$start/$project_path/softwares/CPC*/bin/CPC2.py -i $start/$project_path/lncRNA_Identification/CPC/$Name.Intergenic.RC.fa -o $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.RC.Out

#	awk ' $8 == "noncoding"' $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Out.txt > $start/$project_path/lncRNA_Identification/CPC/$Name.Noncoding.txt
#	awk ' $8 == "noncoding"' $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Out.txt > $start/$project_path/lncRNA_Identification/CPC/$Name.Noncoding.RC.txt

#	awk '{ print $1}' $start/$project_path/lncRNA_Identification/CPC/$Name.Noncoding.txt > $start/$project_path/lncRNA_Identification/CPC/comb.NC.txt
#	awk '{ print $1}' $start/$project_path/lncRNA_Identification/CPC/$Name.Noncoding.RC.txt > $start/$project_path/lncRNA_Identification/CPC/comb.NC.RC.txt
#	sort -u $start/$project_path/lncRNA_Identification/CPC/comb.*.txt >> $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.txt
#	rm $start/$project_path/lncRNA_Identification/CPC/comb.*.txt

#	sed -i "s,:,\t,g" $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.txt
#	sed -i "s,-,\t,g" $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.txt
#	mv $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.txt $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.bed


#	grep 'chrUn*' $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.bed > $start/$project_path/lncRNA_Identification/Pfam/$Name.UnknownChr.bed
#	grep -v 'chrUn*' $start/$project_path/lncRNA_Identification/CPC/$Name.CPC.Final.bed > $start/$project_path/lncRNA_Identification/Pfam/1.$Name.KnownChr.bed
#	tail -n +2 $start/$project_path/lncRNA_Identification/Pfam/1.$Name.KnownChr.bed > $start/$project_path/lncRNA_Identification/Pfam/$Name.KnownChr.bed
#	rm $start/$project_path/lncRNA_Identification/Pfam/1.$Name.KnownChr.bed
#	$start/$project_path/softwares/bedtools getfasta -fi $reference_fa -bed $start/$project_path/lncRNA_Identification/Pfam/$Name.KnownChr.bed -fo $start/$project_path/lncRNA_Identification/Pfam/$Name.KnownChr.fa
#	$start/$project_path/softwares/bedtools getfasta -fi $unscaffolded -bed $start/$project_path/lncRNA_Identification/Pfam/$Name.UnknownChr.bed -fo $start/$project_path/lncRNA_Identification/Pfam/$Name.UnknownChr.fa
#	cat $start/$project_path/lncRNA_Identification/Pfam/$Name.KnownChr.fa $start/$project_path/lncRNA_Identification/Pfam/$Name.UnknownChr.fa > $start/$project_path/lncRNA_Identification/Pfam/$Name.CPC.fa
	

#	filein1="$start/$project_path/lncRNA_Identification/Pfam/$Name.CPC.fa"
#	filein1_1=`echo $filein1`

#	cd $start/$project_path/softwares/lnc*/
#	cp sixFrameTranslate.pl sixFrameTranslate_1.pl	

#	sed -i "s,file_in_replace,${filein1_1},g" $start/$project_path/softwares/lnc*/sixFrameTranslate_1.pl


#	cp $start/$project_path/softwares/lnc*/codon.txt $start/$project_path/lncRNA_Identification/Pfam
#	perl $start/$project_path/softwares/lnc*/sixFrameTranslate.pl -i $start/$project_path/lncRNA_Identification/Pfam/$Name.CPC.fa -o $start/$project_path/lncRNA_Identification/Pfam/$Name.Pfam.faa -l 6
#	awk 'BEGIN {RS = ">" ; ORS = ""} length($2) <= 100000 {print ">"$0}' $start/$project_path/lncRNA_Identification/Pfam/$Name.Pfam.faa > $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamReady.faa
#	split --additional-suffix=.faa $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamReady.faa $start/$project_path/lncRNA_Identification/Pfam/part


#	for files in $start/$project_path/lncRNA_Identification/Pfam/part*.faa
#	do
#	$start/$project_path/softwares/hmmer-3.3.2/src/hmmsearch --tblout $files.txt $start/$project_path/lncRNA_Identification/Pfam/Pfam-A.hmm $files
#	done

#	cat $start/$project_path/lncRNA_Identification/Pfam/part*.txt > $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamTable.txt	
#	rm $start/$project_path/lncRNA_Identification/Pfam/part*


	awk ' (NR==1) || ($5 > .1 ) ' $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamTable.txt > $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamFiltered.txt
	awk '{print $1 }' $start/$project_path/lncRNA_Identification/Pfam/$Name.PfamFiltered.txt > $start/$project_path/lncRNA_Identification/Pfam/list.txt
	grep -v "# *" $start/$project_path/lncRNA_Identification/Pfam/list.txt > $start/$project_path/lncRNA_Identification/Pfam/tmp.txt
	mv $start/$project_path/lncRNA_Identification/Pfam/tmp.txt $start/$project_path/lncRNA_Identification/Pfam/list.txt
	sed -i 's/|/\t/g' $start/$project_path/lncRNA_Identification/Pfam/list.txt
	awk '{print $1 }' $start/$project_path/lncRNA_Identification/Pfam/list.txt > $start/$project_path/lncRNA_Identification/Pfam/list2.txt
	sort -u $start/$project_path/lncRNA_Identification/Pfam/list2.txt > $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.txt
	rm $start/$project_path/lncRNA_Identification/Pfam/list* 

	echo Check Filtering
	read status_of_filtering

	sed -i 's/:/\t/g' $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.txt
	sed -i 's/-/\t/g' $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.txt
	cp $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.txt $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.bed
	grep 'chrUn*' $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.bed > $start/$project_path/lncRNA_Identification/BLAST/$Name.UnknownChr.bed
	grep -v 'chrUn*' $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.Final.List.bed > $start/$project_path/lncRNA_Identification/BLAST/1.$Name.KnownChr.bed
	tail -n +2 $start/$project_path/lncRNA_Identification/BLAST/1.$Name.KnownChr.bed > $start/$project_path/lncRNA_Identification/BLAST/$Name.KnownChr.bed
	rm $start/$project_path/lncRNA_Identification/BLAST/1.$Name.KnownChr.bed
	$start/$project_path/softwares/bedtools getfasta -fi $reference_fa -bed $start/$project_path/lncRNA_Identification/BLAST/$Name.KnownChr.bed -fo $start/$project_path/lncRNA_Identification/BLAST/$Name.KnownChr.fa
	$start/$project_path/softwares/bedtools getfasta -fi $unscaffolded -bed $start/$project_path/lncRNA_Identification/BLAST/$Name.UnknownChr.bed -fo $start/$project_path/lncRNA_Identification/BLAST/$Name.UnknownChr.fa
	cat $start/$project_path/lncRNA_Identification/BLAST/$Name.KnownChr.fa $start/$project_path/lncRNA_Identification/BLAST/$Name.UnknownChr.fa > $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.fa
	
	echo End here
	read BLAST_ready

	$start/$project_path/softwares/ncbi*/bin/blastn –db $start/$project_path/lncRNA_Identification/BLAST/nt -query $start/$project_path/lncRNA_Identification/BLAST/$Name.Pfam.fa -out $start/$project_path/lncRNA_Identification/BLAST/$Name.BLAST.txt
done


	##Remove those with hits
	#awk -F, '$11 > "10\r"' $Name_2.BLAST.out > $Name_2.BLAST.removed.out
	


