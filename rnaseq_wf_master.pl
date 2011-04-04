#!/usr/bin/perl -w

use strict;
use warnings;

#Master scripts for analysing rnaseq data from Illumina pipleine using tophat. The program performs QC in FASTQC, aligns reads to human genome using Tophat program and generates read counts using htseqCount.
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
rnaseq_wf_master.pl  -id [indir...n] -a [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata]

=head2 COMMANDS AND OPTIONS

-id/--infiledir Infile dir(s)

-a/--projectid The project ID

-s/--sampleid The sample ID(s)

-em/--email

-odf/--outdirdata The data files output directory (defaults to data)

-ods/--outdirscript The script files output directory (defaults to wf_scripts)

-pFQC/--fastqc Flag for running FASTQC (defaults to yes (=1))

-pTH/--tophat Flag for running Tophat (defaults to yes (=1))

-pSTV/--samtoolsview Flag for running samtools view (defaults to yes (=1))

-pHTSC/--htseqcount Flag for running htseqCount (defaults to yes (=1))

-pRPKM/--rpkmforgenes Flag for running rpkmForGenes (defaults to yes (=1))

-pCHTSC/--concat_htseqCounts Flag for running concat_htseqCounts_wf.1.0.pl (defaults to yes (=1))

-schtsc/--sampleid_concathtseqcounts sampleid for header columns in concatenated htseqCount data

-ochtsc/--outfile Output file (Defaults to concat_htseqCount_wf.txt)

-pDES/--deseq Flag for running DESeq (defaults to yes (=1))

-sidDES/--deseqsampleid DESeq sample id (defaults to deseq)

-nrcond1DES/--nrcondone DESeq number of conditions one (Must be supplied to run DESeq)

-nrcond2DES/--nrcondtwo DESeq number of conditions two (Must be supplied to run DESeq)

=head3 I/O

Input format ( dir/infiles.fastq )

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{rnaseq_wf_master.pl -id [indir...n] -a [projectid] -s [sampleid] -em [e-mail] -ods [outdirscripts] -odf [outdirdata]
	       -id/--infiledir Infile dir(s), comma sep
	       -a/--projectid The project ID
	       -s/--sampleid The sample ID(s),comma sep 
	       -em/--email e-mail
	       -odf/--outdirdata The data files output directory (defaults to data)
	       -ods/--outdirscript The script files output directory (defaults to wf_scripts)
	       -pFQC/--fastqc Flag for running FASTQC (defaults to yes (=1))
	       -pTH/--tophat Flag for running Tophat (defaults to yes (=1))
	       -pSTV/--samtoolsview Flag for running samtools view (defaults to yes (=1))
	       -pHTSC/--htseqcount Flag for running htseqCount (defaults to yes (=1))
	       -pRPKM/--rpkforgenes Flag for running rpkmForGenes (defaults to yes (=1))
	       -pCHTSC/--concat_htseqCounts Flag for running concat_htseqCounts_wf.1.0.pl (defaults to yes (=1))
	       -schtsc/--sampleid_concathtseqcounts Sampleid for header columns in concatenated htseqCount data
	       -ochtsc/--outfile Output file (Defaults to concat_htseqCount_wf.txt)
	       -pDES/--deseq Flag for running DESeq (defaults to yes (=1))
	       -sidDES/--deseqsampleid DESeq sample id (defaults to deseq)
	       -nrcond1DES/--nrcondone DESeq number of conditions one (Must be supplied to run DESeq)
	       -nrcond2DES/--nrcondtwo DESeq number of conditions two (Must be supplied to run DESeq)	  
	   };
    
}

my ($aid,$em,$odf,$ods,$fnend,$ochtsc, $sidDES, $nrcond1DES, $nrcond2DES, $schtsc, $filename, $filename2, $fnt, $fnt2,$help) = (0,0,"data", "wf_scripts", ".sh", "concat_htseqCount_wf.txt", "deseq"); #Arguments for project
my ($pFQC, $pTH, $pSTV, $pHTSC, $pRPKM, $pCHTSC, $pDES) = (1, 1, 1, 1, 1, 1, 1); #Arguments for running programs
my (@inid,@sid, @refchr, @sidhtseq); #infiledir, sampleid (dir) ,reference annotaion, sampleid (ind)
my (%infiles);

GetOptions('id|infiledir:s'  => \@inid, #Comma separeted list
	   'a|projectid:s'  => \$aid,
	   's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
	   'em|email:s'  => \$em,
	   'odf|outdirdata:s'  => \$odf, #One above sample id
	   'ods|outdirscript:s'  => \$ods,
	   'pFQC|fastqc:n' => \$pFQC,
	   'pTH|tophat:n' => \$pTH,
	   'pSTV|samtoolsview:n' => \$pSTV,
	   'pHTSC|htseqcount:n' => \$pHTSC,
	   'pCHTSC|concathtseqcount:n' => \$pCHTSC,
	   'schtsc|sampleidconcathtseqcounts:s'  => \@sidhtseq, #Comma separeted list, ind sample
	   'ochtsc|outfileconcathtseqcounts:s'  => \$ochtsc, #Concatenated outout file of htseqCounts data
	   'pRPKM|rpkmforgenes:n' => \$pRPKM,
	   'pDES|deseq:n' => \$pDES,
	   'sidDES|sampleiddes:s' => \$sidDES, #sample id for sample vs sample
	   'nrcond1DES|nrcondonedes:s' => \$nrcond1DES, #sample id for sample vs sample
	   'nrcond2DES|nrcondtwodes:s' => \$nrcond2DES, #sample id for sample vs sample
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if (@inid == 0) {
   my $verbosity = 2;
 print"\n";
 pod2usage({-message => "Must supply an infile directory as comma separeted list.\n",
     -verbose => $verbosity
   });
}
if ($aid eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a project ID", "\n\n";
    die $USAGE;
}
if ( scalar(@sid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
    die $USAGE;
}
if ($pDES eq 1 && $nrcond1DES eq "" || $nrcond2DES eq "") { #if run Deseq nr cond must be specified
    
    print STDERR "\n";
    print STDERR "Must supply nr of conditions to contrast in DESeq", "\n\n";
    die $USAGE;
}

@inid = split(/,/,join(',',@inid)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs
@sidhtseq = split(/,/,join(',',@sidhtseq)); #Enables comma separeted list of sample IDs for concat_htseqCounts_wf.1.0.pl


for (my $inputdir=0;$inputdir<scalar(@inid);$inputdir++) { #Collects inputfiles
    
    my @infiles = `cd $inid[ $inputdir ]/fastq;ls *.fastq;`; #cd to input dir and collect .fastq files
    
    print STDERR "\nSample ID", "\t", $sid[$inputdir],"\n";
    print STDERR "Inputfiles", "\n", @ { $infiles{ $sid[$inputdir] }  =[@infiles] }; #hash with sample id as key and inputfiles in dir as array 
 
}



#########################
###Run program part######
#########################

if ($pFQC eq 1) { #Run FASTQC

    print STDERR "\nFASTQC", "\n";
    print STDERR "Creating sbatch script FASTQC and writing script file(s) to: ", $ods,"/sampleid/fastqc/fastqc_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script FASTQC data files will be written to: ", $odf,"/sampleid/fastqc/", "\n";

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Fastqc($sid[$sampleid]);
	
    }

}

if ($pTH eq 1) { #Run Tophat

    print STDERR "\nTophat", "\n";
    print STDERR "Creating sbatch script Tophat and writing script file(s) to: ", $ods,"/sampleid/tophat/tophat_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script Tophat data files will be written to: ", $odf,"/sampleid/tophat/filename", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Tophat($sid[$sampleid]);
	
    }

}

if ($pSTV eq 1) { #Run Samtools view 
    
    print STDERR "\nSamtools view", "\n";
    print STDERR "Creating sbatch script samtools view and writing script file(s) to: ", $ods,"/sampleid/samtools/samtools_view_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script rpkmforgenes data files will be written to: ", $odf,"/sampleid/tophat/filename", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	
	Samtoolsview($sid[$sampleid]);
	
    }
}

if ($pTH eq 1) { #Run HTSeqCount 

    print STDERR "\nHTSeqCount", "\n";
    print STDERR "Creating sbatch script HTSeqCount and writing script file(s) to: ", $ods,"/sampleid/htseq/htseqCount_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script HTSeqCount data files will be written to: ", $odf,"/sampleid/htseq/filename", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Htseqcount($sid[$sampleid]);
	
    }

}

if ($pRPKM eq 1) { #Run Rpkmforgenes and cut_off.1.0.R 

    print STDERR "\nRpkmforgenes", "\n";
    print STDERR "Creating sbatch script rpkmforgenes and writing script file(s) to: ", $ods,"/sampleid/rpkmforgenes/rpkmForGenes_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script rpkmforgenes data files will be written to: ", $odf,"/sampleid/rpkmforgenes/filename", "\n";


    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Rpkmforgenes($sid[$sampleid]);
	
    }
}

if ($pCHTSC eq 1) { #Run concat_htseqCounts_wf.1.0.pl 

    print STDERR "\nConcat_htseqCounts_wf.1.0.pl", "\n";
    print STDERR "Creating sbatch script concat_htseqCounts and writing script file(s) to: ", $ods,"/concat_htseqc/concat_htseqc_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script concat_htseqCounts data files will be written to: ", $odf,"/concat_htseqc/filename", "\n";

    
	Concat_htseqCD();
}

if ($pDES eq 1) { #Run Deseq_script_1_0.R 

    print STDERR "\nDeseq_script_1_0.R", "\n";
    print STDERR "Creating sbatch script deseq_wf and writing script file(s) to: ", $ods,"/deseq/deseq_wf.sh", "\n";
    print STDERR "Sbatch script deseq_wf data files will be written to: ", $odf,"/deseq/filename", "\n";

    
	DESeq(); #NOTE: package DESeq and gplots are not installed on uppmax as yet
}

######################
###Sub Routines#######
######################

sub DESeq { 
    
#Generates a sbatch script and to run DESEQ
    
    `mkdir -p $odf/deseq/info;`; #Creates the htseq and info data file directory
    `mkdir -p $ods/deseq;`; #Creates the concat_htseqCounts_wf.1.0.pl script directory
    
    $filename = "/bubo/proj/$aid/private/$ods/deseq/deseq_wf.";
    Checkfnexists($filename, $fnend);

    open (DES, ">$filename") or die "Can't write to $filename: $!\n";
    
    print DES "#! /bin/bash -l", "\n";
    print DES "#SBATCH -A ", $aid, "\n";
    print DES "#SBATCH -t 00:10:00", "\n"; 
    print DES "#SBATCH -J DES", "\n";
    print DES "#SBATCH -n 1 ", "\n";
    print DES "#SBATCH -e $odf/deseq/info/deseq_wf.", $fnt ,".stderr.txt", "\n";
    print DES "#SBATCH -o $odf/deseq/info/deseq_wf.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print DES "#SBATCH --mail-type=All", "\n";
	print DES "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print DES 'echo "Running on: $(hostname)"',"\n\n";
    print DES "module add bioinfo-tools", "\n";
    print DES "module load R/2.12.2", "\n";
    print DES "#Samples", "\n";
    print DES 'inSampleDir="', "/bubo/proj/$aid/private/$odf/concat_htseqc", '"', "\n\n";
    print DES 'outSampleDir="', "/bubo/proj/$aid/private/$odf", '"', "\n\n"; #Deseq dir is created in r script
    
    print DES "Rscript /bubo/proj/$aid/private/$ods/Deseq_script.R ",'${inSampleDir}', "/$ochtsc", " $sidDES $nrcond1DES $nrcond2DES 0.05 ", '${outSampleDir}', "/$ochtsc";
    
    print DES "wait", "\n";
    close(DES);
    return;
    
}

sub Concat_htseqCD {
    
#Generates a sbatch script and runs concat_htseqCounts_wf.1.0.pl
    
    `mkdir -p $odf/concat_htseqc/info;`; #Creates the htseq and info data file directory
    `mkdir -p $ods/concat_htseqc;`; #Creates the concat_htseqCounts_wf.1.0.pl script directory
    
    $filename = "/bubo/proj/$aid/private/$ods/concat_htseqc/concat_htseqc_wf.";
    Checkfnexists($filename, $fnend);

    open (CHTSC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CHTSC "#! /bin/bash -l", "\n";
    print CHTSC "#SBATCH -A ", $aid, "\n";
    print CHTSC "#SBATCH -t 00:10:00", "\n"; 
    print CHTSC "#SBATCH -J CHTSC", "\n";
    print CHTSC "#SBATCH -n 1 ", "\n";
    print CHTSC "#SBATCH -e $odf/concat_htseqc/info/concat_htseqc_wf.", $fnt ,".stderr.txt", "\n";
    print CHTSC "#SBATCH -o $odf/concat_htseqc/info/concat_htseqc_wf.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CHTSC "#SBATCH --mail-type=All", "\n";
	print CHTSC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CHTSC 'echo "Running on: $(hostname)"',"\n\n";
    print CHTSC "#Samples", "\n";
    print CHTSC 'outSampleDir="', "/bubo/proj/$aid/private/$odf/concat_htseqc", '"', "\n\n";
    
    print CHTSC "perl /bubo/proj/$aid/private/$ods/concat_htseqCounts_wf.1.0.pl -i ";
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { 
	
	$_[0] = $sid[$sampleid];
	
	for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	    
	    my $tempinfile = $infiles{$_[0]}[$infile];
	    chomp $tempinfile; #Removing "\n";
	   
	    if ($infile eq scalar( @{ $infiles{$_[0]} })-1 && $sampleid eq scalar(@sid)-1) { #Last file and last sampleid
		
		print CHTSC "/bubo/proj/$aid/private/$odf/$_[0]/htseq/htseqCounts_$tempinfile.txt ";	
	    }
	    else {
		
		print CHTSC "/bubo/proj/$aid/private/$odf/$_[0]/htseq/htseqCounts_$tempinfile.txt,";
	    }
	}
    }
    print CHTSC "-s ";
    for (my $sidhtseqc=0;$sidhtseqc<scalar(@sidhtseq);$sidhtseqc++) {
	
	if ($sidhtseqc eq scalar(@sidhtseq)-1) {
	    
	    print CHTSC "$sidhtseq[$sidhtseqc] ";
	}
	else {
	    
	    print CHTSC "$sidhtseq[$sidhtseqc],";
	}
    }
    print CHTSC "-o ", '${outSampleDir}', "/$ochtsc", "\n\n";
    print CHTSC "wait", "\n";
    close(CHTSC);
    return;
    
}


sub Rpkmforgenes { 
    
#Generates sbatch scripts and runs Rpkmforgenes
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/rpkmforgenes/info;`; #Creates the rpkmforgenes and info data file directory
    `mkdir -p $ods/$_[0]/rpkmforgenes;`; #Creates the rpkmforgenes script directory
    
    $filename = "/bubo/proj/$aid/private/$ods/$_[0]/rpkmforgenes/rpkmForGenes.2.0_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    
    my $t = ceil(2.5*scalar( @{ $infiles{$_[0]} }) ); #Roughly 2.5 h for indexed reads going through samtools view and rpkmforgenes
    open (RPKM, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RPKM "#! /bin/bash -l", "\n";
    print RPKM "#SBATCH -A ", $aid, "\n";
    print RPKM "#SBATCH -t 1:00:00", "\n";
    print RPKM "#SBATCH -J RPKM", $_[0], "\n";
    print RPKM "#SBATCH -p node -n 8 ", "\n";
    print RPKM "#SBATCH -e $odf/$_[0]/rpkmforgenes/info/rpkmForGenes2.0_wf_$_[0].", $fnt, ".stderr.txt", "\n";
    print RPKM "#SBATCH -o $odf/$_[0]/rpkmforgenes/info/rpkmForGenes2.0_wf_$_[0].", $fnt, ".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print RPKM "#SBATCH --mail-type=All", "\n";
	print RPKM "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RPKM 'echo "Running on: $(hostname)"',"\n\n";
    print RPKM "module add bioinfo-tools", "\n";
    print RPKM "module load python/2.6.6", "\n";
    print RPKM "module load samtools/0.1.8", "\n\n";
    print RPKM "module load R/2.12.2", "\n\n";
    print RPKM "#Annotationfile","\n";
    print RPKM 'annotation_backg="', "/bubo/proj/$aid/private/$ods/references/background_hg19_noest.gtf",'"', "\n\n";
    print RPKM "#Samples", "\n";
    print RPKM 'outSampleDir="', "/bubo/proj/$aid/private/$odf/$_[0]/rpkmforgenes", '"', "\n";
    
    my $k=1;
    my $allp=0; #Required for correct timing of wait
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} } );$infile++) { #For all files platform
	
	if ($allp eq $k*8) { #Using only 8 cores
	    
	    print RPKM "wait", "\n";
	    $k=$k+1;
	}  
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print RPKM 'inSampleDir="',"/bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", '"', "\n\n";
	
	@refchr = `head -n 1 /bubo/proj/$aid/private/$ods/references/concat.fa.fa`; #Check for chr1 or simply 1 for adjusting reference
	#shift @refchr; #Removes header from array
	my $ref=0;
	for (my $pos=0;$pos<scalar(@refchr);$pos++) { #All headers (should not matter) and positions
	    
	    if($refchr[$pos] =~ /chr/) { #Adjust reference for chr or simply 1
		
		$ref = 1;
		$pos = scalar(@refchr);
	    }
	}
	if ($ref eq 1) {
	    print RPKM 'annotation="', "/bubo/proj/$aid/private/$ods/references/refGene-hg19.gtf",'"', "\n\n";
	}
	else {
	    
	    print RPKM 'annotation="', "/bubo/proj/$aid/private/$ods/references/refGene_no_chr-hg19.gtf",'"', "\n\n";
	}
	print RPKM "python /bubo/proj/$aid/private/$ods/rpkmforgenes.2.0.py -sam -gffann -readcount -i ", '${inSampleDir}', "/accepted_hits.sam -o ", '${outSampleDir}',"/$tempinfile", "_refg_background_rpkmout.txt -a ", '${annotation_backg}', " & \n\n";
	print RPKM "python /bubo/proj/$aid/private/$ods/rpkmforgenes.2.0.py -sam -readcount -i ", '${inSampleDir}', "/accepted_hits.sam -o ", '${outSampleDir}',"/$tempinfile", "_ref_rpkmout.txt -a ", '${annotation}', " & \n\n";
	$allp=$allp+2; #Two commands at once
    }
    
    print RPKM "wait", "\n\n";
    $k=1;
    print RPKM "cd ", '${outSampleDir}', "\n\n"; 
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} } );$infile++) { #For all files platform
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print RPKM "wait", "\n";
	    $k=$k+1;
	}  
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";

	print RPKM "Rscript /bubo/proj/$aid/private/$ods/cut_off.1.0.R ", '${outSampleDir}', "/$tempinfile", "_refg_background_rpkmout.txt ", '${outSampleDir}', "/$tempinfile", "_ref_rpkmout.txt $tempinfile > ", '${outSampleDir}', "/$tempinfile", "_cut_off.1.0_out.txt &", "\n\n";

    }
    print RPKM "wait", "\n";
    close(RPKM);
    return;
    
}


sub Htseqcount {
    
#Generates a sbatch script and runs HTSeqCount
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/htseq/info;`; #Creates the htseq and info data file directory
    `mkdir -p $ods/$_[0]/htseq;`; #Creates the htseq script directory
    
    $filename = "/bubo/proj/$aid/private/$ods/$_[0]/htseq/htseqCount_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(1.5*scalar( @{ $infiles{$_[0]} }) ); #One full lane on Hiseq takes approx. 1.5 h for HTSEQ to process, round up to nearest full hour.
    open (HTSC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print HTSC "#! /bin/bash -l", "\n";
    print HTSC "#SBATCH -A ", $aid, "\n";
    print HTSC "#SBATCH -t 1:00:00", "\n"; ##SBATCH -t $t:00:00", "\n";
    print HTSC "#SBATCH -J HTSC", $_[0], "\n";
    print HTSC "#SBATCH -p node -n 8 ", "\n";
    print HTSC "#SBATCH -e $odf/$_[0]/htseq/info/htseqCount_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print HTSC "#SBATCH -o $odf/$_[0]/htseq/info/htseqCount_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print HTSC "#SBATCH --mail-type=All", "\n";
	print HTSC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print HTSC 'echo "Running on: $(hostname)"',"\n\n";
    print HTSC "module add bioinfo-tools", "\n";
    print HTSC "module load python/2.6.6", "\n";
    print HTSC "module load htseq/0.4.6", "\n";
    print HTSC "module load samtools/0.1.8", "\n\n";
    print HTSC "#Annotationfile","\n";
    print HTSC 'annotation="', "/bubo/proj/$aid/private/$ods/references/Homo_sapiens.GRCh37.59.gtf",'"', "\n\n";
    print HTSC "#Samples", "\n";
    print HTSC 'outSampleDir_htsc="', "/bubo/proj/$aid/private/$odf/$_[0]/htseq", '"', "\n\n";
    
    my $k=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print HTSC "wait", "\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print HTSC 'inSampleDir="',"/bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", '"', "\n";
	print HTSC "cd /bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", "\n\n";
	print HTSC "python -m HTSeq.scripts.count -m intersection-strict -s no -t exon -i gene_id ", '${inSampleDir}', "/accepted_hits.sam ", '${annotation}', " > ", '${outSampleDir_htsc}', "/htseqCounts_$tempinfile.txt &", "\n\n";	
    }
    print HTSC "wait", "\n";
    close(HTSC);
    return;
    
}

sub Samtoolsview {
    
#Generates a sbatch script and runs Samtools view
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/tophat/info;`; #Creates the tophat and info data file directory
    `mkdir -p $ods/$_[0]/samtools;`; #Creates the samtools script directory
    
    $filename = "/bubo/proj/$aid/private/$ods/$_[0]/samtools/samtools_view_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(1.5*scalar( @{ $infiles{$_[0]} }) ); #One full lane on Hiseq takes approx. 1.5 h for samtools to process, round up to nearest full hour.
    open (STV, ">$filename") or die "Can't write to $filename: $!\n";
    
    print STV "#! /bin/bash -l", "\n";
    print STV "#SBATCH -A ", $aid, "\n";
    print STV "#SBATCH -t 1:00:00", "\n"; ##SBATCH -t $t:00:00", "\n";
    print STV "#SBATCH -J STV", $_[0], "\n";
    print STV "#SBATCH -p node -n 8 ", "\n";
    print STV "#SBATCH -e $odf/$_[0]/tophat/info/samtools_view_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print STV "#SBATCH -o $odf/$_[0]/tophat/info/samtools_view_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print STV "#SBATCH --mail-type=All", "\n";
	print STV "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print STV 'echo "Running on: $(hostname)"',"\n\n";
    print STV "module add bioinfo-tools", "\n";
    print STV "module load samtools/0.1.8", "\n\n";
    print STV "#Annotationfile","\n";
    print STV "#Samples", "\n";
    
    my $k=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print STV "wait", "\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print STV 'inSampleDir="',"/bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", '"', "\n";
	print STV 'outSampleDir="', "/bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", '"', "\n\n";
	print STV "cd /bubo/proj/$aid/private/$odf/$_[0]/tophat/$infile.$fnt", "\n\n";
	print STV "samtools view -h -o ", '${outSampleDir}', "/accepted_hits.sam ", '${inSampleDir}', "/accepted_hits.bam &", "\n\n";	
    }
    print STV "wait", "\n\n";

    print STV "cd /bubo/proj/$aid/private/", "\n\n"; #Required for next steps

    if ($pHTSC eq 1) { #If run HTSeqCount
	
	$filename2 = "/bubo/proj/$aid/private/$ods/$_[0]/htseq/htseqCount_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print STV "sbatch $filename2","\n\n";
	print STV "wait", "\n\n";
    }
    if ($pRPKM eq 1) { #If run RPKMforgenes
	
	$filename2 = "/bubo/proj/$aid/private/$ods/$_[0]/rpkmforgenes/rpkmForGenes.2.0_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print STV "sbatch $filename2","\n\n";
	print STV "wait", "\n\n";
    }
    close(STV);
    return;
    
}


sub Tophat {
#Generates sbatch scripts and runs Tophat on files frpm platform. 
#Note! Does not use SLURM module but locally installed top binaries. As well as reference taken from uppmax and exchanged chr names from chr1 to 1 etc.
#Last sbatch script will run two inputfiles and start HTSeqCount and RPKMforgenes when finished. The following assumptions is made to ensure that this script will finish last.
#1. It will be added to SLURM as the last script. 
#2. It will allocate twice the amount of time since it is running two input files, hence will find it harde to grasp a node. 
#3. It contains two input files and will therefore take longer to run, given the same nr of nodes. 
#$_[0] = sampleid
    
    `mkdir -p /bubo/proj/$aid/private/$odf/$_[0]/tophat/info;`; #Creates the tophat folder and info data file directory
    `mkdir -p /bubo/proj/$aid/private/$ods/$_[0]/tophat;`; #Creates the tophat script directory
    
    my $k=0;
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} } );$infile++) { #For all files platform 
	
	$filename = "/bubo/proj/$aid/private/$ods/$_[0]/tophat/tophat_wf_$_[0]_$k.";
	Checkfnexists($filename, $fnend);

	open (TH, ">$filename") or die "Can't write to $filename: $!\n";
	
	print TH "#! /bin/bash -l", "\n";
	print TH "#SBATCH -A ", $aid, "\n";
	print TH "#SBATCH -p node -n 8 ", "\n";
	print TH "#SBATCH -C thin", "\n";
	
	if ($infile eq ( scalar( @{ $infiles{$_[0]} }) -2) ) {
	    print TH "#SBATCH -t 60:00:00", "\n"; ##SBATCH -t 60:00:00", "\n";
	}
	else {
	    print TH "#SBATCH -t 30:00:00", "\n"; ##SBATCH -t 30:00:00", "\n";
	}
	
	print TH "#SBATCH -J TH", "$_[0]_",$k, "\n";
	print TH "#SBATCH -e /bubo/proj/$aid/private/$odf/$_[0]/tophat/info/tophat_wf_$_[0]_", $k,".$fnt.stderr.txt", "\n";
	print TH "#SBATCH -o /bubo/proj/$aid/private/$odf/$_[0]/tophat/info/tophat_wf_$_[0]_", $k,".$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print TH "#SBATCH --mail-type=All", "\n";
	    print TH "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print TH 'cd $HOME', "\n\n";
	print TH 'echo "Running on: $(hostname)"',"\n\n";
	print TH "module load bioinfo-tools", "\n\n"; 
	print TH 'module="', "/programs/tophat-1.1.4.Linux_x86_64",'"', "\n\n";
	print TH "module load samtools/0.1.8", "\n\n";

	print TH "#Reference", "\n";
	print TH 'reference="', "/bubo/proj/$aid/private/$ods/references/concat.fa",'"', "\n\n";

	print TH "#Annotationfile","\n";
	print TH 'annotation="', "/bubo/proj/$aid/private/$ods/references/Homo_sapiens.GRCh37.59.gtf",'"', "\n\n"; 

	print TH "#Samples", "\n";
	print TH 'inSampleDir="',"/bubo/proj/$aid/private/$odf/$_[0]/fastq", '"', "\n";
	
	
	if ($infile eq ( scalar( @{ $infiles{$_[0]} }) -2) ) { #To print the last two file in the same sbatch 

	    my $tempinfile = $infiles{$_[0]}[$infile];
	    chomp($tempinfile);
	    print TH "mkdir -p /bubo/proj/$aid/private/$odf/$_[0]/tophat/$k.$fnt", "\n"; #Creates the tophat folder and info data file directory
	    print TH 'outSampleDir="', "/bubo/proj/$aid/private/$odf/$_[0]/tophat/$k.$fnt", '"', "\n\n";
	    print TH '.${module}', "/tophat -o ", '${outSampleDir}', " --solexa1.3-quals -p 8 --rg-platform illumina --rg-platform-unit 1 --rg-center SciLife --GTF ", '${annotation} ${reference} ${inSampleDir}', "/$tempinfile", "\n\n";
	    print TH "wait", "\n\n";
	    
	    my $tempinfile2 = $infiles{$_[0]}[($infile+1)];
	    chomp($tempinfile2);
	    my $infile2 = $infile+1;
	    print TH "mkdir -p /bubo/proj/$aid/private/$odf/$_[0]/tophat/",$k+1, ".$fnt", "\n"; #Creates the tophat folder and info data file directory
	    print TH 'outSampleDir="', "/bubo/proj/$aid/private/$odf/$_[0]/tophat/", $k+1,".$fnt", '"', "\n\n";
	    print TH '.${module}', "/tophat -o ", '${outSampleDir}', " --solexa1.3-quals -p 8 --rg-platform illumina --rg-platform-unit 1 --rg-center SciLife --GTF ", '${annotation} ${reference} ${inSampleDir}', "/$tempinfile2", "\n\n";
	    print TH "wait", "\n\n";
	    
	    $infile=$infile+1; #Compensate for printing two files in last sbatch script
	    
	    if ($pSTV eq 1) { #If run HTSeqCount
		
		$filename2 = "/bubo/proj/$aid/private/$ods/$_[0]/samtools/samtools_view_wf_$_[0].";
		Checkfn2exists($filename2, $fnend);
		print TH "sbatch $filename2","\n\n";
		print TH "wait", "\n\n";
	    }
	    
	}
	else { #If not the last two files
	    
	    my $tempinfile = $infiles{$_[0]}[$infile];
	    chomp($tempinfile);
	    print TH "mkdir -p /bubo/proj/$aid/private/$odf/$_[0]/tophat/$k.$fnt", "\n"; #Creates the tophat folder and info data file directory
	    print TH 'outSampleDir="', "/bubo/proj/$aid/private/$odf/$_[0]/tophat/$k.$fnt", '"', "\n\n";
	    print TH '.${module}', "/tophat -o ", '${outSampleDir}', " --solexa1.3-quals -p 8 --rg-platform illumina --rg-platform-unit 1 --rg-center SciLife --GTF ", '${annotation} ${reference} ${inSampleDir}', "/$tempinfile", "\n\n";
	    print TH "wait", "\n\n";
	    
	    if ($pSTV eq 1 && scalar( @{ $infiles{$_[0]} } ) < 2 ) { #If run samtools view and only 1 infile
		
		$filename2 = "/bubo/proj/$aid/private/$ods/$_[0]/samtools/samtools_view_wf_$_[0].";
		Checkfn2exists($filename2, $fnend);
		print TH "sbatch $filename2","\n\n";
		print TH "wait", "\n\n";
	    }
	    
	    $k=$k+1; #Tracks nr of sbatch scripts
	}
	close(TH);
	#my $ret = `sbatch $filename`;
	#my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	#print STDERR "Sbatch script submitted, job id: $jobID\n";
	#print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
	#print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
    }
    
    return;
}

sub Fastqc {
    
#Generates a sbatch script and runs FASTQC
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/fastqc/info;`; #Creates the fastqc and info data file directory
    `mkdir -p $ods/$_[0]/fastqc;`; #Creates the fastqc script directory
    
    $filename = "$ods/$_[0]/fastqc/fastqc_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(0.5*scalar( @{ $infiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.
    open (FASTQC, ">$filename") or die "Can't write to $filename: $!\n";

    print FASTQC "#! /bin/bash -l", "\n";
    print FASTQC "#SBATCH -A ", $aid, "\n";
    print FASTQC "#SBATCH -t $t:00:00", "\n";
    print FASTQC "#SBATCH -J FQC", $_[0], "\n";
    print FASTQC "#SBATCH -p node -n 8 ", "\n";
    print FASTQC "#SBATCH -e $odf/$_[0]/fastqc/info/fastqc_$_[0].", $fnt ,".stderr.txt", "\n";
    print FASTQC "#SBATCH -o $odf/$_[0]/fastqc/info/fastqc_$_[0].", $fnt ,".stdout.txt", "\n";

    unless ($em eq 0) {
	
	print FASTQC "#SBATCH --mail-type=All", "\n";
	print FASTQC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print FASTQC 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQC "cd /bubo/proj/$aid/private/$odf/$_[0]/fastq", "\n\n";
    print FASTQC "module add bioinfo-tools", "\n\n";
    print FASTQC "module add FastQC/0.7.2", "\n\n";
    print FASTQC "#Samples", "\n";
    print FASTQC 'outsampleDir="',"/bubo/proj/$aid/private/$odf/$_[0]/fastqc", '"',"\n\n";
 
    my $k=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print FASTQC "wait", "\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print FASTQC "fastqc ", $tempinfile, ' -o ${outsampleDir}';
	print FASTQC " &", "\n\n";

    }
    
    print FASTQC "wait", "\n";
     
    close(FASTQC);
    my $ret = `sbatch $filename`;
    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
    print STDERR "Sbatch script submitted, job id: $jobID\n";
    print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
    print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
	
}

sub Checkfnexists {
    
#$_[0] = complete filepath
#$_[1] = file ending

    my $fn;
    $fnt = 0; #Nr of sbatch with identical filenames
    for (my $i=0;$i<999;$i++) { #Number of possible files with the same name
	
	$fn = "$_[0]$i$_[1]"; #filename, filenr and fileending
	$fnt = $i; #Nr of sbatch with identical filenames, global variable
	if (-e $fn) { #if file exists 
	}
	else {
	    $i=999; #Exit loop
	}
	
    }
    $filename = $fn; #Transfer to global variable
    return;
}

sub Checkfn2exists {
    
#$_[0] = complete filepath
#$_[1] = file ending

    my $fn;
    for (my $i=0;$i<999;$i++) { #Number of possible files with the same name
	
	$fn = "$_[0]$i$_[1]"; #filename, filenr and fileending
	if (-e $fn) { #if file exists 
	}
	else {
	    $i=999; #Exit loop
	}
	
    }
    $filename2 = $fn; #Transfer to global variable
    return;
}
