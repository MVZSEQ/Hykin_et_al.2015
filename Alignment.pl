#!/usr/bin/perl
#use Your::Power;
use strict;
use Getopt::Std;
use File::Basename;
use warnings;
#Ke Bi (kebi@berkeley.edu)

die(qq/
Usage: Alignment.pl [options] 

external dependencies: novoalign, SAMtools

Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    FILE    Reference fasta sequence file 
              (for example: ABC_123.fasta)
-r    DIR     Path to the directory of cleaned reads
-o    DIR     Path to the results directory
-i    INT     Avg. insert size for the libraries
-v    INT     STD insert size for the libraries (usually 0.1*i)
-l    INT     read length [100] 
-t    INT     Maximum alignment score acceptable for the best 
              alignment. Use 90 for population samples (roughly 
              3 mismatches are allowed per pair). If use 150, 
              roughly 5 mismatches are allowed per read in PE. 
              For very divergent genomes, use default value 
              by not defining -t [null]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

\n\n/) unless (@ARGV);


my %opts = (f=>undef, r=>undef, v=>undef, o=>undef, i=>undef, t=>undef, l=>100);
getopts('f:r:o:i:t:v:l:', \%opts);

my $length = $opts{l} + $opts{l};

my $data_dir;
if ($opts{r} =~ m/\/$/ ){
  $data_dir = $opts{r}; 
}
else {
  $data_dir = $opts{r} . "/";
}

my $res_dir;
if ($opts{o} =~ m/\/$/ ){
  $res_dir = $opts{o}; 
}
else {
  $res_dir = $opts{o} . "/";
}


my $indexed_assemblies_in_target =  substr ($opts{f}, 0, -5) . "nix";
system ("novoindex $indexed_assemblies_in_target $opts{f}");


my @reads1 = <$data_dir*_1_final.txt>;

foreach my $read1 (@reads1) {
  my $read2 = $read1; $read2 =~ s/_1_final/_2_final/;
  my $readSolo = $read1; $readSolo =~ s/_1_final/_u_final/;
  my $lib = $1 if basename($read1) =~ m/(\S+)_1_final/;
  my $sorted_in_target_bams = $res_dir . $lib . "_sorted";
  if ($opts{t}) {
    my $alignment_paired_to_in_target_assemblies= system("novoalign -R 30 -t $opts{t} -n $length  -d $indexed_assemblies_in_target -f $read1  $read2  -i PE $opts{i}, $opts{v}  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outPairedSam1");
    
    my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -t $opts{t}   -d $indexed_assemblies_in_target -f $readSolo  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM > outSoloSam1");
    
  }
  
  if (!$opts{t}) {
    
    my $alignment_paired_to_in_target_assemblies= system("novoalign -R 30 -n $length   -d $indexed_assemblies_in_target -f $read1  $read2  -i PE $opts{i}, $opts{v}  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outPairedSam1");
    
    my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -d $indexed_assemblies_in_target -f $readSolo  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM > outSoloSam1");
    
    
  }
  
  system("grep -v ZS:Z:NM outPairedSam1 > paired_in_target.sam");
  system("grep -v ZS:Z:NM outSoloSam1 > solo_in_target.sam");
  
  system("samtools view -bS paired_in_target.sam >  paired_in_target.bam");
  system("samtools view -bS solo_in_target.sam > solo_in_target.bam");
  
  
  system("samtools merge raw.bam paired_in_target.bam solo_in_target.bam");
  
  system("samtools sort raw.bam $sorted_in_target_bams");
  
  my $sorted_in_target_bams2 = $sorted_in_target_bams. '.bam';
  
  system("samtools index $sorted_in_target_bams2");  

  system("rm outPairedSam1 outSoloSam1 paired_in_target.sam solo_in_target.sam paired_in_target.bam solo_in_target.bam raw.bam"); 
  
}

