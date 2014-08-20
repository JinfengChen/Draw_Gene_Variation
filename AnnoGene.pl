#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"gff:s","project:s","help");


my $help=<<USAGE;
perl $0 --gff
Read GFF file, annotate gene in gff with SNP, INDEL, SV effect, as well as  *.anno and *.inf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "QTL.mPing";
my $refanno=readanno("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.anno");
my $refinf =readinf("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/RiceGene/MSU7.gene.inf");
my $snp =readSNP("../input/HEG4.SNP.annotation.gene.anno");
my $indel=readSNP("../input/HEG4.INDEL.annotation.gene.anno");
my $del  =readSV("../input/HEG4.Deletion.final.valid.eff.sum");
my $int  =readSV("../input/HEG4.insertion.final.valid.eff.sum");
#my $fpkm =readcufflink("../input/cufflink.table");
my $fpkm=readfpkm("/rhome/cjinfeng/BigData/02.Transcription/Transcriptome/bin/DiffExpression/NB_HEG4_EG4.table");
readgff($opt{gff},$refanno,$refinf,$snp,$indel,$del,$int,$fpkm);

#Chr1    MSU_osa1r7      gene    4234514 4244531 .       -       .       ID=LOC_Os01g08540;Name=LOC_Os01g08540;Note=DNA%20mismatch%20repair
sub readgff
{
my ($file,$refanno,$refinf,$snp,$indel,$del,$int,$fpkm)=@_;
print "Gene\tSNPeff\tINDELeff\tDELeff\tINTeff\tNB_fpkm\tEG4_fpkm\tFDR\tNB_fpkm\tHEG4_fpkm\tFDR\tAnnotation\tOryzabase\n";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2] eq "gene" and $unit[8]=~/ID=(.*?);/){
       my $snpeff=$snp->{$1} ? join(";",@{$snp->{$1}}) : "NA";
       my $indeleff=$indel->{$1} ? join(";",@{$indel->{$1}}) : "NA";
       my $deleff=$del->{$1} ? join(";",@{$del->{$1}}) : "NA";
       my $inteff=$int->{$1} ? join(";",@{$int->{$1}}) : "NA";
       my $expr  =$fpkm->{$1} ? $fpkm->{$1} : "NA\tNA\tNA";
       print "$1\t$snpeff\t$indeleff\t$deleff\t$inteff\t$expr\t$refanno->{$1}\t$refinf->{$1}\n";
    }
}
close IN;
}
 


sub readanno
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0] = $1 if ($unit[0]=~/(.*?)\.\d+/);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub readinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[2]}=$_;
}
close IN;
return \%hash;
}



#LOC_Os01g01070  DOWNSTREAM      expressed protein       
#LOC_Os01g01070  INTRON  expressed protein       
#LOC_Os01g01070  NON_SYNONYMOUS_CODING   expressed protein       
sub readSNP
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}}, $unit[1];
}
close IN;
return \%hash;
}

#Insertion_1     LOC_Os01g01369  Intron  Intron  3-beta-hydroxysteroid-Delta-isomerase,putative,expressed        NA
#Insertion_1     LOC_Os01g01380  Upstream        Upstream        expressed protein       NA
#Insertion_10    LOC_Os01g05020  Upstream        Upstream        expressed protein       NA
#Insertion_1101  LOC_Os03g53340  Partial-gene    Utr3prime;Intron        HSF-type DNA-binding domain containing protein,expressed        8548    Os03g0745000
sub readSV
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[1]}}, $unit[3];
}
close IN;
return \%hash;
}


#Gene    EG4     HEG4    NB
#LOC_Os01g01010.1        8.99969 12.3017 12.5327
sub readcufflink
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $gene=$1 if ($unit[0]=~/(LOC_\w+)\.\d+/);
    $hash{$gene}="$unit[1]\t$unit[2]\t$unit[3]";
}
close IN;
return \%hash;
}



#Gene    NB      EG4     FoldChange      P-value FDR     NB      HEG4    FoldChange      P-value FDR
#LOC_Os01g01010  1810.32 1390.17 0.77    0.02    0.39    1318.05 1414.15 1.07    0.65    1.00    Up
sub readfpkm
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $gene=$unit[0];
    $hash{$gene}="$unit[1]\t$unit[2]\t$unit[5]\t$unit[6]\t$unit[7]\t$unit[10]";
}
close IN;
return \%hash;
}

