#!/usr/bin/perl
use strict;
use warnings;
my $dir_to_process = '/database and method code/MDGM database/Function dataset/Annotation/NCBI functional annotation/';
open OUT,">/dataset and method code/MDGM database/Functional dataset/Annotation/gene.gff" or die "$!";    
opendir my $dh,$dir_to_process or die "Cannot open $dir_to_process:$!";
   foreach my $file(readdir $dh) {
open MZ,"</dataset and method code/MDGM database/Functional dataset/Annotation/NCBI functional annotation/$file" or die "$!";
my @data2; 
while(<MZ>) {
        my $tag = $_;
        chomp $tag;
      if($tag=~/NC_*/) {
                    my $mm = $_;
                    chomp $mm;
                    my @data2 = split(/\s+/,$mm);
                    #print OUT "$data2[2]\n";
                    if ($data2[2]=~"gene")
                    {print OUT "$mm\n";}
                        }
             }
close MZ;
      }
close OUT;
closedir $dh;


