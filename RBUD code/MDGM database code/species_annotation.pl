

    #!/usr/bin/perl
    use strict;
    use warnings;
my $dir_to_process = '/database and method code/MDGM database/Species databset/Annotation/NCBI species annotation/all.gbk';
open OUT,">/database and method code/MDGM database/Species dataset/Annotation/Species annotation" or die "$!";    
opendir my $dh,$dir_to_process or die "Cannot open $dir_to_process:$!";
   foreach my $file(readdir $dh) {
open MZ,"</database and method code/MDGM database/Species dataset/Annotation/NCBI species annotation/all.gbk/$file" or die "$!"; 
while(<MZ>) {
	my $tag = $_;
	chomp $tag;
      if($tag=~/LOCUS/) {
                       my $mm = $_;
                       chomp $mm;
                       my @data2 = split(/\s+/,$mm);
                       print OUT "$data2[1]\n";
                        }
        
      if($tag =~/ORGANISM/) {
			 my $nm = $_;
			 chomp $nm;
		       substr($nm,1,9,"organism\;");
                         print OUT "$nm\;";
                         my @data = <MZ>;       
                         foreach my $seq (@data) {
                                                chomp $seq;
                                                unless($seq=~/REFERENCE/) {
                                                                     print OUT "$seq";
                                                                          }
                                                  else{@data=();}
                                                 }
              
                              
                           }
       else {next;}
       print OUT "\n";
                    
           }
   close MZ;
                                 } 
   close OUT; 
   closedir $dh;                    
