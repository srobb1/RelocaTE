#!/usr/bin/perl -w
use strict;
#use *confident.nonref.txt 

#TE	TSD	Exper	chromosome	insertion_site	strand	left_flanking_read_count	right_flanking_read_count	left_flanking_seq	right_flanking_seq
#1360	CACATAG	ab_neurons	2R:1372610..1372616	+	2	25	AAGCGACATATTTAAGAGCACTCATCAACAGAATAGTGCTCTTAAGATGCTATAACGTATTTAGTATTCTAAATCACAATCAAAAGTGGTAACATGCACA	TAGATTAAATACAAACCGAATAGATCCCAGCTTGCGTCATGAACACAGGGATAAGTACTGTAAAGCATTCAAATTATATTTATTGGCTAGGCATCATATG
#1360	GTACTAC	ab_neurons	2R:5086262..5086268	+	30	43	

#strain  TE      TSD   chromosome.pos    strand  avg_flankers    spanners        status
#EG4_2   mping   TTA   Chr1:2925..2927   +  32   0       homozygous

my %inserts;
my @insert_files = @ARGV;
foreach my $file (@insert_files){
  open IN, "$file" or die "Can't open $file $!\n";
  while (my $line = <IN>){
   chomp $line;
   next if $line =~ /chromosome.pos.+avg_flan/;
      my ($te, $tsd, $strain, $pos, $strand, $left_flankers, $right_flankers) = split /\t/ , $line;
      #$class =~ s/\?//;
      $inserts{$pos}{$strain}=$te;
  }
}
my %uniq;
foreach my $pos (keys %inserts){
  my @strains = keys %{$inserts{$pos}};
  my $strains = join ( ',' , (sort @strains));
  ## keep count of how many positions this combination of strains have in common 
  $uniq{$strains}{count}++;
  my @TEs;
  foreach my $strain (sort @strains){
    ## sort classes in same order as strains
    push @TEs , $inserts{$pos}{$strain};
  }
  my $TEs = join ( ',' , (@TEs));
  ## keep track of the classifications for this combination of strains at this position
  $uniq{$strains}{pos}{$pos}=$TEs;
}
my %TE_count;
print "----each position and the stains that share it and the TE\n";
foreach my $strains (keys %uniq){
  foreach my $pos ( keys %{$uniq{$strains}{pos}}){
    #my $class_count;
    my $TEs  = $uniq{$strains}{pos}{$pos};
    $TE_count{$strains}{$TEs}++;;
    print "$strains\t$pos\t$TEs\n";
  }
}
print "\n\n-----count of how many times this combination of strains share this same set of TEs\n";
foreach my $strains (keys %TE_count){
  foreach my $TEs (keys %{$TE_count{$strains}}){
    my $count = $TE_count{$strains}{$TEs};
    print "$strains\t$TEs\t$count\n";
  }
}
print "\n\n-----count of the times this combination of strains share a position\n";
foreach my $strains (keys %uniq){
  my $count = $uniq{$strains}{count};
  print "$strains\t$count\n";
}

