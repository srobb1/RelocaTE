#!/usr/bin/perl -w
use strict;
#use inserts.characterixed.txt

#strain  TE      TSD   chromosome.pos    strand  avg_flankers    spanners        status
#EG4_2   mping   TTA   Chr1:2925..2927   +  32   0       homozygous

my %inserts;
my @insert_files = @ARGV;
foreach my $file (@insert_files){
  open IN, "$file" or die "Can't open $file $!\n";
  while (my $line = <IN>){
   chomp $line;
   next if $line =~ /chromosome.pos.+avg_flan/;
      my ($strain, $te, $tsd, $pos, $strand, $spanners, $flankers, $class) = split /\t/ , $line;
      $class =~ s/\?//;
      $inserts{$pos}{$strain}=$class;
  }
}
my %uniq;
foreach my $pos (keys %inserts){
  my @strains = keys %{$inserts{$pos}};
  my $strains = join ( ',' , (sort @strains));
  ## keep count of how many positions this combination of strains have in common 
  $uniq{$strains}{count}++;
  my @classes;
  foreach my $strain (sort @strains){
    ## sort classes in same order as strains
    push @classes , $inserts{$pos}{$strain};
  }
  my $classes = join ( ',' , (@classes));
  ## keep track of the classifications for this combination of strains at this position
  $uniq{$strains}{pos}{$pos}=$classes;
}
my %class_count;
print "----each position and the stains that share it and the classifications\n";
foreach my $strains (keys %uniq){
  foreach my $pos ( keys %{$uniq{$strains}{pos}}){
    #my $class_count;
    my $classes  = $uniq{$strains}{pos}{$pos};
    $class_count{$strains}{$classes}++;;
    print "$strains\t$pos\t$classes\n";
  }
}
print "\n\n-----count of how many times this combination of strains share this same set of classifications\n";
foreach my $strains (keys %class_count){
  foreach my $classes (keys %{$class_count{$strains}}){
    my $count = $class_count{$strains}{$classes};
    print "$strains\t$classes\t$count\n";
  }
}
print "\n\n-----count of the times this combination of strains share a position\n";
foreach my $strains (keys %uniq){
  my $count = $uniq{$strains}{count};
  print "$strains\t$count\n";
}

