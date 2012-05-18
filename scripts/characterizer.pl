#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Cwd;
## script take output from find_TE_insertion.pl pipeline : ***.te_insertion.all.txt
## and a bam file of the same reads used to id insertiosn
## aligned to reference genome to look for spanner
## ***.te_insertion.all.txt has the following format:
## mping   A119    Chr1    1448    C:1     R:0     L:1
##
## examples to run this script
## for i in `seq 1 12`;do mping_homozygous_heterozygous.pl Chr$i.mping.te_insertion.all.txt ../bam_files/A119_500bp_Chr$i.sorted.bam > A119.Chr$i.inserts_characterized.txt ; done
## or
## mping_homozygous_heterozygous.pl A119.mping.te_insertion.all.txt ../bam_files/ > A119.inserts_characterized.txt
##
## 031412 added a check for XM and NM tags in sam line, if >0 then it is not 100% Match
## 022712 changed the logic for calling 'homo', 'het', 'new ins/somat' 

my $sites_file   = shift;
my $bam_dir      = shift;
my $genome_fasta = shift;
my $cwd          = getcwd();
my @bam_files

  ;    ## can be a single .bam file or a direcotory containing .bam files
if ( -d $bam_dir ) {

  #remove trailing '/'
  $bam_dir =~ s/\/$//;
  @bam_files = <$bam_dir/*bam>;
}
elsif ( -f $bam_dir or -l $bam_dir ) {
  push @bam_files, $bam_dir;
}
open INSITES, "$sites_file" or die "cannot open $sites_file $!\n";
my @dir_path = split '/', $sites_file;
my $filename = pop @dir_path;
$cwd =~ s/\/$//;    #remove trailing /
open OUTGFF, ">$cwd/$filename.homo_het.gff";
print
  "strain\tTE\tTSD\tchromosome.pos\tavg_flankers\tspanners\tstatus\n"; #\t$Smatch\t$cigar_all\n";
my %matches;

while ( my $line = <INSITES> ) {
  next if $line =~ /=/;
  next if $line =~ /^\s/;
  chomp $line;

  # mping   A119    Chr1    1448    C:1     R:0     L:1
  #my ( $te, $exp, $chromosome, $pos, $total_string, $right_string, $left_string ) = split /\t/, $line;
  #my $TSD = 'TTA';
  my ($te,$TSD,$exp,$chromosome, $pos, $total_string, $right_string, $left_string) = split /\t/, $line;
  my ($total_count) = $total_string =~ /C:(\d+)/;
  my ($left_count)  = $left_string  =~ /L:(\d+)/;
  my ($right_count) = $right_string =~ /R:(\d+)/;

  my $Mmatch = 0;
  my $cigar_all;
  if ( $left_count > 1 and $right_count > 1 ) {
    my @sam_all;
    foreach my $bam_file (@bam_files) {
      ## get any alignments that overlap the insertion site
      my @sam_out = `samtools view $bam_file \'$chromosome:$pos-$pos\'`;

#my @header = `samtools view -HT $genome_fasta $bam_file \'$chromosome:$pos-$pos\'`;
#my $header = join "\n" , @header;
#$matches{"$chromosome.pos"}{header} = $header;
      push @sam_all, @sam_out;
    }

    #remove redundant lines in sam file
    my %sorted_sam;
    my $order;
    foreach my $line (@sam_all) {
      $order++;
      if ( !exists $sorted_sam{$line} ) {
        $sorted_sam{$line} = $order;
      }
    }

    #make new sorted sam array by sorting on the value of the sort hash
    my @sorted_sam =
      sort { $sorted_sam{$a} <=> $sorted_sam{$b} } keys %sorted_sam;

    foreach my $sam_line (@sorted_sam) {
      chomp $sam_line;
      my @sam_line = split /\t/, $sam_line;
      my $cigar    = $sam_line[5];
      my $flag     = $sam_line[1];
      my $seqLen   = length $sam_line[9];
      my $start    = $sam_line[3];
      my $end      = $start + $seqLen - 1;
      next unless $end >= $pos + 5;
      next unless $start <= $pos - 5;
      ## must be a all M match no soft clipping
      #if ( $cigar =~ /^[1-3]?S?\d+M[1-3]?S?$/ ) { ##allow for a couple of soft clippings
      if ( $cigar =~ /^\d+M$/ ) {
        my ($NM) = $sam_line =~ /NM:i:(\d+)/; ## edit distance used
        my ($XM) = $sam_line =~ /XM:i:(\d+)/; ## bwa specific: mismatch count, often the same as NM, have not seen a XM>0 and NM==0
        if (defined $XM and $XM == 0){
          $Mmatch++;
        }elsif(defined $NM and $NM == 0){
          $Mmatch++;
        }elsif (!defined $NM and !defined $XM){
          $Mmatch++;
        }else{
          $matches{"$chromosome.$pos"}{sam}{$sam_line} = 1;
        }
      }
      #elsif ( $cigar !~ /S/ and $cigar =~ /[IND]/ ) {
      elsif ( $cigar =~ /[IND]/ ) {

        #push @{$matches{"$chromosome.$pos"}{sam}} , $sam_line;
        $matches{"$chromosome.$pos"}{sam}{$sam_line} = 1;
      }else {
        $matches{"$chromosome.$pos"}{sam}{$sam_line} = 1;

      }

      #$cigar_all.="$cigar,";
    }
    my $spanners         = $Mmatch;
    my $average_flankers = $total_count / 2;
    my $status           = 0;
    ## this characterization needs some statistics
    #if ( $average_flankers < 5 and ( ( $spanners - $average_flankers ) >= 10 ) )
    #{
    #  $status = 'new_insertion';
    #}
    #elsif ( ( ( $average_flankers - $spanners ) > 10 )
    #  and $spanners < 5 )
    #{
    #  $status = 'somatic_excision';
    #}
    ##avg_flankers already has to be (left>1 right>1 therefore atleast 2+2/2 =>) 2
    if ( $spanners == 0 ) {
      $status = 'homozygous';
    }
    elsif ( ($spanners - $average_flankers) > ($spanners + $average_flankers)/2){
      $status = 'new_insertion';
    }
    elsif ( abs( $average_flankers - $spanners ) <= 5 ) {
      $status = 'heterozygous';
    }
    elsif ( $average_flankers >= 9 and $spanners < 5 ) {
      $status = 'homozygous?';
    }
    elsif ( $average_flankers >= 5 and $spanners < 5 ) {
      $status = 'homozygous?';
    }
    elsif (
      abs( $average_flankers - $spanners ) <=
      ( ( $average_flankers + $spanners ) / 2 ) )
    {
      $status = 'heterozygous?';
    }
    elsif ($spanners < ($average_flankers * .2) and $spanners <= 10){
      $status = 'homozygous?';
    }
    #elsif ( ( $spanners - $average_flankers ) > 10 ) {
    #  $status = 'new_insertion?';
    #}
    $matches{"$chromosome.$pos"}{status} = $status;
    print "$exp\t$te\t$TSD\t$chromosome.$pos\t$average_flankers\t$spanners\t$status\n"
      ;    #\t$Smatch\t$cigar_all\n";
    print OUTGFF
"$chromosome\t$exp\ttransposable_element_attribute\t$pos\t$pos\t.\t.\t.\tID=$chromosome.$pos.spanners;avg_flankers=$average_flankers;spanners=$spanners;type=$status;TE=$te;TSD=$TSD\n";

#print "$chromosome.$pos\t$total_count\t$left_count\t$right_count\t$Mmatch\t$status\n";#\t$Smatch\t$cigar_all\n";
  }
}
##generate vcf of spanners looking for excision events
__END__
my @unlink_files;
foreach my $pos ( keys %matches ) {
  my ( $target, $loc ) = split /\./, $pos;
  my $range      = "$target:$pos";
  my $sam        = "$cwd/$pos.sam";
  my $bam        = "$cwd/$pos.bam";
  my $sorted_bam = "$cwd/$pos.sorted";
  my @sam_lines  = keys %{ $matches{$pos}{sam} };
  if ( @sam_lines > 1 ) {
    my $pos_sam = join "\n", @sam_lines;
    open POSSAM, ">$sam";
    print POSSAM $pos_sam;
    `samtools view -bT $genome_fasta $sam > $bam`;
    `samtools sort $bam $sorted_bam`;
    `samtools index $sorted_bam.bam`;

`samtools mpileup -C50 -ugf $genome_fasta -r $range  $sorted_bam.bam | bcftools view -bvcg - > $cwd/$pos.var.raw.bcf`;
`bcftools view $cwd/$pos.var.raw.bcf | vcfutils.pl varFilter -D 100 > $cwd/$pos.var.flt.vcf`;

    push @unlink_files, $sam, $bam, "$sorted_bam.bam.bai","$sorted_bam.bam", "$cwd/$pos.var.raw.bcf";
    close POSSAM;
  }
}
unlink @unlink_files;
