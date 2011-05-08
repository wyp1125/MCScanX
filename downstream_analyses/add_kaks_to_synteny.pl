use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::DNAStatistics;
use Getopt::Std;
%options=();
getopts("i:d:o:", \%options);
if(! exists $options{i}||! exists $options{d} ||!exists $options{o})
{
print "Usage:\nperl add_kaks_to_synteny.pl -i synteny_file -d cds_file -o output_file\n";
exit;
}
open(input0,$options{d}) or die "Cannot open cds_file!\n";
open(input1,$options{i}) or die "Cannot open synteny_file!\n";
open(output,">$options{o}") or die "Cannot open output_file!\n";
%s=();
$num=0;
while($line=<input0>)
{
if($line=~/^\>/)
{
if($num>0)
{
$s{$head}=$seq;
}
chomp($line);
$head=substr($line,1);
$seq="";
$num++;
}
else
{
$seq=$seq.$line;
}
}
$s{$head}=$seq;
while($line=<input1>)
{
chomp($line);
if($line eq "" ||$line=~/^\#/)
{
print output "$line\n";
next;
}
@a=split("\t",$line);
if(!exists $s{$a[1]}||!exists $s{$a[2]})
{
print output "$line\t-2\t-2\n";
}
#if(exists $s{$a[1]} && exists $s{$a[2]})
else
{
open(output1,">temp7734.cds");
print output1 "\>$a[1]\n$s{$a[1]}\>$a[2]\n$s{$a[2]}";
close(output1);
$tempcds = Bio::SeqIO -> new(-file => 'temp7734.cds', -format => 'fasta');
%dna_hash;
$cds1=$tempcds->next_seq;
$dna_hash{$cds1 -> display_id} = $cds1;
$cds2=$tempcds->next_seq;
$dna_hash{$cds2 -> display_id} = $cds2;
$os_prot = Bio::SeqIO -> new(-file=> '>temp7734.pro', -format=>'fasta');
$os_prot -> write_seq($cds1 -> translate());
$os_prot -> write_seq($cds2 -> translate());
system("clustalw -infile='temp7734.pro'");
$get_prot_aln = Bio::AlignIO -> new(-file=>"temp7734.aln", -format=>"CLUSTALW");
$prot_aln = $get_prot_aln -> next_aln();
$dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);
eval{
#print "here\n";
$stats = Bio::Align::DNAStatistics->new();
   $result = $stats->calc_all_KaKs_pairs($dna_aln);

   my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
   for my $an (@$result)
   {
     for (sort keys %$an )
     {
          next if /Seq/;
          if($_ eq "D_n"){$Dn = $an->{$_}};
          if($_ eq "D_s"){$Ds = $an->{$_}};
          if($_ eq "S_d"){$S_d = $an->{$_};}
          if($_ eq "N_d"){$N_d = $an->{$_};}
          if($_ eq "S"){$S = $an->{$_};}
          if($_ eq "N"){$N = $an->{$_};}
     }
   }
   if($Dn !~ /\d/){$Dn = -2;}
   if($Ds !~ /\d/){$Ds = -2;}
print output "$line\t$Dn\t$Ds\n";
};
if($@)
{
print output "$line\t-2\t-2\n";
next;
}
}
}
system("rm temp7734.*");
