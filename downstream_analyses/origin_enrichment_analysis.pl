use Getopt::Std;
%options=();
getopts("i:d:o:", \%options);
if(! exists $options{i}||! exists $options{d})
{
print "Usage:\nperl origin_enrichment_analysis.pl -i gene_family_file -d gene_type_file\n";
print "optional: -o output_file\n";
exit;
}
open(input0,$options{d}) or die "Cannot open gene_type_file!\n";
open(input1,$options{i}) or die "Cannot open gene_family_file!\n";
$write_file=0;
if(exists $options{o})
{
open(output,">$options{o}") or die "Cannot open output_file!\n";
$write_file=1;
}
%h=();
@num=(0,0,0,0,0);
while($line=<input0>)
{
chomp($line);
@a=split("\t",$line);
if($a[1]!=0)
{
$h{$a[0]}=$a[1];
$num[$a[1]]++;
}
}
$num[0]=scalar keys %h,"\n";

if($write_file==1)
{
print output "ID\tDipsersed\tProximal\tTandem\tWGD or segmental\n";
}
else
{
print "ID\tDipsersed\tProximal\tTandem\tWGD or segmental\n";
}

while($line=<input1>)
{
chomp($line);
@a="";
@a=split("\t",$line);
if($#a<=1)
{
next;
}
@dn=(0,0,0,0,0);
for($j=1;$j<=$#a;$j++)
{
if(exists $h{$a[$j]})
{
$dn[$h{$a[$j]}]++;
}
}
$dn[0]=$dn[1]+$dn[2]+$dn[3]+$dn[4];
for($j=1;$j<=4;$j++)
{
$pval[$j]=min(1,4*pvalue($dn[$j],$dn[0],$num[$j],$num[0]));
}
if($write_file==1)
{
print output "$a[0]\t$pval[1]\t$pval[2]\t$pval[3]\t$pval[4]\n";
}
else
{
print "$a[0]\t$pval[1]\t$pval[2]\t$pval[3]\t$pval[4]\n";
}
#print "$dn[0] $dn[1] $dn[2] $dn[3] $dn[4]\n";
}
#print pvalue(2,9,10,19);
sub pvalue
{
my $k=$_[0];
my $n=$_[1];
my $C=$_[2];
my $G=$_[3];
my $um=min($n,$C);
my $lm=max(0,$n+$C-$G);
if($um==$lm)
{
return 1.0;
}
my $cutoff = hypergeometric_probability($k, $n, $C, $G);
my $right_tail = 0;
for(my $i=$lm;$i<$um+1;$i++)
{
my $p = hypergeometric_probability($i, $n, $C, $G);
if ($i>=$k)
{
$right_tail += $p;
}
}
$right_tail = min($right_tail, 1);
return $right_tail;
}
sub min
{
return $_[0]<$_[1]?$_[0]:$_[1];
}
sub max
{
return $_[0]>$_[1]?$_[0]:$_[1];
}
sub hypergeometric_probability 
{
my $i=$_[0];
my $n=$_[1];
my $C=$_[2];
my $G=$_[3];
return exp(lncombination($C,$i)+lncombination($G-$C,$n-$i)-lncombination($G,$n));
}
sub lncombination
{
my $n=$_[0];
my $p=$_[1];
return lnfactorial($n) - lnfactorial($p) - lnfactorial($n - $p);
}

sub lnfactorial
{
my $n=$_[0];
if($n <= 1)
{return 0;}
else
{return lngamma($n + 1);}
}

sub lngamma
{
my $z=$_[0];
my $x=0;
$x += 0.1659470187408462e-06 / ($z + 7);
$x += 0.9934937113930748e-05 / ($z + 6);
$x -= 0.1385710331296526 / ($z + 5);
$x += 12.50734324009056 / ($z + 4);
$x -= 176.6150291498386 / ($z + 3);
$x += 771.3234287757674 / ($z + 2);
$x -= 1259.139216722289 / ($z + 1);
$x += 676.5203681218835 / ($z);
$x += 0.9999999999995183;

return log($x) - 5.58106146679532777 - $z + ($z - 0.5) * log($z + 6.5);
}
