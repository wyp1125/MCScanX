use Getopt::Std;
%options=();
getopts("i:d:o:", \%options);
if(! exists $options{i}||! exists $options{d})
{
print "Usage:\nperl detect_collinearity_within_gene_families.pl -i gene_family_file -d collinearity_file\n";
print "optional: -o output_file\n";
exit;
}
open(input0,$options{d}) or die "Cannot open collinearity_file!\n";
open(input1,$options{i}) or die "Cannot open gene_family_file!\n";
$write_file=0;
if(exists $options{o})
{
open(output,">$options{o}") or die "Cannot open output_file!\n";
$write_file=1;
}
%g=();
%p=();
while($line=<input0>)
{
chomp($line);
if($line eq "" ||$line=~/^\#/)
{
next;
}
@a=split("\t",$line);
$g{$a[1]}=1;
$g{$a[2]}=1;
$gid=$a[1]."\:".$a[2];
if($a[1] gt $a[2])
{
$gid=$a[2]."\:".$a[1];
}
$p{$gid}=1;
}
############################################
while($line=<input1>)
{
#print $line;
chomp($line);
@a="";
@a=split("\t",$line);
if($#a<=1)
{
next;
}
$num=0;
for($j=1;$j<=$#a;$j++)
{
if(exists $g{$a[$j]})
{
$gene[$num]=$a[$j];
$num++;
}
}
$oline=$a[0];
for($i=0;$i<$num-1;$i++)
{
for($j=$i+1;$j<$num;$j++)
{
$gid=$gene[$i]."\:".$gene[$j];
if($gene[$i] gt $gene[$j])
{
$gid=$gene[$j]."\:".$gene[$i];
}
if(exists $p{$gid})
{
$oline=$oline."\t".$gid;
}
}
}
if($write_file==1)
{
print output "$oline\n";
}
else
{
print "$oline\n";
}
}
