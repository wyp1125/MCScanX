use Getopt::Std;
%options=();
getopts("i:o:", \%options);
if(! exists $options{i}||! exists $options{o})
{
print "Usage:\nperl group_syntenic_genes.pl -i synteny_file -o output_file\n";
exit;
}
open(input,$options{i}) or die "Cannot open synteny_file!\n";
open(output,">$options{o}") or die "Cannot open output_file!\n";
%h=();
$id=0;
while($line=<input>)
{
chomp($line);
if($line eq "" ||$line=~/^\#/)
{
next;
}
@a=split("\t",$line);
if(!exists $h{$a[1]} && !exists $h{$a[2]})
{
$h{$a[1]}=$id;
$h{$a[2]}=$id;
$id++;
}
elsif(!exists $h{$a[1]} && exists $h{$a[2]})
{
$h{$a[1]}=$h{$a[2]};
}
elsif(exists $h{$a[1]} && !exists $h{$a[2]})
{
$h{$a[2]}=$h{$a[1]};
}
else
{
$minv=$h{$a[0]};
$maxv=$h{$a[1]};
if($h{$a[1]}<$h{$a[0]})
{
$minv=$h{$a[1]};
$maxv=$h{$a[0]};
}
for $key (keys %h)
{
if($h{$key}==$maxv)
{
$h{$key}=$minv;
}
}
}
}
%g=();
for $key (keys %h)
{
$g{$h{$key}}=$g{$h{$key}}."\t".$key;
}
%f=();
for $key (keys %g)
{
@a="";
@a=split("\t",$g{$key});
$f{$g{$key}}=$#a;
}
$id=0;
foreach $value (sort {$f{$b} <=> $f{$a} }keys %f)
{
print output "$id$value\n";
$id++;
}




