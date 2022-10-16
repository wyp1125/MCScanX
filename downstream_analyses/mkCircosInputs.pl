#!/usr/bin/perl -w
$extend = 10000;

if( @ARGV != 2 ){
	print "USAGE: $0 gffFile collinearityFile\n";
	exit;
}
chomp $ARGV[0];
chomp $ARGV[1];

( $output = $ARGV[0] ) =~ s/gff$/karyotype/;
open( IN, "$ARGV[0]" ) || die( "Unable to open $ARGV[0]\n" );
while( <IN> ){
	chomp $_;
	@data = split( /\t/, $_ );
	$chr = substr( $data[0], 2 );
	next if( $chr !~ /^\d+$/ & ( $chr !~ /X|Y/i ) );
	$species = substr( $data[0], 0, 2 );
	if( !defined( $max -> { $species } -> { $chr } ) ){
		$max -> { $species } -> { $chr } = $data[3];
	}else{
		next if( $max -> { $species } -> { $chr } >= $data[3] );
		$max -> { $species } -> { $chr } = $data[3];
	}
	if( !defined( $min -> {	$species } -> {	$chr } ) ){
		$min -> { $species } -> { $chr } = $data[2];
	}else{
		next if( $min -> { $species } -> { $chr	} <= $data[3] );
		$min -> { $species } -> { $chr } = $data[3];
        }
	
}
close( IN );

open( OUT, ">$output" );
foreach $s ( sort keys %$min ){
	foreach $l ( sort { $a <=> $b } keys %{ $min -> { $s } } ){
		#print "$s\t$l\t", $min -> { $s } -> { $l }, "\t", $max -> { $s } -> { $l }, "\n";
		#print OUT "chr\ \-\ $s$l\ $l\ 0\ ",( $max -> { $s } -> { $l } + $extend ), "\ chr$l\n";
		print OUT "chr\ \-\ $s$l\ $s$l\ 0\ ",( $max -> { $s } -> { $l } + $extend ), "\ chr$l\n";
	}
}
close( OUT );

open( IN, "<$ARGV[0]" ) || die( "Unable to open $ARGV[0]\n" );
while( <IN> ){
	chomp $_;
	@data = split( /\t/, $_ );
#######
# at1	AT1G01280.1	112263	113947
#######
	$link{ $data[1] } = $data[0]." ".$data[2]." ".$data[3];
}
close( IN );

( $output = $ARGV[1] ) =~ s/collinearity/links/;
open( IN, "$ARGV[1]" ) || die( "Unable to open $ARGV[1]\n" );
open( OUT, ">$output" );
while( <IN> ){
	next if( $_ =~ /^\#/ );
	#next if( $_ !~ /\-\ \ 0\:\t/ );
	chomp $_;
	@data = split( /\t/, $_ );
	#print "$data[1]\t$data[2]\n";
	if( !defined( $link{ $data[1] } ) || !defined( $link{ $data[2] } ) ){
		print "Undefined $_\n";
	}else{
		print OUT $link{ $data[1] }, " ", $link{ $data[2] }, "\n";
	}
}
close( IN );
close( OUT );
