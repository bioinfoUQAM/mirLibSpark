use strict;
use warnings;

my $infile = $ARGV[0];
open my $IN, $infile or die $!;

print "\n";
while (my $line = <$IN>){
	print $line;
}

close $IN;