use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;
#use Config::Tiny;

my ($fastq,$name,$read_len,$outdir);

GetOptions(
	"fq:s"        => \$fastq,        # can be .gz or .fq or .fastq
	"n:s"         => \$name,         # Need
	"rlen:i"      => \$read_len,     # Default: 200bp SE
	"od:s"        => \$outdir,       # Need
	) or die;


# check args
if (not defined $fastq || not defined $name || not defined $outdir){
	die "please check you args!\n";
}

# default value
if (not defined $read_len){
	$read_len = 200;
}

# outfile
my $out_fa = "$outdir/$name/$name.SE$read_len\.fasta";

if (!-d "$outdir/$name"){
	`mkdir -p $outdir/$name`;
}


open FA, ">$out_fa" or die;

if ($fastq =~ /\.fq$/ || $fastq =~ /\.fastq$/){
	open IN, "$fastq" or die;
}else{
	open IN, "gunzip -dc $fastq |" or die;
}


while (<IN>){
	chomp;
	my $header = $_;
	$header =~ s/^\@//;
	my $seq = <IN>;
	chomp $seq;
	<IN>;
	<IN>; # qual

	my $seq_len = length($seq);
	my $bin_count = int($seq_len/$read_len);
	
	
	my $new_header;
	for my $idx (1..$bin_count){
		my $start_pos = $read_len * ($idx - 1);
		my $subseq = substr($seq,$start_pos,$read_len);
		$new_header = "$header\_$idx";
		print FA "\>$new_header\n";
		print FA "$subseq\n";
	}
}
close IN;
close FA;
