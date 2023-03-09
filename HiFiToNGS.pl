use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;
#use Config::Tiny;



my ($fastq,$name,$ins_frag_len,$read_len,$outdir);

GetOptions(
	"fq:s"        => \$fastq,        # can be .gz or .fq or .fastq
	"n:s"         => \$name,         # Need
	"ilen:i"      => \$ins_frag_len, # Default: 500
	"rlen:i"      => \$read_len,     # Default: PE150
	"od:s"        => \$outdir,       # Need
	) or die;


# check args
if (not defined $fastq || not defined $name || not defined $outdir){
	die "please check you args!\n";
}

# default value
if (not defined $ins_frag_len){
	$ins_frag_len = 500;
}

if (not defined $read_len){
	$read_len = 150;
}


my $out_fq1 = "$outdir/$name/$name\.PE$read_len\.R1.fastq";
my $out_fq2 = "$outdir/$name/$name\.PE$read_len\.R2.fastq";

if (!-d "$outdir/$name"){
	`mkdir -p $outdir/$name`;
}

open FQ1, ">$out_fq1" or die;
open FQ2, ">$out_fq2" or die;


if ($fastq =~ /\.fq$/ || $fastq =~ /\.fastq$/){
	open IN, "$fastq" or die;
}else{
	open IN, "gunzip -dc $fastq |" or die;
}


while (<IN>){
	chomp;
	my $header = $_;
	my $seq = <IN>;
	chomp $seq;
	
	<IN>;
	
	my $qual = <IN>;
	chomp $qual;

	#print "$seq\n";
	#print "$qual\n";
	# overlap is 200bp
	# if ccs len is 10kb, the ins_frag len is 500, read len is 150
	# read1: ----- 500bp
	# read2:    ----- 500bp

	my $overlap_len = int($ins_frag_len * 0.8); # if ins frag len is 500, then overlap len is 400
	my $skip_len = $ins_frag_len - $overlap_len;
	my $seq_len = length($seq);
	my $qual_len = length($qual);

	my $idx = 0;
	my $flag = 0;
	my $max_idx = $seq_len - $ins_frag_len; # 10 - 5 (6/7/8/9/10)
	while ($idx <= $max_idx){
		my $frag = substr($seq,$idx,$ins_frag_len);
		my $fwd_r1 = substr($frag,0,$read_len);
		my $frag_rev_comp = reverse($frag);
		$frag_rev_comp =~ tr/ATCG/TAGC/;
		my $rev_r2 = substr($frag_rev_comp,0,$read_len);

		# get qual info
		my $frag_qual = substr($qual,$idx,$ins_frag_len);
		my $fwd_r1_qual = substr($frag_qual,0,$read_len);
		my $frag_qual_rev = reverse($frag_qual);
		my $rev_r2_qual = substr($frag_qual_rev,0,$read_len);


		# <read>:<is filtered>:<control number>:<sample number>
		$flag += 1;
		print FQ1 "$header\_$flag 1\:N\:0\:1\n"; # https://help.basespace.illumina.com/files-used-by-basespace/fastq-files
		print FQ1 "$fwd_r1\n";
		print FQ1 "\+\n";
		print FQ1 "$fwd_r1_qual\n";

		print FQ2 "$header\_$flag 2\:N\:0\:1\n";
		print FQ2 "$rev_r2\n";
		print FQ2 "\+\n";
		print FQ2 "$rev_r2_qual\n";

		$idx += $skip_len;
	}
}
close IN;

close FQ1;
close FQ2;







