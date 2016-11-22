use strict;

if (!(defined($ARGV[0]) && defined($ARGV[1]) && defined($ARGV[2])))
{
	print "Usage: perl check_circular.pl (input fasta file) (reads pair1) (reads pair2) (output fasta with circular seqs) (output fasta without circular seqs)\n";
	exit;
}

my $input_f = $ARGV[0];
my $reads1_f = $ARGV[1];
my $reads2_f = $ARGV[2];
my $out_f = $ARGV[3];
my $out_notc_f = $ARGV[4];

my $CHECK_LEN = 100;
my $PAIR_LEN = 500;
my $IDE_CUTOFF = 97;
my $MIN_LEN = 2000;
my $READS_CUTOFF = 2;

my $terminal_file;
my $cmd;
my %hash1;
my %hash2;
my $tmp;
my $line;
my $i;

$terminal_file = get_terminal($input_f, $CHECK_LEN);

print "Running BLAST...\n";
$cmd = "formatdb -p F -i $terminal_file";
system($cmd);
$cmd = "blastn -query $terminal_file -db $terminal_file -evalue 1e-5 -outfmt 6 -out $terminal_file.m8 -num_threads 8";
system($cmd);
print "BLAST done\n";

unlink("$terminal_file.nin");
unlink("$terminal_file.nhr");
unlink("$terminal_file.nsq");

$tmp = filter("$terminal_file.m8");
%hash1 = %$tmp;

$terminal_file = get_terminal($input_f, $PAIR_LEN);
#split_terminal($terminal_file);

print "Running Bowtie2\n";
$cmd = "bowtie2-build $terminal_file $terminal_file";
system($cmd);
$cmd = "bowtie2 -x $terminal_file -1 $reads1_f -2 $reads2_f -S $terminal_file.sam -p 8";
system($cmd);
print "Bowtie2 done\n";

unlink("$terminal_file.1.bt2");
unlink("$terminal_file.2.bt2");
unlink("$terminal_file.3.bt2");
unlink("$terminal_file.4.bt2");
unlink("$terminal_file.rev.1.bt2");
unlink("$terminal_file.rev.2.bt2");

$tmp = filter2("$terminal_file.sam");
%hash2 = %$tmp;

open(FILE, "<$input_f");
open(OUT1, ">$out_f");
open(OUT2, ">$out_notc_f");
while(defined($line = <FILE>))
{
	if ($line =~ /^>/)
	{
		$tmp = substr($line, 1);
		chomp($tmp);
		if ((exists $hash1{$tmp}) && (exists $hash2{$tmp}))
		{
			$i = 1;
		}
		else
		{
			$i = 0;
		}
	}
	if ($i == 1)
	{
		print OUT1 $line;
	}
	elsif ($i == 0)
	{
		print OUT2 $line;
	}
}
close(FILE);
close(OUT1);
close(OUT2);

unlink($terminal_file);
unlink("$terminal_file.m8");
unlink("$terminal_file.sam");
unlink("formatdb.log");


sub filter
{
	my $m8_s = $_[0];
	my $line;
	my @arr;
	my $tmp1;
	my $tmp2;
	my %hash;
	my $d1;
	my $d2;

	open(FILTER, "<$m8_s");
	while(defined($line = <FILTER>))
	{
		@arr = split(/\t/, $line);
		$tmp1 = substr($arr[0], 0, rindex($arr[0], "_"));
		$tmp2 = substr($arr[1], 0, rindex($arr[1], "_"));
		$d1 = $arr[7] - $arr[6];
		$d2 = $arr[9] - $arr[8];
		$d1 = $d1 * $d2;
		if ($arr[0] ne $arr[1] && $tmp1 eq $tmp2 && $d1 > 0 && $arr[2] >= $IDE_CUTOFF)
		{
			if (!(exists $hash{$tmp1}))
			{
				$hash{$tmp1} = 0;
			}
		}
	}
	close(FILTER);

	return \%hash;
}

sub filter2
{
	my $sam = $_[0];
	my $l1;
	my $l2;
	my @arr1;
	my @arr2;
	my %hash;
	my $tmp1;
	my $tmp2;
	my $n1;
	my $n2;

	open(F1, "<$sam");
	while(!(eof(F1)))
	{
		$l1 = <F1>;
		if ($l1 !~ /^\@/)
		{
			$l2 = <F1>;
			@arr1 = split(/\t/, $l1);
			@arr2 = split(/\t/, $l2);
			if ($arr1[2] eq "*" || $arr2[2] eq "*")
			{
				next;
			}
			$tmp1 = substr($arr1[2], 0, rindex($arr1[2], "_"));
			$n1 = substr($arr1[2], rindex($arr1[2], "_") + 1);
			$tmp2 = substr($arr2[2], 0, rindex($arr2[2], "_"));
			$n2 = substr($arr2[2], rindex($arr2[2], "_") + 1);
			if ($tmp1 eq $tmp2 && $n1 != $n2)
			{
				if (defined($hash{$tmp1}))
				{
					$hash{$tmp1}++;
				}
				else
				{
					$hash{$tmp1} = 1;
				}
			}
		}
	}
	close(F1);
	close(F2);

	foreach $tmp1 (keys %hash)
	{
		if ($hash{$tmp1} < $READS_CUTOFF)
		{
			delete($hash{$tmp1});
		}
	}
	return \%hash;
}

sub get_terminal
{
	my $seq_s = $_[0];
	my $cutlen = $_[1];
	my $out_s = $seq_s . ".term";
	my $tmp_s = $seq_s . ".tmp";
	my $line;
	my $len;
	my $header;
	my $tmp;
	my $i;
	my %lenhash;

	open(SEQFILE, "<$seq_s") || die "Cannot open sequence file $seq_s\n";
	while(defined($line = <SEQFILE>))
	{
		chomp($line);
		if ($line =~ /^>/)
		{
			$header = $line;
			$lenhash{$header} = 0;
		}
		else
		{
			$lenhash{$header} += length($line);
		}
	}

	seek(SEQFILE, 0, 0);
	open(TERMFILE, ">$out_s");
	while(defined($line = <SEQFILE>))
	{
		chomp($line);
		if ($line =~ /^>/)
		{
			if (-e $tmp_s && $lenhash{$header} >= $MIN_LEN)
			{
				close(TMPFILE);
				open(TMPFILE, "<$tmp_s");
				read TMPFILE, $tmp, $cutlen;
				print TERMFILE "$header\_1\n$tmp\n";
				seek(TMPFILE, -1 * $cutlen, 2);
				read TMPFILE, $tmp, $cutlen;
				print TERMFILE "$header\_2\n$tmp\n";
				close(TMPFILE);
			}
			else
			{
				close(TMPFILE);
			}
			$header = $line;
			$len = 0;
			open(TMPFILE, ">$tmp_s");
		}
		else
		{
			print TMPFILE $line;
			$len = $len + length($line);
		}
	}
	if (-e $tmp_s && $lenhash{$header} >= $MIN_LEN)
	{
		close(TMPFILE);
		open(TMPFILE, "<$tmp_s");
		read TMPFILE, $tmp, $cutlen;
		print TERMFILE "$header\_1\n$tmp\n";
		seek(TMPFILE, -1 * $cutlen, 2);
		read TMPFILE, $tmp, $cutlen;
		print TERMFILE "$header\_2\n$tmp\n";
	}
	close(TMPFILE);
	close(SEQFILE);
	close(TERMFILE);

	unlink($tmp_s);

	return($out_s);
}

sub split_terminal
{
	my $term_f = $_[0];
	my $line;
	my $i;
	my $line;

	$i = 0;
	open(TMPFILE, "<$_[0]");
	open(TMPOUT1, ">$_[0].1");
	open(TMPOUT2, ">$_[0].2");
	while(defined($line = <TMPFILE>))
	{
		if ($line =~ /^>[A-Za-z0-9_\-=]+_([12])$/)
		{
			$i = $1;
		}
		if ($i == 1)
		{
			print TMPOUT1 $line;
		}
		elsif ($i == 2)
		{
			print TMPOUT2 $line;
		}
	}
	close(TMPFILE);
	close(TMPOUT1);
	close(TMPOUT2);
}

