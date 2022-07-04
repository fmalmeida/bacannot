#!/usr/bin/perl

# splitgb.pl
#
# 160804
# v.0.1
# ksaitoh
#
# Split individual sequences of concatenated GenBank format DNA sequence file
# Save split sequences as their accession number plus ".gb"
#
# usage > perl splitgb.pl infile

$infile = shift @ARGV;

$begin = "LOCUS";
$end = "//";
@entry = ();
$ext = ".gbk";

open(IN, "<$infile") or die "Failed to open $infile\n";

$no = 0;

while ($line = <IN>)  {
    chomp ($line);
    $line =~ s/\r//;
    if ($line =~ /^$begin/) {
        ++$no;
        @entry = split /\s+/, $line;
        $outfile = $entry[1].$ext;
        open(OUT, ">$outfile") or die "Failed to open $outfile\n";
        print OUT "$line\n";
        do {
                $line = <IN> or die "End of file\n";
                chomp ($line);
                $line =~ s/\r//;
                print OUT "$line\n";
        } while ($line ne $end);
        close OUT;
        @entry = ();
    }
}

print "$no seqs saved separately.\n";

close (IN);