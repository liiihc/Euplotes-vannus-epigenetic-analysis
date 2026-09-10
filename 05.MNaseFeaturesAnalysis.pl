#!/usr/bin/perl
# Usage: perl script.pl MNase.sam dyad.txt genome.fasta features.bed output_prefix

use strict;
use warnings;

open (IN1, "<", "$ARGV[0]") or die "Cannot open MNase file: $!";
open (IN2, "<", "$ARGV[1]") or die "Cannot open Dyad file: $!";
open (IN3, "<", "$ARGV[2]") or die "Cannot open Genome file: $!";
open (IN4, "<", "$ARGV[3]") or die "Cannot open BED file: $!";

my $prefix = $ARGV[4] || "result";
open (OUT1, ">", "$prefix.density.txt") or die "Cannot write to $prefix.density.txt: $!";
open (OUT_LEN, ">", "$prefix.at_len.txt") or die "Cannot write to $prefix.at_len.txt: $!";
open (OUT_META, ">", "$prefix.at_meta.txt") or die "Cannot write to $prefix.at_meta.txt: $!";
open (OUT_BIN, ">", "$prefix.at_bin.txt") or die "Cannot write to $prefix.at_bin.txt: $!";

# Global hash declarations
my (%nuc, %fra, %seq, %ori_map);
my %at_tract_len_dist; 
my %at_meta_nuc;       
my %at_meta_count;     
my %at_bin_count;      
my (%NUC_plot, %GC_plot, %AT_plot);

# 1. Process MNase-seq data and calculate nucleosome occupancy
while (<IN1>){
    next if $_ =~ /^@/; # Skip SAM header lines
    if ($_ =~ /^.+?\s+.+?\s+(.+?)\s+(\d+)\s+.+?\s+=\s+\d+\s+(\d+)/){
        my ($chr, $start, $len) = ($1, $2, $3);
        if ($len >= 130 && $len <= 170) {
            my $end = $start + $len;
            foreach my $pos ($start..$end){
                $nuc{"$chr\t$pos"}++;
            }
            $fra{$chr}++;
        }
    }
}

# 2. Read Genome Fasta sequences
my $name;
while (<IN3>){
    chomp;
    if ($_ =~ />(\S+)/){
        $name = $1;
    } elsif (defined $name && length($_) > 0) {
        $seq{$name} .= uc($_); 
    }
}

# 3. Read orientation/strand information
while (<IN4>){
    if ($_ =~ /^(.+?)\s+(\d+)\s+(\d+)\s+.+?\s+.+?\s+(.)/){
        $ori_map{"$1\t$2"} = $4; 
    }
}

# --- AT Tract Analysis ---
print "Analyzing AT tracts and nucleosome correlation...\n";
foreach my $chr (keys %seq) {
    next unless (defined $fra{$chr} && $fra{$chr} > 0);
    my $genome_seq = $seq{$chr};
    my $chr_len = length($genome_seq);
    next if $chr_len == 0;
    
    # Identify poly(A) or poly(T) tracts with length >= 5
    while ($genome_seq =~ /([A]{5,}|[T]{5,})/g) {
        my $tract = $1;
        my $len = length($tract);
        my $end_pos = pos($genome_seq);
        my $start_pos = $end_pos - $len + 1;
        my $mid_pos = int(($start_pos + $end_pos) / 2);
        
        # A. Length statistics
        $at_tract_len_dist{$len}++;
        
        # B. Chromosome 20-bin distribution analysis
        my $bin = int(($mid_pos - 1) / $chr_len * 20) + 1;
        $bin = 20 if $bin > 20;
        $bin = 1 if $bin < 1;
        $at_bin_count{$bin}++;

        # C. Meta-analysis around AT tract center
        for my $offset (-200 .. 200) {
            my $abs_pos = $mid_pos + $offset;
            if (exists $nuc{"$chr\t$abs_pos"}) {
                $at_meta_nuc{$offset} += ($nuc{"$chr\t$abs_pos"} / $fra{$chr});
            }
            $at_meta_count{$offset}++;
        }
    }
}

# 4. Analysis based on input dyad/feature sites
seek(IN2, 0, 0); 
while (<IN2>){
    if ($_ =~ /^(.+?)\s+(\d+)\s+(\d+)/){
        my ($chr, $s, $e) = ($1, $2, $3);
        next unless (defined $seq{$chr} && length($seq{$chr}) > 0);
        
        my $mid = int(($s + $e) / 2);
        my $ori = $ori_map{"$chr\t$s"} || "+";
        
        my @chars = split("", $seq{$chr});
        foreach my $a (0..600){
            my $b = ($ori eq "+") ? ($mid + $a) : ($mid - $a);
            my $pos_key = "$chr\t$b";
            
            if (defined $fra{$chr} && $fra{$chr} > 0 && exists $nuc{$pos_key}){
                $NUC_plot{$a} += ($nuc{$pos_key} / $fra{$chr});
            }
            
            my $idx = $b - 1;
            if ($idx >= 0 && $idx < scalar @chars) {
                my $base = $chars[$idx];
                if ($base =~ /[GC]/) { $GC_plot{$a}++; }
                elsif ($base =~ /[AT]/) { $AT_plot{$a}++; }
            }
        }
    }
}

# --- Final Output ---

# Output 1: Occupancy and GC around target sites
print OUT1 "pos\tnuc\tGC\n";
foreach my $a (0..600){
    my $gc = $GC_plot{$a} || 0;
    my $at = $AT_plot{$a} || 0;
    my $gc_ratio = ($gc + $at > 0) ? ($gc / ($gc + $at)) : 0;
    my $nuc_val = $NUC_plot{$a} || 0;
    print OUT1 "$a\t$nuc_val\t$gc_ratio\n";
}

# Output 2: AT tract length distribution
print OUT_LEN "Length\tFrequency\n";
foreach my $l (sort {$a <=> $b} keys %at_tract_len_dist) {
    print OUT_LEN "$l\t$at_tract_len_dist{$l}\n";
}

# Output 3: Nucleosome depletion around AT tracts (Meta-analysis)
print OUT_META "Relative_Pos\tMean_Nuc_Occupancy\n";
foreach my $off (sort {$a <=> $b} keys %at_meta_count) {
    my $meta_nuc = $at_meta_nuc{$off} || 0;
    my $meta_cnt = $at_meta_count{$off} || 0;
    my $mean_nuc = ($meta_cnt > 0) ? ($meta_nuc / $meta_cnt) : 0;
    print OUT_META "$off\t$mean_nuc\n";
}

# Output 4: AT tract distribution across 20 chromosome bins
print OUT_BIN "Bin\tAT_Tract_Count\n";
foreach my $b (1..20) {
    my $count = $at_bin_count{$b} || 0;
    print OUT_BIN "$b\t$count\n";
}

close IN1; close IN2; close IN3; close IN4;
close OUT1; close OUT_LEN; close OUT_META; close OUT_BIN;
