open (IN1,"<","$ARGV[0]");#SRR5230786.1987516      147     AMCR01000001.1  82      42      75M     =       31      -126    TGAAA
open (IN2,"<","$ARGV[1]");#chr	pos	penetrance
open (IN3,"<","$ARGV[2]");#genome.gtf AMCR01000001.1  Genbank gene    116     1895    .       +       .       gene_id "OXYTRI_00001"; gbkey "Gene"; gene_biotype "protein_coding"; locus_tag "OXYTRI_00001"; 
open (IN4,"<","$ARGV[3]");#fai
open (OUT1,">","$ARGV[4]");
while (defined ($_=<IN1>)){
	if ($_=~/\s+(AMCR.+?)\s+(\d+)\s+.+?=\s+\d+\s+(\d+)/){
		$chr=$1;
		$a=$2;
		$pos=$a+int($3/2);
		$tmp="$chr\t$pos";
		$MNase{$tmp}++;
		$fra{$chr}++}}
while (defined ($_=<IN2>)){
	if ($_=~/^(AMCR.+?\s+\d+)\s+(.+)/){
		$m6A{$1}=$2}}
while (defined ($_=<IN4>)){
	if ($_=~/^(AMCR.+?)\s+(\d+)/){
		$len{$1}=$2}}
$bs=100;
while (defined ($_=<IN3>)){
	if ($_=~/^(AMCR.+?)\s+.+?\s+gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+/ && $len{$1}>=1000 && $fra{$1}>0){
		$chr=$1;
		$a=$2;
		$b=$3;
		$ori=$4;
		$len=$len{$chr};
		$G="other";
		if ($ori eq "+"){
			$start=$a;
			$TSS=$start;
			foreach $num(1..4){
				$tmp1=$num*$bs;
				$num1=$num-1;
				$tmp2=$num1*$bs;
				if ($TSS<=$tmp1 && $TSS>$tmp2){
					$G="($tmp1,$tmp2]";
					$G{$G}=1}}
			foreach $pos(1..1000){
				$tmp="$chr\t$pos";
				$MNase{$pos}{$G}+=10000*$MNase{$tmp}/$fra{$chr};
				$m6A{$pos}{$G}+=10000*$m6A{$tmp};
				if ($pos==$TSS){
					$TSS{$pos}{$G}++;$sum_TSS++}}}
		elsif ($ori eq "-"){
                        $start=$b;
                        $TSS=$len-$start+1;
                        foreach $num(1..4){
                                $tmp1=$num*$bs;
                                $num1=$num-1;
                                $tmp2=$num1*$bs;
                                if ($TSS<=$tmp1 && $TSS>$tmp2){
                                        $G="($tmp2,$tmp1]";
					$G{$G}=1}}
                        foreach $pos(1..1000){
				$pos=$len-$pos+1;
                                $tmp="$chr\t$pos";
                                $MNase{$pos}{$G}+=10000*$MNase{$tmp}/$fra{$chr};
                                $m6A{$pos}{$G}+=10000*$m6A{$tmp};
                                if ($pos==$TSS){
                                        $TSS{$pos}{$G}++;$sum_TSS++}}}}}
print OUT1 "pos\tgroup\tMNase\tm6A\tTSS\n";
foreach $G(sort {$a cmp $b} keys %G){
	foreach $pos(1..1000){
		if (!exists $MNase{$pos}{$G}){$MNase{$pos}{$G}=0}
		if (!exists $m6A{$pos}{$G}){$m6A{$pos}{$G}=0}
		if (!exists $TSS{$pos}{$G}){$TSS{$pos}{$G}=0}
		$sum_m6A{$G}+=$m6A{$pos}{$G};$sum_MNase{$G}+=$MNase{$pos}{$G}}}
foreach $G(sort {$a cmp $b} keys %G){
        foreach $pos(1..1000){
		$MNase{$pos}{$G}=$MNase{$pos}{$G}/$sum_MNase{$G};
		$TSS{$pos}{$G}=$TSS{$pos}{$G}/$sum_TSS;
		$m6A{$pos}{$G}=$m6A{$pos}{$G}/$sum_m6A{$G};
		print OUT1 "$pos\t$G\t$MNase{$pos}{$G}\t$m6A{$pos}{$G}\t$TSS{$pos}{$G}\n"}}

