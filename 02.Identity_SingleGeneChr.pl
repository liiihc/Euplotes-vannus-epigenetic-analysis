open (IN1,"<","$ARGV[0]");#fai
open (IN2,"<","$ARGV[1]");#gtf
open (OUT1,">","$ARGV[2]");#SG.bed
open (OUT2,">","$ARGV[3]");#MG.bed
open (OUT3,">","$ARGV[4]");#NG.bed
while (defined ($_=<IN2>)){
	if ($_=~/^(.+?)\s+.+?\s+gene\s+\d+\s+\d+\s+.+?\s+(.)\s+/){
		$hash{$1}++;
		$ori{$1}=$2}}
while (defined ($_=<IN1>)){
	if ($_=~/^(.+?)\s+(\d+)\s+/){
		$chr=$1;
		$len=$2;
		if ($hash{$chr}==1){
			print OUT1 "$chr\t1\t$len\t.\t$chr\t$ori{$chr}\n"}
		elsif ($hash{$chr}>1){
                        print OUT2 "$chr\t1\t$len\t.\t$chr\t$ori{$chr}\n"}
		elsif(!exists $hash{$chr}){
                        print OUT3 "$chr\t1\t$len\t.\t$chr\t+\n"}}}
