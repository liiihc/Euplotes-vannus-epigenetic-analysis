open (IN1,"<","$ARGV[0]");#SRR11906179.39995232    163     chr_001 45      1       34M     =       63      52      CCCAACCCCAAACCCCCACCCCAACCCCAACCCC      FFFFF:FFF,F:FFFF,:,FFF,:,FFF:,:FFF      AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:16A17      YS:i:-4 YT:Z:CP
open (IN2,"<","$ARGV[1]");#chr  len
open (IN3,"<","$ARGV[1]");
open (OUT1,">","$ARGV[2]");
while (defined ($_=<IN1>)){
	if ($_=~/^.+?\s+\d+\s+(.+?)\s+/){
		$count{$1}++;
		$total++}}
while (defined ($_=<IN2>)){
	if ($_=~/^(.+?)\s+(\d+)/){
		$chr=$1;
		$len=$2;
		$ratio=(1000000*$count{$chr})/($total*$len);
		$ratio{$chr}=$ratio;
		$sum+=$ratio;
		$num++}}
$ave=$sum/$num;
while (defined ($_=<IN3>)){
	if ($_=~/^(.+?)\s+(\d+)/){
		$chr=$1;
		$val=$ratio{$chr}/$ave;
		print OUT1 "$chr\t$val\n"}}
