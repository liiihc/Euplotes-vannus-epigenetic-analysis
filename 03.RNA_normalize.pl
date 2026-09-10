open (IN1,"<","$ARGV[0]");#Geneid  Chr     Start   End     Strand  Length  Eup.RNA_rmdup.bam
open (IN2,"<","$ARGV[1]");#chr copynumber 
open (IN3,"<","$ARGV[2]");
open (OUT1,">","$ARGV[3]");
while (defined ($_=<IN3>)){
	if ($_=~/^(.+)$/){
		$hash{$1}=1}}
while (defined ($_=<IN2>)){
	if ($_=~/^(.+?)\s+(.+)/ && exists $hash{$1}){
		$CN{$1}=$2}}
while (defined ($_=<IN1>)){
	if ($_=~/^EVAN/){
		@arr=split("\t",$_);
		$gene=$arr[0];
		$chr=$arr[1];
		$start=$arr[2];
		$end=$arr[3];
		$ori=$arr[4];
		$len=$arr[5];
		$count=$arr[6];
		$CN=$CN{$chr};
		if ($CN>0 ){$exp=1000000*$count/($len*$CN);
		print OUT1 "$chr\t$start\t$end\t$gene\t$exp\t$ori\n"}}}
