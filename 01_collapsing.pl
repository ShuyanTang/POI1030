#!/bin/perl
## collapsing variant alleles into cases in gene level 

use List::Util;

my $outputfile = shift;
my $freq_cut = 0.001;

my $case_file= "POI1030.snv";
my $ctrl_file="HB5000.snv";

my (%sample,%final,%an);

open INFO,"POI1030.info";
while(<INFO>){
	chomp();
	$sample{'Case'}{$_}=();
}
close INFO;

open INFO,"HB5000.info";
while(<INFO>){
	chomp();
	$sample{'Control'}{$_}=();
}
close INFO;

mut_type(\%final,$case_file,'Case');
mut_type(\%final,$ctrl_file,'Control');

foreach my $id(keys %final){
	$an{$id}{'Control'}=(exists $final{$id}{'Control_AN'})? sum0(keys $final{$id}{'Control_AN'})/(scalar keys $final{$id}{'Control_AN'} )/ : 0;
	$an{$id}{'Case'}=(exists $final{$id}{'Case_AN'})? sum0(keys $final{$id}{'Case_AN'})/(scalar keys $final{$id}{'Case_AN'}) : 0;
}

my @title;
open SAVE,">$outputfile";
print SAVE "ID\tGene\t";
foreach my $type ('AN','syn','LoF','Cadd_20','Cadd_10','SPM','REVEL','M-CAP','MetaSVM'){
	foreach my $c ('Case','Control'){
		push @title,"$type\_$c";
	}
}
print SAVE (join "\t",@title)."\n";
foreach my $id (keys %final){
	print SAVE "$id\t$final{$id}{'Gene'}\t";
	
	my @content;
	print SAVE "$an{$id}{'Case'}\t$an{$id}{'Control'}\t";
	foreach my $type ('syn','LoF','Cadd_20','Cadd_10','SPM','REVEL','M-CAP','MetaSVM'){
		foreach my $c ('Case','Control'){
			my $count=0;
			foreach my $case (keys $sample{$c}){
				my $i;
				if (not exists $final{$id}{$c}{$type}{$case}){
					$i=0;
				}elsif($final{$id}{$c}{$type}{$case} >=2){
					$i=2;
				}else{
					$i=1;
				}
				$count+=$i;
			}
			push @content,$count;
		}
	}
	print SAVE (join "\t",@content)."\n";
}
close SAVE;


sub mut_type{
	my ($hash,$file,$gp)=@_;
	my (%final,%type,%content,%name);

	my %mcap;
	open MCAP,"~/Annotation/hg19_mcap_v1_4_D.txt";
	while (<MCAP>){
		chomp();
		my ($chr,$start,$end,$ref,$alt,$score,$sensitivity) = split /\t/;
		$mcap{"$chr\t$start\t$end\t$ref\t$alt"}='D';
	}
	close MCAP;

	my %cadd;

	open CADD,"~/CADD1.6/snv-cadd1.6.annovar";
	while (<CADD>){
		chomp();
		my ($chr,$start,$end,$ref,$alt,$raw,$phred) = split /\t/;
		$cadd{"$chr\t$start\t$end\t$ref\t$alt"}=$phred;
	}
	close CADD;
	
	open FILE,"$file" or die "$file\t$!";
	my $head=<FILE>;chomp($head);
	my $i=0;
	my %heads;
	foreach my $h(split /\t/,$head){
		$heads{$h}=$i++;
	}
	while(my $content = <FILE>){
		$content =~ s/[\r\n]//g;
		my @a = split /\t/,$content;
		my $id=$a[$heads{'Gene ID'}];

		my $function = $a[$heads{'Function'}];
		my $chr=$a[$heads{'Chrs'}];
		my $start=$a[$heads{'Start'}];
		my $end=$a[$heads{'End'}];
		my $ref=$a[$heads{'Ref Allele'}];
		my $alt=$a[$heads{'Alt Allele'}];
		my $snv = join "\t",($chr,$start,$end,$ref,$alt);

		my $an = $a[$heads{'Case_AN'}];	
		$hash->{$id}{"$gp\_AN"}{$an}=();

		my $aa=($a[$heads{'AA_Change_Position'}] =~ /\d+/)?$a[$heads{'AA_Change_Position'}]:0 ;
		my $length=($a[$heads{'Protein_Length'}] =~ /\d+/)?$a[$heads{'Protein_Length'}]:1 ;
		my $change=$aa/$length;

		foreach my $case( split ",",$a[$heads{'Case_het'}]){
			next if $case =~ /^\s*$/;
			$hash->{$id}{$gp}{'LoF'}{$case}++ if ($function =~  /stop_gain|frameshift|splice site|start_lost/ and $change <0.98 );
			$hash->{$id}{$gp}{'syn'}{$case}++ if ($function =~ /synonymous/);
			if($function =~ /missense/){
				$hash->{$id}{$gp}{'SPM'}{$case}++ if ($a[$heads{'SIFT_Pred'}] =~ /D/ and  $a[$heads{'PolyPhen-2_Pred'}] =~ /D|P/ and  $a[$heads{'MutationTaster_Pred'}] =~ /D/ );
				$hash->{$id}{$gp}{'REVEL'}{$case}++ if ($a[$heads{'REVEL'}]!~ /^\s*$/ and $a[$heads{'REVEL'}] >0.75);
				$hash->{$id}{$gp}{'M-CAP'}{$case}++ if (exists $mcap{$snv});
				$hash->{$id}{$gp}{'MetaSVM'}{$case}++ if ($a[$heads{'MetaSVM_pred'}] =~ /D/ );
				$hash->{$id}{$gp}{'Cadd_20'}{$case}++ if (exists $cadd{$snv} and $cadd{$snv} >= 20 ); 
				$hash->{$id}{$gp}{'Cadd_10'}{$case}++ if (exists $cadd{$snv} and $cadd{$snv} >= 10 );
			}
		}
		foreach my $case( split ",",$a[$heads{'Case_hom'}]){
			next if $case =~ /^\s*$/;
			$hash->{$id}{$gp}{'LoF'}{$case}+=2 if ($function =~  /stop_gain|frameshift|splice site|start_lost/ and $change <0.98 );
			$hash->{$id}{$gp}{'syn'}{$case}+=2 if ($function =~ /synonymous/);
			if($function =~ /missense/){
				$hash->{$id}{$gp}{'SPM'}{$case}+=2 if ($a[$heads{'SIFT_Pred'}] =~ /D/ and  $a[$heads{'PolyPhen-2_Pred'}] =~ /D|P/ and  $a[$heads{'MutationTaster_Pred'}] =~ /D/ );
				$hash->{$id}{$gp}{'REVEL'}{$case}+=2 if ($a[$heads{'REVEL'}]!~ /^\s*$/ and$a[$heads{'REVEL'}] >0.75);
				$hash->{$id}{$gp}{'M-CAP'}{$case}+=2 if (exists $mcap{$snv});
				$hash->{$id}{$gp}{'MetaSVM'}{$case}+=2 if ($a[$heads{'MetaSVM_pred'}] =~ /D/ );
				$hash->{$id}{$gp}{'Cadd_20'}{$case}+=2 if (exists $cadd{$snv} and $cadd{$snv} >= 20 ); 
				$hash->{$id}{$gp}{'Cadd_10'}{$case}+=2 if (exists $cadd{$snv} and $cadd{$snv} >= 10 );
			}
		}
	}
	close FILE;
}