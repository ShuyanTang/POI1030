#!/bin/perl
## generate geneset files for further wilcoxon test
use warnings;
use strict;
use feature qw(say);


our $BurdenFile="~/POI1030/all.burden_result.tsv";
our $GenesetFile="~/POI1030/Geneset.tsv";
our $OutputDir="~/POI1030/GeneSet_wilcoxon";mkdir $OutputDir;
our (%geneset,%all,%pathway,%ispath);

## get genes in gene sets
open GENESET, $GenesetFile;
while(<GENESET>){
	chomp();
	my ($set,$g,$ids)=split /\t/;
	map{$geneset{$_}{$set}=()}(split /,/,$ids);
}
close GENESET;

## get gene-level P-values
open BURDEN,"$BurdenFile" or die "$BurdenFile $!";
my $head=<BURDEN>;chomp($head);
my $i=0;
my %heads;
foreach my $h(split /\t/,$head){
	$heads{$h}=$i++;
}
while(<BURDEN>){
	chomp();
	my @a=split /\t/;
	my $id=$a[$heads{'ID'}];
	my $g=$a[$heads{'Gene'}];
	my $p= $a[$heads{'LoF.Control.P'}];
	if (exists $geneset{$id}){
		foreach my $p(keys $geneset{$id}){
			$pathway{$p}{$id}=();
		}
	}
	next if $p eq 'NA';
	$all{$id}{'LoF.p'}=$p;
	$all{$id}{'Gene'}=$g;
}
close BURDEN;

## get backgroud genes
open GENEPOS,"~/analysis/gene_pos_sorted.txt";
<GENEPOS>;
my $j=0;
my (%rank,%bg,%pos);
while(<GENEPOS>){
	chomp();
	my ($chr,$start,$end,$length,$gene,$id) = split /\t/;
	next if not exists $all{$id};
	$pos{$id}=join ",",($chr,$start,$end,$j);
	$rank{$j}=$id;
	$j++;
}
close GENEPOS;

foreach my $id(keys %all){
	next if (not exists $pos{$id});
	my ($chr,$start,$end,$n) = split /,/,$pos{$id};
	foreach my $m($n-25..$n+25){
		my $id2 = $rank{$m};
		next if (not defined $id2 or $id eq $id2); 
		my ($chr2,$start2,$end2) = split /,/,$pos{$id2};
		$bg{$id}{$id2}=() if ($end2 < $start or $end < $start2);
	}
}


my @title =('Gene','ID','Group','LoF.p');

foreach my $path (keys %pathway){
	say "Wilcoxon\tGeneSet:$path";
	my %final;
	foreach my $id(keys $pathway{$path}){
		foreach my $id2(keys $bg{$id}){
			$final{$id2}{'Group'}='Background';
			$final{$id2}{'Gene'}=$all{$id2}{'Gene'};
			$final{$id2}{'ID'}=$id2;
			$final{$id2}{'LoF'}= $all{$id2}{'LoF.p'};
		}
	}
	foreach my $id(keys $pathway{$path}){
		$final{$id}{'Group'}='Gene set';
		$final{$id}{'Gene'}=$all{$id}{'Gene'};
		$final{$id}{'ID'}=$id;
		foreach my $key('LoF'){
			$final{$id}{$key}= $all{$id}{$key};
		}
	}
	open WIL,">$dir/geneset.$path.wilcoxon";
	print WIL (join "\t",@title)."\n";
	foreach my $id (keys %final){
		my @content = map{$final{$id}{$_}}(@title);
		print WIL (join "\t",@content)."\n";
	}
	close WIL;
}

