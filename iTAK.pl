#!/usr/bin/perl

=head1
[Sep-08-2014][v1.5]: combination rules for plantTFDB and plnTFDB, new classification system for future rule update
[Jan-19-2014][v1.4]: new category system for plant protein kinase, from RKD and HMM build by Shiu et al. 2012 
		     update the hmmscan to version 3.1, 2x faster than hmm 3.0
[Jun-22-2013][v1.3]: update some small bugs, Pfam V27, WAK and WAKL family, a switch for transponsase filter
[Aug-26-2011][v1.2]: report unusual sequences
[Jun-03-2011][v1.1]: remove unsignificiant domain using GA score
[Dec-14-2010][v1.0]: first stable version 
=cut

use strict;
use warnings;
use Cwd;
use FindBin;
use File::Basename;
use Bio::SeqIO;
use IO::File;
use Bio::SearchIO;
use Getopt::Long;

my $usage = q/
UPDATE: Jan-19-2014
VERSION: v1.5
USAGE: 
	perl iTAK.pl [parameters]

	-i  [String]	Name of input sequence file in FASTA format (required)
	-s  [String]	Type of input sequences ('p' for protein | 'n' for 
			nucleotide). (default = p)
	-m  [String]	Type of analysis ('t' for TF identification | 'p' for 
			protein kinase identification |'b' for both). (default 
			= b)
	-a  [Integer]	number of CPUs used for hmmscan. (default = 1)
	-o  [String]	Name of the output directory. ( default = 'input file 
			name' + '_output')

/;

#################################################################
# Input Parameters                                              #
#################################################################
my ($help, $input_seq, $seq_format, $mode, $cpus, $output_fd, $debug);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$input_seq,
	"s=s"	=> \$seq_format,
	"m=s"	=> \$mode,
	"a=i"	=> \$cpus,
	"o=s"	=> \$output_fd,
	"d"	=> \$debug
);

die $usage if $help;
die $usage unless $input_seq;
$seq_format ||= "p";
$mode ||= "b";
$cpus ||= 1;
$output_fd ||= $input_seq."_output";
$debug ||= 0;

#################################################################
# Check input files and parameters 				#
#################################################################

#################################################################
# check mode parameters						#
# 1. mode t, just identify transcription factors		#
# 2. mode p, just identify protein kinases			#
# 3. mode b, identify both					#
#################################################################
unless ($mode =~ m/^t|p|b$/i) { die $usage."\nPlease input correct mode\n"; }

#################################################################
# set transposase						#
#################################################################
my %transposase = (
	#"PF01498" => 23.4,	# HTH_Tnp_Tc3_2 	Transposase
	#"PF11427" => 20.8,	# HTH_Tnp_Tc3_1 	Tc3 transposase
	#"PF03004" => 21.6,	# Transposase_24 	Plant transposase (Ptta/En/Spm family)
);

#################################################################
# check input_seq files						#
# 1. seq_format eq n means input seqs are DNA			#
# 2. seq_format eq p (default) means input seqs are Proteins	#
# if you input DNA, it can be converted to protein, if your	#
# input seq and seq_format are not match, will die		#
#################################################################
my %seq_hash;

my ($input_seq_hash, $converted_seq_hash) = seq_to_hash($input_seq);

if ($seq_format eq "N" || $seq_format eq "n")
{
	if ( scalar(keys(%$input_seq_hash)) > 0 ) 
	{ 
		print "\nYou input contains protein or unusual sequences!\n\n";

		my %unusual_list;
		foreach my $id (sort keys %$input_seq_hash) {
			print ">".$id."\n".$$input_seq_hash{$id}."\n";
			my $pid = $id; $pid =~ s/_\d$//ig;
			$unusual_list{$pid} = 1;
                }

		print "\nUnusual sequence IDs:\n\n";
		foreach my $id (sort keys %unusual_list) {
			print $id."\t";
		}
		print "\n\n";
		exit;
	}
	else { %seq_hash = %$converted_seq_hash; }
}
elsif ($seq_format eq "P" || $seq_format eq "p")
{
	if ( scalar(keys(%$converted_seq_hash)) > 0) 
	{
		print "\nYour input contains DNA or unusual sequences!\n\n";
		
		my %unusual_list;
		foreach my $id (sort keys %$converted_seq_hash) { 
			print ">".$id."\n".$$converted_seq_hash{$id}."\n";
			my $pid = $id; $pid =~ s/_\d$//ig;
			$unusual_list{$pid} = 1;
		}

		print "\nUnusual sequence IDs:\n\n";
		foreach my $id (sort keys %unusual_list) {
			print $id."\t";
		}
		print "\n\n";
		exit;
	}
	else { %seq_hash = %$input_seq_hash; }
}
else
{
	print "[ERR]paramter -s $seq_format\n$usage"; 
	exit;
}

#################################################################
# create temp folder and then store seq to temp			#
#################################################################
my $temp_dir = $input_seq."_itak_temp";
if (-e $temp_dir) { die "[ERR]Temp folder $temp_dir exist, please delete it before running iTAK\n";  }
else { mkdir($temp_dir); }

my $protein_seq = $temp_dir."/protein_seq";
my $pfh = IO::File->new(">".$protein_seq) || die $!;
foreach my $pid (sort keys %seq_hash)
{
	print $pfh ">".$pid."\n".$seq_hash{$pid}."\n";
}
$pfh->close;

#################################################################
# check and create output folder				#
#################################################################
my $output_dir = $output_fd;
if (-e $output_dir ) { die "[ERR]output folder $output_dir exist, please delete it before running iTAK\n"; }
else { mkdir($output_dir); }

#################################################################
#------------------Varialbes that Default Set-------------------#
#################################################################
my $bin_dir = ${FindBin::RealBin}."/bin";
my $dbs_dir = ${FindBin::RealBin}."/database";
unless (-e $bin_dir) { die "[ERR]bin directory does not exist.\n$bin_dir\n"; }
unless (-e $dbs_dir) { die "[ERR]database directory does not exist.\n $dbs_dir\n"; }

#################################################################
#--------Default Variables for Transcription Factors------------#
#################################################################

#################################################################
# Database for identify transcription factors			#
# This database include Pfam-A and customized HMM Profiles	#
#################################################################
my $hmm3_db = $dbs_dir."/TFHMM_3.hmm";

#################################################################
# Rules for identify Transcription Factors			#
#################################################################
my $tf_rule = $dbs_dir."/TF_Rule";
unless (-s $tf_rule) { die "[ERR]Rules for Transcription Factors domains do not exist.\n"; }

# put TFs and TRs family to hash
my $tf_cat = $dbs_dir."/TF_family";
my %tf_family_cat = tf_family_cat_to_hash($tf_cat);

#################################################################
# GA score table for Transcription Factors			#
#################################################################
my $ga_table = $dbs_dir."/GA_table";	# GA score and Desc from PfamA, dependend one rules, prepared file
unless (-s $ga_table) { die "GA score and Desc do not exist.\n"; }
my ($ga_hash, $desc_hash) = ga_to_hash($ga_table);	

#################################################################
#-------------Default Variables for Protein Kinases-------------#
#################################################################

#################################################################
# Self-build Protein Kinases HMM profile for classification 	#
# base on PlantsP sequences, choose below HMM profile		#
# PlantsPHMM3.hmm ; Hmmer3 format				#
#################################################################
my $plantp_hmm_3  = $dbs_dir."/PlantsPHMM3_89.hmm";
my $rkd_hmm_3 = $dbs_dir."/PlantsPHMM3_89.hmm";
my $shiu_hmm_3 = $dbs_dir."/Plant_Pkinase_fam.hmm";

#################################################################
# protein kinases Cat ID and It's desc to hash.			#
#################################################################
my $pk_desc = $dbs_dir."/PK_class_desc";
unless (-s $pk_desc) { die "protein kinase descriptions do not exist.\n"; }
my $pkid_des = pk_to_hash($pk_desc);

# ======== parameters end ========

#################################################################
# main								#
#################################################################

# Step 1 
# using all input seqs to perform hmmscan against PfamA + Self-build HMM
# parse the hmmscan results 
my $tmp_hmm_result = $temp_dir."/temp_hmm3_result"; 
my $hmmscan_command = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_hmm_result $hmm3_db $protein_seq";
print $hmmscan_command if $debug == 1;
unless (-s $tmp_hmm_result){ system($hmmscan_command) && die "Error at hmmscan command: $hmmscan_command\n"; }
my ($all_hits, $all_detail) = parse_hmmscan_result($tmp_hmm_result);

# Step 2
# the parsed hmmscan result was used in TFs prediction
my $tmp_rule = $tf_rule;
my ($tmp_align, $tmp_family) = identify_domain($all_hits, $all_detail, $tmp_rule, \%transposase);

# Step3
# output the TFs prediction result 
if ($mode =~ m/^t|b$/i)
{
	if ($tmp_family)
    	{
		my $tf_family = $output_dir."/".$input_seq."_tf_family";
		my $tf_align  = $output_dir."/".$input_seq."_tf_align";
		my $tf_seq    = $output_dir."/".$input_seq."_tf_seq";

		my %pid; my %tf;

		my @num = split(/\n/, $tmp_family);
		foreach my $line (@num)
		{
			my @a = split(/\t/, $line);
			$pid{$a[0]} = 1;
			$tf{$a[1]} = 1;
		}

		my $aafh = IO::File->new(">".$tf_family) || die "Can not open transcription factor classification file: $tf_family $!\n";
		foreach my $tf_line (@num)
		{
			my @tf_m = split(/\t/, $tf_line);
			print $aafh $tf_m[0]."\t".$tf_m[1]."\t".$tf_family_cat{$tf_m[1]}."\n";
		}
		$aafh->close;

		my $bbfh = IO::File->new(">".$tf_align) || die "can not open transcription factor alignment file: $tf_align $!\n";
		print $bbfh $tmp_align;
		$bbfh->close;

		my $ccfh = IO::File->new(">".$tf_seq) || die "can not open transcription factor alignment file: $tf_seq $!\n";
		my @s = split(/\n/, $tmp_family);
		foreach my $ccc (@s)
		{
			my @ss = split(/\t/, $ccc, 2);
			print $ccfh ">".$ss[0]."\n".$seq_hash{$ss[0]}."\n";
		}
		$ccfh->close;

		print scalar(keys(%pid))." transcription factors were identified.\n";
	}
	else
	{
		print "No transcription factor was identified.\n";
    	}
}

#########################################################
# After Protein Kinase Prediction, Two Char Produced	#
# 1. temp_all_hmmscan_family				#
# 2. temp_all_hmmscan_domain				#
#########################################################

if ($mode =~ m/^p|b$/i)
{
	# Step 3.1 produce protein kinase sequence
	# pkinase_id has the protein ID with high GA score in Pfam Kinase domain
	# key: gene_id, value: 1; 
	my %pkinase_id = cutoff($all_hits, \%transposase);

	# pkinase_aln
	# key: gene_id, value: align detail for gene
	my %pkinase_aln = aln_to_hash($all_detail);

	my $protein_kinase_seq = $output_dir."/".$input_seq."_pkseq";

	my $pk_seq_num = scalar(keys(%pkinase_id));

	if ($pk_seq_num > 0)
	{
		# save protein kinase sequence to fasta file (with protein kinases domain)
		my $pk_fh = IO::File->new(">".$protein_kinase_seq) || die "Can not open protein kinase sequence file: $protein_kinase_seq\n";
		foreach my $pid (sort keys %pkinase_id)
		{
			if (defined $seq_hash{$pid}) 
			{
				print $pk_fh ">".$pid."\n".$seq_hash{$pid}."\n";
			}
			else
			{
				print "Error! no sequences match to this id $pid\n";
			}
		}
		$pk_fh->close;

		# Step 3.3 Get Protein Kinases Classification using hmmscan
		# Step 3.3.1 hmmscan
		my $tmp_plantsp_hmm_result = "$temp_dir/temp_plantsp_hmm_result";
		my $tmp_rkd_hmm_result 	   = "$temp_dir/temp_rkd_hmm_result";
		my $tmp_shiu_hmm_result    = "$temp_dir/temp_shiu_hmm_result";

		my $plantsp_hmmscan_command = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_plantsp_hmm_result $plantp_hmm_3 $protein_kinase_seq";
		my $rkd_hmmscan_command     = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_rkd_hmm_result     $rkd_hmm_3    $protein_kinase_seq";
		my $shiu_hmmscan_command    = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_shiu_hmm_result    $shiu_hmm_3   $protein_kinase_seq";

		print $plantsp_hmmscan_command if $debug == 1;
		print $rkd_hmmscan_command if $debug == 1;
		print $shiu_hmmscan_command if $debug == 1;

		unless (-s $tmp_plantsp_hmm_result) { system($plantsp_hmmscan_command) && die "Error at hmmscan command: $plantsp_hmmscan_command\n"; }
		#unless (-s $tmp_rkd_hmm_result)  { system($rkd_hmmscan_command) && die "Error at hmmscan command: $rkd_hmmscan_command\n"; }	
		unless (-s $tmp_shiu_hmm_result) { system($shiu_hmmscan_command) && die "Error at hmmscan command: $shiu_hmmscan_command\n"; }
		
		# Step 3.3.2 parse hmmscan result
		my ($hmm3_simple_result, $hmm3_simple_align) = parse_hmmscan_result($tmp_plantsp_hmm_result);
		#my ($hmm3_rkd_result, $hmm3_rkd_align) = parse_hmmscan_result($tmp_rkd_hmm_result);
		my ($hmm3_shiu_result, $hmm3_shiu_align) = parse_hmmscan_result($tmp_shiu_hmm_result);

		# Step 3.3.3 get classification info base simple hmminfo, means get best hits of simple result, then add annotation and output alignment file
		my %pkinases_cat = get_classification($hmm3_simple_result);
		my %shiu_cat = get_classification($hmm3_shiu_result);
		#my %rkd_cat = get_classification($hmm3_rkd_result);

		#################################################
		# output PlantsP alignment and classification 	#
		#################################################
		my $protein_kinase_aln1 = $output_dir."/".$input_seq."_pkaln1";
		my $protein_kinase_cat1 = $output_dir."/".$input_seq."_pkcat1";
	
		my $ca_fh1 = IO::File->new(">".$protein_kinase_cat1) || die "Can not open protein kinase sequence file: $protein_kinase_cat1 $!\n";
		my $al_fh1 = IO::File->new(">".$protein_kinase_aln1) || die "Can not open protein kinase sequence file: $protein_kinase_aln1 $!\n";

		foreach my $pid (sort keys %pkinases_cat)
		{
			if (defined $pkinase_aln{$pid})
			{
				print $al_fh1 $pkinase_aln{$pid};
			}
			else
			{
				die "Error! Do not have alignments in hmm3 parsed result\n";
			}
			delete $pkinase_id{$pid};	
		}

		foreach my $pid (sort keys %pkinase_id)
		{
			print $ca_fh1 $pid."\tPPC:1.Other\n";

			if (defined $pkinase_aln{$pid})
			{
				print $al_fh1 $pkinase_aln{$pid};
			}
			else
			{
				die "Error! Do not have alignments in hmm3 parsed result\n";
			}
		}
		$al_fh1->close;

		# Step 3.3.4 find WNK1 domain in 4.1.5 cat;
		# $input, $hmm, $score, $cat_id, $new_cat_id
		%pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/wnk1_hmm_domain/WNK1_hmm" , "30" ,"PPC:4.1.5", "PPC:4.1.5.1");
		$$pkid_des{"PPC:4.1.5.1"} = "WNK like kinase - with no lysine kinase";

		%pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/mak_hmm_domain/MAK_hmm" , "460.15" ,"PPC:4.5.1", "PPC:4.5.1.1");
		$$pkid_des{"PPC:4.5.1.1"} = "Male grem cell-associated kinase (mak)";
	
		foreach my $pid (sort keys %pkinases_cat) {
			print $ca_fh1 $pid."\t".$pkinases_cat{$pid}."\t".$$pkid_des{$pkinases_cat{$pid}}."\n";
		}
		$ca_fh1->close;


		#################################################
                # output Shiu alignment and classification   	#
                #################################################
                my $protein_kinase_aln2 = $output_dir."/".$input_seq."_pkaln2";
                my $protein_kinase_cat2 = $output_dir."/".$input_seq."_pkcat2";

                my $ca_fh2 = IO::File->new(">".$protein_kinase_cat2) || die "Can not open protein kinase sequence file: $protein_kinase_cat2 $!\n";
                my $al_fh2 = IO::File->new(">".$protein_kinase_aln2) || die "Can not open protein kinase sequence file: $protein_kinase_aln2 $!\n";
		foreach my $pid (sort keys %shiu_cat) { print $ca_fh2 $pid."\t".$shiu_cat{$pid}."\n"; }
		print $al_fh2 $hmm3_shiu_align;
		$ca_fh2->close;
		$al_fh2->close;

		#################################################
		# output rkd alignment and classification       #
		#################################################
		my $protein_kinase_aln3 = $output_dir."/".$input_seq."_pkaln3";
		my $protein_kinase_cat3 = $output_dir."/".$input_seq."_pkcat3";

		#my $ca_fh3 = IO::File->new(">".$protein_kinase_cat3) || die "Can not open protein kinase sequence file: $protein_kinase_cat3 $!\n";
		#my $al_fh3 = IO::File->new(">".$protein_kinase_aln3) || die "Can not open protein kinase sequence file: $protein_kinase_aln3 $!\n";
		#$ca_fh3->close;
		#$al_fh3->close;

		# report protein kinase number
		print $pk_seq_num." protein kinases were identified.\n";
    	}
	else
	{
		print "No protein kinase was identified.\n";
	}
}

#################################################################
# delete temp folder						#
#################################################################
my $remove_temp_cmd = "rm -rf $temp_dir";
print "remove temp folder: ".$remove_temp_cmd."\n";
#system($remove_temp_cmd) && die "Error at $remove_temp_cmd\n";

#################################################################
# finished							#
#################################################################


#################################################################
# kentnf : subroutines						#
#################################################################

=head2 ga_to_hash

 Function: reading GA score table to hash. Because of PfamA file is too large(about 2GB) and parse this file will spend a lot of time, so put GA score to file first, and then take this info directly can save lot of time.

 Input: GA score table file name, this file in bin folder. Format is : DomainID \t GA Score \n

 Output: GA score hash, key is domain id, value is GA score.
=cut
sub ga_to_hash
{
	my $table = shift; 
	my @m; 
	my %ga;
	my $tfh = IO::File->new($table) || die "Can not open GA score table file: $table $!\n";
	while(<$tfh>)
	{
		chomp;
		my @m = split(/\t/, $_);
		$m[0] =~ s/\..*//;
		$ga{$m[0]} = $m[1];
	}
	$tfh->close;
	return (\%ga);
}

=head2 parse_hmmscan_result
 
 Function: parse hmmscan result 

 Input: hmmscan result file name

 Output1: detail hits info, format is below:
          GeneID	PfamID		GA	Evalue
	  AT2G34620.1	PF02536.7	242.2	4.5e-72

 Output2: alignment of hits into, format is below:
	  1. GeneID      -- AT1G01140.1
	  2. PfamID      -- PF00069.18
	  3. Query Start -- 19
          4. Query End   -- 274
	  5. Hit Start   -- 1
          6. Hit End     -- 260
          7. Query Seq   -- YEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVF....
	  8. Alignment   -- ye++++lG+Gsf+kV  ak+  tg++ A+Kil++e+  ....
	  9. Hit Seq     -- yelleklGsGsfGkVykakkkktgkkvAvKilkkeeek....
          10. GA Score   -- 241.4
          11. Evalue     -- 6.9e-72 
          12. Description-- Protein kinase domain
          13. Qeury Len  -- 448
=cut
sub parse_hmmscan_result
{
	my $hmm_result = shift;

	my ($result_out1, $result_out2, $one_result, $align_detail, $hits);

	my $rfh = IO::File->new($hmm_result) || die "Can not open hmmscan result file : $hmm_result $! \n";
	while(<$rfh>)
	{
		unless (/^#/)
		{
			if (/^\/\//)
			{
				#################################
				# parse domain info		#
				#################################

				############################################################
				# short info for every hsp				   #
				# $hits = query_id \t domain_id \t evalue \t score desc \n #
				############################################################
	
				#$hits = parse_Hits($one_result);
				my @hits_content = split(/>>/, $one_result);

				#########################################################
				# parse hit head content like below			#
				#########################################################
				#
				# Match a line like this
				# E-value  score  bias    E-value  score  bias    exp  N   Sequence Description
				# -------  -----  -----   -------  ------ -----   ---- --  -------- -----------
				#  4e-83   285.8  10.0    5.3e-83  285.5  7.0     1.1  1   Q14SN3.1 Q14SN3_9HEPC Polyprotein (Fragment).
				#######################################################################################################
				# previous hmmer3 parse function need this part for best one domain.
				# New version just need Query name, length and no hit information in this part.
				#######################################################################################################
				my @hit_head = split(/\n/, $hits_content[0]);

				my ($query_name, $query_length);

				my $jumper = 0;
				foreach my $hit_head_line (@hit_head)
				{
					if ($hit_head_line =~ m/^Query:\s+(\S+)\s+\[M\=(\d+)\]/) 
					{
						$query_name = $1; $query_length = $2;
					}
					elsif ($hit_head_line =~ m/^Query:\s+(\S+)\s+\[L\=(\d+)\]/) 
					{
						$query_name = $1; $query_length = $2;
					}
					elsif ($hit_head_line =~ m/No hits detected that satisfy reporting thresholds/)
					{
						$jumper = 1;
						#die "$hit_head_line\n ";
					}
					else { next; }
				}

				#########################################################
				# parse hsp part of hits				#
				#########################################################

				my $hsp_detail = ""; my $hsp_info = "";

				unless( $jumper == 1)
				{
					for(my $ih=1; $ih<@hits_content; $ih++)
					{
						my $one_hit = ">>".$hits_content[$ih];
						($hsp_info, $hsp_detail) = parse_align($one_hit, $query_name, $query_length);
						$result_out1.= $hsp_info;
						$result_out2.= $hsp_detail;
					}
				}

				#########################################################
				# init 							#
				#########################################################
			    	$one_result = "";
			}
			else 
			{
				#store all one protein hmmscan info to this char;
				$one_result.=$_;
			}
		}
	}
	$rfh->close;

	return ($result_out1, $result_out2);
}

sub parse_align
{
	my ($hsp_info, $query_name, $query_length) = @_;

	my $output1 = ""; my $output2 = "";
	
	#########################################################
	# get hit id, hit desc and hsp info form one hit	#
	#########################################################
	my @hsp_line = split(/\n/, $hsp_info);

	my ($hit_id, $hit_desc);
	my %info1 = ();

	for(my $i=0; $i<@hsp_line; $i++)
	{
		if ($hsp_line[$i] =~ m/^>>\s+/)
		{
			my @aa = split(/\s+/, $hsp_line[$i], 3);
			$hit_id = $aa[1];
			$hit_desc = $aa[2];
		}
		elsif ( $hsp_line[$i] =~ m/^\s+(\d+)\s+\W\s+/)
		{
			$info1{$1} = $hsp_line[$i];
		}
		else
		{
			next;	
		}
	}

	#########################################################
	# get query string, hit string, match string of HSP	#
	#########################################################
	my ($query_string, $hit_string, $match_string, $hsp_length, $align_pos, $match_start);

	my @domain = split(/== domain/, $hsp_info);

	for(my $j=1; $j<@domain; $j++)
	{
		my @domain_line = split(/\n/, $domain[$j]);

		my @info = split(/\s+/, $info1{$j});

		for(my $k=1; $k<@domain_line; $k++)
		{
			#################################################
			# get hit string of HSP				#
			#################################################
			if ($domain_line[$k] =~ m/^\s+\Q$hit_id\E\s+\d+\s+(\S+)\s+\d+/ )
			{
				my @fff = split(/\s+/, $domain_line[$k]);

				$hit_string = $fff[3];
	
				$hsp_length = length($hit_string);

				$align_pos = index($domain_line[$k], $fff[3]);

				$match_start = 1;

				if ($align_pos < 0)
				{
					die "Error! Align start position is negative: $align_pos\n$domain_line[2]\n$hit_string\n";
				}
			}

			#################################################
			# get match string of HSP			#
			#################################################
			elsif (defined $match_start && $match_start == 1 )
			{
				$match_string = substr($domain_line[$k], $align_pos, $hsp_length);
				$match_start = 0;
			}

			#################################################
			# get query string of HSP			#
			#################################################
			elsif ($domain_line[$k] =~ m/^\s+\Q$query_name\E/ ) 
			{
				my @qqq = split(/\s+/, $domain_line[$k]);
				$query_string = $qqq[3];
			}
			else
			{
				next;
			}

		}
		
		#where these code come from? 
		#$score   = $dMatch[3];
		#$evalue  = $dMatch[6];
		#$hmmfrom = $dMatch[7];
		#$hmmto	  = $dMatch[8];
		#$seqfrom = $dMatch[10];
		#$seqto   = $dMatch[11];

		$output1.="$query_name\t$hit_id\t$info[3]\t$info[6]\n";
		$output2.="$query_name\t$hit_id\t$info[10]\t$info[11]\t$info[7]\t$info[8]\t$query_string\t$match_string\t$hit_string\t$info[3]\t$info[6]\t$hit_desc\t$query_length\n";
	}

	return ($output1,$output2);
}

=head2 parse_format_result

 Function: filter hmmscan parsed result, if socre >= GA score or evalue <= 1e-3, it will be seleted.

 Input: 1. formated result of hmmscan: query name; hit name;score; evalue; 
	2. gathering score hash, PF id ; GA score; 
                                 key      vaule

 Return: hash1 -- key: seq_id,    value: domain_id1 \t domain_id2 \t ... \t domain_idn
 Return: hash2 -- key: domain_id, value: PfamID
 Retrun: hash3 -- key: domain_id, value: SeqID \t PfamID \t 
=cut
sub parse_format_result
{
	my ($in_file, $ga_score_hash) = @_;

	my %out_hash; my %hsp_hit_id; my %hsp_detail;

	my @in_file = split(/\n/, $in_file);

	my $len = length(scalar(@in_file)); 

	my $uid = 0;
	
	for(my $in=0; $in<@in_file; $in++)
	{
		#################################################
		# filter HSP detail result with GA or e-value	#
		#################################################
		chomp($in_file[$in]);
		my @fmm = split(/\t/, $in_file[$in]);

		$fmm[1] =~ s/\..*//; #this is Pfam or self-build domain ID

		my $valued = 0; my $gaScore;

		if (defined $$ga_score_hash{$fmm[1]})
		{
			$gaScore = $$ga_score_hash{$fmm[1]};

			if ( $fmm[9] >= $gaScore ) { $valued = 1; }
		}
		else
		{ 
			if ($fmm[10] <= 1e-3) { $valued = 1; }
		}

		#################################################
		# if valued means we select this domain info	#
		# then creat three hashes base on this		#
		#################################################
		if ($valued == 1)
		{
			$uid++; my $zero = "";
			my $rlen = $len-length($uid);
			for(my $l=0; $l<$rlen; $l++)
			{
				$zero.="0";
			}

			$hsp_detail{$zero.$uid} = $in_file[$in];
			$hsp_hit_id{$zero.$uid} = $fmm[1];
			
			if (defined $out_hash{$fmm[0]})
			{
				$out_hash{$fmm[0]}.= "\t".$zero.$uid;
			}
			else
			{
				$out_hash{$fmm[0]}.= $zero.$uid;
			}
		}
	}	
	return (\%out_hash, \%hsp_hit_id, \%hsp_detail);
}

=head2


=head2 parse_rule

 Function: parse rule file, return three hash.

 Input: rule_list (file)

 Return:
	return hash1 -- key: family_name, value: domain_id1 # domain_id2 # ... # domain_idn
	return hash2 -- key: family_name, value: domain_id1 # domain_id2 # ... # domain_idn
	return hash3 -- key: family_name, value: mode
=cut

sub parse_rule
{
	my $rule_file = shift;
	
	# put rules to hash: key: rid (order), name, required, required num, forbidden 
	my %rules;
	my $fh = IO::File->new($rule_file) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/:/, $_, 2)
	}
	$fh->close;
}


sub parse_rule
{
	my $rule_file = shift;
	my %required; my %forbidden; my %mode;
	my ($family_name, $required_d, $forbidden_d, $required_m);
 
	my $fh = IO::File->new($rule_file) || die "Can not open rule file $rule_file $! \n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/:/, $_, 2);

		if ($a[0] eq 'Family Name') 	{ $family_name = $a[1]; }
		if ($a[0] eq 'Required Domain')	{ $required_d = $a[1];  }
		if ($a[0] eq 'Forbidden Domain'){ $forbidden_d = $a[1]; }
		if ($a[0] eq 'Required Mode')	{ $required_m = $a[1];  }
		if ($_ eq "//")
		{
			$required{$family_name} = $required_d;
			$forbidden{$family_name}= $forbidden_d;
			$mode{$family_name} = $required_m;
			$family_name = ""; $required_d = ""; $forbidden_d = ""; $required_m = "";
		}
	}
	$fh->close;

	return (\%required, \%forbidden, \%mode);



}

=head2 get_domain_id
 
 Function: Get domain ID from rule_list file; return domain id hash

 Input: rule_list_file
	mode

 Return: mode 1 : return all unique domain ids
	 mode 2 : return all required domain ids 

=cut
sub get_domain_id
{
	my ($rule_list_file, $mode) = @_;

	my %out_hash; my $family_name;

	my $fh = IO::File->new($rule_list_file) || die "can not open rule list file $!\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/:/, $_, 2);
		if ($a[0] eq 'Family Name') { $family_name = $a[1]; }

		if ($a[0] eq 'Required Domain' || $a[0] eq 'Forbidden Domain')
		{
			my @b;
			if ($a[1]) { @b = split(/#/, $a[1]); }

			if ($mode == 1)
			{
				foreach my $out_domain_id (@b)
				{
				    if ($out_domain_id =~ m/\w+/)
				    {
					$out_hash{$out_domain_id} = 1;
				    }
				}
			}
			elsif( $mode == 2 && $a[0] eq 'Required Domain')
			{
				foreach my $out_domain_id (@b)
				{
				    if ($out_domain_id =~ m/\w+/)
				    {
					if (defined $out_hash{$out_domain_id}) { $out_hash{$out_domain_id}.= "\t".$family_name; }
					else { $out_hash{$out_domain_id} = $family_name; }
				    }
				}
			}
			else
			{

			}
		}	
	}
	$fh->close;

	return %out_hash;
}

=head2 check_family

 Function: checking TFs families for every seqs base seq scan result;

 Input: all_hmmscan_domains (array);
	required_domains (array);
	forbidden_domains (array);
	mode (char);

 Return: is_family (char) : 1, is family. 2, not family

=cut
sub check_family
{
	my $domain = shift;	my $required = shift;	    my $forbidden = shift;
	my @domain = @$domain;  my @required = @$required;  my @forbidden = @$forbidden;

	my $r_mode = shift;

	my $is_family = 0; my $off = 0;

	my %required; my $required_i = 0;
	foreach my $re_id (@required)
	{
		$required_i++;
		$required{$required_i} = $re_id;
	}

	my %domain; my $domain_i = 0;
	foreach my $domain_id (@domain)
	{
		$domain_i++;
		$domain{$domain_i} = $domain_id;
	}

	#########################################################
	# for all mode, delete every required                   #
	#########################################################

	if ($r_mode eq "all")
	{
		foreach my $domain_order (sort keys %domain)
		{
			my $domain_id = $domain{$domain_order};
		
			foreach my $required_order (sort keys %required)
			{
				my $required_id = $required{$required_order};
				
				if ($domain_id eq $required_id)
				{
					delete $required{$required_order};
					delete $domain{$domain_order};
					last;
				}
			}
		}

		if ( scalar(keys(%required)) > 0 )
		{
			$off = 1;
		}
	}

	#########################################################
        # If find one forbidden domain, turn off it.            #
        # After turn off, it is not identify as TF              #
        #########################################################
	foreach (my $iii =0; $iii< @domain; $iii++)
	{
		my $d_id = $domain[$iii];
		
		foreach my $f_id (@forbidden)
		{
			if ($d_id eq $f_id) { $off = 1; }
		}
	}

	#########################################################
        # using off value to decide family value                #
        #########################################################
	
	if ( $off == 0) { $is_family = 1; }

	return $is_family;
}

=head2 identify_domain
 identify Transcription Factors conde + finde Protein Kinases domains 
=cut
sub identify_domain
{
	my ($all_hits, $all_detail, $rule, $transposase) = @_;

	my $alignment = "";
	my $family = "";

	# 1. parse hmmscan result
	# my ($all_hits, $all_detail) = parse_hmmscan_result($hmm_result);

	# 2. Using GA score and e-value filter the all hits
	# $pid_did    -- key: protein sequence id;  value: domain order id
	# $hsp_hit    -- key: domain order id;      value: PfamID
	# $hsp_detail -- key: domain order id;      value: aligment_detail
	my ($pid_did, $hsp_hit, $hsp_detail) = parse_format_result($all_detail, $ga_hash);

	# 3. Parse rule list file to produce rules
	my ($required_pack, $forbidden_pack, $rule_mode) = parse_rule($rule);
	my %required = %$required_pack; my %forbidden = %$forbidden_pack; my %rule_mode = %$rule_mode;

	my %required_domain = get_domain_id($rule, 2);
	print "\nThere are ".scalar(keys(%required_domain))." required domains for TF classification\n" if $debug == 1;

	# 4. Classify the proteins and output the results to files;
    	foreach my $protein_id (sort keys %$pid_did)
    	{

		foreach my $rid (sort keys %rules)
		{

		}
		




		my $is_family = 0;

		my @did = split(/\t/, $$pid_did{$protein_id});

		my $convert_domain = "";	# convert uid to domain id for one hit;
		my $hit_alignment = "";		# get alignment for one hit;

		for(my $ci=0; $ci<@did; $ci++)
		{
			my $c_domain_id = $$hsp_hit{$did[$ci]};
			$convert_domain = $convert_domain."\t".$c_domain_id;

			my $hsp_alignment = $$hsp_detail{$did[$ci]};
			$hit_alignment.=$hsp_alignment."\n";
		}
		$convert_domain =~ s/^\t//;

		my @convert_domain = split(/\t/, $convert_domain);

		# checking transposase
		my $has_transposase = 0;
		foreach my $domain_id ( @convert_domain )
		{
			if (defined $$transposase{$domain_id}) { $has_transposase = 1; }
		}

		if ($has_transposase == 1) 
		{
			print "Protein $protein_id was removed for including transposase: $convert_domain\n";
			next; 
		}




		# checking family
		for(my $di=0; $di<@did; $di++)
		{
			my $domain_id = $$hsp_hit{$did[$di]};  # convert uid to domain id

			if (defined $required_domain{$domain_id} && $is_family == 0 )
			{
				my @familys = split(/\t/, $required_domain{$domain_id});

				foreach my $fn (@familys)
				{
					my @this_required = split(/#/, $required{$fn});

					my @this_forbidden;

					if ( $forbidden{$fn} )
					{
						@this_forbidden = split(/#/, $forbidden{$fn});
					}

					my $mode_r = $rule_mode{$fn};	
	
					$is_family = check_family(\@convert_domain, \@this_required, \@this_forbidden, $mode_r);

					if ($is_family == 1) 
					{
						if ($fn eq "ARR-B_A") { $fn = "ARR-B"; }
				 
						$family.=$protein_id."\t".$fn."\n";
						$alignment.= $hit_alignment;
						last; 
					}
				}
			}

			if ($is_family == 1) { last; }
		}	

	    	unless ($is_family == 1) 
	    	{
			#################################################
			# to check myb-related				#
			#################################################
			my $is_a = 0; my $not_is_a = 0;
			for (my $dii=0; $dii<@did; $dii++)
			{
				if (    $$hsp_hit{$did[$dii]} eq "PF00249" ) { $is_a = 1; }
				if (    $$hsp_hit{$did[$dii]} eq "PF01388" ||
					$$hsp_hit{$did[$dii]} eq "PF00072" ||
					$$hsp_hit{$did[$dii]} eq "PF00176" ||
					$$hsp_hit{$did[$dii]} eq "G2-like" ||
					$$hsp_hit{$did[$dii]} eq "Trihelix" )
				{ $not_is_a = 1; }
			}

			if ($is_a == 1 && $not_is_a == 0)
			{
				$family.=$protein_id."\t"."MYB\n";
				$alignment.= $hit_alignment;
			}

			#################################################
			# to check orphans				#
			#################################################
			my $is_orphans = 0;
			for(my $dii=0; $dii<@did; $dii++)
			{
				if (
					$$hsp_hit{$did[$dii]} eq "PF06203" || 
					$$hsp_hit{$did[$dii]} eq "PF00643" || 
					$$hsp_hit{$did[$dii]} eq "PF00072" || 
					$$hsp_hit{$did[$dii]} eq "PF00412" ||
					$$hsp_hit{$did[$dii]} eq "PF02671" || 
					$$hsp_hit{$did[$dii]} eq "PF03925" || 
					$$hsp_hit{$did[$dii]} eq "PF09133" || 
	                                $$hsp_hit{$did[$dii]} eq "PF09425" 
				   )
				{
					$is_orphans = 1;
				}
			}

			if ($is_orphans == 1)
			{
				$family.=$protein_id."\t"."Orphans\n";
				$alignment.= $hit_alignment;
			}
	 	}
    	}
	return ($alignment, $family);
}

=head2 compare_a_b


=cut
sub compare_a_b
{
        my ($ha, $hb) = @_;

        my %ha = %$ha; my %hb = %$hb;

        my %hash1; my %hash2; my %match;

        foreach my $ida (sort keys %ha)
        {
                my $aa = $ida."\t".$ha{$ida};
                $hash1{$aa} = 1;
        }

        foreach my $idb (sort keys %hb)
        {
                my $bb = $idb."\t".$hb{$idb};
                $hash2{$bb} = 1;
        }

        foreach my $key1 (sort keys %hash1)
        {
                if (defined $hash2{$key1})
                {
                        $match{$key1} = 1;
                        delete $hash1{$key1};
                        delete $hash2{$key1};
                }
        }

        return (\%match, \%hash1, \%hash2);
}

=head2 seq_to_hash

 Function: split nucleotide/protein sequence, convert nucleotide to protien
 Input:    fasta sequence
 Output:   1) protein seq hash
           2) RNA/DNA translated seq hash
=cut
sub seq_to_hash
{
	my $seq = shift;
	my %seq_hash = (); my %converted_seq_hash = ();

	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$seq);
	while(my $inseq = $in->next_seq)
	{
		if ( $inseq->alphabet eq "rna" || $inseq->alphabet eq "dna" )
		{
			my @prots = Bio::SeqUtils->translate_6frames($inseq);
                        for (my $i = 0; $i < @prots; $i++)
                        {
                                my $tseq = $prots[$i]->seq;
				my $tid  = $inseq->id."_$i";
				$converted_seq_hash{$tid} = $tseq;
                        }

		}
		else
		{
			$seq_hash{$inseq->id} = $inseq->seq;
		}
	}
	return (\%seq_hash,\%converted_seq_hash);
}

=head2 pk_to_hash

=cut
sub pk_to_hash
{
	my $file = shift;
	my %hash;
	my $pfh = IO::File->new($file) || die "Can not open protein kinase description file: $file $!\n";
	while(<$pfh>)
	{
		chomp;
		my @pm = split(/\t/, $_, 2);
		$hash{$pm[0]} = $pm[1];
	}
	$pfh->close;
	return (\%hash);
}


=head2 cutoff

parse hmm3 parsed result then select pkinase seq id from parsed result

input: hmm3 parsed result
output: selected pkinase id

=cut
sub cutoff
{
	my ($parsed_result, $transposase) = @_;
	my @parsed_result = split(/\n/, $parsed_result);

	my %pk_id; my %transposase_id;

	# check the input sequence
	for(my $i=0; $i<@parsed_result; $i++)
	{
		my @a = split(/\t/, $parsed_result[$i]);

		# check protein kinase ID and put it to hash
		if ($a[1] eq "PF00069.20" && $a[2] >= 20.4 ) { $pk_id{$a[0]} = 1; }
		if ($a[1] eq "PF07714.12" && $a[2] >= 20.3 ) { $pk_id{$a[0]} = 1; }

		# check transposase ID and put it to hash
		my $pfam_id = $a[1]; $pfam_id =~ s/\..*//ig;
		if (defined $$transposase{$pfam_id} && $a[2] >= $$transposase{$pfam_id}) 
		{
			$transposase_id{$a[0]} = 1;
		}
	}

	# remove the transposase from the protein kinase
	foreach my $protein (sort keys %pk_id)
	{
		if (defined $transposase_id{$protein})
		{
			delete $pk_id{$protein};
			print "Protein $protein was removed for including transposase, PKs identify process\n";
		}
	}

	return %pk_id;
}

=head2 get_classification

 output the hit with highest score

=cut
sub get_classification
{
	my $simple_hmm2_result = shift;
	my @array = split(/\n/, $simple_hmm2_result);

	my %hit; my %score;
	for(my $i=0; $i<@array; $i++)
	{
		my @a = split(/\t/, $array[$i]);

		if (defined $hit{$a[0]})
		{
			if ($a[2] > $score{$a[0]})
			{
 				$hit{$a[0]} = $a[1];
				$score{$a[0]} = $a[2];
			}
		}
		else
		{
			$hit{$a[0]} = $a[1];
			$score{$a[0]} = $a[2];
 		}
	}
	return %hit;
}

=head2 aln_to_hash
=cut
sub aln_to_hash
{
	my $aln_detail = shift;

	my %aln_hash;

	my @line = split(/\n/, $aln_detail);

	for(my $i=0; $i<@line; $i++)
	{
		my @a = split(/\t/, $line[$i]);

		# filte the aligment by GA score
		my $pfam_id = $a[1];
		$pfam_id =~ s/\..*//;

		if ( $a[9] >= $$ga_hash{$pfam_id} )
		{
			#print $a[9]."\t$pfam_id\t".$$ga_hash{$pfam_id}."\n";
			if (defined $aln_hash{$a[0]})
			{
				$aln_hash{$a[0]}.=$line[$i]."\n";
			}
			else
			{
				$aln_hash{$a[0]} = $line[$i]."\n";
			}
		}
	}
	return %aln_hash;
}

=head2 get_wnk1
 easy to use, after hmm2 classification scan. put ppc id and hmm profile to this function, then the function 
 will compare the seq belong to this id and hmm file, next re-assign new classification to it.
=cut
sub get_wnk1
{
	my ($input, $hmm, $score, $cat_id, $new_cat_id) = @_;

	my %pk_cat = %$input;

	#########################################################
	# find seq base on ppc cat id				#
	#########################################################
	my $ppcseq = "$temp_dir/temp_ppc_seq";

	my $ppfh = IO::File->new(">".$ppcseq) || die "Can not open ppc sequence file: $ppcseq\n";

	my $s_num = 0;

	foreach my $seq_id (sort keys %pk_cat)
	{
		if ($pk_cat{$seq_id} eq $cat_id)
		{
			if (defined $seq_hash{$seq_id})
			{
				print $ppfh ">".$seq_id."\n".$seq_hash{$seq_id}."\n";
				$s_num++;
			}
			else
			{
				die "Error: $seq_id do not have sequence in sequence file\n";
			}
		}
	}
	$ppfh->close;

	#########################################################
	# do hmmscan , compare ppcseq and hmm			#
	#########################################################
    if ($s_num >=1)
    {
	my $ppc_hmm_result = $temp_dir."/temp_ppc_hmm3_result"; 

	my $hmm_cmd = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $ppc_hmm_result $hmm $ppcseq";

	print $hmmscan_command if $debug == 1;

	system($hmm_cmd) && die "Error at hmmscan command: $hmm_cmd\n";

	#########################################################
	# parse hmmscan result	get simple result and new cat	#
	# base on score						#
	#########################################################

	my ($ppc_hits, $ppc_detail) = parse_hmmscan_result($ppc_hmm_result);

	my @hit = split(/\n/, $ppc_hits);

	for (my $i=0; $i<@hit; $i++)
	{
		my @a = split(/\t/, $hit[$i]);

		if ($a[2] >= $score)
		{
			$pk_cat{$a[0]} = $new_cat_id;
		}
	}	
    }
	#########################################################
	# output the result					#
	#########################################################

	return %pk_cat;
}

=head2

=cut
sub tf_family_cat_to_hash
{
	my $input_file = shift;

	my %hash_cat = ();

	my $zfh = IO::File->new($input_file) || die "Can not open category file: $input_file\n";

	while(<$zfh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$hash_cat{$a[0]} = $a[1];
	}
	$zfh->close;

	return %hash_cat;
}
