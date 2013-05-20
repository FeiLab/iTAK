#!/usr/bin/perl

=head
v 1.2 Aug-26-2011
      report unusual sequences

v 1.1 Jun-03-2011
      remove unsignificiant domain using GA score

v 1.0 Dec-14-2010
      first stable version 

=cut

#use strict;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use IO::File;
use Getopt::Std;
use Bio::SearchIO;

my $usage = q/
UPDATE: Aug-26-2011
VERSION: v1.2
USAGE: 
	Perl iTAK.pl [parameters]

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
my %options;

getopts('i:m:a:s:o:d', \%options) || die "$usage\nError getting options!";

my $input_seq = 	$options{i} || die $usage."\nYou must provide input file (-i)!\n\n";

my $seq_format = 	$options{s} || "p";

my $mode = 		$options{m} || "b";

my $cpus = 		$options{a} || 1;

my $output_fd = 	$options{o} || $input_seq."_output";

my $debug = 		$options{d} || 0;

#################################################################
#--------------Check input files and parameters-----------------#
#################################################################

#################################################################
# check mode parameters						#
# 1. mode t, just identify transcription factors		#
# 2. mode p, just identify protein kinases			#
# 3. mode b, identify both					#
#################################################################
unless ($mode =~ m/^t|p|b$/i) { die $usage."\nPlease input right parameter for '-m'\n"; }

#################################################################
# check input_seq files						#
# 1. seq_format eq 1 means input seqs are DNA			#
# 2. seq_format eq 0 (default) means input seqs are Proteins	#
# if you input DNA, it can be converted to protein, if your	#
# input seq and seq_format are not match, will die		#
#################################################################
my %seq_hash;

my ($input_seq_hash, $converted_seq_hash) = seq_to_hash($input_seq);

if ($seq_format eq "N" || $seq_format eq "n")
{
	if ( scalar(keys(%$input_seq_hash)) > 0 ) { 
		print "\nYou input contains protein or unusual sequences!\n\n";

		my %unusual_list;

		foreach my $aaa (sort keys %$input_seq_hash) {
			print ">".$aaa."\n".$$input_seq_hash{$aaa}."\n";
			my $bbb = $aaa; $bbb =~ s/_\d$//ig;
			$unusual_list{$bbb} = 1;
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
	if ( scalar(keys(%$converted_seq_hash)) > 0) {

		print "\nYour input contains DNA or unusual sequences!\n\n";
		
		my %unusual_list;

		foreach my $aaa (sort keys %$converted_seq_hash) { 
			print ">".$aaa."\n".$$converted_seq_hash{$aaa}."\n";
			my $bbb = $aaa; $bbb =~ s/_\d$//ig;
			$unusual_list{$bbb} = 1;
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
	die "Error at paramter -s $seq_format\n$usage";
}

#################################################################
#------------------Default and output Floders-------------------#
#################################################################

#################################################################
# create temp folder and then store seq to temp			#
#################################################################

my $current_dir = getcwd;
my $temp_dir = $current_dir."/".$input_seq."_temp";

if (-e $output_dir)
{
	system "rm -r -f $temp_dir";
	mkdir($temp_dir);
}
else
{
	mkdir($temp_dir);
}

my $protein_seq = $temp_dir."/protein_seq";
my $pfh = IO::File->new(">".$protein_seq) || die "can not open protein sequence file: $protein_seq $!\n";
foreach my $pid (sort keys %seq_hash)
{
	print $pfh ">".$pid."\n".$seq_hash{$pid}."\n";
}
$pfh->close;

#################################################################
# check and create output folder				#
#################################################################

my $output_dir = $current_dir."/".$output_fd;

if (-e $output_dir)
{
	system "rm -rf $output_dir";
	mkdir($output_dir);
}
else
{
	mkdir($output_dir);
}

#################################################################
#------------------Varialbes that Default Set-------------------#
#################################################################

my $program_dir;
if ($0 =~ m{^/}) { $program_dir = dirname($0);}
else { $program_dir = dirname("$current_dir/$0"); }
my $bin_dir =  $program_dir."/bin";
my $dbs_dir =  $program_dir."/database";

unless (-e $bin_dir) { die "bin directory does not exist.\n"; }
unless (-e $dbs_dir) { die "database directory does not exist.\n"; }

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
unless (-s $tf_rule) { die "Rules for Transcription Factors domains do not exist.\n"; }

#################################################################
# category of transcription factor and regulator		#
#################################################################
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

#################################################################
# protein kinases Cat ID and It's desc to hash.			#
#################################################################
my $pk_desc = $dbs_dir."/PK_class_desc";
unless (-s $pk_desc) { die "protein kinase descriptions do not exist.\n"; }
my $pkid_des = pk_to_hash($pk_desc);

#################################################################
# parameters end                                                #
#################################################################

#################################################################
#--------------------------Main Part----------------------------#
#################################################################

# Step1. using all seqs do HMM search against PfamA+Self-build HMM or plantPKhmm
my $tmp_hmm_result = $temp_dir."/temp_hmm3_result"; 

my $hmmscan_command = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_hmm_result $hmm3_db $protein_seq";

print $hmmscan_command if $debug == 1;

unless (-s $tmp_hmm_result)
{	
	system($hmmscan_command) && die "Error at hmmscan command: $hmmscan_command\n";
}

my ($all_hits, $all_detail) = parse_hmmscan_result($tmp_hmm_result);

# Step2 parse hmmscan result using TFs prediction subroutine
my $tmp_rule = $tf_rule;

my ($tmp_align, $tmp_family) = identify_domain($all_hits, $all_detail, $tmp_rule);

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
	my %pkinase_id = cutoff($all_hits);

	my %pkinase_aln = aln_to_hash($all_detail);

	my $protein_kinase_seq = $output_dir."/".$input_seq."_pkseq";

	my $pk_seq_num = scalar(keys(%pkinase_id));

    if ($pk_seq_num > 0)
    {
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
	my $tmp_pk_hmmpfam_3 = "$temp_dir/temp_hmm3_result2";

	my $hmmscan_command3 = $bin_dir."/hmmscan --acc --notextw --cpu $cpus -o $tmp_pk_hmmpfam_3 $plantp_hmm_3 $protein_kinase_seq";

	unless (-s $tmp_pk_hmmpfam_3)
	{
		system($hmmscan_command3) && die "Error at hmmscan command: $hmmscan_command3\n";
	}
	
	print $hmmscan_command3 if $debug == 1;
	
	# Step 3.3.2 parse hmmscan result
	my ($hmm3_simple_result, $hmm3_simple_align) = parse_hmmscan_result($tmp_pk_hmmpfam_3);

	# Step 3.3.3 get classification info base simple hmminfo, means get best hits of simple result, then add annotation and output alignment file
	my %pkinases_cat = get_classification($hmm3_simple_result);

	my $protein_kinase_aln = $output_dir."/".$input_seq."_pkaln";
	my $protein_kinase_cat = $output_dir."/".$input_seq."_pkcat";
	
	my $ca_fh = IO::File->new(">".$protein_kinase_cat) || die "Can not open protein kinase sequence file: $protein_kinase_cat";
	my $al_fh = IO::File->new(">".$protein_kinase_aln) || die "Can not open protein kinase sequence file: $protein_kinase_aln";

	foreach my $pid (sort keys %pkinases_cat)
	{
		if (defined $pkinase_aln{$pid})
		{
			print $al_fh $pkinase_aln{$pid};
		}
		else
		{
			die "Error! Do not have alignments in hmm3 parsed result\n";
		}

		delete $pkinase_id{$pid};	
	}

	foreach my $ppid (sort keys %pkinase_id)
	{
		print $ca_fh $ppid."\t1.Other\n";

		if (defined $pkinase_aln{$ppid})
		{
			print $al_fh $pkinase_aln{$pid};
		}
		else
		{
			die "Error! Do not have alignments in hmm3 parsed result\n";
		}
	}
	$al_fh->close;

	# Step 3.3.4 find WNK1 domain in 4.1.5 cat;
	# $input, $hmm, $score, $cat_id, $new_cat_id
	%pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/wnk1_hmm_domain/WNK1_hmm" , "30" ,"PPC:4.1.5", "PPC:4.1.5.1");
	$$pkid_des{"PPC:4.1.5.1"} = "WNK like kinase - with no lysine kinase";

	%pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/mak_hmm_domain/MAK_hmm" , "460.15" ,"PPC:4.5.1", "PPC:4.5.1.1");
	$$pkid_des{"PPC:4.5.1.1"} = "Male grem cell-associated kinase (mak)";
	
	foreach my $pid (sort keys %pkinases_cat)
	{
		print $ca_fh  $pid."\t".$pkinases_cat{$pid}."\t".$$pkid_des{$pkinases_cat{$pid}}."\n";
	}
	$ca_fh->close;

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
my $cmd = "rm -rf $temp_dir";

system($cmd) && die "Error at $cmd\n";

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

 Output1: best hits info
 Output2: alignment of best hits info
 
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
						$one_hit = ">>".$hits_content[$ih];
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
	my ($query_string, $hit_string, $match_string);

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
			elsif ( $match_start == 1 )
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
			$score   = $dMatch[3];
				$evalue  = $dMatch[6];

				$hmmfrom = $dMatch[7];
				$hmmto	 = $dMatch[8];
				$seqfrom = $dMatch[10];
				$seqto   = $dMatch[11];

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

 Return: hash1 : seq_id   domain_id1 \t domain_id2 \t ... \t domain_idn
 Return: hash2 : 

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
			$uid++; $zero = "";
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

=head2 parse_rule

 Function: parse rule file, return two hash : family_name  required_domains && family_name  forbidden_domains

 Input: rule_list (file)

 Return:
	returen array: (\%hash1, \%hash2) ; hash1 is required, hash2 is fobidden
	                key             value
	return hash1 : family_name     domain_id1 \t domain_id2 \t ... \t domain_idn
	                key             value
	return hash2 : family_name     domain_id1 \t domain_id2 \t ... \t domain_idn

=cut
sub parse_rule
{
	my $rule_file = shift;

	my %required; my %forbidden; my %mode;

	my $family_name; my $required_d; my $forbidden_d; my $required_m;
 
	my $rufh = IO::File->new($rule_file) || die "Can not open rule file $rule_file $! \n";
	while(my $rufh_line = <$rufh>)
	{
		chomp($rufh_line);

		my @rmm = split(/:/, $rufh_line);

		if ($rufh_line =~ m/^Family Name/)
		{
			$family_name = $rmm[1];
		}

		if ($rufh_line =~ m/^Required Domain/)
		{
			$required_d = $rmm[1];
		}

		if ($rufh_line =~ m/^Forbidden Domain/)
		{
			$forbidden_d = $rmm[1];
		}

		if ($rufh_line =~ m/^Required Mode/)
		{
			$required_m = $rmm[1];
		}

		if ($rufh_line eq "//")
		{
			$required{$family_name} = $required_d;
			$forbidden{$family_name}= $forbidden_d;
			$mode{$family_name} = $required_m;
			$family_name = ""; $required_d = ""; $forbidden_d = ""; $required_m = "";
		}
	}
	$rufh->close;

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

	my $rlfh = IO::File->new($rule_list_file) || die "can not open rule list file $!\n";

	my %out_hash;

	my $family_name;

	while(my $rlfh_line = <$rlfh>)
	{
		chomp($rlfh_line);

		if ($rlfh_line =~ m/^Family Name:/) 
		{
			$family_name = $rlfh_line;
			$family_name =~ s/^Family Name://;
		}

		if ($rlfh_line =~ m/^Required Domain:/ || $rlfh_line =~ m/^Forbidden Domain:/)
		{
			my @mma = split(/:/, $rlfh_line);
			my @mmb;
			if ($mma[1])
			{
				@mmb = split(/#/, $mma[1]);
			}

			if ($mode == 1 )
			{
				foreach my $out_domain_id (@mmb)
				{
				    if ($out_domain_id =~ m/\w+/)
				    {
					unless (defined $out_hash{$out_domain_id}) { $out_hash{$out_domain_id} = 1; }
				    }
				}
			}
			elsif( $mode == 2 && $rlfh_line =~ m/^Required Domain:/)
			{
				foreach my $out_domain_id (@mmb)
				{
				    if($out_domain_id =~ m/\w+/)
				    {
					if (defined $out_hash{$out_domain_id}) { $out_hash{$out_domain_id} .= "\t".$family_name; }
					else { $out_hash{$out_domain_id} = $family_name; }
				    }
				}
			}
			else
			{

			}
		}	
	}
	$rlfh->close;

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
	my ($all_hits, $all_detail, $rule) = @_;

	my $alignment = "";
	my $family = "";

	# 1. parse hmmscan result
	#my ($all_hits, $all_detail) = parse_hmmscan_result($hmm_result);

	# 2. Using GA score and e-value filter the all hits
	my ($pid_did, $hsp_hit, $hsp_detail) = parse_format_result($all_detail, $ga_hash);

	# 3. Parse rule list file to produce rules
	my ($required_pack, $forbidden_pack, $rule_mode) = parse_rule($rule);

	my %required = %$required_pack; my %forbidden = %$forbidden_pack; my %rule_mode = %$rule_mode;

	my %required_domain = get_domain_id($rule, 2);
	print "\nThere are ".scalar(keys(%required_domain))." required domains for TF classification\n" if $debug == 1;

	# 4. Classify the proteins and output the results to files;

    foreach my $protein_id (sort keys %$pid_did)
    {
		#########################################################
		# for transcription factors				#
		#########################################################
		my $is_family = 0;

		my @did = split(/\t/, $$pid_did{$protein_id});

		my $convert_domain = ""; # convert uid do domain id for one hit;
		my $hit_alignment = ""; # get alignment for one hit;
		for(my $ci=0; $ci<@did; $ci++)
		{
			my $c_domain_id = $$hsp_hit{$did[$ci]};
			$convert_domain = $convert_domain."\t".$c_domain_id;

			my $hsp_alignment = $$hsp_detail{$did[$ci]};
			$hit_alignment .= $hsp_alignment."\n"
		}
		$convert_domain =~ s/^\t//;

		my @convert_domain = split(/\t/, $convert_domain);

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

				if ($forbidden{$fn} )
				{
					@this_forbidden = split(/#/, $forbidden{$fn});
				}

				my $mode_r = $rule_mode{$fn};	
	
				$is_family = check_family(\@convert_domain, \@this_required, \@this_forbidden, $mode_r);
				#$is_family =1;
				#print $is_family."\n";
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
    #######
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

 Function:

 Input:

 Output:

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

=head2 parse_hmmscan_result
=cut
sub parse_parsed_hmmscan_result
{
	my $hmm_scan_file = shift;

	my %hmm_result;

	my $all_line;
	my $hfh = IO::File->new($hmm_scan_file) || die "Can not open hmm scan result: $hmm_scan_file \n";
	while(<$hfh>)
	{
		if ($_ =~ m/Query:/ && $j == 12)
		{
			$all_line.= "//\n".$_;
		}
		else
		{
			$all_line.= $_;
		}
	}
	$hfh->close;
	
	my @result = split(/\/\/\n/, $all_line);

	my ($qid, $hr);

	foreach my $res (@result)
	{
		my @line = split(/\n/, $res);
		
		for(my $i=0; $i<@line; $i++)
		{
			if ($i == 0)
			{
				my @a = split(/\s+/, $line[$i]);
				$qid = $a[1];
			}
			elsif ($i == 5)
			{
				my $hr = $line[$i];
				$hr =~ s/\s+/\t/;
				if ($hr =~ m/PPC:/) { $hmm_result{$qid} = $hr; }
			}
			else
			{
				next;
			}
		}
	}
	return %hmm_result;
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
	my $parsed_result = shift;
	my @parsed_result = split(/\n/, $parsed_result);

	my %pk_id;

	for(my $i=0; $i<@parsed_result; $i++)
	{
		my @a = split(/\t/, $parsed_result[$i]);
		if ($a[1] eq "PF00069.18" && $a[2] >= 20.4 ) { $pk_id{$a[0]} = 1; }
		if ($a[1] eq "PF07714.10" && $a[2] >= 20.3 ) { $pk_id{$a[0]} = 1; }
	}
	return %pk_id;
}

=head2 get_classification
=cut
sub get_classification
{
	my $simple_hmm2_result = shift;
	my @array = split(/\n/, $simple_hmm2_result);

	my %hit;

	for(my $i=0; $i<@array; $i++)
	{
		my @a = split(/\t/, $array[$i]);

		if (defined $hit{$a[0]})
		{
                	my $high_score = $score{$a[0]};

			if ($a[2] > $high_score)
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
