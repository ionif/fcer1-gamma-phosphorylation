#!/usr/bin/perl -w

# Copyright 2014 Translational Genomics Research Institute

# This program is free software: you can redistribute it and/or modifyd
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# You can contact the developer by email at brandonthomas@nau.edu

use strict;

my $version = "1.0";

# Automatically reap forks that aren't done so with waitpid().
# Currently this in place to kill fork killer forks. Ideally
# they would be killed with a stopping waitpid(), however.
$SIG{CHLD} = 'IGNORE';

use Data::Dumper;
use Cwd qw(getcwd chdir abs_path);
use Config;
use File::Path;
use File::Copy qw(copy move);
use File::Path qw/mkpath rmtree/;
use File::Basename qw(fileparse);
use POSIX qw(floor ceil);
use POSIX qw( WNOHANG );    # In case we want to use a non-blocking waitpid()
use Safe;
use List::Util qw(sum max min);
use Scalar::Util qw(looks_like_number);

# Check to see if local machine has the libraries needed to produce graphical output. The conditional is further down.
eval {
	require GD::Graph::area;
	require GD::Graph::lines;
	require GD::Graph::mixed;
	GD::Graph::area->import();
	GD::Graph::lines->import();
	GD::Graph::mixed->import();
};

sub print_usage();
sub parse_file(@);
sub update_values($);
sub create_pbs_command($$$$);
sub verbose_mkdir($);
sub calculate_var_sets($$);
sub run_models_pbs($$$);
sub run_models_fork($$);
sub generate_model_set_files($$$$$);
sub load_exp($$);
sub load_data($);
sub store_factors($);
sub store_newcolumns($);
sub create_columns($$$$$$);
sub get_model_params($$);
sub pick_weighted($$$);
sub store_mutate($);
sub read_summary($);
sub run_generation($$$$$);
sub run_analyze($$$$);
sub run_genetic($$$$$);
sub create_new_config($$$$$$);
sub load_all_summaries($);
sub run_consolidate($$$);
sub consolidate($$@);
sub run_submit($$$);
sub run_resume($$);
sub zero_rand($);
sub rlnorm($$);
sub empty_sort(@);
sub write_all_summary($$$);
sub numalpha_sort;
sub combine_columns($$$@);
sub welcome_message($);
sub finished($$$);
sub graph($);
sub walltime_watch($$);
sub smooth_runs($$);
sub get_arch_info($);
sub bootstrap_run($$$);
sub show_output;
sub generate_bootstrap_file($);
sub run_monitor($);
sub check_queue($$);

# Get our root directory
my ( $filename_junk, $medir ) = fileparse( abs_path($0) );

# Path to our PBS submission script
my $pbs_generic_script = "${medir}BioNetFit.pbs";

my $verbosity;									# Terminal output verbosity
my $tolerance                   = 1e-6;          # Tolerance to round off in time comparison
my $max_retry_different_parents = 100;			# Maximum number of times we will try to select different parents during breeding
my $username                    = $ENV{USER};    # For default email address.
my $qsub                        = 0;			# Not using qsub by default

# Our main program commands.  These are submitted on the command line as arguments to the script
my $submit_command      = 'submit';
my $resume_command      = 'resume';
my $generation_command  = 'generation';
my $genetic_command     = 'genetic';
my $analyze_command     = 'analyze';
my $consolidate_command = 'results';
my $rlnorm_test         = 'rlnorm';
my $monitor_command     = 'monitor';
my $finish_command      = 'results';

# Make sure we always know where we started
my $top_directory = getcwd;

# Used later on to store file and path names
my @filenames = qw(bng_command model exp_file net_file);
my @pathnames = qw(output_dir);

# Parameters that are used to assign model variables.
# These parameters will not be included when outputting new config file
# after running genetic algorithm because they are only used to generate the 1st
# generation's config files
my @var_variables = qw(substitute_var fixed_var list_var linear_var log_var lognormrandom_var lograndom_var random_var static_list_var loguniform_var);

# Parameters that are assigned a list of values.
# These parameters will not be included when outputting new config file
# after running genetic algorithm because variables are being set specifically.
my @list_variables = qw(newcolumn factor mutate exp_file scan_parameter);

# var_variables and sensitivity_variables are list variables as well.
push @list_variables, @var_variables;

# Parameters that should not be passed to pbs
my @skip_variables = qw(pbs_mail_options email job_name walltime analyze mail_when_analyzed job_sleep generation max_generations default_mutate skip_genetic cluster_parallel);

# Required variables.
my @required_variables = qw(job_name model bng_command max_generations permutations exp_file output_dir );

# Variable defaults.
my %var_defaults = (
	parallel_count          => 2,
	job_sleep               => 1,
	force_different_parents => 1,
	mail_when_aborted       => "no",
	mail_when_begins        => "no",
	mail_when_terminates    => "no",
	show_welcome_message    => 1,
	ask_create              => 1,
	ask_overwrite           => 1,
	skip_genetic            => 0,
	stop_when_stalled       => 1,
	verbosity               => 1,
	make_plots              => 0,
	max_retries             => 3,
	retry_count             => 0,
	swap_rate               => 0.5,
	use_cluster             => 0,
	min_objfunc_value       => 0,
	smoothing               => 1,
	extra_weight            => 0,
	keep_parents            => 0,
	run_job_from_worknode   => 1,
	objfunc                 => 1,
	delete_old_files        => 1,
	bootstrap               => 0,
	bootstrap_num           => 1,
	save_cluster_output     => 0,
	bootstrap_retries       => 3,
	bootstrap_retry         => 0,
	first_gen_permutations  => 0,
	multisim                => 1,
	job_limit               => 5000,
	job_counter             => 0,
	max_walltime            => "08:00:00",
	walltime_limit			=> "240:00:00",
	cluster_parallel        => 4,
);

#############################################################
### Start of main program
#############################################################
{
	# Get our subcommand from the command line
	my $subcommand;
	my $new_max_generations;
	
	if ( 2 == @ARGV ) {
		$subcommand = shift @ARGV;
	}
	elsif ( 1 == @ARGV ) {
		$subcommand = $submit_command;
	}
	elsif ( 3 == @ARGV ) {
		$subcommand = shift @ARGV;
		$new_max_generations = shift @ARGV;
	}
	else {
		print_usage();
	}

	# Store path to our original .conf file
	my ($orig_config_file) = @ARGV;

	# Make sure we're using an absolute path for the config file
	$orig_config_file = abs_path($orig_config_file);

	unless ( -f $orig_config_file ) {
		die "BioNetFit couldn't find the .conf file you specified. Please check that the path and filename to the .conf file is correct\n";
	}

	# Read settings from configuration file(s).
	my @config_filenames;
	push @config_filenames, $orig_config_file;

	# Need to change to the BioNetFit directory in order to parse relative paths in .conf file.
	# Save CWD so we can get back after .conf parsing
	chdir $medir;

	# Create our $values_ref from the config file(s).  $values_ref is used throughout the script.
	my $values_ref = parse_file(@config_filenames);

	# Create defaults and modify values where appropriate.
	update_values($values_ref);

	# Chang back to cwd
	chdir $top_directory;

	# $qsub can be used anywhere in the script to know if we are on a cluster or not. It is used primarily to store the output filename for script output when
	# running BioNetFit on a cluster.
	if ( $values_ref->{use_cluster} ) {
		$qsub = "$medir/$values_ref->{job_name}";
	}

	# If we are bootstrapping, let's save our original config for later use
	$values_ref->{orig_config} = $orig_config_file
	  if ( $values_ref->{bootstrap} && !exists $values_ref->{generation} );

	# Set global verbosity
	$verbosity = $values_ref->{verbosity};

	# If we did not find graphing modules ($@ is set) but user wants to make graphs, tell them they can't.
	if ( $@ && $values_ref->{make_plots} ) {
		print "You've set make_plots to true but we can't find the required Perl modules.  Please see the manual for instructions on installing these modules.\n"
		  if ( $subcommand eq $submit_command );
		$values_ref->{make_plots} = 0;
	}

	# Initialize random number seed if seed is given.
	srand( $values_ref->{seed} )
	  if ( $values_ref->{seed} );

	# Be sure our model(s) end in .txt or .bngl
	my ($root) = $values_ref->{model} =~ /(.*)\.txt/;
	($root) = $values_ref->{model} =~ /(.*)\.bngl/ unless ($root);
	die "Error: Unrecognized model name $values_ref->{model}, needs .txt or .bngl extension.\n"
	  unless ($root);

	# If $values_ref->{generation} exists, it means we're PAST the first generation
	unless ( exists $values_ref->{generation} ) {
		if ( $values_ref->{max_generations} > 0 ) {

			# If multiple generations, start with 1.
			$values_ref->{generation} = 1;
		}
		else {
			# Not doing multiple generations, so start with 0.
			$values_ref->{generation} = 0;
		}
	}

	print "Getting architecture information...\n"
	  if ( $verbosity >= 3 );

	# Get architecture information
	get_arch_info($values_ref)
	  if ( !exists $values_ref->{command_success} && $subcommand ne $monitor_command );

	if ( $values_ref->{generation} == 1 || $subcommand eq $consolidate_command ) {
		my $model_full = $values_ref->{model};

		# Check too see if model file exists
		die "Error: Unable to find a model file. Please be sure all 'model' files in your .conf file exist.\n\n"
		  unless ( -f $model_full );
		my ( $model, $path, $suffix );

		( $model, $path, $suffix ) = fileparse( $model_full, qr/\.[^\.]*$/ );
		$model = "$path$model";

		my @exp_names;

		# Here we convert the .exp to unix-style line endings.
		# Read lines into array, first converting line endings
		foreach my $exp_file ( @{ $values_ref->{exp_file} } ) {
			my ( $basename, $path, $ext ) = fileparse( $exp_file, qr/\.[^\.]*$/ );
			push @exp_names, $basename;

			my $changes_made = 0;

			die "Error: Unable to find the .exp file: $exp_file\n"
			  unless ( -f $exp_file );
			open EXPFILE, '<', $exp_file
			  or show_output "die", "Error: Couldn't open .exp file for unix conversion.";
			my @exp_old;
			while ( my $line = <EXPFILE> ) {
				my $oldline = $line;
				$line =~ s/[\x0A\x0D]+/\n/g;
				$changes_made = 1
				  if ( $oldline ne $line );
				push @exp_old, $line;
			}
			close EXPFILE;

			# Only output new .exp if we changed anything
			if ($changes_made) {

				# Now output the converted file
				open EXPFILE, '>', $exp_file
				  or show_output "die", "Error: Couldn't open .exp file for unix conversion.";

				foreach (@exp_old) {
					print EXPFILE $_;
				}

				close EXPFILE;
			}
		}

		open MODELFILE, '<', $model_full
		  or show_output "die", "Error: Couldn't open model file for consistency checks.";
		my @orig_model = <MODELFILE>;
		close MODELFILE;

		my @new_model;

		# Replace line break delimiters, putting split lines onto the same line
		for ( my $i = 0 ; $i < @orig_model ; $i++ ) {
			if ( $orig_model[$i] =~ /\\/ ) {
				$orig_model[$i] =~ s/\\\W*//;
				$orig_model[$i] = $orig_model[$i] . $orig_model[ $i + 1 ];
				push @new_model, $orig_model[$i] unless $orig_model[$i] =~ /\\/;
				$i++;
				if ( $orig_model[ $i - 1 ] =~ /\\/ ) {
					$orig_model[ $i - 1 ] =~ s/\\\W*//;
					$orig_model[ $i - 1 ] = $orig_model[ $i - 1 ] . $orig_model[ $i + 1 ];
					push @new_model, $orig_model[ $i - 1 ] unless $orig_model[ $i - 1 ] =~ /\\/;
					$i++;
					if ( $orig_model[ $i - 2 ] =~ /\\/ ) {
						$orig_model[ $i - 2 ] =~ s/\\\W*//;
						$orig_model[ $i - 2 ] = $orig_model[ $i - 2 ] . $orig_model[ $i + 1 ];
						push @new_model, $orig_model[ $i - 2 ] unless $orig_model[ $i - 2 ] =~ /\\/;
						$i++;
						if ( $orig_model[ $i - 3 ] =~ /\\/ ) {
							$orig_model[ $i - 3 ] =~ s/\\\W*//;
							$orig_model[ $i - 3 ] = $orig_model[ $i - 3 ] . $orig_model[ $i + 1 ];
							push @new_model, $orig_model[ $i - 3 ] unless $orig_model[ $i - 3 ] =~ /\\/;
							$i++;
							if ( $orig_model[ $i - 4 ] =~ /\\/ ) {
								$orig_model[ $i - 4 ] =~ s/\\\W*//;
								$orig_model[ $i - 4 ] = $orig_model[ $i - 4 ] . $orig_model[ $i + 1 ];
								push @new_model, $orig_model[ $i - 4 ];
								$i++;
							}
						}
					}
				}
			}
			else {
				push @new_model, $orig_model[$i];
			}
		}

		my @suffix_lines;
		my @scan_parameter_lines;

		# Store any model lines containing suffixes or scan parameters. Also store whether
		# we are using an ode solver
		for (@new_model) {
			next if /^#/;
			#print $_;
		       if (/suffix=>(“|’|”|’|")\w+(“|’|”|’|")/i) {
				push @suffix_lines, $_;
				#print 2;
			}
			  if (/parameter_scan/i) {
				push @scan_parameter_lines, $_;
				#print $_;
			}
			$values_ref->{ode} = 1
			  if ( /^simulate.+method=>(“|’|”|’)ode(“|’|”|’)/i || /^simulate_ode/i );
			die "Error: It looks like you have a prefix specified in your model action command. This causes BioNetFit to get confused. Please use only suffixes, which must correspond to .exp filenames.\n\n"
			  if (/prefix=>(“|’|”|’|")\w+(“|’|”|’|")/i);
		}
		$values_ref->{ode} = 0 unless $values_ref->{ode};

		my $num_name_matches = 0;
		foreach my $suffix_line (@suffix_lines) {
			for (@exp_names) {
				#print $suffix_line;
				  if ( $suffix_line =~ /suffix=>(“|’|”|’|")$_(“|’|”|’|")/ ){
					$num_name_matches += 1;
				   	#print 1;
				  }
			}
		}

		my $need_to_die = 0;

		# See if our exp filenames have matching suffixes specified in model file
		if ( $num_name_matches != @exp_names ) {
			#print Dumper @exp_names;
			print "\nError: It looks like the suffixes in your model action commands do not match your .exp filenames. Please be sure that they match,  For example, for every .exp file used with name: \"myData_A.exp\" there should be a corresponding suffix in your model command:\n\t\"simulate({suffix=>myData_A\",...});\n\n";
			$need_to_die = 1;
		}

		my @scan_suffixes;
		my @scan_parameters;
		for (@scan_parameter_lines) {
			if (/suffix=>/) {
				my ( $junk, $suffix, $morejunk ) = $_ =~ m/suffix=>(“|’|”|’|")(\w+)(“|’|”|’|")/g;
				push @scan_suffixes, $suffix;
			}
			if (/parameter=>/) {
				my ( $evenmorejunk, $parameter, $wowthatsalotofjunk ) = $_ =~ m/parameter=>(“|’|”|’|")(\w+)(“|’|”|’|")/g;
				push @scan_parameters, $parameter;
			}
		}
		unless ( scalar(@scan_suffixes) == scalar(@scan_parameters) ) {
			print "\nError: It looks like your parameter_scan commands don't contain both a suffix and parameter. Please see the BioNetFit documentation\n\n";
			$need_to_die = 1;
		}

		die
		  if $need_to_die;

		for ( my $i = 0 ; $i < @scan_suffixes ; $i++ ) {
			#print $scan_suffixes[$i];
			$values_ref->{scan_parameter}[$i] = "$scan_suffixes[$i] $scan_parameters[$i]";
		}
	}

	# If user didn't tell us what type of cluster they're using, let's figure it out ourselves
	if ( $values_ref->{use_cluster} && !defined $values_ref->{cluster_software} && $subcommand ne $monitor_command ) {
		print "Figuring out what type of cluster we're on...\n"
		  if ( $verbosity >= 3 );

		my @slurm_tests;
		my @torque_tests;
		my @ge_tests;

		push @slurm_tests, (`which sbatch 2> /dev/null`);

		if (@slurm_tests) {
			$values_ref->{cluster_software}      = "slurm";
			$values_ref->{run_job_from_worknode} = 1;
		}

		# If any of these binaries exist, we're using torque
		push @torque_tests, ( `which maui 2> /dev/null`, `which moab 2> /dev/null` ) unless ( $values_ref->{cluster_software} );

		if (@torque_tests) {
			$values_ref->{cluster_software}      = "torque";
			$values_ref->{run_job_from_worknode} = 1;
		}

		# If any of these binaries exist, we're using grid engine
		push @ge_tests, ( `which sge_execd 2> /dev/null`, `which qconf 2> /dev/null`, `which qmon 2> /dev/null`, `which qhost 2> /dev/null` ) unless ( $values_ref->{cluster_software} );

		if (@ge_tests) {
			$values_ref->{cluster_software}      = "ge";
			$values_ref->{run_job_from_worknode} = 0;

			die "Error: You since you're using an SGE-based cluster, you must specify a parallel environment using the 'pe_name' option in your .conf file\n"
			  unless $values_ref->{pe_name};
		}

		# If we can't figure it out, ask the user
		if ( !defined $values_ref->{cluster_software} ) {
			my $answer;
			do {
				print "BioNetFit was unable to determine if your cluster runs SLURM, TORQUE or Grid Engine software. Please enter (S) for SLURM, (T) for a torque-based cluster, or (G) for a grid-engine based cluster. If you are unsure, try (T). You can set this option permanently using 'cluster_software={torque,ge} in your .conf file. ";
				$answer = <STDIN>;
				chomp $answer;
			} until ( ( $answer eq 'T' ) || ( $answer eq 't' ) || ( $answer eq 'G' ) || ( $answer eq 'g' ) || ( $answer eq 'S' ) || ( $answer eq 's' ) );

			if ( $answer eq "s" || $answer eq "S" ) {
				$values_ref->{cluster_software}      = "slurm";
				$values_ref->{run_job_from_worknode} = 1;
			}
			elsif ( $answer eq "T" || $answer eq "t" ) {
				$values_ref->{cluster_software}      = "torque";
				$values_ref->{run_job_from_worknode} = 1;
			}
			else {
				$values_ref->{cluster_software}      = "ge";
				$values_ref->{run_job_from_worknode} = 0;
			}
		}

		# If we're using torque, check to see if we're using maui or moab.  This is necessary because our CPU requirement directives
		# are different for maui.
		if ( $values_ref->{cluster_software} eq "torque" ) {
			my $maui_test = `which maui &> /dev/null`;
			$values_ref->{maui} = 1
			  unless $maui_test =~ / no maui /;
		}
	}

	# Choose the correct submission command based on cluster software
	if ( $values_ref->{use_cluster} && $subcommand ne $monitor_command ) {
		if ( $values_ref->{cluster_software} eq "slurm" ) {
			$values_ref->{cluster_command} = "sbatch";
		}
		else {
			$values_ref->{cluster_command} = "qsub";
		}
	}

	# Add ending slashes to output directory if user failed to do it
	if ( substr( $values_ref->{output_dir}, -1, 1 ) ne "/" ) {
		$values_ref->{output_dir} .= "\/";
	}

	my $current_directory  = $values_ref->{output_dir};
	my $previous_directory = $values_ref->{output_dir};
	my $next_directory;

	$current_directory  .= $values_ref->{job_name};
	$previous_directory .= $values_ref->{job_name};

	$next_directory = $current_directory;
	$current_directory .= "/$values_ref->{generation}";
	my $gennum = $values_ref->{generation};
	$next_directory .= "/" . ( $values_ref->{generation} + 1 );

	# Get our directory vars right for the current generation
	if ( $values_ref->{generation} && ( $values_ref->{generation} > 1 ) ) {
		$previous_directory .= "/" . ( $values_ref->{generation} - 1 );
	}
	else {
		$previous_directory = "";
	}

	# Remove path from config file if it is present.
	my $stripped_config_file = $orig_config_file;
	my $outdir               = $values_ref->{output_dir};
	chop($outdir);
	$stripped_config_file =~ s/$values_ref->{job_name}(\/[\d]+)?\///;
	$stripped_config_file =~ s/$outdir(\/[\d]+)?\///;

	if ($subcommand) {

		# If we're submitting a job (the default command)
		if ( $subcommand eq $submit_command ) {
			print "Doing configuration consistency checks and making sure all required files are in place...\n"
			  if ( $verbosity >= 3 );
			  
			my @exp_files = @{ $values_ref->{exp_file} };

			# Load exp file so we can check for formatting errors
			load_exp( \@exp_files, $values_ref );

			# If this is truly our first run, and not a bootstrap re-run
			if ( -f ".lock_$values_ref->{job_name}" && $values_ref->{ask_overwrite} ) {
				my $answer;
				do {
					print "", "You're trying to overwrite a job which might still be in progress. Do you really want to do this? ";
					$answer = <STDIN>;
					chomp $answer;
				} until ( ( $answer eq 'y' ) || ( $answer eq 'n' ) || ( $answer eq 'Y' ) || ( $answer eq 'N' ) );

				die "Please change your job name and/or starting directory.\n"
				  if ( ( $answer eq 'N' ) || ( $answer eq 'n' ) );
				unlink(".lock_$values_ref->{job_name}");
			}

			if ( $values_ref->{first_gen_permutations} ) {
				if ( $values_ref->{first_gen_permutations} <= $values_ref->{permutations} ) {
					die "Error: You set 'first_gen_permutations' lower than 'permutations' in your .conf file but 'first_gen_permutations' must be higher than 'permutations'. Please comment or increase 'first_gen_permutations' in your .conf file.\n";
				}
			}

			if ( $values_ref->{use_cluster} ) {
				my $perms = $values_ref->{permutations};
				$perms = $values_ref->{first_gen_permutations} if ( $values_ref->{first_gen_permutations} );
				my $estimated_job_count = ( ( $perms * $values_ref->{smoothing} ) / ( $values_ref->{cluster_parallel} * $values_ref->{multisim} ) + 2 );
				if ( $estimated_job_count >= $values_ref->{job_limit} ) {
					die "Error: The current values of 'job_limit' is lower than the anticipated number of jobs needed ($estimated_job_count). Please adjust 'job_limit', 'permutations', 'smoothing', 'cluster_parallel', and/or 'multisim' to fix this\n";
				}
			}

			# Be sure bootstrap_chi is set if we are bootstrapping
			die "Error: You told BioNetFit to bootstrap but did not set 'bootstrap_chi' in your .conf file. Please turn off bootstrapping or set bootstrap_chi.\n"
			  if ( $values_ref->{bootstrap} && !defined $values_ref->{bootstrap_chi} );

			die "Error: You told BioNetFit to email about job status, but didn't supply an email address in your .conf file. Please do this using the 'email' option.\n"
			  if ( ( $values_ref->{email_after_generation} || $values_ref->{email_after_run} || $values_ref->{email_on_abort} ) && !defined $values_ref->{email} );

			# Make sure all the required files are in place.
			warn "Can't find BioNetFit.pbs. This will prevent BioNetFit from running on a cluster.\n"
			  unless ( -f $pbs_generic_script );

			# Check to see if BNG2.pl exists
			die "Error: Unable to find BioNetGen executable. Please fix 'bng_command' in your .conf file\n\n"
			  unless ( -f $values_ref->{bng_command} );

			# Check to see if qsub exists, but only if the user provided a full path
			if ( $values_ref->{use_cluster} && ( fileparse( $values_ref->{cluster_command} ) ne $values_ref->{cluster_command} ) ) {
				die "Error: Unable to find Qsub executable. Please fix 'cluster_command' in your .conf file\n\n"
				  unless ( -f $values_ref->{cluster_command} );
			}

			# If our output directory already exists, see if we should overwrite it.
			if ( -d "$values_ref->{output_dir}$values_ref->{job_name}" ) {
				if ( $values_ref->{ask_overwrite} ) {
					my $answer;
					do {
						print "", "Directory $values_ref->{output_dir}$values_ref->{job_name} already exists. Overwrite? ";
						$answer = <STDIN>;
						chomp $answer;
					} until ( ( $answer eq 'y' ) || ( $answer eq 'n' ) || ( $answer eq 'Y' ) || ( $answer eq 'N' ) );

					die "Please change your job name and/or starting directory.\n"
					  if ( ( $answer eq 'N' ) || ( $answer eq 'n' ) );
				}

				if ( $verbosity >= 2 && $values_ref->{bootstrap_num} < 2 ) {
					print "Deleting old files. This may take a minute.\n";
				}
				elsif ( $verbosity >= 2 && $values_ref->{bootstrap_num} >= 2 ) {
					show_output "Deleting old files. This may take a minute.\n";
				}

				# Deleting a directory full of many files can take a long time. It's quicker to use a system command to move the directory,
				# then fork a process that will stay to delete everything in it.
				verbose_mkdir("$values_ref->{output_dir}/to_be_deleted")
				  unless ( -d "$values_ref->{output_dir}/to_be_deleted" );

				my $command = "mv $values_ref->{output_dir}/$values_ref->{job_name} $values_ref->{output_dir}/to_be_deleted";
				$command .= "&& mv values_ref->{output_dir}/$values_ref->{job_name}_cluster_output $values_ref->{output_dir}/to_be_deleted"
				  if ( $values_ref->{use_cluster} && $values_ref->{bootstrap_num} == 1 && -d "$values_ref->{output_dir}$values_ref->{job_name}_cluster_output" && $values_ref->{save_cluster_output} );

				if ( -d "$values_ref->{output_dir}/to_be_deleted/$values_ref->{job_name}" ) {
					my $pid = fork();
					if ( $pid == 0 ) {

						# Remove the old output directory
						rmtree("$values_ref->{output_dir}/to_be_deleted")
						  or warn "Warning: Could not delete $values_ref->{output_dir}$values_ref->{job_name}. This may result in errors when running BioNetFit. Please be sure none of your old output files are currently in use or opened. If using qsub, check qstat -u [username] to see if you have any jobs still running.\n";
						exit;
					}
					waitpid( $pid, WNOHANG );
				}
				else {
					rmtree("$values_ref->{output_dir}/$values_ref->{job_name}");
					rmtree("$values_ref->{output_dir}/$values_ref->{job_name}_cluster_output")
					  if ( $values_ref->{use_cluster} && $values_ref->{bootstrap_num} == 1 && -d "$values_ref->{output_dir}$values_ref->{job_name}_cluster_output" && $values_ref->{save_cluster_output} );
				}

				# Remove old bootstrap directory ( If this truly is our first run, and not a bootstrap run or redo)
				if ( $values_ref->{bootstrap} && $values_ref->{bootstrap_num} == 1 && $orig_config_file !~ /bootstrap/ ) {
					rmtree("$values_ref->{output_dir}$values_ref->{job_name}_bootstrap")
					  or warn "Warning: Could not delete $values_ref->{output_dir}$values_ref->{job_name}_bootstrap. This may result in errors when running BioNetFit. Please be sure none of your old output files are currently in use or opened. If using qsub, check qstat -u [username] to see if you have any jobs still running.\n";
				}

				# Make sure old output directory was actually deleted
				if ( -d "$values_ref->{output_dir}$values_ref->{job_name}" ) {
					warn "Warning: Could not delete $values_ref->{output_dir}$values_ref->{job_name}. This may result in errors when running BioNetFit. Please be sure none of your old output files are currently in use or opened. If using qsub, check qstat -u [username] to see if you have any jobs still running.\n";
				}
			}

			# Create a new qsub output directory
			verbose_mkdir("$values_ref->{output_dir}$values_ref->{job_name}_cluster_output")
			  if ( $values_ref->{use_cluster} && $values_ref->{bootstrap_num} == 1 && $values_ref->{save_cluster_output} && !-d "$values_ref->{output_dir}$values_ref->{job_name}_cluster_output" );

			# Create a new bootstrap output directory
			verbose_mkdir("$values_ref->{output_dir}$values_ref->{job_name}_bootstrap")
			  if ( $values_ref->{bootstrap} && !-d "$values_ref->{output_dir}$values_ref->{job_name}_bootstrap" );

			print "Generating new bootstrap map\n"
			  if ( $verbosity >= 3 );

			# Generate our bootstrap map and save to file
			generate_bootstrap_file($values_ref)
			  if ( $values_ref->{bootstrap} );

			print "Creating .lock file\n"
			  if ( $verbosity >= 3 );

			open( TMP, ">>.lock_$values_ref->{job_name}" );
			close TMP;

			show_output "Submitting job...\n"
			  if ( $verbosity >= 2 );

			# Finally, run the submission function
			run_submit( $values_ref, $current_directory, $orig_config_file );
		}
		elsif ( $subcommand eq $resume_command ) {
			print "Attemping to resume fitting...\n"
			  if ( $verbosity >= 1 );

			verbose_mkdir( $values_ref->{output_dir} )
			  unless ( -d $values_ref->{output_dir} );

			verbose_mkdir("$values_ref->{output_dir}/$values_ref->job_name}")
			  unless ( -d "$values_ref->{output_dir}/$values_ref->{job_name}" );

			opendir( CUR, "$values_ref->{output_dir}/$values_ref->{job_name}/" )
			  or die "Can't open dir $values_ref->{output_dir}/$values_ref->{job_name}/: $!\n";

			my @files = readdir(CUR);
			closedir(CUR);

			my @generations;
			my $current_config;
			my $resume_generation;

			for (@files) {
				push @generations, $_
				  if ( looks_like_number($_) );
			}

			@generations = sort { $b <=> $a } @generations;

			my $orig_config;
			if ( $values_ref->{bootstrap} ) {
				$orig_config = fileparse( $values_ref->{orig_config} );
			}
			else {
				$orig_config = fileparse($orig_config_file);
			}

			if ( !defined $orig_config ) {
				die "Error: The run you're trying to resume doesn't have the data necessary to resume the run. You will need to start over.\n";
			}

			if (@generations) {
				for ( my $i = 0 ; $i < @generations ; $i++ ) {
					if ( -f "$values_ref->{output_dir}/$values_ref->{job_name}/$generations[$i]/$orig_config" || -f "$values_ref->{output_dir}/$values_ref->{job_name}/$generations[$i]/bootstrap.conf" ) {
						$current_config = "$values_ref->{output_dir}/$values_ref->{job_name}/$generations[$i]/$orig_config";
						$current_config = "$values_ref->{output_dir}/$values_ref->{job_name}/$generations[$i]/bootstrap.conf"
						  unless ( -f $current_config );
						$resume_generation = $generations[$i];
						last;
					}
				}
			}
			if ( $values_ref->{bootstrap} && !defined $current_config ) {
				$current_config = "$values_ref->{output_dir}/$values_ref->{job_name}_bootstrap/bootstrap.conf"
				  if ( -f "$values_ref->{output_dir}/$values_ref->{job_name}_bootstrap/bootstrap.conf" );
				$resume_generation = 1;
			}

			if ( $resume_generation == 1 ) {
				$values_ref->{generation} = 1;
			}

			if ( !defined $current_config ) {
				die "Error: The run you're trying to resume doesn't have the data necessary to resume the run. You will need to start over.\n";
			}

			my $resume_directory = "$values_ref->{output_dir}$values_ref->{job_name}/$resume_generation";

			verbose_mkdir($resume_directory)
			  unless ( -d $current_directory );
			
			my $output = "Found config file. Resuming run from generation $resume_generation.\n";
			$output = "Found config file. Resuming run from generation $resume_generation, bootstrap #$values_ref->{bootstrap_num}.\n"
			  if ( $values_ref->{bootstrap} );
			
			my @config_files   = ($current_config);
			my $values_ref_new = parse_file(@config_files);
			update_values($values_ref_new);

			generate_bootstrap_file($values_ref)
			  if ( $values_ref_new->{bootstrap} && $resume_generation == 1 );

			$values_ref_new->{generation}            = $resume_generation;
			$values_ref_new->{cluster_software}      = $values_ref->{cluster_software};
			$values_ref_new->{run_job_from_worknode} = $values_ref->{run_job_from_worknode};
			$values_ref_new->{max_generations}		 = $new_max_generations if $new_max_generations;
			 
			if ($new_max_generations) {
				open CURR_CONF, '<', $current_config
				  or show_output "Warning: Can't open file: $current_config to set new max_generations value\n";
				my @lines = <CURR_CONF>;
				close CURR_CONF;
				
				open CURR_CONF, ">", $current_config
				  or show_output "die", "Warning: Can't open file $current_config to set new max_generations value\n";

				for (@lines) {
					s/max_generations \d+/max_generations $new_max_generations/;
					print CURR_CONF $_;
				}
				close CURR_CONF;
				
				$output .= "Setting the new maximum number of generations to $new_max_generations\n";
			}
			
			print $output
			  if ( $verbosity >= 1 );
			
			# Delete files in old target directory (except the config file)
			opendir( CUR, "$resume_directory/" ) or show_output "Can't open dir $resume_directory/ to delete files from interrupted run.\n";
			@files = readdir(CUR);
			closedir(CUR);

			print "Deleting old simulation files in $resume_directory\n"
			  if ( $verbosity >= 2 );

			foreach (@files) {
				if ( -f "$resume_directory/$_" && $_ !~ /.conf/ ) {

					#show_output "Deleting $current_directory/$_...\n";
					unlink "$resume_directory/$_" or show_output "Can't delete $resume_directory/$_: $!\n";
				}
			}

			# Delete all generation directories past the generation we are resuming into
			opendir( CUR, "$values_ref->{output_dir}$values_ref->{job_name}/" ) or show_output "Can't open $values_ref->{output_dir}$values_ref->{job_name}/ to delete data from interrupted run.\n";
			@files = readdir(CUR);
			closedir(CUR);

			foreach (@files) {
				next unless ( looks_like_number($_) && $_ > $resume_generation );
				rmtree "$values_ref->{output_dir}$values_ref->{job_name}/$_"
				  or print "Can't remove $values_ref->{output_dir}/$_. This may cause issues with your resumed run. Be sure none of the files are being used and try again.\n";
			}

			# Finally, run the submission function
			run_resume( $values_ref_new, $current_config );
		}
		elsif ( $subcommand eq $monitor_command ) {
			run_monitor($values_ref);
		}

		# If we're running the generation command.. (this is run to begin each new generation)
		elsif ( $subcommand eq $generation_command ) {

			# Create the directory that this generation's output will fall into
			verbose_mkdir($current_directory)
				unless (-d $current_directory);

			# Give the user option to run extra commands before we start the new generation
			if ( $values_ref->{command_extra} ) {
				system( $values_ref->{command_extra} )
				  or show_output "Couldn't execute command_extra.  Be sure your path and permissions are correct\n";
			}

			# If user specifies, delete unnecessary files as we go so as to preserve disk space
			if ( $values_ref->{delete_old_files} && exists $values_ref->{generation} && $values_ref->{generation} > 1 && ( -d $previous_directory ) ) {
				opendir( CUR, "$previous_directory/" )
				  or show_output "Can't open dir $previous_directory/ to delete previous generation's files.\n";
				my @files = readdir(CUR);
				closedir(CUR);

				for (@files) {
					if (/\.BNG_OUT$|\.xml$|\.cdat$|\.finished$|\.species$/) {
						unlink "$previous_directory/$_"
						  or show_output "Couldn't delete $_ from previous directory.\n";
					}
					elsif (/._\d+\.gdat$/) {
						unlink "$previous_directory/$_"
						  or show_output "Couldn't delete $_ from previous directory.\n";
					}
					elsif ( -d "$previous_directory/$_" && $_ ne ".." and $_ ne "." ) {
						rmtree "$previous_directory/$_";
					}
				}
				if ( $values_ref->{smoothing} > 1 ) {
					for ( 0 .. ($values_ref->{smoothing} - 1) ) {
						rmtree("$previous_directory/$_")
						  if ( -d "$previous_directory/$_" );
					}
				}
			}

			# And start the generation function
			run_generation( $values_ref, $current_directory, $next_directory, $orig_config_file, 0 );
		}

		# If we're running the genetic command, let's get everything prepped.
		elsif ( $subcommand eq $genetic_command ) {
			verbose_mkdir($current_directory);
			$values_ref->{generation}++;

			# Create new configuration file name.
			my $old_config      = fileparse( ${stripped_config_file} );
			my $new_config_file = "$values_ref->{output_dir}$values_ref->{job_name}/" . "$values_ref->{generation}/$old_config";

			# And run the genetic function
			run_genetic( $values_ref, $current_directory, $previous_directory, $new_config_file, $orig_config_file );
		}

		# If none of the above commands, let's prepare for output analysis/consolidation
		else {
			if ( $values_ref->{first_gen_permutations} && $values_ref->{generation} == 1 ) {
				$values_ref->{permutations} = $values_ref->{first_gen_permutations};
			}

			# If we have the analyze command, start analysis function
			if ( $subcommand eq $analyze_command ) {
				run_analyze( $values_ref, $current_directory, $values_ref->{permutations}, $orig_config_file );
			}

			elsif ( $subcommand eq $consolidate_command ) {
				my $message = "Preliminary results have been placed in $values_ref->{output_dir}$values_ref->{job_name}/Results\n";
				run_consolidate( $values_ref, $orig_config_file, $message );
			}
		}
	}
	else {
		show_output "die", "No subcommand given for $0.\n";
	}
	exit 0;
}

sub welcome_message($) {
	my ($values_ref) = @_;

	my $parallel = ( $values_ref->{use_cluster} ) ? $values_ref->{cluster_parallel} : $values_ref->{parallel_count};

	# Fill output variables with relevant info
	my $max_generations = $values_ref->{max_generations};
	my $max_parents     = $values_ref->{max_parents};
	my $model           = $values_ref->{model};
	my $walltime        = $values_ref->{max_walltime}
	  if ( exists $values_ref->{max_walltime} );
	my $jobpath    = "$values_ref->{output_dir}$values_ref->{job_name}";
	my $make_plots = ( $values_ref->{make_plots} ) ? "Yes" : "No";
	my $smoothing  = ( $values_ref->{smoothing} > 1 ) ? "$values_ref->{smoothing} runs" : "No";

	#my $bootstrapping = ($values_ref->{bootstrap}) ? "$values_ref->{bootstrap} times" : "No";

	# And output
	show_output "\n===========================================================\n";
	show_output "Running BioNetFit $version with the following options:\n\n";
	show_output "Job path: $jobpath\n";
	show_output "Model: $model\n";
	for ( @{ $values_ref->{exp_file} } ) {
		show_output "Data: $_\n";
	}
	show_output "Max walltime: $walltime\n"
	  if ($walltime);

	show_output "Parallel count: $parallel\n";
	show_output "Generations: $max_generations\n";
	show_output "First Gen Permutations: $values_ref->{first_gen_permutations}\n"
	  if ( $values_ref->{first_gen_permutations} );
	show_output "Permutations: $values_ref->{permutations}\n";
	show_output "Multisim: $values_ref->{multisim}\n"
	  if ( $values_ref->{multisim} > 1 && $values_ref->{use_cluster} );
	show_output "Max parents: $max_parents\n";
	show_output "Keep $values_ref->{keep_parents} parents\n"
	  if ( $values_ref->{keep_parents} );
	show_output "Smoothing: $smoothing\n";

	#show_output "Bootstrapping? $bootstrapping\n";
	show_output "Objective function: $values_ref->{objfunc}\n";
	show_output "Creating graphs?: $make_plots\n";
	my $timestring = "Start timestamp: " . time() . "\n";
	show_output $timestring;
	show_output "============================================================\n\n";
}

sub run_submit($$$) {
	show_output "\n*** Function: run_submit() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $current_directory, $orig_config_file ) = @_;

	# Just submitted the job, so let's show the welcome message
	welcome_message($values_ref)
	  if ( $values_ref->{show_welcome_message} );

	# Create the command to run the next generation
	my $command = "perl " . abs_path($0) . " $generation_command $orig_config_file";

	# If we're chaining qsub, create the command to run the first generation
	if ( $values_ref->{use_cluster} ) {
		if ( $values_ref->{run_job_from_worknode} ) {
			verbose_mkdir($current_directory);
			$command = create_pbs_command( $pbs_generic_script, "G1", $values_ref, { command => $command } );
			my $job = submit_qsub_job( $values_ref, $command );

			exit 0
			  if ( $values_ref->{bootstrap} && $values_ref->{bootstrap_num} > 1 );
			run_monitor($values_ref);
		}
		else {
			if ( fork() == 0 ) {
				( $values_ref->{command_success} == system($command) )
				  or show_output "die", "Error: Unable to run $command\n";
			}
			if ( fork() == 0 ) {
				print "Forking process to monitor and display job output...\n"
				  if ( $verbosity >= 3 );
				run_monitor($values_ref);
			}
			exit 0;
		}

	}
	else {
		show_output "\nStarting new generation with command: $command\n"
		  if ( $verbosity >= 3 );

		if ( $values_ref->{use_cluster} && $values_ref->{bootstrap_num} == 1 ) {
			if ( fork() == 0 ) {
				print "Forking process to monitor and display job output...\n"
				  if ( $verbosity >= 3 );
				run_monitor($values_ref);
			}
		}
		( $values_ref->{command_success} == system($command) )
		  or show_output "die", "Error: Unable to run $command\n";
	}
}

sub run_resume($$) {
	show_output "\n*** Function: run_resume() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $orig_config_file ) = @_;

	# Create the command to run the next generation
	my $command = "perl $0 $generation_command $orig_config_file";

	# If we're chaining qsub, create the command to run the first generation
	if ( $values_ref->{use_cluster} && $values_ref->{run_job_from_worknode} ) {
		$command = create_pbs_command( $pbs_generic_script, "G$values_ref->{generation}", $values_ref, { command => $command } );
	}

	# Run the command.  If we're chaining qsub, this will submit a new qsub job.
	#if ($values_ref->{use_cluster}) {
	#	my $pid = fork();
	#	if ($pid == 0) {
	#		run_monitor($values_ref);
	#		exit 0;
	#	}
	#}

	# Run the command.  If we're chaining qsub, this will submit a new qsub job.
	if ( $values_ref->{use_cluster} && $values_ref->{run_job_from_worknode} ) {
		my $job = submit_qsub_job( $values_ref, $command );
	}
	else {
		( $values_ref->{command_success} == system($command) )
		  or show_output "die", "Error: Unable to run $command\n";
	}
}

sub run_analyze($$$$) {
	show_output "\n*** Function: run_analyze() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $current_directory, $perm_count, $orig_config_file ) = @_;

	my $newcolumns_ref;

	if ( defined $values_ref->{newcolumn} ) {
		$newcolumns_ref = store_newcolumns( $values_ref->{newcolumn} );
	}

	my @exp_files       = @{ $values_ref->{exp_file} };
	my %model_diff      = ();
	my @perm_diff       = ();
	my @perm_model_diff = ();
	my @var_names       = ();
	my @var_sets        = ();

	# 	@perm_diff structure
	#	[
	#		{
	#			value 	=> 4513.53
	#			perm 	=> 1
	#		}
	#		{
	#			value 	=> 4513.53
	#			perm 	=> 2
	#		}
	#	]
	#
	# 	%model_diff
	#	{
	#		model.bngl =>
	#			[
	#				[				#timepoint 0
	#					1234.53			#perm 1
	#					5678			#perm 2
	#					123				#perm 3
	#				]
	#				[				#timepoint 1
	#					1234.53
	#					5678
	#					123
	#				]
	#			]
	#	}
	#
	#	@perm_model_diff
	#	[
	#		{								#perm 1
	#			model.bngl	=>	45688897	#total chi
	#		}
	#		{								#perm 2
	#			model.bngl	=>	98746652	#total chi
	#		}
	#		{								#perm 3
	#			model.bngl	=>	45688897	#total chi
	#		}
	#	]

	show_output "Loading experiment data...\n"
	  if ( $verbosity >= 3 );

	# Load all experimental data.
	my ( $expdata, $expcolumns ) = load_exp( \@exp_files, $values_ref );
	my @exp_roots = keys %$expdata;

	#print Dumper @model_roots;
	# Create new columns in experimental data if needed.
	# Not sure why new column would be needed.  I've never had any reason to use it
	# and don't know exactly what it does.
	#foreach my $root_filename (keys %$expdata) {
	#my $expfile = $values_ref->{exp_file};

	#if (exists $newcolumns_ref->{$expfile}) {
	#	create_columns($expfile, $root_filename, \@model_roots, $expdata, $expcolumns, $newcolumns_ref->{$expfile});
	#}
	#}

	show_output "Loading model data and calculating fits...\n"
	  if ( $verbosity >= 3 );

	# Load all model data, apply appropriate factors and calculate differences.
	for ( my $perm = 0 ; $perm < $perm_count ; $perm++ ) {
		my $perm_p1 = $perm + 1;
		my %data    = ();
		my %columns = ();

		my $tempfile = fileparse( $values_ref->{model} );

		$tempfile =~ s/\.bngl$/_perm${perm_p1}\.bngl/;
		$tempfile =~ s/\.txt$/_perm${perm_p1}\.txt/;

		# Load parameter values from one of the models.
		my $param_filename = "${current_directory}/$tempfile";

		my $failed_file = $param_filename;
		$failed_file =~ s/\.bngl$/\.failed/;
		$failed_file =~ s/\.txt$/\.failed/;

		# If a run failed, give it a realllly high chi^2 so it doesn't get used in breeding.
		if ( -f $failed_file ) {
			show_output "Skipping analysis of $param_filename because run failed.\n"
			  if ( $verbosity >= 1 );

			$perm_diff[$perm]{perm} = $perm;

			$perm_diff[$perm]{value} = 2**53;

			foreach my $root_filename (@exp_roots) {
				$perm_model_diff[$perm]{$root_filename} = 2**53;
				$model_diff{$root_filename}[1][$perm] = 2**53;
			}
			push @var_sets, [@var_names];
			next;
		}

		show_output "die", "Error: Cannot find $param_filename with .txt or .bngl extension!\n"
		  unless ( -f $param_filename );

		my ( $names_ref, $vars_ref ) = get_model_params( $param_filename, $values_ref );
		# $names_ref and $vars_ref contain the variables and values from the given model file (only one model is being looked at)

		if (@var_names) {
			my $string = "Error: Number of variables does not match for $param_filename.\n var_names:\n" . Dumper(@var_names) . "names_ref:\n" . Dumper(@$names_ref);
			show_output "die", $string
			  unless ( @var_names == @$names_ref );
			for ( my $k = 0 ; $k < @var_names ; $k++ ) {
				show_output "die", "Error: Variable names do not match.\n"
				  unless ( $var_names[$k] eq $names_ref->[$k] );
			}
		}
		else {
			# Var_names not initialized yet.
			@var_names = @$names_ref;
		}
		push @var_sets, $vars_ref;

		# Load output for current permutation.
		# This loops once for each of the models being used
		foreach my $root_filename (@exp_roots) {
			my $filetemp = fileparse( ${root_filename} );
			$filetemp .= "_perm";

			my $perm_filename = "${current_directory}/$filetemp${perm_p1}.gdat";

			my ( $data_ref, @columns ) = load_data($perm_filename);
			$data{$root_filename} = $data_ref;
		}

		# Apply appropriate factors and calculate differences.
		# This loop only runs once if we are using a single model

		my $model_name = fileparse( $values_ref->{model} );
		$model_name =~ s/.\w+$//;

		foreach my $root_filename (@exp_roots) {
			my $perm_filename = "${current_directory}/${root_filename}_perm${perm_p1}.gdat";

			# Check to see if model has enough points to compute value
			my @exp = @{ $expdata->{$root_filename} };
			my @sim = @{ $data{$root_filename} };

			my $filename = fileparse($root_filename);
			$filename =~ s/${model_name}_//;

			# Control column will not be 'time' in dose-response data
			my $control_col = "time";
			if ( $values_ref->{scan_parameter} ) {
				foreach ( @{ $values_ref->{scan_parameter} } ) {
					if ( ( split( ' ', $_ ) )[0] eq $filename ) {
						$control_col = ( split( ' ', $_ ) )[1];
					}
				}
			}

			# Fill our random map for bootstrapping
			my @bootstrap = ( {} );
			if ( $values_ref->{bootstrap} ) {
				open EXP_DATA, '<', "$values_ref->{output_dir}$values_ref->{job_name}_bootstrap/$filename.dat"
				  or show_output "die", "Error: Can't open file: $values_ref->{output_dir}$values_ref->{job_name}_bootstrap/$filename.dat\n";
				my @lines = <EXP_DATA>;
				close EXP_DATA;
				for ( my $i = 0 ; $i < @exp ; $i++ ) {
					my @vars = split( ' ', $lines[$i] );
					my $j = 0;
					foreach my $col ( sort numalpha_sort @{ $expcolumns->{$root_filename} } ) {
						next if ( $col eq $control_col || $col =~ /_SD$/ );
						$bootstrap[$i]{$col} = $vars[$j];
						$j++;
					}
				}
			}

			# This checks to make sure our simulation has at least as many timepoints as our experimental data,
			# and skips the permutation (via high chi^2) if it doesn't.
			if ( ( $exp[ scalar(@exp) - 1 ]{$control_col} - $tolerance ) > $sim[ scalar(@sim) - 1 ]{$control_col} ) {
				my $string = "Final time in ${perm_filename}[$perm] (" . $sim[ scalar(@sim) - 1 ]{$control_col} . ") is less than final time in experimental data (" . $exp[ scalar(@sim) - 1 ]{$control_col} . "). [skipping $perm_filename]\n";
				show_output $string;

				# 2^53 is maximum integer size on 32-bit systems.
				$perm_diff[$perm]{perm}                 = $perm;
				$perm_diff[$perm]{value}                = 2 ^ 53;
				$perm_model_diff[$perm]{$root_filename} = 2 ^ 53;
				$model_diff{$root_filename}[0][$perm]   = 2 ^ 53;
				next;
			}

			# Create new columns before scaling factors are applied to any model.
			if ( exists $newcolumns_ref->{$root_filename} ) {
				create_columns( $perm_filename, $root_filename, \@exp_roots, \%data, \%columns, $newcolumns_ref->{$root_filename} );
			}

			# Check for divide by zero if we are dividing by equilibrium value
			#if ( $values_ref->{divide_by_init} ) {
			#	foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
			#		next if ($col eq $control_col);
			#		show_output "die", "Error: You told BioNetFit to divide by equilibrium, but the equilbrium value for $col in $root_filename.gdat is 0. We cannot divide by zero.\n"
			#		  if ( $sim[0]{$col} == 0 );
			#	}
			#}

			my %col_means;
			my %col_pt_sd;
			my %col_sd;

			# Calculate our column averages for use in SOS calculation
			# TODO: Having the calculation done here leaves us with a lot of redundancy, and is unnecessarily resource heavy.
			# Would be nice to move it out of the perm loop.
			if ( $values_ref->{objfunc} == 4) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					next
					  if ( $col eq $control_col || $col =~ /_SD$/ );
					my $col_sum;
					my $num_points;
					for ( my $i = 0 ; $i < @exp ; $i++ ) {
						next
						  if ( $exp[$i]{$col} eq "NaN" );
						$num_points += 1;
						$col_sum += $exp[$i]{$col};
					}
					$col_means{$col} = $col_sum / $num_points;
					show_output "die", "Error: You are using objfunc=4 (divide by column average), but the average for column $col is 0. We cannot divide by zero.\n"
					  if ( $col_means{$col} == 0 && $values_ref->{objfunc} == 4);
				}

				#if ($values_ref->{objfunc} == 5) {
				#foreach my $col (@{$expcolumns->{$root_filename}}) {
				#next
				#if ($col eq $control_col || $col =~ /_SD$/);

				#my $stdev;
				#my $average = $col_means{$col};
				#my $sq_diff;
				#my $num_points;
				#for (my $i = 0; $i < @exp; $i++) {
				#next
				#if ( $exp[$i]{$col} eq "NaN");

				#$sq_diff += ($exp[$i] - $average) ** 2;
				#$num_points++;
				#}
				#$stdev = sqrt($sq_diff/$num_points);
				#$col_sd{$col} = $stdev;
				#}
				#}
			}
			if ( $values_ref->{objfunc} == 2 ) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					next unless ( $col =~ /_SD$/ );
					my $sd_arr = [];
					for ( my $i = 0; $i < @exp; $i++ ) {
						show_output "die", "You are using standard deviations in your .exp file, but we found a missing/non-numeric/zero value. All standard deviation values must be present\n"
							unless ( looks_like_number($exp[$i]{$col} && abs($exp[$i]{$col}) != 0) );
						push @$sd_arr, $exp[$i]{$col};
					}
					$col_pt_sd{$col} = $sd_arr;
				}
			}
			
			#print "SDs: \n";
			#print Dumper %col_pt_sd;
			
			#print "Unaltered sim: \n";
			#print Dumper @sim;
			
			#$skipsim = 0;
			
			# Divide all values by the values at t=0
			if ($values_ref->{divide_by_init}) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					next if ( $col =~ /_SD$/ || $col eq $control_col);
					if ($sim[0]{$col} == 0) {
						$sim[0]{$col} = 1e-6;
						#show_output "die", "Error: You told BioNetFit to divide by equilibrium, but the equilbrium value for $col is $sim[0]{$col} in $perm_filename. We cannot divide by zero. Quitting.\n";
						#skipsim = 1;
						#$sim[0]{$col} = 1;
						
					}
					for (my $i = 1; $i < @_; $i++) {
						$sim[$i]{$col} = $sim[$i]{$col}/$sim[0]{$col};
					}
					$sim[0]{$col} = 1;
				}
			}
			
			#next if ($skipsim);
			
			#print "After div_by_init: \n";
			#print Dumper @sim;
			
			# log2 transform if desired
			if ($values_ref->{log_transform_sim_data}) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					next if ( $col =~ /_SD$/ || $col eq $control_col);
					for (my $i = 0; $i < @sim; $i++) {
						if ($sim[$i]{$col} == 0) {
							$sim[$i]{$col} = 1e-6;
						}
						$sim[$i]{$col} = log($sim[$i]{$col})/log($values_ref->{log_transform_sim_data});
					}
				}
			}
			
			#next if ($skipsim);
			
			#print "After log transform: \n";
			#print Dumper @sim;
			
			# Standardize simulation output such that a column mean is equal to 0 and its SD is equal to 1
			if ($values_ref->{standardize_sim_data}) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					next if ( $col =~ /_SD$/ || $col eq $control_col);
					
					# Calculate mean of column
					my $sum;
					
					for (my $i = 0; $i < @sim; $i++) {
						$sum += $sim[$i]{$col};
					}
					
					my $average = $sum/@sim;

					#show_output "die", "You chose to standardize simulation output to mean of 0 and std of 1, but we can't do this because the mean of the column $col is equal to $average in $perm_filename. This might happen if you are dividing by equilibrium then log transforming when the simulation output values are all the exact same."
					next
						if ($average == 0);
					
					# Calculate SD of column
					my $sqtotal = 0;
					
					for (my $i = 0; $i < @sim; $i++) {
						$sqtotal += ($average - $sim[$i]{$col}) ** 2;
					}
					
					my $std = ($sqtotal / (@sim-1)) ** 0.5;
					
					for (my $i = 0; $i < @sim; $i++) {
						# Standardize
						$sim[$i]{$col} = ($sim[$i]{$col} - $average)/$std
					}
				}
			}
			
			#print "After standardization: \n";
			#print Dumper @sim;
			
			# This needs moved out of loop
			my @exp_std = ();
			
			if ($values_ref->{standardize_exp_data}) {
				foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
					
					my $sum;
					my $num_pts = 0;
					for (my $i = 0; $i < @exp; $i++) {
						next if ($exp[$i]{$col} eq "NaN");
						$sum += $exp[$i]{$col};
						$num_pts += 1;
					}
					
					my $average = $sum/$num_pts;
					
					#Calculate SD of column
					my $sqtotal = 0;
					
					#print Dumper @exp
					#	if ($average == 0);
					
					#show_output "die", "You chose to standardize .exp data to mean of 0 and std of 1, but we can't do this because the mean of the column $col is equal to $average. Quitting."
					#	if ($average == 0);
					
					for (my $i = 0; $i < @exp; $i++) {
						next if ($exp[$i]{$col} eq "NaN");

						$sqtotal += ($average - $exp[$i]{$col}) ** 2;
					}
					
					my $std = ($sqtotal / ($num_pts - 1)) ** 0.5;
					
					for (my $i = 0; $i < @exp; $i++) {
						if ($exp[$i]{$col} eq "NaN") {
							$exp_std[$i]{$col} = "NaN";
							next;
						}
						elsif ($col eq $control_col) {
							$exp_std[$i]{$control_col} = $exp[$i]{$control_col};
							next;
						}
						
						elsif ($col =~ /_SD$/) {
							$exp_std[$i]{$col} = $exp[$i]{$col};
							next;
						}
						
						if ($average == 0) {
							$exp_std[$i]{$col} = $exp[$i]{$col};
						}
						else {
							$exp_std[$i]{$col} = ($exp[$i]{$col} - $average) / $std;
						}
					}
				}
			}
			
			#if ($perm_p1 == 1) {
			#	print "Exp after standardization:\n";
			#	print Dumper @exp_std;
			#	
			#	print "Sim after standardization:\n";
			#	print Dumper @sim;
			#}
			
			#print Dumper @{ $expcolumns->{$root_filename} };
			#print Dumper @exp;

			## Here we calculuate the window average
			#elsif ($values_ref->{objfunc} == 4) {
			#foreach my $col (@{$expcolumns->{$root_filename}}) {
			#next
			#if ($col eq $control_col || $col =~ /_SD$/);
			#my $col_sum;
			#for (my $i = 0; $i < @exp; $i++) {
			#next
			#if ( $exp[$i]{$col} eq "NaN");

			#if ($i == 0) {
			#$col_means{$col}[$i] = ($exp[$i]{$col} + $exp[$i+1]{$col}) / 2;
			#}
			#elsif ($i == $#exp) {
			#$col_means{$col}[$i] = ($exp[$i]{$col} + $exp[$i-1]{$col}) / 2;
			#}
			#else {
			#$col_means{$col}[$i] = ($exp[$i]{$col} + $exp[$i-1]{$col} +  $exp[$i+1]{$col}) / 3;
			#}

			#if ($col_means{$col}[$i] == 0) {
			#my $string = "Error: You are using objfunc=4 (divide by window average), but the window average for column $col at timepoint " . ($i + 1) . " is 0. We cannot divide by zero.\n";
			#show_output "die", $string;
			#}
			#}

			#}
			#}
			
			my @exp_fit;
			if (@exp_std) {
				@exp_fit = @exp_std;
			}
			else {
				@exp_fit = @exp;
			}

			# Calculate difference sums for each time step.
			my $j = 0;    # $j loops through experimental time points
			              # $i will loop through simulation time points
			for ( my $i = 0 ; $i < @sim ; $i++ ) {

				my $line_sum = 0;

				# Stop if at end of data file.
				last
				  unless ( $exp_fit[$j] );

				# Ignore timepoints that are not included in experimental data.
				next
				  if abs( $exp_fit[$j]{$control_col} - $sim[$i]{$control_col} ) > $tolerance;

				# Ignore any columns that are not listed in the experiment data.
				foreach my $col ( sort numalpha_sort @{ $expcolumns->{$root_filename} } ) {

					# Don't include time in calculation.
					next
					  if ( $col eq "time" || $col eq $control_col || $col =~ /_SD$/ || $exp_fit[$j]{$col} eq "NaN" );

					# And show_output "die", if we're missing information
					show_output "die", "Error: Column '$col' in $perm_filename is missing data that is present in your .exp file. Check to be sure that all data points in your .exp file have a corresponding data point in your simulation output."
					  unless defined $sim[$i]{$col};

					my ( $sim_val, $exp_val );
					#if ( $values_ref->{divide_by_init} ) {
					#	$sim_val = $sim[$i]{$col} / $sim[0]{$col};
					#}
					#else {
						$exp_val = $exp_fit[$j]{$col};
						$sim_val = $sim[$i]{$col};
					#}

					# Calculate squared residuals
					# Absolute
					if ( $values_ref->{objfunc} == 1 ) {
						$line_sum += ( $exp_val - $sim_val )**2;
					}

					# Relative
					elsif ( $values_ref->{objfunc} == 3 ) {
						if ( $exp_val == 0 ) {
							show_output "die", "Error: You chose objfunc=3 which divides by .exp value, but one of your .exp values is a 0. We cannot divide by 0. Please choose a different objfunc in your .conf file\n";
						}
						$line_sum += ( ( $exp_val - $sim_val ) / $exp_val )**2;
					}

					# Div by column mean
					elsif ( $values_ref->{objfunc} == 4 ) {
						$line_sum += ( ( $exp_val - $sim_val ) / $col_means{$col} )**2;
					}

					# Div by column window mean
					#elsif ($values_ref->{objfunc} == 4) {
					#	$line_sum += (($exp_val - $sim_val) / $col_means{$col}[$j]) ** 2;
					#}
					# Div by SD
					elsif ( $values_ref->{objfunc} == 2 ) {							
						my $sd_col_name = "${col}_SD";
						show_output "die", "Can't divide by 0 in objfunc=2 (col: $sd_col_name, timepoint: $j)\n"
							if (abs($col_pt_sd{$sd_col_name}[$j]) == 0);
						$line_sum += ( ( $exp_val - $sim_val ) / $col_pt_sd{$sd_col_name}[$j] )**2;
						
						if ($perm_p1 == 1) {
							#print "col: $col\n";
							#print "sd_col: $sd_col_name\n";
							#print "sos: $exp_val - $sim_val / $col_pt_sd{$sd_col_name}[$j] ** 2 = $line_sum\n";
						}
					}

					#elsif ($values_ref->{objfunc} == 5) {
					#	$line_sum += (($sim_val - $col_means{$col}) / $col_sd{$col}) ** 2;
					#}
					else {
						# Default to absolute
						$line_sum += ( $exp_val - $sim_val )**2;
					}

					# Multiply sum by bootstrap map value if we are bootstrapping
					$line_sum = $line_sum * $bootstrap[$j]{$col}
					  if ( $values_ref->{bootstrap} );
				}
				
				if ($perm_p1 == 1) {
					#print "line sum: $line_sum\n";
				}

				# Store Results
				$model_diff{$root_filename}[$j][$perm] = $line_sum;
				$perm_diff[$perm]{perm} = $perm;
				$perm_diff[$perm]{value} += $line_sum;
				$perm_model_diff[$perm]{$root_filename} += $line_sum;
				
				# Go to next experimental timepoint.
				$j++;
			}
			if ($perm_p1 == 1) {
				#print "all sum: $perm_model_diff[$perm]{$root_filename}\n";
			}
		}
	}

	#print "perm diff:\n";
	#print Dumper @perm_diff;
	
	show_output "die", "Error: No SOS calculations could be made between .exp file and simulation data. Please be sure that simulation and experimental data have matching timepoints and matching column names\n."
	  unless @perm_diff;

	show_output "Outputting generation summary...\n"
	  if ( $verbosity >= 3 );

	# Output sorted summary of all permutations.
	my $summarydiff_filename = "${current_directory}_summary_diff.txt";

	open SUMMARY_DIFF, '>', $summarydiff_filename
	  or show_output "die", "Error: Unable to open $summarydiff_filename for writing.\n";

	# Sort to find best matches.
	my @sorted_perm_diff = sort { $a->{value} <=> $b->{value} } @perm_diff;

	show_output "Outputting best fit combined with experimental data...\n"
	  if ( $verbosity >= 3 );

	foreach my $root_filename (@exp_roots) {
		
		my ( $exp_name, $path, $suffix ) = fileparse($root_filename);
		my $best_perm  = $sorted_perm_diff[0]{perm};
		my $model_name = fileparse( $values_ref->{model} );
		$model_name =~ s/.\w+$//;

		show_output "die", "Error: Nothing found to sort, summary_diff empty?\n"
		  unless defined $best_perm;

		show_output "die", "Error: Unable to extract model from $root_filename\n"
		  unless $exp_name;

		my ($model_exp_filename) = $exp_name;
		$model_exp_filename =~ s/${model_name}_//;

		my $control_col = "time";
		if ( $values_ref->{scan_parameter} ) {
			foreach ( @{ $values_ref->{scan_parameter} } ) {
				if ( ( split( ' ', $_ ) )[0] eq $model_exp_filename ) {
					$control_col = ( split( ' ', $_ ) )[1];
				}
			}
		}

		$model_exp_filename = $path . $model_exp_filename . ".exp";

		my $gdat_filename = $exp_name . "_perm";
		$gdat_filename = "${current_directory}/$gdat_filename" . ( $best_perm + 1 ) . ".gdat";

		my $curr_gen = fileparse($current_directory);

		my $combined_filename = "${current_directory}/${exp_name}_gen${curr_gen}_perm" . ( $best_perm + 1 ) . ".txt";

		#combine_columns( $control_col, $combined_filename, $model_exp_filename, $gdat_filename );
	}

	my %name_lookup;
	for ( my $j = 0 ; $j < @var_names ; $j++ ) {
		$name_lookup{ $var_names[$j] } = $j;
	}

	# Output header
	#my @sorted_var_names = sort numalpha_sort @var_names;
	my @sorted_var_names = @var_names;

	#printf SUMMARY_DIFF "%-14s", "Permutation";
	#printf SUMMARY_DIFF " %-14s", "Chi-Sq";
	printf SUMMARY_DIFF "%-1s",  "Permutation";
	printf SUMMARY_DIFF " %-1s", "Chi-Sq";
	foreach (@sorted_var_names) {

		#printf SUMMARY_DIFF " %-14s", $_;
		printf SUMMARY_DIFF " %-1s", $_;
	}
	printf SUMMARY_DIFF "\n";

	foreach my $permref (@sorted_perm_diff) {

		#printf SUMMARY_DIFF "%-14s", "gen" . $values_ref->{generation} .
		#			 "perm" . ($permref->{perm} + 1);
		printf SUMMARY_DIFF "%-1s", "gen" . $values_ref->{generation} . "perm" . ( $permref->{perm} + 1 );

		printf SUMMARY_DIFF " %.16f", sqrt( $permref->{value} );

		# Print out parameters appropriate for current permutation.
		unless ( $var_sets[ $permref->{perm} ] ) {
			show_output "Undefined: $permref->{perm}\n";
			print Dumper @sorted_perm_diff;
			print Dumper $permref;
			print Dumper @var_sets;
		}

		my @list = @sorted_var_names;
		foreach (@list) {
			if ( looks_like_number( $var_sets[ $permref->{perm} ]->[ $name_lookup{$_} ] ) ) {
				printf SUMMARY_DIFF " %.8e", $var_sets[ $permref->{perm} ]->[ $name_lookup{$_} ];
			}
			else {
				printf SUMMARY_DIFF " %.8e", 0;
			}
		}
		print SUMMARY_DIFF "\n";
	}
	close SUMMARY_DIFF;

	# Output sorted list of each permutation with each model diff.
	my $summarymodeldiff_filename = "${current_directory}/perm_model_diff.txt";
	open SUMMARYMODEL_DIFF, '>', $summarymodeldiff_filename
	  or show_output "die", "Error: Unable to open $summarymodeldiff_filename for writing.\n";

	# Print header.
	print SUMMARYMODEL_DIFF "Permutation";
	foreach my $model ( sort numalpha_sort keys %model_diff ) {
		print SUMMARYMODEL_DIFF "\t$model";
	}
	print SUMMARYMODEL_DIFF "\n";

	foreach my $sorted_permref (@sorted_perm_diff) {
		my $permref = $perm_model_diff[ $sorted_permref->{perm} ];
		print SUMMARYMODEL_DIFF "perm", $sorted_permref->{perm} + 1;
		foreach my $model ( sort numalpha_sort keys %$permref ) {
			print SUMMARYMODEL_DIFF "\t", $permref->{$model};
		}
		print SUMMARYMODEL_DIFF "\n";
	}
	close SUMMARYMODEL_DIFF;

	show_output "Finished analysis. Cleaning up and getting ready for next step...\n"
	  if ( $verbosity >= 3 );

	# Finish up if we've reached max_generation, do boostrapping stuff if necessary
	if ( $values_ref->{generation} >= $values_ref->{max_generations} ) {

		my $exit_message;
		my $best_run_filename;
		my $chi;

		if ( ( $values_ref->{use_cluster} && $values_ref->{run_job_from_worknode} ) || !$values_ref->{use_cluster} ) {
			$exit_message = "\nReached generation $values_ref->{generation}, all done.\n";
			( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
		}
		else {
			verbose_mkdir("$values_ref->{output_dir}$values_ref->{job_name}/Results/");
		}

		if ( $values_ref->{bootstrap} && $values_ref->{run_job_from_worknode} ) {
			bootstrap_run( $values_ref, $chi, $best_run_filename );
		}
		exit 0;
	}
	else {
		my $command = "perl $0 $genetic_command $orig_config_file";

		show_output "Running genetic algorithm...\n"
		  if ( $verbosity >= 2 );
		exec($command);
	}
}

sub run_generation($$$$$) {
	show_output "\n*** Function: run_generation() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $current_directory, $next_directory, $orig_config_file, $run_sensitivity ) = @_;

	if ( $values_ref->{first_gen_permutations} && $values_ref->{generation} == 1 ) {
		$values_ref->{perms}        = $values_ref->{permutations};
		$values_ref->{permutations} = $values_ref->{first_gen_permutations};
	}

	show_output "Generating list of variable substitutions for next generation's models...\n"
	  if ( $verbosity >= 3 );

	# Generate permutations of all variable substitutions.
	my @names = ();
	my @var_sets = calculate_var_sets( \@names, $values_ref );
	
	#print "names before sorting:\n";
	#print Dumper @names;
	
	#print "values before sorting\n";
	#print Dumper @var_sets;
	
	my @idx = sort { $names[$a] cmp $names[$b] } 0 .. $#names;
	@names = @names[@idx];
	for (@var_sets) {
		@$_ = @$_[@idx];
	}

	if ( $values_ref->{first_gen_permutations} && $values_ref->{generation} == 2 ) {
		my $difference = $values_ref->{first_gen_permutations} - $values_ref->{permutations};
		splice( @var_sets, -$difference );
	}

	show_output "Generating model files for next generation...\n"
	  if ( $verbosity >= 3 );

	# Generate model files for each var_set.
	my @models     = $values_ref->{model};
	my @model_sets = ();

	foreach (@models) {
		if (@var_sets) {
			push @model_sets, generate_model_set_files( $current_directory, $_, \@names, \@var_sets, $values_ref );
		}
		else {
			show_output "die", "Error: No permutations to run. Did you set any free parameters in your .conf file?\n";
		}
	}

	show_output "die", "Error: You need more than two permutations.\n\n"
	  if ( scalar(@model_sets) <= 2 );

	# Run each model.
	my $string = "\nGeneration $values_ref->{generation}, running " . scalar(@model_sets) . " models...\n";
	show_output $string
	  if ( $verbosity >= 1 );

	my $my_analyze_command = $analyze_command;

	if ( $values_ref->{use_cluster} ) {

		my @jobs = run_models_pbs( $current_directory, \@model_sets, $values_ref );
		my $perm_count = scalar(@model_sets);

		$perm_count = $perm_count * ( $values_ref->{smoothing} );

		my $ticker = 0;
		my $failed_runs;
		my @directory_contents;

		show_output "Waiting for simulations to finish...\n"
		  if ( $verbosity >= 2 );

		if ( defined $values_ref->{max_walltime} ) {
			my @walltimesplit = split( /:/, $values_ref->{max_walltime} );
			my $walltime = ( $walltimesplit[0] * 3600 ) + ( $walltimesplit[1] * 60 ) + ( $walltimesplit[2] );

			my $numstarted  = 0;
			my $numfinished = 0;

			my %babysitter = ();

			while (1) {
				my @status_list = qw(Q qw PD);
				my $queued = check_queue( $values_ref, \@status_list );

				# If we have any queued runs...
				# We don't keep track of each sim and the chunk it's running in. So if we have any queued runs, the best thing we can do is
				# "pause" walltime of all runs by adding 5 seconds and looping again.
				if ($queued) {
					while ( my ( $key, $value ) = each %{ $babysitter{STARTED} } ) {
						$babysitter{STARTED}{$key}{LIMIT} += 5;
					}
					sleep 5;
					next;
				}

				opendir CURDIR, "${current_directory}"
				  or show_output "die", "Error: Walltime loop unable to open ${current_directory}\n";

				@directory_contents = readdir(CURDIR);
				closedir(CURDIR);

				foreach (@directory_contents) {
					my ( $basename, $path, $ext ) = fileparse( $_, qr/\.[^\.]*$/ );

					# Gather our started runs
					if ( defined $ext && $ext eq ".BNG_OUT" && $numstarted < $perm_count ) {
						if ( !exists $babysitter{STARTED}{$basename} ) {
							$numstarted++;
							$babysitter{STARTED}{$basename} = { LIMIT => ( $walltime + time() ), };
						}
					}

					# Gather our finished runs
					elsif ( defined $ext && ( $ext eq ".finished" || $ext eq ".failed" ) && $numfinished < $perm_count ) {
						if ( !exists $babysitter{FINISHED}{$basename} ) {
							$numfinished++;
							$babysitter{FINISHED}{$basename} = { END => time(), };
						}
					}
				}

				while ( my ( $key, $value ) = each %{ $babysitter{STARTED} } ) {

					# If our STARTED entry doesn't yet have a FINISHED entry, let's check how long it's been running
					if ( !exists $babysitter{FINISHED}{$key} ) {

						# Check to see if it's gone over walltime
						if ( time() > $babysitter{STARTED}{$key}{LIMIT} ) {
							$numfinished++;
							$babysitter{FINISHED}{$key} = { END => time(), };
							open( TMP, ">>$current_directory/$key.failed" );
							close TMP;
						}
					}
				}

				last
				  if ( $numfinished >= $perm_count );

				$ticker++;
				sleep(1);
			}
		}
		else {
			while (1) {
				opendir CURDIR, "${current_directory}"
				  or next;

				@directory_contents = readdir(CURDIR);
				closedir(CURDIR);

				my $finished_count = grep( /.finished$|.failed$/, @directory_contents );

				last
				  if ( $finished_count >= $perm_count );

				$ticker++;
				sleep(5);
			}
		}

		show_output "Simulations finished. Moving and renaming files...\n"
		  if ( $verbosity >= 3 );

		my @exp_files = @{ $values_ref->{exp_file} };
		foreach (@exp_files) {
			fileparse $_;
			s/.exp$//;
		}
		my $model_name = fileparse( $values_ref->{model} );
		$model_name =~ s/.\w+$//;

		# Move and rename output data
		unless ( $values_ref->{smoothing} > 1 ) {
			opendir CURDIR, "$current_directory"
			  or show_output "Unable to open $current_directory to move .gdat files\n";

			my @directory_contents = readdir(CURDIR);
			closedir(CURDIR);

			for (@directory_contents) {
				if ( /\.scan$/ | /\.gdat$/ ) {
					my $filename_1 = $_;
					$filename_1 =~ s/^${model_name}_//;
					$filename_1 =~ s/^perm\d+_//;
					$filename_1 =~ s/\.\w+$//;
					if ( grep /$filename_1/, @exp_files ) {
						my $filename = $_;
						$filename =~ s/scan$/gdat/;
						my ($perm_string) = $_ =~ /(_perm\d+)/;
						$filename =~ s/_perm\d+//;
						$filename =~ s/\.gdat$/$perm_string\.gdat/;
						move "$current_directory/$_", "$current_directory/$filename";
					}
				}
			}
		}
		else {
			for ( my $i = 0 ; $i < $values_ref->{smoothing} ; $i++ ) {
				opendir CURDIR, "$current_directory/$i"
				  or show_output "die", "Unable to open $current_directory/$i to move .gdat files\n";

				my @directory_contents = readdir(CURDIR);
				closedir(CURDIR);

				for (@directory_contents) {
					if ( /\.scan$/ | /\.gdat$/ ) {
						my $filename_1 = $_;
						$filename_1 =~ s/^${model_name}_//;
						$filename_1 =~ s/^perm\d+_//;
						$filename_1 =~ s/\.\w+$//;
						if ( grep /$filename_1/, @exp_files ) {
							my $filename = $_;
							$filename =~ s/scan$/gdat/;
							my ($perm_string) = $_ =~ /(_perm\d+)/;
							$filename =~ s/_perm\d+//;
							$filename =~ s/\.gdat$/$perm_string\.gdat/;
							$filename =~ s/\.gdat$/_$i.gdat/;

							move "$current_directory/$i/$_", "$current_directory/$filename";
						}
					}
				}
			}
		}

		show_output "Averaging simulation outputs...\n"
		  if ( $verbosity >= 2 );

		sleep 1;

		smooth_runs( $current_directory, $values_ref )
		  if ( $values_ref->{smoothing} > 1 );

		show_output "Making sure simulation outputs look okay...\n"
		  if ( $verbosity >= 3 );

		$failed_runs = grep( /.failed/, @directory_contents );

		if ( $failed_runs > ( $values_ref->{permutations} - 3 ) ) {
			if ( $values_ref->{retry_count} < $values_ref->{max_retries} ) {
				show_output "The previous generation has too many failed runs, and must be re-run. Be sure your runs aren't going over walltime, and that NFsim is completing simulations successfully (check your .BNG_OUT and .NF_OUT files for BioNetGen and NFsim output)\n";

				#($values_ref->{command_success} == system("rm -r $current_directory"))
				#	   or show_output "die", "Couldn't delete old files before re-running generation\n";

				rmtree $current_directory
				  or show_output "Warning: Couldn't delete current generation's output before re-running generation.\n";

				$values_ref->{retry_count}++;
				mkpath($current_directory);
				run_generation( $values_ref, $current_directory, $next_directory, $orig_config_file, 0 );
			}
			else {
				# If we have enough simulations (at least 3 generations) let's try to compile results from them
				if ( $values_ref->{generation} >= 3 ) {
					my $exit_message = "Tried to re-run the same generation too many times. Attempting to compile results anyway\n";
					rmtree $current_directory;
					my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
					if ( $values_ref->{bootstrap} ) {
						bootstrap_run( $values_ref, $chi, $best_run_filename );
					}
				}

				# If not, we can't do anything else and must show_output "die", .
				else {
					mkpath("$values_ref->{output_dir}$values_ref->{job_name}/Results");
					show_output "die", "Tried to re-run the same generation too many times. We don't have enough generations to compile reults. Exiting.";
				}
			}
		}

		my $analyze_name = "A$values_ref->{generation}";
		$analyze_name .= "-G" . ( $values_ref->{generation} + 1 )
		  if ( $values_ref->{generation} != $values_ref->{max_generations} );

		my $command = create_pbs_command( $pbs_generic_script, $analyze_name, $values_ref, { command => "$medir/" . fileparse($0) . " $my_analyze_command $orig_config_file" } );

		show_output "Analyzing simulation data...\n"
		  if ( $verbosity >= 2 );

		my $job = submit_qsub_job( $values_ref, $command );

		if ( !$values_ref->{run_job_from_worknode} ) {
			my $config_base = fileparse($orig_config_file);
			show_output "Waiting for analysis to finish...\n"
			  if ( $verbosity >= 3 );

			while (1) {
				if ( -d "$values_ref->{output_dir}$values_ref->{job_name}/Results/" ) {

					my $exit_message = "\nReached generation $values_ref->{generation}, all done.\n";
					my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );

					if ( $values_ref->{bootstrap} ) {

						bootstrap_run( $values_ref, $chi, $best_run_filename );
					}

					exit 0;
				}

				last
				  if ( -f "$next_directory/$config_base" );
				sleep(1);
			}

			$orig_config_file = fileparse($orig_config_file);
			$command          = "perl $0 $generation_command $next_directory/$orig_config_file";
			show_output "Starting new generation...\n"
			  if ( $verbosity >= 2 );
			( $values_ref->{command_success} == system($command) )
			  or show_output "die", "Error: Unable to run command: $command\n";
		}
	}
	else {
		my $child_status = run_models_fork( \@model_sets, $values_ref );

		if ( $child_status < 3 ) {
			if ( $values_ref->{retry_count} < $values_ref->{max_retries} ) {
				show_output "The previous generation has too many failed runs, and must be re-run. Be sure your runs aren't going over walltime, and that NFsim is completing simulations successfully (check your .BNG_OUT and .NF_OUT files for BioNetGen and NFsim output)\n\n";

				#($values_ref->{command_success} == system("rm -r $current_directory"))
				#	or show_output "die", "Couldn't delete old files before re-running generation\n";

				rmtree $current_directory
				  or show_output "Warning: Couldn't delete current generation's output before re-running generation.\n";

				$values_ref->{retry_count}++;
				mkpath($current_directory);
				run_generation( $values_ref, $current_directory, $next_directory, $orig_config_file, 0 );
			}
			else {
				# If we're re-run our generation too many times we should stop.
				# If we have enough simulations (at least 3 generations) let's try to compile results from them
				if ( $values_ref->{generation} >= 3 ) {
					my $exit_message = "Tried to re-run the same generation too many times. Attempting to compile results anyway\n";
					rmtree $current_directory;
					my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
					if ( $values_ref->{bootstrap} ) {
						bootstrap_run( $values_ref, $chi, $best_run_filename );
					}
				}

				# If not, we can't do anything else and must show_output "die", .
				else {
					mkpath("$values_ref->{output_dir}$values_ref->{job_name}/Results");
					show_output "die", "Tried to re-run the same generation too many times. We don't have enough generations to compile reults. Exiting.";
				}
			}
		}

		my $command = "perl $0 $my_analyze_command $orig_config_file";

		show_output "Analyzing simulation data...\n"
		  if ( $verbosity >= 2 );
		exec($command);
	}
}

sub load_all_summaries($) {
	show_output "\n*** Function: load_all_summaries() ***\n"
	  if ( $verbosity >= 4 );
	my ($values_ref) = @_;

	# Read in all sorted summaries.
	my $job_dir = "$values_ref->{output_dir}$values_ref->{job_name}";
	opendir JOB_DIR, $job_dir
	  or show_output "die", "Error: Unable to open the directory " . "\n";
	my @summary_files =
	  grep( /^\d+_summary_diff\.txt$/, readdir(JOB_DIR) );
	closedir(JOB_DIR);

	show_output "die", "Error: No summary files to consolidate.\n"
	  unless (@summary_files);
	my $var_names_ref;
	my @all_summaries = ();
	foreach my $file (@summary_files) {
		my $var_data_ref;
		( $var_names_ref, $var_data_ref ) = read_summary("$job_dir/$file");
		$file =~ /^(\d+)_/;
		if ($1) {    #TODO What is this?
			         #map {$_->[0] .= "_gen$1"} @$var_data_ref;
		}
		push @all_summaries, @$var_data_ref;
	}

	# Sort to get lowest chi-squared values.
	@all_summaries = sort { $a->[1] <=> $b->[1] } @all_summaries;

	return $var_names_ref, \@all_summaries;
}

sub run_consolidate($$$) {
	show_output "\n*** Function: run_consolidate() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $orig_config_file, $message ) = @_;

	my ( $var_names_ref, $all_summaries_ref ) = load_all_summaries($values_ref);

	my $consolidated_directory = "$values_ref->{output_dir}$values_ref->{job_name}/Results/";
	my $consolidated_filename  = "${consolidated_directory}sorted_params.txt";

	verbose_mkdir($consolidated_directory)
	  unless ( -d $consolidated_directory );

	write_all_summary( $consolidated_filename, $var_names_ref, $all_summaries_ref );

	my @list = @$all_summaries_ref;

	my $summaries = shift @list;
	my $best_run  = @$summaries[0];
	my $best_chi  = @$summaries[1];

	my ($best_gen) = $best_run =~ /(\d+)/;
	my $best_perm = ( split( /perm/, $best_run ) )[1];

	my $model_root = $values_ref->{model};
	my ($suffix) = $model_root =~ /(\.[^.]*)$/;
	$model_root =~ s/(\.[^.]*)$//;
	$model_root = fileparse($model_root);

	my $best_gdat_filename_full  = "$values_ref->{output_dir}$values_ref->{job_name}/$best_gen/${model_root}_perm$best_perm";
	my $best_model_filename_full = "$values_ref->{output_dir}$values_ref->{job_name}/$best_gen/${model_root}_perm${best_perm}$suffix";

	foreach ( @{ $values_ref->{exp_file} } ) {
		my $exp_name = fileparse($_);
		$exp_name =~ s/.exp$//;

		my $best_gdat_filename = "${model_root}_${exp_name}_perm$best_perm.gdat";

		copy( "$values_ref->{output_dir}$values_ref->{job_name}/$best_gen/$best_gdat_filename", "$values_ref->{output_dir}$values_ref->{job_name}/Results/${exp_name}_bestfit.gdat" )
		  or show_output "Error: Couldn't copy output from best run into the Results directory\n";

		copy( "$values_ref->{output_dir}$values_ref->{job_name}/$best_gen/${model_root}_${exp_name}_gen${best_gen}_perm$best_perm.txt", "$values_ref->{output_dir}$values_ref->{job_name}/Results/${exp_name}_combined.txt" )
		  or show_output "Error: Couldn't copy output from best run into the Results directory\n";
	}

	copy( $orig_config_file, $consolidated_directory );

	copy( $best_model_filename_full, "$values_ref->{output_dir}$values_ref->{job_name}/Results/" )
	  if ( -f $best_model_filename_full );

	if ( $values_ref->{ode} ) {
		$best_model_filename_full =~ s/$suffix/\.net/;
		copy( $best_model_filename_full, "$values_ref->{output_dir}$values_ref->{job_name}/Results/" )
		  if ( -f $best_model_filename_full );
	}

	if ( $values_ref->{make_plots} ) {
		show_output "\nGenerating plots...\n"
		  if ( $verbosity >= 3 );

		graph_chi($values_ref);
		graph_outs($values_ref);
	}

	show_output $message;

	return ( $best_gdat_filename_full, $best_chi, $best_perm );
}

sub write_all_summary($$$) {
	show_output "\n*** Function: write_all_summary() ***\n"
	  if ( $verbosity >= 4 );

	my ( $consolidated_filename, $var_names_ref, $all_summaries_ref ) = @_;
	open CONS, '>', $consolidated_filename
	  or show_output "die", "Error: Unable to open $consolidated_filename for writing.\n";

	#printf CONS "%-16s", "Permutation";
	#printf CONS " %-12s", "Chi-Sq";

	printf CONS "%-1s",  "Permutation";
	printf CONS " %-1s", "Chi-Sq";

	foreach (@$var_names_ref) {
		$_ =~ s/\s+//g;

		#printf CONS " %-14s", $_;
		printf CONS " %-1s", $_;
	}
	printf CONS "\n";

	foreach (@$all_summaries_ref) {
		my @list   = @$_;
		my $perm   = shift @list;
		my $chi_sq = shift @list;

		#printf CONS "%-16s", $perm;
		printf CONS "%-1s", $perm;

		#printf CONS " %-12f", $chi_sq;
		printf CONS " %.16f", $chi_sq;

		foreach (@list) {
			printf CONS " %.8e", $_;
		}
		print CONS "\n";
	}
	close CONS;
}

sub read_summary($) {
	show_output "\n*** Function: read_summary() ***\n"
	  if ( $verbosity >= 4 );

	my ($summarydiff_filename) = @_;

	# Read in sorted summary of all permutations.
	unless ( open SUMMARY_DIFF, $summarydiff_filename ) {
		$summarydiff_filename =~ s/_1_((parent_)?summary)/_$1/
		  or show_output "die", "Error: Unable to open summary $summarydiff_filename.\n";

		open SUMMARY_DIFF, $summarydiff_filename
		  or show_output "die", "Error: Unable to open summary $summarydiff_filename.\n";
	}

	# Read in header.
	$_ = <SUMMARY_DIFF>;
	my (@summary_header) = split;
	my (@var_names)      = @summary_header[ 2 .. $#summary_header ];
	show_output "die", "Error: Unable to read variable names from $summarydiff_filename"
	  unless @var_names;
	my @var_data = ();

	# Read in data
	while (<SUMMARY_DIFF>) {
		my @line = split;
		push @var_data, \@line;
	}

	close SUMMARY_DIFF;

	return \@var_names, \@var_data;
}

sub run_genetic($$$$$) {
	show_output "\n*** Function: run_genetic() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $current_directory, $previous_directory, $new_config_file, $orig_config_file ) = @_;

	# Parses mutation rates and returns them to $mutate_ref
	my $mutate_ref;

	if ( $values_ref->{mutate} ) {
		show_output "Parsing mutation rates...\n"
		  if ( $verbosity >= 3 );
		$mutate_ref = store_mutate( $values_ref->{mutate} );
	}

	show_output "Parsing last results from generation...\n"
	  if ( $verbosity >= 3 );

	# Read in sorted summary of all permutations.
	# var_names_ref contains only parameter names in a single dimensional array
	#current_var_data_ref contains all subsequent rows, including permutation names, chi^2 values, and parameter values
	my ( $var_names_ref, $current_var_data_ref ) = read_summary("${current_directory}_summary_diff.txt");    # outputroot/2_summary_diff.txt

	# Read in sorted summary of previous generations' permutations.
	my ( $prev_var_names_ref, $prev_var_data_ref );
	if ($previous_directory) {
		( $prev_var_names_ref, $prev_var_data_ref ) = read_summary("${previous_directory}/parent_summary_diff.txt");
		# Note that parent_summary_diff can contains results from more than 1 previous generation
	}
	else {
		$prev_var_data_ref = [];
	}

	if ( ($previous_directory) && ( $values_ref->{stop_when_stalled} ) ) {

		# Check to see if parameters in current and previous generation are the same
		#print Dumper($current_var_data_ref, $prev_var_data_ref);
		my $difference = 0;

		for ( my $i = 1 ; $i < ( scalar @$current_var_data_ref ) ; $i++ ) {
			for ( my $j = 2 ; $j < ( scalar @{ $$current_var_data_ref[0] } ) ; $j++ ) {

				#show_output "$$current_var_data_ref[$i][$j], $$prev_var_data_ref[$i][$j]\n";
				if ( $$current_var_data_ref[$i][$j] ne $$prev_var_data_ref[$i][$j] ) {
					$difference = 1;
				}
			}
		}

		if ( !$difference ) {
			my $exit_message = "\n***Current parameters are the same as previous generation's parameters! Stalled at generation $values_ref->{generation}.***\n";
			my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
			if ( $values_ref->{bootstrap} ) {
				bootstrap_run( $values_ref, $chi, $best_run_filename );
			}
			exit 0;
		}
	}

	# Put previous and current generation data into a combined ordered list.
	# @var_data first dimension corresponds to a particular permutation
	# @var_data second dimension corresponds to parameter values in that permutation
	my @var_data = sort { $a->[1] <=> $b->[1] } @$current_var_data_ref, @$prev_var_data_ref;

	# Only keep the top (n) values where (n) represents number of permutations in a generation
	my $var_data_ref = [ @var_data[ 0 .. ( @$current_var_data_ref - 1 ) ] ];

	show_output "Outputting parent_summary_diff...\n"
	  if ( $verbosity >= 3 );

	# Save current parents for next generation's use.
	write_all_summary( "${current_directory}/parent_summary_diff.txt", $var_names_ref, $var_data_ref );

	my $reached_min_chi = 0;

	# Does our list contain a chi^2 value that is below our desired chi^2 value? If so, clean up and exit.
	if ( $var_data[1][1] <= $values_ref->{min_objfunc_value} ) {
		$reached_min_chi = 1;
	}

	if ($reached_min_chi) {
		my $exit_message = "Reached our minimum chi^2 of $values_ref->{min_objfunc_value}, all done.\n";
		my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
		if ( $values_ref->{bootstrap} ) {
			bootstrap_run( $values_ref, $chi, $best_run_filename );
		}
		exit 0;
	}

	my @new_var_data = ();

	if ( $values_ref->{skip_genetic} ) {

		# Skip genetic step, just simulate with best data sets.
		@new_var_data = @$var_data_ref;

		# Remove setting so that it is not passed on to next generation.
		delete( $values_ref->{skip_genetic} );
	}

	else {
		
		#print Dumper $var_data_ref;
		
		show_output "Breed and mutate...\n"
		  if ( $verbosity >= 3 );

		my %all_chosen_parents;

		# Perform genetic algorithm.
		my @survived_var_data;

		# Select all permutations with chi^2 less than the max threshold
		if ( $values_ref->{max_objfunc_value} ) {

			# Select survivors.
			@survived_var_data = grep { $_->[1] <= $values_ref->{max_objfunc_value} } @$var_data_ref;
		}
		else {
			# No one jumps off a cliff
			@survived_var_data = @$var_data_ref;
		}

		@survived_var_data = sort { ( $a ? $a->[1] : 0 ) <=> ( $b ? $b->[1] : 0 ) } @survived_var_data;

		my @parents_to_keep = @survived_var_data[ 0 .. ( $values_ref->{keep_parents} - 1 ) ]
		  if $values_ref->{keep_parents};

		$values_ref->{max_parents} = $values_ref->{first_gen_permutations}
		  if ( $values_ref->{first_gen_permutations} && $values_ref->{generation} == 1 );

		# Select only top (n) permutations where (n) is equal to max_parents
		if ( $values_ref->{max_parents} ) {

			# Only keep the best [max_parents] matches.
			splice @survived_var_data, ( $values_ref->{max_parents} )
			  if ( @survived_var_data > $values_ref->{max_parents} );
		}

		# Need at least 3 survivors because the lowest is always ignored.
		if ( @survived_var_data < 3 ) {
			my $exit_message = "After culling possible parents that were either over the chi^2 threshold or removed as a result of max_parents, we were left with only " . scalar(@survived_var_data) . " breeders. This is too few to continue.\n";
			my ( $best_run_filename, $chi ) = finished( $values_ref, $orig_config_file, $exit_message );
			if ( $values_ref->{bootstrap} ) {
				bootstrap_run( $values_ref, $chi, $best_run_filename );
			}
			exit 0;
		}

		# Create the same number of permutations that we had before.
		my @weights = map { $_->[1] } @survived_var_data;
		
		#print "weights:\n";
		#print Dumper @weights;

		my $max_weight = max(@weights);

		#show_output "max weight: $max_weight\n";

		# Subtract chi^2's from value of highest chi^2, resulting in high values for low chi^2s.
		@weights = map { $max_weight - $_ } @weights;
		
		#print "new weights:\n";
		#print Dumper @weights;

		my $weight_sum = sum(@weights);
		#print "weight sum: " . $weight_sum . "\n";

		# Each pair of parents creates two offspring.
		my $parent_pairs = @$var_data_ref / 2;

		#TODO
		# If there are 10 permutations in a generation and none were removed in the culling, this should return 5.
		# Is there a reason we're dividing generation size by 2 when we could be dividing survivors by 2 instead?
		for ( my $i = 0 ; $i < $parent_pairs ; $i++ ) {
			
			# pick_weighted returns an array element from the list. The element
			# corresponds to an index in $survived_var_data
			my $p1 = pick_weighted( $weight_sum, \@weights, $values_ref );
			my $p2 = pick_weighted( $weight_sum, \@weights, $values_ref );
			
			# If we want to force different parents, make sure p1 and p2 are different
			if ( $values_ref->{force_different_parents} ) {
				my $retry_count = 0;
				while ( $p1 == $p2 ) {
					$p2 = pick_weighted( $weight_sum, \@weights, $values_ref );
					$retry_count++;
					if ( $retry_count > $max_retry_different_parents ) {
						show_output "Tried too many times to select different parents for breeding, giving up. We will just select the top two.\n"
						  if ( $verbosity >= 2 );

						$p1 = 0;
						$p2 = 1;
						last;
					}
				}
			}
			
			#print "p1: " . $survived_var_data[$p1][1] . "\n";
			#print "p2: " . $survived_var_data[$p2][1] . "\n";

			#print Dumper %all_chosen_parents;

			# Breed!
			
			# Create the children arrays
			my @child1 = ( "dummy_perm", 0 );
			my @child2 = ( "dummy_perm", 0 );

			# Loop as (n) times where (n) is the number of parameters
			for ( my $j = 0 ; $j < @$var_names_ref ; $j++ ) {
				my $new_value1;
				my $new_value2;

				if ( rand(100) < ( $values_ref->{swap_rate} * 100 ) ) {

					# Do not swap.
					$new_value1 = $survived_var_data[$p1][ $j + 2 ];
					$new_value2 = $survived_var_data[$p2][ $j + 2 ];
				}
				else {
					# Swap.
					$new_value2 = $survived_var_data[$p1][ $j + 2 ];
					$new_value1 = $survived_var_data[$p2][ $j + 2 ];
				}

				# Mutate.
				my $mut;

				if ( exists $mutate_ref->{ $var_names_ref->[$j] } ) {
					$mut = $mutate_ref->{ $var_names_ref->[$j] };
				}
				elsif ( exists $mutate_ref->{default} ) {
					$mut = $mutate_ref->{default};
				}

				# If the random number is less than our mutation probability
				if ( ($mut) && ( rand() < $mut->{prob} ) ) {
					my $max_change = $new_value1 * $mut->{percent_change};

					# Allow for +/-.
					my $change = rand( $max_change * 2 ) - $max_change;
					$new_value1 += $change;
				}
				if ( ($mut) && ( rand() < $mut->{prob} ) ) {
					my $max_change = $new_value2 * $mut->{percent_change};

					# Allow for +/-.
					my $change = rand( $max_change * 2 ) - $max_change;
					$new_value2 += $change;
				}

				push @child1, $new_value1;
				push @child2, $new_value2;
			}

			push @new_var_data, \@child1;
			push @new_var_data, \@child2
			  unless ( $parent_pairs * 2 ) % 2 == 1 && $parent_pairs - $i == 0.5;
		}

		if (@parents_to_keep) {

			# Make our saved parents looks like the rest of our data
			foreach (@parents_to_keep) {
				$_->[0] = "dummy_perm";
				$_->[1] = 0;
			}

			# Put unchanged parent data back into our data set
			unshift @new_var_data, @parents_to_keep;

			# And compensate for the extra data points by removing an equal number of the last (worst) fits
			for ( my $i = 0 ; $i < scalar(@parents_to_keep) ; $i++ ) {
				my $junk = pop @new_var_data;
			}
		}

		#my @parents_list;
		#for (keys %all_chosen_parents) {
		#	push @parents_list, "${_}.$all_chosen_parents{$_}";
		#}
		#@parents_list = sort { $a <=> $b } @parents_list;
		#my $message;
		#for (@parents_list) {
		#	$message .= "$_\n";
		#}
		#show_output $message;
	}
	
	show_output "Creating .conf file for next generation...\n"
	  if ( $verbosity >= 3 );

	create_new_config( $values_ref, \@new_var_data, $var_names_ref, \@var_variables, \@list_variables, $new_config_file );
	my $command = "perl $0 $generation_command $new_config_file";

	# If we're chaining qsub or not using qsub at all we need to get our next generation going
	#if (!$values_ref->{use_cluster}) {
	if ( ( $values_ref->{use_cluster} && $values_ref->{run_job_from_worknode} ) || !$values_ref->{use_cluster} ) {
		show_output "Starting new generation...\n"
		  if ( $verbosity >= 2 );
		exec($command);
	}

	# But if we're not chaining, we let this run drop off and the parent (on submission node) will keep us going
	else {
		exit 0;
	}
}

sub create_new_config($$$$$$) {
	show_output "\n*** Function: create_new_config() ***\n"
	  if ( $verbosity >= 4 );

	my ( $values_ref, $new_var_data_ref, $var_names_ref, $var_variables_ref, $list_variables_ref, $new_config_file ) = @_;

	$new_config_file =~ /($values_ref->{job_name}\/\d+\/)/;

	my $dir_tree = "$values_ref->{output_dir}/$1";
	verbose_mkdir($dir_tree);

	open NEW_CONFIG, '>', $new_config_file
	  or show_output "die", "Error: Unable to open $new_config_file for writing.\n";

	print NEW_CONFIG "# Automatically generated, job:$values_ref->{job_name}, ", "generation:$values_ref->{generation}\n";

	# Skip perm name and chi-squared value.
	for ( my $i = 0 ; $i < @$var_names_ref ; $i++ ) {
		print NEW_CONFIG "static_list_var $var_names_ref->[$i]";
		foreach (@$new_var_data_ref) {
			print NEW_CONFIG "\t$_->[$i + 2]";
		}
		print NEW_CONFIG "\n";
	}

	foreach my $name ( sort keys %$values_ref ) {
		# Don't output variable definitions because we are redefining them.
		next if ( grep( /^$name$/, @$var_variables_ref ) );

		if ( grep( /^$name$/, @$list_variables_ref ) ) {
			for ( my $i = 0 ; $i < @{ $values_ref->{$name} } ; $i++ ) {
				print NEW_CONFIG "$name $values_ref->{$name}[$i]\n";
			}
		}
		else {
			print NEW_CONFIG "$name $values_ref->{$name}\n";
		}
	}
	close NEW_CONFIG;
}

sub store_mutate($) {
	my ($list_ref) = @_;

	my %mutate;
	foreach (@$list_ref) {
		my ( $column_name, $prob, $percent_change ) = split;
		show_output "die", "Error: Syntax error for mutate $_\n"
		  unless ( defined $percent_change );
		$mutate{$column_name} = { prob => $prob, percent_change => $percent_change };
	}

	return \%mutate;
}

sub pick_weighted($$$) {
	my ( $weight_sum, $weights_ref, $values_ref ) = @_;

	# Pick random number between 0 and $weight_sum. Multiply it by extra_weight factor.
	my $chosen = ( rand($weight_sum) * ( 1 - ( $values_ref->{extra_weight} ) / 10 ) );
	#print "chosen: " . $chosen . "\n";
	my $current_sum = 0;
	for ( my $i = 0 ; $i < @$weights_ref ; $i++ ) {
		$current_sum += $weights_ref->[$i];
		if ( $current_sum >= $chosen ) {
			return $i;
		}
	}

	# Fell off the end, simply return the last value.
	return $#{$weights_ref};
}

sub get_model_params($$) {
	show_output "\n*** Function: get_model_params() ***\n"
	  if ( $verbosity >= 4 );

	my ( $file, $values_ref ) = @_;

	$file =~ s/\.bngl$/\.net/ if $values_ref->{ode};
	$file =~ s/\.txt$/\.net/  if $values_ref->{ode};

	open FD, $file
	  or show_output "die", "Error: Unable to open model file $file to read parameters.\n";

	my @names  = ();
	my @values = ();
	while (<FD>) {
		my ( $name, $value ) = /^# (\w+) changed to ([\d\-+e.]+$)/;
		last unless ($name);

		push @names,  $name;
		push @values, $value;
	}

	return \@names, \@values;
}

sub create_columns($$$$$$) {
	show_output "\n*** Function: create_columns() ***\n"
	  if ( $verbosity >= 4 );

	my ( $model_file, $root_filename, $model_roots, $data_ref, $column_ref, $newcolumns ) = @_;

	foreach my $col ( keys %$newcolumns ) {
		my $factor = $newcolumns->{$col};
		if ( defined $factor ) {

			# Replace any column references with variable references.
			foreach my $model_root (@$model_roots) {
				my $status = $factor =~ s/$model_root:([\w]+)/\$data_ref->{$model_root}[\$i]{$1}/g;
				if ($status) {
					if ( !exists $data_ref->{$model_root}[0]{$1} ) {
						show_output "die", "Error in newcolumn $newcolumns->{$col}:\n Column $1 not present in model $model_root\n";
					}
					if ( @{ $data_ref->{$model_root} } != @{ $data_ref->{$root_filename} } ) {
						show_output "die", "Error in newcolumn $col $root_filename $newcolumns->{$col}:\n $model_root and $root_filename have different numbers of timesteps.\n";
					}
				}
			}
			foreach my $col2 ( keys %{ $data_ref->{$root_filename}[0] } ) {

				#show_output "root=$root_filename, col2=$col2 ";
				#show_output "factor=$factor\n";
				#show_output "first: $col2(?!})";
				#show_output "data:\n";
				#print Dumper $data_ref;
				$factor =~ s/$col2(?!})/\$data_ref->{$root_filename}[\$i]{$col2}/g;
			}

			# Apply factor to create new column.
			my $expression = "\$data_ref->{$root_filename}[\$i]{$col} = " . $factor;

			my $string = "Expression for new column ${model_file}[\$i]{$col}:\n" . "$expression\n";

			#show_output $string
			#	if ($verbosity >= 3);

			my $current_line;
			my $safe = new Safe;
			$safe->permit_only(qw(:default));
			$safe->share('$current_line');

			for ( my $i = 0 ; $i < @{ $data_ref->{$root_filename} } ; $i++ ) {

				#$safe->reval($expression);
				eval($expression);
				if ($@) {
					if ( $@ =~ /("[^"]*")/ ) {
						show_output "die", "Error: Unexpected string $1 in newcolumn expression:\n $newcolumns->{$col}\n";
					}
					elsif ( $@ =~ /Illegal division by zero/ ) {

						# Ignore divide by zero.
						$data_ref->{$root_filename}[$i]{$col} = 0;
						undef($@);
					}
					else {
						show_output "die", "Error applying '$factor' to ${model_file}[$col]:\n $expression\n $@\n";
					}
				}
			}

			push @{ $column_ref->{$root_filename} }, $col;
		}
	}
}

sub load_exp($$) {
	show_output "\n*** Function: load_exp() ***\n"
	  if ( $verbosity >= 4 );

	my ( $exp_files, $values_ref ) = @_;

	my %expdata_hash;
	my %expcolumns_hash;

	my ($model_root) = $values_ref->{model} =~ /(.*)\.txt/;
	($model_root) = $values_ref->{model} =~ /(.*)\.bngl/ unless ($model_root);
	$model_root = fileparse($model_root);

	my $total_sd_columns;
	foreach (@$exp_files) {
		my ($exp_root) = /(.*)\.exp/;
		my ($prefix,$path,$suffix) = fileparse($exp_root);
		my $num_sd_columns;

		# Load experimental data.
		my $exp_file = $_;
		my ( $expdata_ref, @expcolumns ) = load_data("${exp_file}");

		my $control_col = "time";
		if ( $values_ref->{scan_parameter} ) {
			foreach ( @{ $values_ref->{scan_parameter} } ) {
				if ( ( split( ' ', $_ ) )[0] eq $prefix ) {
					$control_col = ( split( ' ', $_ ) )[1];
				}
			}
		}

		if ( !defined $$expdata_ref[0]{$control_col} ) {
			my $message;

			if ( $control_col ne "time" ) {
				$message = "Error: Your .exp file $_ doesn't contain the column header specified by suffix in your .bngl file. Please fix this and re-run\n";
			}
			else {
				$message = "Error: Your .exp file $_ doesn't contain a 'time' column.\n";
			}

			die $message;
		}

		if ( $values_ref->{objfunc} == 2 ) {
			for (@expcolumns) {
				$num_sd_columns += 1 if (/_SD/gi);
			}

			if ( ( ( scalar(@expcolumns) - 1 ) / 2 ) != $num_sd_columns ) {
				show_output "die", "\nIt looks like you are using objfunc=4 which requires a column name ending in '_SD' for every column containing experimental data, but $_ doesn't have the correct number of _SD columns. Please see the documentation for instructions on using objfunc=4\n";
			}
		}

		$expdata_hash{"${path}${model_root}_$prefix"}    = $expdata_ref;
		$expcolumns_hash{"${path}${model_root}_$prefix"} = \@expcolumns;

	}
	return ( \%expdata_hash, \%expcolumns_hash );
}

sub store_factors($) {
	my ($list_ref) = @_;

	my @list = @$list_ref;

	my %factors;
	foreach (@list) {
		my ( $column_name, $file_name, $factor ) = /^\s*(\S+)\s+(\S+)\s+(.*)/
		  or show_output "die", "Syntax error for factor $_\n";
		$column_name =~ /^[a-zA-Z\d_]/
		  or show_output "die", "Error with factor $_:\n Column name must begin with a letter, number or underscore.\n";
		$column_name =~ /^[\w]+$/
		  or show_output "die", "Error with factor $_:\n Column name must only have letters, numbers or underscores.\n";
		$factor =~ /^[+\-*\/%]/
		  or show_output "die", "Error with factor $_:\n Factor must begin with '+', '-', '*', '/', '**' or '%'\n";
		$factors{$column_name} = () unless ( $factors{$column_name} );
		$factors{$column_name}{$file_name} = $factor;
	}

	return \%factors;
}

sub store_newcolumns($) {
	show_output "\n*** Function: store_newcolumns() ***\n"
	  if ( $verbosity >= 4 );

	my ($list_ref) = @_;

	my @list = @$list_ref;

	my %newcolumns;
	foreach (@list) {
		my ( $column_name, $file_name, $factor ) = /^\s*(\S+)\s+(\S+)\s+(.*)/
		  or show_output "die", "Syntax error for newcolumn $_\n";
		$column_name =~ /^[a-zA-Z\d_]/
		  or show_output "die", "Error with newcolumn $_:\n Column name must begin with a letter, number or underscore.\n";
		$column_name =~ /^[\w]+$/
		  or show_output "die", "Error with newcolumn $_:\n Column name must only have letters, numbers or underscores.\n";
		$factor =~ /^[^+\-*\/%]/
		  or show_output "die", "Error with newcolumn $_:\n Factor must not begin with '+', '-', '*', '/', '**' or '%'\n";
		$newcolumns{$file_name} = () unless ( $newcolumns{$file_name} );
		$newcolumns{$file_name}{$column_name} = $factor;
	}

	return \%newcolumns;
}

sub load_data($) {
	my ($model_data) = @_;

	open FH, $model_data
	  or show_output "die", "Error: Unable to open $model_data";

	# Undefines the line delimiter so whole file can be read in at once
	local $/ = undef;

	my $wholeFile = <FH>;
	my @fileLines = split /\r\n|\r|\n/, $wholeFile;

	my ( $hash, @column_names ) = split /\s+/, $fileLines[0];

	show_output "die", "Error: Expected hash as first column name in $model_data.\n"
	  unless $hash eq "#";

	show_output "die", "Error: No columns found in $model_data"
	  unless scalar(@column_names);

	my %column_names = ();

	for ( my $i = 0 ; $i < @column_names ; $i++ ) {
		if ( $column_names[$i] =~ /^\d/ ) {
			show_output "die", "Error: Column name '$column_names[$i]' cannot begin with a number.\n";
		}
		elsif ( $column_names[$i] =~ /\(\)/ ) {
			$column_names[$i] =~ s/\(\)//;
		}

		#unless ($column_names[$i] =~ /^[\w]+$/) {
		#show_output "die", "Column name '$column_names[$i]' can only have letters, numbers ",
		#"and underscores.\n";
		#}
		$column_names{ $column_names[$i] } = $i;
	}

	my $lineno = 0;
	my $safe   = new Safe;
	$safe->permit_only(qw(:base_core));
	my @data;
	foreach (@fileLines) {
		my %current_line;
		$safe->share('%current_line');
		$lineno++;

		# Remove line endings no matter what the source OS was.
		tr/\n\r//d;

		# Ignore blank lines.
		next unless $_;

		my @current_line = split;

		#show_output "die", "Number of columns on line $lineno of $model_data do not match header\n"
		#	unless (@column_names == @current_line) || ($lineno == 1);

		if ( $lineno != 1 ) {
			%current_line = ();
			for ( my $i = 0 ; $i < @current_line ; $i++ ) {

				#my $data = $current_line[$i];
				#$data = sprintf("%.12f", $data);
				$current_line{ $column_names[$i] } = $current_line[$i];
			}
			push @data, \%current_line;
		}
	}

	close FH;

	return \@data, @column_names;
}

sub print_usage() {
	print "BioNetFit v${version}\nUsage:\n\n$0 config_file.conf\n$0 results original_config_file.conf\n$0 resume original_config_file.conf\n$0 resume max_gens original_config_file.conf\n";
	exit 0;
}

sub parse_file(@) {
	my (@config_filenames) = @_;

	my %values = ();

	foreach my $filename (@config_filenames) {
		open CONFFILE, $filename
		  or show_output "die", "Error: Unable to open config file $filename for parsing\n";

		my $lineno = 0;

		while (<CONFFILE>) {
			chomp;
			$lineno++;

			# Skip any line beginning with '#'.
			next if (/^\s*#/);

			# Skip any blank lines.
			next if (/^\s*$/);

			# Match variable assignment.
			/^\s*([a-zA-Z_]+)\s*=\s*(.*)$/;
			my $name  = $1;
			my $value = $2;

			unless ( defined $value ) {
				/^\s*([a-zA-Z_]+)\s+(.*)$/;
				$name  = $1;
				$value = $2;
			}

			show_output "die", "Error: Syntax error on $filename:$lineno: $_\n" unless defined $value;

			# Remove leading and trailing white space.
			$value =~ s/^\s+//;
			$value =~ s/\s+$//;

			# Remove first level quotes if present.
			$value =~ s/"(.*)"/$1/;

			# Save in hash.
			if ( grep /^$name$/, @list_variables ) {

				# Expect multiple values with the same name, so create list.
				@{ $values{$name} } = () unless $values{$name};

				# Add to end of list.
				unshift @{ $values{$name} }, $value;
			}
			else {
				# Only one value expected, so overwrite any existing values.
				$values{$name} = $value;
			}
		}
		close CONFFILE;
		
		#if ($values{"static_list_var"}) {
		#	@{ $values{"static_list_var"} } = empty_sort (@{ $values{"static_list_var"} });
		#}
		#print "values_ref after sorting:\n";
		#print Dumper %values;
	}
	
	return \%values;
}

sub update_values($) {
	my ($values_ref) = @_;

	my $status = 1;

	# Verify required variables are present.
	foreach (@required_variables) {
		unless ( $values_ref->{$_} ) {
			printf "$_ must be defined in configuration file.\n";
			$status = 0;
		}
	}

	$values_ref->{max_parents} = $values_ref->{permutations} unless $values_ref->{max_parents};

	# Set any missing variables to default values.
	foreach ( keys %var_defaults ) {
		$values_ref->{$_} = $var_defaults{$_} unless defined( $values_ref->{$_} );
	}

	# Make filenames absolute.
	foreach my $option_name (@filenames) {
		next unless ( defined $values_ref->{$option_name} );

		my @toparse;
		my $isarray = 0;
		if ( ref( $values_ref->{$option_name} ) eq "ARRAY" ) {
			push @toparse, @{ $values_ref->{$option_name} };
			@{ $values_ref->{$option_name} } = ();
			$isarray = 1;
		}
		else {
			push @toparse, $values_ref->{$option_name};
		}

		foreach (@toparse) {

			my $file     = $_;
			my $abs_file = abs_path $file;

			if ( !$abs_file ) {
				print qq|The file: "$file" specified in your .conf could not be found. Please be sure the file exists and the path is correct.\n|;
				$status = 0;
			}
			elsif ( !-e $abs_file ) {
				print qq|The file: "$file" specified in your .conf could not be found. Please be sure the file exists and the path is correct.\n|;
				$status = 0;
			}
			if ($isarray) {
				push @{ $values_ref->{$option_name} }, $abs_file;
			}
			else {
				$values_ref->{$option_name} = $abs_file;
			}
		}
	}

	# Make pathnames absolute, creating if necessary.
	foreach (@pathnames) {
		next unless ( defined $values_ref->{$_} );

		my $path = $values_ref->{$_};

		verbose_mkdir("$path") if $values_ref->{initialize};

		my $abs_path = abs_path $path;

		if ( !$abs_path || !-e $abs_path) {
			if ($values_ref->{ask_create}) {
				my $answer;
				do
				{
					print "The directory $path specified in your .conf file does not exist. Would you like to create it now? ";
					$answer = <STDIN>;
					chomp $answer;
				} until (($answer eq 'y') || ($answer eq 'n') || ($answer eq 'Y') || ($answer eq 'N'));
							
				die "Error: Please check to be sure your directory exists and the path in your .conf file is correct.\n" 
					if (($answer eq 'N') || ($answer eq 'n'));
			}
			verbose_mkdir($values_ref->{output_dir});
		}

		$values_ref->{$_} = $abs_path;
	}

	unless ($status) {
		show_output "Correct any errors listed above and try again.\n";
		exit 1;
	}
}

sub create_pbs_command($$$$) {
	show_output "\n*** Function: create_pbs_command() ***\n"
	  if ( $verbosity >= 4 );

	my ( $pbs_script, $job_name, $values_ref, $extra_values ) = @_;

	my @flags = qw(R W qw r PD);
	$values_ref->{counter} += 1;
	if ( check_queue( $values_ref, \@flags ) >= $values_ref->{job_limit} || $values_ref->{counter} >= $values_ref->{job_limit} ) {
		show_output "die", "Error: Trying to create a PBS command but the nubmer of queued or running jobs is greater than the maximum allowed by 'job_limit'. Please use either a) alter 'cluster_parallel', 'multisim', or 'permutations' to ensure that BioNetFit doesn't send so many jobs at once; or b) Change 'job_limit' to a higher number.\n";
	}

	my $num_commands = grep( /^command/, keys(%$extra_values) );

	$num_commands = 1
	  unless $num_commands;

	my $pbs_mail_options;

	my $prefix = ( substr( $values_ref->{job_name}, 0, 4 ) );
	$job_name = $prefix . "-" . $job_name;
	$job_name .= "_B$values_ref->{bootstrap_num}"
	  if $values_ref->{bootstrap};

	if ( $values_ref->{email_after_generation} || $values_ref->{email_after_run} || $values_ref->{email_on_abort} ) {
		my $final_gen = "${prefix}-A$values_ref->{max_generations}";
		$final_gen = "${prefix}-G$values_ref->{max_generations}"
		  if ( $job_name =~ /${prefix}-G$values_ref->{max_generations}/ );

		if ( $values_ref->{email_after_run} && $job_name =~ /$final_gen/ ) {
			if ( $values_ref->{cluster_software} eq "slurm" ) {
				$pbs_mail_options = "END";
			}
			else {
				$pbs_mail_options .= "e";
			}
		}
		elsif ( $values_ref->{email_after_generation} && ( $job_name =~ /${prefix}-A\d+.+/ || $job_name =~ /${prefix}-G\d+/ ) ) {
			if ( $values_ref->{cluster_software} eq "slurm" ) {
				$pbs_mail_options = "END";
			}
			else {
				$pbs_mail_options .= "e";
			}
		}

		if ( $values_ref->{email_on_abort} ) {
			if ( $values_ref->{cluster_software} eq "slurm" ) {
				$pbs_mail_options .= "," if $pbs_mail_options;
				$pbs_mail_options .= "FAIL";
			}
			else {
				$pbs_mail_options .= "a";
			}
		}
	}

	my $email = $values_ref->{email};

	my @pbs_l_options;
	my @pbs_other_options;

	# PBS -l options for torque
	if ( $values_ref->{cluster_software} eq "torque" ) {

		# If we're creating a command to run simulations let's specify our CPU requirements
		if ( $job_name =~ /${prefix}-S\d+/ ) {
			push @pbs_l_options, "nodes=1:ppn=$values_ref->{cluster_parallel}";
			push @pbs_l_options, "walltime=$values_ref->{max_walltime}"
			  if ( $values_ref->{max_walltime} );
		}
		else {
			my @limitsplit = split( /:/, $values_ref->{walltime_limit} );
			my @walltimesplit = split( /:/, $values_ref->{max_walltime} );
			
			my $max = ( $limitsplit[0] * 3600 ) + ( $limitsplit[1] * 60 ) + ( $limitsplit[2] );
			my $walltime = ( $walltimesplit[0] * 3600 ) + ( $walltimesplit[1] * 60 ) + ( $walltimesplit[2] );
				
			my $seconds = ( ( $values_ref->{permutations} * $values_ref->{smoothing} ) / ( $values_ref->{cluster_parallel} * $values_ref->{multisim} ) ) * $walltime / 2 ;
					
			if ($seconds > $max) {
				 push @pbs_l_options, "walltime=" . $values_ref->{walltime_limit};
			}
			else {
				# Analyze and Generation should have a walltime that is equal to that of combined walltimes of all individual chunks.
				push @pbs_l_options, "walltime=00:" . int($seconds/60) . ":00";
			}
		}
	}
	elsif ( $values_ref->{cluster_software} eq "ge" ) {
		if ( $job_name =~ /${prefix}-S\d+/ ) {
			push @pbs_l_options, "h_rt=$values_ref->{max_walltime}";
		}
	}

	my @vars = ();
	foreach ( keys %$values_ref ) {

		# Skip values that should be skipped
		next if ( grep /^$_$/, @skip_variables, @list_variables );

		show_output "die", "Error: Problem with $_\n" unless defined $values_ref->{$_};

		push @vars, qq|$_="$values_ref->{$_}"|;
	}

	foreach ( keys %$extra_values ) {
		show_output "die", "Error: Problem with $_\n"
		  unless defined $extra_values->{$_};

		push @vars, qq|$_="$extra_values->{$_}"|;
	}

	my $vars = "";
	if (@vars) {
		if ( $values_ref->{cluster_software} eq "slurm" ) {
			$vars = "--export " . join ",", @vars;
		}
		else {
			$vars = "-v " . join ",", @vars;
		}
	}

	my $command = "$values_ref->{cluster_command} ";

	if ($pbs_mail_options) {
		if ( $values_ref->{cluster_software} eq "slurm" ) {
			$command .= "--mail-user $email --mail-type $pbs_mail_options ";
		}
		else {
			$command .= "-M $email -m $pbs_mail_options ";
		}
	}

	if ( $values_ref->{cluster_software} eq "slurm" ) {
		$command .= "-J $job_name $vars --share ";

		if ( $job_name =~ /${prefix}-S\d+/ ) {
			$command .= "--mincpus $values_ref->{cluster_parallel} ";
		}
	}
	else {
		$command .= "-N $job_name $vars -S /bin/bash ";
	}

	if ( $values_ref->{save_cluster_output} ) {
		$command .= "-o $values_ref->{output_dir}/$values_ref->{job_name}_cluster_output/$job_name ";
	}
	else {
		$command .= "-o /dev/null ";
	}

	if ( $values_ref->{cluster_software} eq "torque" ) {
		$command .= "-d $values_ref->{output_dir}$values_ref->{job_name}/$values_ref->{generation} -j oe ";
	}
	elsif ( $values_ref->{cluster_software} eq "ge" ) {
		$command .= "-wd $values_ref->{output_dir}$values_ref->{job_name}/$values_ref->{generation} -j y ";
		$command .= "-pe $values_ref->{pe_name} $values_ref->{cluster_parallel} ";
	}

	$command .= "-A $values_ref->{account_name} "
	  if ( $values_ref->{account_name} );

	if ( $values_ref->{queue_name} ) {
		if ( $values_ref->{cluster_software} eq "slurm" ) {
			$command .= " -p $values_ref->{queue_name} ";
		}
		else {
			$command .= " -q $values_ref->{queue_name} ";
		}
	}

	if ( $values_ref->{max_walltime} && $values_ref->{cluster_software} eq "slurm" ) {
		if ( $job_name =~ /${prefix}-S\d+/ ) {
			$command .= "-t $values_ref->{max_walltime} ";
		}
		else {
			my @limitsplit = split( /:/, $values_ref->{walltime_limit} );
			my @walltimesplit = split( /:/, $values_ref->{max_walltime} );
		
			my $max = ( $limitsplit[0] * 3600 ) + ( $limitsplit[1] * 60 ) + ( $limitsplit[2] );
			my $walltime = ( $walltimesplit[0] * 3600 ) + ( $walltimesplit[1] * 60 ) + ( $walltimesplit[2] );
				
			my $seconds = ( ( $values_ref->{permutations} * $values_ref->{smoothing} ) / ( $values_ref->{cluster_parallel} * $values_ref->{multisim} ) ) * $walltime / 2 ;
					
			if ($seconds > $max) {
				$command .= "-t $values_ref->{walltime_limit} ";
			}
			else {
				# Analyze and Generation should have a walltime that is equal to that of combined walltimes of all individual chunks.
				$command .= "-t 00:" . int( $seconds / 60 ) . ":00 ";
			}
		}
	}

	if (@pbs_l_options) {
		$command .= "-l ";
		$command .= shift @pbs_l_options;
		foreach (@pbs_l_options) {
			$command .= ",$_";
		}
		$command .= " ";
	}

	$command .= qq|$pbs_script|;

	show_output "Ran create_pbs_command and came up with this command:\n$command\n"
	  if ( $verbosity >= 3 );

	return $command;
}

sub verbose_mkdir($) {
	my ($dir) = @_;

	my $status = 1;

	unless ( -e $dir ) {
		my $err;

		mkpath($dir);
		if ( $err && @$err ) {
			$status = 0;
			for my $diag (@$err) {
				my ( $file, $message ) = %$diag;
				if ( $file eq '' ) {
					show_output "general error: $message\n";
				}
				else {
					show_output "problem creating $file: $message\n";
				}
			}
		}
	}

	return $status;
}

sub calculate_var_sets($$) {
	show_output "\n*** Function: calculate_var_sets() ***\n"
	  if ( $verbosity >= 4 );

	my ( $names_ref, $values_ref ) = @_;

	my @var_sets = ();
	
	# Add substitution variables.
	foreach ( empty_sort @{ $values_ref->{substitute_var} } ) {
		my ( $name, $new_name, $junk ) = split;
		show_output "die", "Error parsing substitute_var '$_'\n" unless $new_name;
		show_output "die", "Error: Too many parameters for substitute_var '$_'\n" if $junk;

		push @$names_ref, $name;

		my @new_var_sets = ();
		if ( 0 == @var_sets ) {
			push @new_var_sets, [$new_name];
		}
		else {
			foreach my $var_set (@var_sets) {
				push @new_var_sets, [ @$var_set, $new_name ];
			}
		}

		@var_sets = @new_var_sets;
	}

	# Add static list variables.
	foreach ( empty_sort @{ $values_ref->{static_list_var} } ) {
		my ( $name, @current_values ) = split;
		show_output "die", "Error parsing list_var '$_'\n" unless @current_values;

		push @$names_ref, $name;

		# Do not permute.
		my @new_var_sets = ();

		if ( 0 == @var_sets ) {
			foreach my $current_value (@current_values) {
				push @var_sets, [$current_value];
			}
		}
		else {
			show_output "die", "Error: Number of elements in static list $_ doesn't match.\n"
			  unless ( @var_sets == @current_values );

			for ( my $i = 0 ; $i < @var_sets ; $i++ ) {
				push @{ $var_sets[$i] }, $current_values[$i];
			}
		}
	}

	# Add list variables.
	foreach ( empty_sort @{ $values_ref->{list_var} } ) {
		my ( $name, @current_values ) = split;
		show_output "die", "Error parsing list_var '$_'\n" unless @current_values;

		push @$names_ref, $name;

		# Permute current_values with existing var_sets.
		my @new_var_sets = ();
		foreach my $current_value (@current_values) {
			if ( 0 == @var_sets ) {
				#show_output "cv: $current_value\n";
				push @new_var_sets, [$current_value];
			}
			else {
				foreach my $var_set (@var_sets) {
					push @new_var_sets, [ @$var_set, $current_value ];
				}
			}
		}

		@var_sets = @new_var_sets;
	}

	# Add linear variables.
	foreach ( empty_sort @{ $values_ref->{linear_var} } ) {
		my ( $name, $first, $last, $steps, $junk ) = split;
		show_output "die", "Error parsing linear_var '$_'\n" unless $steps;
		show_output "die", "Error: Too many parameters for linear_var '$_'\n" if $junk;

		my @current_values = ();

		push @$names_ref, $name;

		my $step_size = 0;
		$step_size = ( $last - $first ) / ( $steps - 1 ) if ( $steps > 1 );

		# Generate all values except last.
		for ( my $i = 0, my $current = $first ; $i < ( $steps - 1 ) ; $i++, $current += $step_size ) {
			push @current_values, $current;
		}

		# Make sure last value is exactly $last regardless of roundoff error.
		if ( $steps > 1 ) {
			push @current_values, $last;
		}
		else {
			push @current_values, $first;
		}

		# Permute current_values with existing var_sets.
		my @new_var_sets = ();

		foreach my $current_value (@current_values) {
			if ( 0 == @var_sets ) {
				push @new_var_sets, [$current_value];
			}
			else {
				foreach my $var_set (@var_sets) {
					push @new_var_sets, [ @$var_set, $current_value ];
				}
			}
		}

		@var_sets = @new_var_sets;
	}

	# Add log variables.
	foreach ( empty_sort @{ $values_ref->{log_var} } ) {
		my ( $name, $first, $last, $steps, $mantissa, $base, $junk ) = split;
		show_output "die", "Error parsing log_var '$_'\n" unless $base;
		show_output "die", "Too many parameters for log_var '$_'\n" if $junk;

		my @current_values = ();

		push @$names_ref, $name;

		my $step_size = 0;
		$step_size = ( $last - $first ) / ( $steps - 1 ) if ( $steps > 1 );

		# Generate all values except last.
		for ( my $i = 0, my $current = $first ; $i < ( $steps - 1 ) ; $i++, $current += $step_size ) {
			if ( 10 == $base ) {
				push @current_values, "${mantissa}e$current";
			}
			else {
				push @current_values, $mantissa * ( $base**$current );
			}
		}

		# Make sure last value is exactly $last regardless of roundoff error.
		if ( $steps > 1 ) {
			if ( 10 == $base ) {
				push @current_values, "${mantissa}e$last";
			}
			else {
				push @current_values, $mantissa * ( $base**$last );
			}
		}
		else {
			if ( 10 == $base ) {
				push @current_values, "${mantissa}e$first";
			}
			else {
				push @current_values, $mantissa * ( $base**$first );
			}
		}

		# Permute current_values with existing var_sets.
		my @new_var_sets = ();

		foreach my $current_value (@current_values) {
			if ( 0 == @var_sets ) {
				push @new_var_sets, [$current_value];
			}
			else {
				foreach my $var_set (@var_sets) {
					push @new_var_sets, [ @$var_set, $current_value ];
				}
			}
		}

		@var_sets = @new_var_sets;
	}

	# Add log normal random variables.
	foreach ( empty_sort @{ $values_ref->{lognormrandom_var} } ) {
		my ( $name, $mean, $stddev, $junk ) = split;
		show_output "die", "Error parsing lnormrandom_var '$_'\n"
		  unless $stddev;

		show_output "die", "Error: Too many parameters for lnormrandom_var '$_'\n"
		  if $junk;

		push @$names_ref, $name;

		if ( @var_sets > 1 ) {

			# Ignore $steps, just add log-normal random value to each set.
			for ( my $i = 0 ; $i < @var_sets ; $i++ ) {
				push @{ $var_sets[$i] }, rlnorm( $mean, $stddev );
			}
		}
		elsif ( 1 == @var_sets ) {

			# Only one var_set, so generate $steps permutations with random numbers.
			my @orig_var_set = @{ $var_sets[0] };
			@var_sets = ();

			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ @orig_var_set, rlnorm( $mean, $stddev ) ];
			}
		}
		else {
			# No steps yet, so generate $steps random numbers.
			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ rlnorm( $mean, $stddev ) ];
			}
		}
	}

	# Add log random variables.
	foreach ( empty_sort @{ $values_ref->{lograndom_var} } ) {
		my ( $name, $first, $last, $mantissa, $base, $junk ) = split;
		show_output "die", "Error parsing lograndom_var '$_'\n" unless $base;
		show_output "die", "Error: Too many parameters for lograndom_var '$_'\n" if $junk;

		my $range = $last - $first;

		push @$names_ref, $name;

		if ( @var_sets > 1 ) {

			# Ignore $steps, just add random value to each set.
			for ( my $i = 0 ; $i < @var_sets ; $i++ ) {
				my $current = int( zero_rand($range) + 0.5 ) + $first;

				# randomize mantissa between 1 and base.
				$mantissa = rand( $base - 1 ) + 1;
				if ( 10 == $base ) {
					push @{ $var_sets[$i] }, "${mantissa}e$current";
				}
				else {
					push @{ $var_sets[$i] }, $mantissa * ( $base**$current );
				}
			}
		}

		elsif ( 1 == @var_sets ) {

			# Only one var_set, so generate $steps permutations with random numbers.
			my @orig_var_set = @{ $var_sets[0] };
			@var_sets = ();

			#for (my $i = 0; $i < $steps; $i++)
			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				my $current = int( zero_rand($range) + 0.5 ) + $first;

				# randomize mantissa between 1 and base
				$mantissa = rand( $base - 1 ) + 1;
				if ( 10 == $base ) {
					push @var_sets, [ @orig_var_set, "${mantissa}e$current" ];
				}
				else {
					push @var_sets, [ @orig_var_set, $mantissa * ( $base**$current ) ];
				}
			}
		}
		else {
			# No steps yet, so generate $steps random numbers.
			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				my $current = floor( zero_rand($range) + 0.5 ) + $first;

				# randomize mantissa between 1 and base
				$mantissa = rand( $base - 1 ) + 1;
				if ( 10 == $base ) {
					push @var_sets, ["${mantissa}e$current"];
				}
				else {
					push @var_sets, [ $mantissa * ( $base**$current ) ];
				}
			}
		}
	}

	# Add Bill's distribution
	foreach ( empty_sort @{ $values_ref->{loguniform_var} } ) {

		#print Dumper @var_sets;
		my ( $name, $min, $max, $junk ) = split;
		show_output "die", "Error parsing random_var '$_'\n" unless $max;
		show_output "die", "Error: Too many parameters for random_var '$_'\n" if $junk;

		push @$names_ref, $name;

		if ( @var_sets > 1 ) {

			# Ignore $steps, just add random value to each set.
			for ( my $i = 0 ; $i < @var_sets ; $i++ ) {
				push @{ $var_sets[$i] }, ( 10**( ( log($min) / log(10) ) + rand() * ( ( log($max) / log(10) ) - ( log($min) / log(10) ) ) ) );
			}
		}
		elsif ( 1 == @var_sets ) {

			# Only one var_set, so generate $steps permutations with random numbers.
			my @orig_var_set = @{ $var_sets[0] };
			@var_sets = ();

			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ @orig_var_set, ( 10**( ( log($min) / log(10) ) + rand() * ( ( log($max) / log(10) ) - ( log($min) / log(10) ) ) ) ) ];
			}
		}
		else {
			# No steps yet, so generate $steps random numbers.
			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ ( 10**( ( log($min) / log(10) ) + rand() * ( ( log($max) / log(10) ) - ( log($min) / log(10) ) ) ) ) ];
			}
		}
	}

	# Add random variables.
	foreach ( empty_sort @{ $values_ref->{random_var} } ) {
		my ( $name, $first, $last, $junk ) = split;
		show_output "die", "Error parsing random_var '$_'\n" unless $last;
		show_output "die", "Error: Too many parameters for random_var '$_'\n" if $junk;

		my $range = $last - $first;
		show_output "die", "Error: Invalid range in $_\n" if ( $range < 0 );

		push @$names_ref, $name;

		if ( @var_sets > 1 ) {

			# Ignore $steps, just add random value to each set.
			for ( my $i = 0 ; $i < @var_sets ; $i++ ) {
				push @{ $var_sets[$i] }, ( zero_rand($range) + $first );
			}
		}
		elsif ( 1 == @var_sets ) {

			# Only one var_set, so generate $steps permutations with random numbers.
			my @orig_var_set = @{ $var_sets[0] };
			@var_sets = ();

			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ @orig_var_set, zero_rand($range) + $first ];
			}
		}
		else {
			# No steps yet, so generate $steps random numbers.
			for ( my $i = 0 ; $i < $values_ref->{permutations} ; $i++ ) {
				push @var_sets, [ zero_rand($range) + $first ];
			}
		}
	}

	return @var_sets;
}

sub generate_model_set_files($$$$$) {
	show_output "\n*** Function: generate_model_set_files() ***\n"
	  if ( $verbosity >= 4 );

	my ( $current_directory, $model, $names_ref, $var_setref, $values_ref ) = @_;

	my ( $root_filename, $path ) = fileparse($model);
	$root_filename =~ s/(\.[^.]*)$//;

	my $extension = $1;
	if ( ( $extension ne '.txt' ) && ( $extension ne '.bngl' ) ) {
		show_output "die", "Error: Model extension needs to be .txt or .bngl\n";
	}

	my @model_sets = ();
	my $file_count = 1;
	my $lim        = 0;
	my @generate_network_lines;
	my @run_network_lines;
	my @model_with_free_params;

	if ( $values_ref->{generation} == 1) {
		open FILE, "<", $model
		  or show_output "die", "Error: Cannot open model $model";
		my @orig_model = <FILE>;
		close FILE;

		my $in_parameter_block = 0;
		my %name_matches;
		my $params_inserted_bool = 0;
		for (@orig_model) {
			$in_parameter_block = 1 if (/begin parameters/);
			$in_parameter_block = 0 if (/end parameters/);
			push @model_with_free_params, $_;

			if ($in_parameter_block) {
				for ( my $m = 0 ; $m < @$names_ref ; $m++ ) {
					my $name;

					$name = @$names_ref[$m];

					if ( !$params_inserted_bool ) {
						my $insertion = "$name=1 #Free parameter definition inserted by BioNetFit\n";
						push @model_with_free_params, $insertion;
					}
					if ( $_ =~ /$name/ ) {
						$name_matches{$name} = 1;
					}
				}
				$params_inserted_bool = 1;
			}
		}

		unless ( scalar( keys %name_matches ) >= scalar(@$names_ref) ) {
			show_output "die", "Error: The free parameters specified in your .conf file weren't found in your .bngl file. Check to be sure that the free parameter names specified in your .conf file have matching parameters in your .bngl file\n";
		}

		open OUT, ">", "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension"
		  or show_output "die", "Error: Unable to open $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension for initial model output";
		print OUT @model_with_free_params;
		close OUT;
	}

	if ( $values_ref->{ode} ) {
		my $gen_net_cmd_found = 0;
		my $gen_net_complete  = 0;

		if ( $values_ref->{generation} > 1 ) {
			open FILE, "<", "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension"
			  or show_output "die", "Error: Cannot open model $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension";
			@model_with_free_params = <FILE>;
			close FILE;
		}

		my $has_begin_actions = 0;
		for (@model_with_free_params) {
			if ($gen_net_cmd_found) {
				push @run_network_lines, $_;
			}
			elsif ( $values_ref->{generation} == 1 ) {
				push @generate_network_lines, $_;
			}

			if (/^begin actions/) {
				$has_begin_actions = 1;
			}

			if (/^generate_network\(/) {
				$gen_net_cmd_found = 1;
			}
			elsif ($gen_net_cmd_found) {
				$gen_net_complete = 1 if /\)/;
			}
		}

		if ( $values_ref->{generation} == 1 && !$values_ref->{net_file} ) {
			my $net_maker_filename = "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.bngl";

			open OUT, ">", $net_maker_filename
			  or show_output "die", "Error: Unable to open $net_maker_filename for initial model output";
			print OUT @generate_network_lines;
			print OUT "end actions\n" if $has_begin_actions;
			close OUT;

			my $command         = "$values_ref->{bng_command} --outdir $values_ref->{output_dir}$values_ref->{job_name} $net_maker_filename $values_ref->{pipe_char} $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.txt && touch $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.finished";
			my $command_success = $values_ref->{command_success};

			##$command_success = 0
			##	if $Config{osname} =~ /darwin/;

			show_output "Generating initial network...\n"
			  if ( $verbosity >= 2 );

			unless ( $command_success == system($command) ) {
				show_output "die", "Error generating network.  Check BNG output at $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.txt";
			}

			while (1) {
				last if ( -f "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.finished" );
				sleep 1;
			}
			$values_ref->{net_file} = "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.net";
		}

		$values_ref->{net_file} = "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}.net"
		  unless ( $values_ref->{net_file} );

		open IN, "<", "$values_ref->{net_file}"
		  or show_output "die", "Error: Unable to open $values_ref->{net_file} for input";

		my @orig_net_file = <IN>;
		close IN;

		foreach (@$var_setref) {
			my @var_values = @$_;

			show_output "die", "Error: Count of names and values do not match.\n"
			  unless ( scalar(@var_values) == scalar(@$names_ref) );

			my $new_net_filename = "${current_directory}/${root_filename}_perm${file_count}.net";
			copy "$values_ref->{net_file}", $new_net_filename;

			## Output the lines used to run the model
			my $new_run_filename = "${current_directory}/${root_filename}_perm${file_count}${extension}";
			open OUT, ">", $new_run_filename
			  or show_output "die", "Error: Unable to open $new_run_filename for output";
			print OUT "begin actions\n" if $has_begin_actions;
			print OUT "readFile({file=>\"${new_net_filename}\"});";
			print OUT @run_network_lines;
			close OUT;

			my @changelog;
			## Replace each instance of $name with $value throughout the file.
			for ( my $m = 0 ; $m < @$names_ref ; $m++ ) {
				my $name;
				my $value;

				$name  = @$names_ref[$m];
				$value = ( $var_values[$m] );

				my $in_parameter_block = 0;
				my $index              = 0;
				for (@orig_net_file) {
					$in_parameter_block = 1 if (/begin parameters/);
					$in_parameter_block = 0 if (/end parameters/);

					if ($in_parameter_block) {
						my $replaced = s/\s$name\s.*/ $name $value/;
						if ($replaced) {
							push @changelog, "# $name changed to $value\n";
							$orig_net_file[$index] = $_;
							last;
						}
					}
					$index++;
				}
			}

			show_output "die", "Error: BioNetFit couldn't find/replace some of the free parameters specified in your .conf file. Check to be sure that parameter names in the .conf file match parameter names in the .bngl file\n"
			  if ( scalar(@changelog) != scalar(@$names_ref) );

			push @changelog, "# End of permute change log\n";
			my @new_net_file = ( @changelog, @orig_net_file );

			## Change values in the new .net file
			open OUT, ">", $new_net_filename
			  or show_output "die", "Error: Unable to open $new_net_filename for output";

			print OUT @new_net_file;

			close OUT;
			push @model_sets, $new_run_filename;
			$file_count++;
		}
	}
	else {
		if ( $values_ref->{generation} > 1 ) {
			open FILE, "<", "$values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension"
			  or show_output "die", "Error: Cannot open model $values_ref->{output_dir}$values_ref->{job_name}/${root_filename}_with_defs$extension";
			@model_with_free_params = <FILE>;
			close FILE;
		}

		foreach (@$var_setref) {
			my @var_values = @$_;
			show_output "die", "Error: Count of names and values do not match.\n"
			  unless ( scalar(@var_values) == scalar(@$names_ref) );

			my $new_filename = "${current_directory}/${root_filename}_perm${file_count}${extension}";

			my @new_model = @model_with_free_params;

			# Replace each instance of $name with $value throughout the parameters block.
			my @changelog;
			for ( my $m = 0 ; $m < @$names_ref ; $m++ ) {
				my $name;
				my $value;

				$name  = @$names_ref[$m];
				$value = $var_values[$m];

				my $in_parameter_block = 0;
				my $index              = 0;
				foreach (@new_model) {
					$in_parameter_block = 1 if (/begin parameters/);
					$in_parameter_block = 0 if (/end parameters/);
					if ($in_parameter_block) {
						if ( /^\s*$name\s+/ || /^\s*$name\s*\=\s*/ ) {
							s/$name.*/$name $value/;
							push @changelog, "# $name changed to $value\n";
							$new_model[$index] = $_;
							last;
						}
					}
					$index++;
				}
			}

			show_output "die", "Error: BioNetFit couldn't find/replace any of the free parameters specified in your .conf file. Check to be sure that parameter names in the .conf file match parameter names in the .bngl file\n"
			  if ( !@changelog );

			push @changelog, "# End of permute change log\n";

			@new_model = ( @changelog, @new_model );

			open OUT, ">", $new_filename
			  or show_output "die", "Error: Unable to open $new_filename for output";

			print OUT @new_model;
			close OUT;

			push @model_sets, $new_filename;
			$file_count++;
		}
	}
	return @model_sets;
}

sub run_models_pbs($$$) {
	show_output "\n*** Function: run_models_pbs() ***\n"
	  if ( $verbosity >= 4 );

	my ( $current_directory, $model_setsref, $values_ref ) = @_;

	my @jobs;
	my $counter = 0;
	chdir $current_directory;

	my @perms = @$model_setsref;

	my $smoothing = $values_ref->{smoothing};

	show_output "Generating and running qsub commands...\n"
	  if ( $verbosity >= 3 );

	# This loop runs at least once, even if we don't need to do smoothing
	for ( my $j = 0 ; $j < $smoothing ; $j++ ) {

		# Refill our model_setsref since we shifted the contents in the last round
		@perms = @$model_setsref
		  if $j > 0;

		for ( ; ; ) {

			# $counter keeps track of our how many "chunks" we're sending our job in.
			# If we have 16 permutations to run, and cluster_parallel is set to 8, then we send
			# the jobs off in 2 separate qsub chunks.
			my $multisims;
			my $parallels;
			my $total_sims = 0;
			my $job_name;
			my %all_commands;
			for ( my $i = 0 ; $i < $values_ref->{cluster_parallel} ; $i++ ) {
				my $full_command;
				my @touch_on_error;
				last unless (@perms);
				for ( my $k = 0 ; $k < $values_ref->{multisim} ; $k++ ) {
					$_ = shift(@perms);
					last unless $_;

					my $perm_name;
					my $perm;
					my ( $name, $path, $suffix ) = fileparse( $_, qr/\.[^\.]*$/ );
					($perm) = /_perm(\d+)/;

					$name = $perm
					  unless $name;

					$perm_name = $perm
					  unless $perm_name;

					unless ($perm) {
						show_output "die", "Error: The file $_ does not have _perm, don't know what to do.\n";
					}

					my $outputdir   = "%%%";
					my $outputdir_f = "";

					if ( $values_ref->{smoothing} > 1 ) {
						$outputdir_f = "/$j";
					}

					verbose_mkdir("$values_ref->{output_dir}$values_ref->{job_name}/$values_ref->{generation}$outputdir_f")
					  unless ( -d "$values_ref->{output_dir}$values_ref->{job_name}/$values_ref->{generation}$outputdir_f" );

					my $command = "$values_ref->{bng_command} --outdir $outputdir$outputdir_f $_ $values_ref->{pipe_char} $outputdir/${name}";
					$command .= ( $values_ref->{smoothing} > 1 ) ? "_$j.BNG_OUT" : ".BNG_OUT";

					if ( $smoothing > 1 ) {
						$command .= " && touch %%%/${name}_${j}.finished ";
						push @touch_on_error, "${name}_$j.failed";
					}
					else {
						$command .= " && touch %%%/${name}.finished ";
						push @touch_on_error, "${name}.failed";
					}
					$full_command .= "&& "
					  if $full_command;
					$full_command .= "$command";
				}
				$all_commands{"command$i"} = $full_command;
				for (@touch_on_error) {
					$all_commands{"touch_on_error$i"} .= "&& "
					  if $all_commands{"touch_on_error$i"};
					$all_commands{"touch_on_error$i"} .= "touch $_ ";
				}
			}

			last unless %all_commands;
			$counter++;

			$job_name = "S$values_ref->{generation}_C${counter}";

			my $command = create_pbs_command( $pbs_generic_script, ( fileparse($job_name) ), $values_ref, \%all_commands );

			show_output "die", "Error: Command not created for job $job_name\n"
			  unless $command;

			last unless $command;

			if ( $values_ref->{job_sleep} ) {
				sleep $values_ref->{job_sleep};
			}

			my $job = submit_qsub_job( $values_ref, $command );
			chomp $job;
			push @jobs, $job;
		}
	}

	# Move back to original directory.
	chdir $medir;

	return @jobs;
}

sub run_models_fork($$) {
	show_output "\n*** Function: run_models_fork() ***\n"
	  if ( $verbosity >= 4 );

	my ( $model_setsref, $values_ref ) = @_;

	my $parallel_count = $values_ref->{parallel_count};
	my $per_child      = floor( @$model_setsref / $parallel_count );

	if ( ( $per_child * $parallel_count ) < @$model_setsref ) {

		# Had remainder, add one so none left over at end.
		$per_child++;
	}

	my $outputdir = "$values_ref->{output_dir}$values_ref->{job_name}/$values_ref->{generation}";

	my $perms = scalar(@$model_setsref);
	my @children;
	my @grandchildren;

	my $smoothing = $values_ref->{smoothing};
	my @model_setsref_backup;

	if ( $smoothing > 1 ) {
		@model_setsref_backup = @$model_setsref;
	}

	# Babysitter is a hash which holds the start time, time limit, PID, and end time for each individual process
	my %babysitter = ();

	my $numstarted  = 0;
	my $numfinished = 0;

	show_output "Starting walltime monitor; generating and running simulation commands...\n"
	  if ( $verbosity >= 3 );

	for ( my $i = 0 ; $i < $smoothing ; $i++ ) {

		# Refill our model_setsref since we shifted the contents in the last round
		@$model_setsref = @model_setsref_backup
		  if $i > 0;

		# If user set a walltime...
		if ( exists $values_ref->{max_walltime} ) {

			# Ticker counts how many times we've looped through.  It's only used if we are debugging and want to output directory
			# contents after a set number of loops
			my $ticker = 0;

			# Store our hours, minutes, and seconds in an array
			my @walltimesplit = split( /:/, $values_ref->{max_walltime} );

			# Then convert to seconds
			my $walltime = ( $walltimesplit[0] * 3600 ) + ( $walltimesplit[1] * 60 ) + ( $walltimesplit[2] );

			# Spawn our killer
			my $kpid = fork();

			if ( $kpid == 0 ) {
				sleep 1;
				my $ppid = getppid;

				while (1) {

					# Get a list of files in the current generation's output directory
					opendir CURDIR, "$outputdir"
					  or show_output "Warning, walltime subrouting is unable to open $outputdir\n";

					my @directory_contents = readdir(CURDIR);
					closedir(CURDIR);

					foreach (@directory_contents) {
						my $basename;

						my $ext = "";
						($ext) = $_ =~ /(\.[^.]+)$/;

						# Count the number of .BNG_OUT files.  Each one is indicative of a permutation being started.
						if ( defined $ext && $ext eq ".BNG_OUT" && $numstarted < ( scalar(@$model_setsref) * ( $i + 1 ) ) ) {
							$basename = ( split( /\./, $_ ) )[0];
							my ($tpid) = $_ =~ /\.(\d+)\./;

							# Store the process info in %babysitter if we haven't already
							if ( !exists $babysitter{STARTED}{$basename} ) {

								$numstarted++;
								$babysitter{STARTED}{$basename} = {
									PID   => $tpid,
									LIMIT => ( $walltime + time() ),
								};

								#show_output "adding started run: $basename $i $tpid\n";
							}
						}

						# Count the number of .finished files.  Each one is indicative of a permutation being finished.
						if ( defined $ext && ( $ext eq ".finished" || $ext eq ".failed" ) && $numfinished < ( scalar(@$model_setsref) * ( $i + 1 ) ) ) {
							$basename = ( split( /\./, $_ ) )[0];

							# Store the process info in %babysitter if we haven't already
							if ( !exists $babysitter{FINISHED}{$basename} ) {
								$numfinished++;
								$babysitter{FINISHED}{$basename} = { END => time(), };

								#show_output "added finished run: $basename $i\n";
							}
						}
					}

					# Loop through and do some checking on each process..
					while ( my ( $key, $value ) = each %{ $babysitter{STARTED} } ) {

						# If our STARTED entry doesn't yet have a FINISHED entry, let's check how long it's been running
						if ( !exists $babysitter{FINISHED}{$key} ) {

							# Check to see if it's gone over walltime
							if ( $babysitter{STARTED}{$key}{LIMIT} < time() ) {
								
								# Create a .failed file so other aspects of BioNetFit can know that we have a failed run
								open( TMP, ">>$outputdir/$key.failed" );
								close TMP;
								
								# Kill the perl process and add to kill list
								my $killstatus = kill( 'TERM', $babysitter{STARTED}{$key}{PID} );

								show_output "Run $key ($babysitter{STARTED}{$key}{PID}) went over walltime and was killed with status: $killstatus\n"
								  if ( $verbosity >= 2 );

								# Kill NFsim processes
								my @nf_pids = `ps aux | grep nfsim | awk '{print \$2}'`;
								chomp @nf_pids;
								$killstatus = kill( 'TERM', @nf_pids ) if @nf_pids;

								$babysitter{FINISHED}{$key} = { END => time(), };

								my $basename = $outputdir . ( split( /_/, $key ) )[0] . "_" . ( split( /_/, $key ) )[1];
								$numfinished++;
							}
						}
					}

					# If all runs are finished we can exit the loop
					last
					  if ( $numfinished == ( scalar(@$model_setsref) * ( $i + 1 ) ) );

					# Be sure to stop if the parent exited prematurely
					last
					  unless ( kill 0, $ppid );

					$ticker++;

					#show_output "\nlooping $i.\nnumstarted: $numstarted\nnumfinished: $numfinished\n";

					select( undef, undef, undef, 0.25 );
				}
				exit 0;
			}

			waitpid( $kpid, WNOHANG );
		}

		for ( 1 .. $parallel_count ) {
			my @child_filenames = ();

			#print Dumper $model_setsref;
			for ( 1 .. $per_child ) {
				last unless @$model_setsref;
				my $temp = shift @$model_setsref;
				push @child_filenames, $temp;
			}

			# Don't fork if nothing for child to do.
			last unless (@child_filenames);

			my $pid;
			$pid = fork();

			if ($pid) {

				# parent
				push @children, $pid;
			}

			elsif ( 0 == $pid ) {

				#show_output "Spawned child $$\n";
				#print Dumper @child_filenames;
				foreach (@child_filenames) {
					my $cpid = fork();
					if ($cpid) {
						push @grandchildren, $cpid;

						#show_output "Spawned grandchild $cpid from $$\n";
					}

					elsif ( 0 == $cpid ) {
						my ( $name, $path, $suffix ) = fileparse( $_, qr/\.[^\.]*$/ );
						my $child_root = /^(.*)\.[^.]+/;
						my $basename   = "$path$name";
						my ( $bng_out, $nf_out );

						if ( $smoothing > 1 ) {
							$bng_out = "${basename}_$i.$$.BNG_OUT";
							$nf_out  = "${basename}_$i.$$.NF_OUT";
						}
						else {
							$bng_out = "${basename}.$$.BNG_OUT";
							$nf_out  = "${basename}.$$.NF_OUT";
						}

						my $command = "";

						if ( $values_ref->{smoothing} > 1 ) {
							$outputdir = "$outputdir/$i";
						}

						verbose_mkdir($outputdir);

						$command .= "$values_ref->{bng_command} --outdir $outputdir $_ $values_ref->{pipe_char} $bng_out";

						my $command_success = $values_ref->{command_success};

						# TODO Figure out osx exit code issues.
						# MAC OS has some weird exit-code issues.  We change to 0 from the default here.
						$command_success = 0
						  if $Config{osname} =~ /darwin/;

						show_output "Running: $command\n\n"
						  if ( $verbosity >= 4 );

						unless ( $command_success == system($command) ) {
							show_output "Error running simulation for $_, results will be ignored. Please check to be sure the corrent NFsim binary is set in your .conf file\n";

							if ( $smoothing > 1 ) {
								( $command_success == system("touch ${basename}_$i.failed") )
								  or show_output "Unable to touch ${basename}_$i.failed.\n";
							}
							else {
								( $command_success == system("touch ${basename}.failed") )
								  or show_output "Unable to touch ${basename}.failed.\n";
							}
						}
						else {
							if ( $smoothing > 1 ) {
								( $command_success == system("touch ${basename}_$i.finished") )
								  or show_output "Unable to touch ${basename}_$i.finished.\n";
							}
							else {
								( $command_success == system("touch ${basename}.finished") )
								  or show_output "Unable to touch ${basename}.finished.\n";
							}
						}
						exit(0);
					}
					waitpid( $cpid, 0 );
				}
				exit(0);
			}
			else {
				show_output "die", "Error: Couldn't fork: $!\n";
			}
		}

		my $child_status;

		foreach (@children) {

			#show_output "waiting on pid $_\n";
			waitpid( $_, 0 );

			#show_output "waitpid status for $_: $child_status\n";
		}
	}

	# Need to pause for a second to make sure .failed and .finished files have been touched
	#sleep(1);

	my @exp_files = @{ $values_ref->{exp_file} };
	foreach (@exp_files) {
		fileparse $_;
		s/.exp$//;
	}
	my $model_name = fileparse( $values_ref->{model} );
	$model_name =~ s/.\w+$//;

	# Move and rename output data
	unless ( $values_ref->{smoothing} > 1 ) {
		opendir CURDIR, "$outputdir"
		  or show_output "Couldn't open $outputdir to rename .gdat files\n";

		my @directory_contents = readdir(CURDIR);
		closedir(CURDIR);

		for (@directory_contents) {
			if ( /\.scan$/ | /\.gdat/ ) {
				my $filename = $_;
				$filename =~ s/^${model_name}_//;
				$filename =~ s/^perm\d+_//;
				$filename =~ s/\.\w+$//;
				if ( grep /$filename/, @exp_files ) {
					my $filename = $_;
					$filename =~ s/scan$/gdat/;
					my ($perm_string) = $_ =~ /(_perm\d+)/;
					$filename =~ s/_perm\d+//;
					$filename =~ s/\.gdat$/$perm_string\.gdat/;
					move "$outputdir/$_", "$outputdir/$filename";
				}
			}
		}
	}
	else {
		for ( my $i = 0 ; $i < $values_ref->{smoothing} ; $i++ ) {
			opendir CURDIR, "$outputdir/$i"
			  or show_output "Couldn't open $outputdir/$i to move .gdat files\n";

			my @directory_contents = readdir(CURDIR);
			closedir(CURDIR);

			for (@directory_contents) {
				if ( /\.scan$/ | /\.gdat$/ ) {
					my $filename = $_;
					$filename =~ s/^${model_name}_//;
					$filename =~ s/^perm\d+_//;
					$filename =~ s/\.\w+$//;
					if ( grep /$filename/, @exp_files ) {
						my $filename = $_;
						$filename =~ s/scan$/gdat/;
						my ($perm_string) = $_ =~ /(_perm\d+)/;
						$filename =~ s/_perm\d+//;
						$filename =~ s/\.gdat$/$perm_string\.gdat/;
						$filename =~ s/\.gdat$/_$i.gdat/;

						move "$outputdir/$i/$_", "$outputdir/$filename";
					}
				}
			}
		}
	}

	show_output "Averaging simulation output...\n"
	  if ( $verbosity >= 3 );

	smooth_runs( $outputdir, $values_ref )
	  if ( $smoothing > 1 );

	opendir CURDIR, "$outputdir"
	  or show_output "die", "Error: Unable to open $outputdir to check for failed runs\n";

	my @directory_contents = readdir(CURDIR);
	closedir(CURDIR);

	my $failed = grep( /.failed/, @directory_contents );

	show_output "$failed runs failed. This could happen because they went over walltime, or because the NFsim executable had a problem. Check the .BNG_OUT and/or .NF_OUT files for more info.\n\n"
	  if ( $failed && $verbosity >= 1 );

	return ( $perms - $failed );
}

sub smooth_runs($$) {
	my ( $current_directory, $values_ref ) = @_;

	my %all_runs = ();

	opendir CURDIR, "$current_directory"
	  or show_output "die", "Error: Unable to open $current_directory for run averaging.\n";
	my @directory_contents = readdir(CURDIR);
	closedir(CURDIR);

	my @exp_files = @{ $values_ref->{exp_file} };
	foreach (@exp_files) {
		fileparse $_;
		s/.exp$//;
	}

	my $model_name = fileparse( $values_ref->{model} );
	$model_name =~ s/.\w+$//;

	my $basename;

	foreach (@directory_contents) {
		my ( $basename, $path, $ext ) = fileparse( $_, qr/\.[^.]*/ );

		if ( defined $ext && $ext eq ".gdat" ) {
			$basename =~ s/_perm\d+.+//;
			$basename =~ s/^${model_name}_//;

			next unless grep /$basename/, @exp_files;
			
			my ($perm) = ($_ =~ /(perm\d+)/);
			$perm =~ s/perm//;
			my ($subperm) = $_ =~ /(_\d+\.gdat$)/;
			$subperm =~ s/_//;
			$subperm =~ s/\.gdat//;
			
			next
				if (-f "${model_name}_perm$perm.failed");
				
			# Skip the .gdat file if it has no sub-perm in the filename (meaning it is incomplete)
			next
			  if ( !defined $subperm );

			my @perm_data = load_data("$current_directory/$_");

			my $increment = 0;
			while (1) {
				if ( !defined $all_runs{$basename}{$perm}[$increment] ) {
					$all_runs{$basename}{$perm}[$increment] = $perm_data[0];
					last;
				}
				$increment++;
			}
		}
	}

	#	%all_runs structure:
	#	1 => {
	#			[							# $all_runs{$permutation}
	#				[						# $all_runs{$permutation}[0]
	#					{	time=>0.00	}	# $all_runs{$permutation}[0][0]
	#					{	time=>1.00	}	# $all_runs{$permutation}[0][1]
	#				],
	#				[
	#					{	time=>0.00	}	# $all_runs{$permutation}[1][0]
	#					{	time=>0.00	}	# $all_runs{$permutation}[1][1]
	#				]
	#			]
	#		}
	#	2 => {
	#			[
	#				[
	#					{	time=>0.00	}
	#				],
	#				[
	#					{	time=>0.00	}
	#				]
	#			]
	#		}

	my $num_iterations = 0;
	my $sum            = 0;

	my @averaged = ( {} );
	my $out_file;

	my $permutation;

	# For each exp file
	foreach $basename ( keys %all_runs ) {
		undef @averaged;
		my @columns_ref;
		for ( keys %{ $all_runs{$basename}{1}[0][0] } ) {
			push @columns_ref, $_;
		}

		my $control_col = "time";
		if ( $values_ref->{scan_parameter} ) {
			foreach ( @{ $values_ref->{scan_parameter} } ) {
				if ( ( split( ' ', $_ ) )[0] eq $basename ) {
					$control_col = ( split( ' ', $_ ) )[1];
				}
			}
		}

		# For each permutation
		foreach $permutation ( keys %{ $all_runs{$basename} } ) {

			# For each column name
			my $column;
			foreach $column ( keys %{ $all_runs{$basename}{$permutation}[0][0] } ) {
				next
				  if ( $column eq "time" || $column eq $control_col );

				# For each time group
				my $time;
				for $time ( 0 .. $#{ $all_runs{$basename}{$permutation}[0] } ) {

					# For each iteration of the current permutation
					my $iteration;
					for $iteration ( 0 .. $#{ $all_runs{$basename}{$permutation} } ) {
						$sum += $all_runs{$basename}{$permutation}[$iteration][$time]{$column};
						$num_iterations++;
					}

					#show_output "total for $column in permutation $permutation at timepoint $time: $sum\n";
					#show_output "average for $column in permutation $permutation at timepoint $time: ($sum / $num_iterations)\n";

					my $average = $sum / $num_iterations;

					$sum            = 0;
					$num_iterations = 0;

					$averaged[$time]{$control_col} = $all_runs{$basename}{$permutation}[0][$time]{$control_col};
					$averaged[$time]{$column}      = $average;
				}
			}

			my $output_file = "${current_directory}/${model_name}_${basename}_perm${permutation}.gdat";

			open( FH, '>', $output_file )
			  or show_output "die", "Error: Unable to open $output_file for writing.\n";

			printf FH "%s", "#";

			if ( $values_ref->{scan_parameter} ) {
				printf FH " %13s", $control_col;
			}
			else {
				printf FH " %13s", "time";
			}

			for my $column_out (@columns_ref) {
				if ( $column_out eq "time" || $column_out eq $control_col ) {
					next;
				}
				else {
					printf FH "   %13s", $column_out;
				}
			}

			printf FH "%s", "\n";

			my $timepoint;
			for $timepoint ( 0 .. $#{averaged} ) {
				printf FH "%-1s", " ";

				if ( $control_col ne "time" ) {
					printf FH "%-1s", $averaged[$timepoint]{$control_col};
				}
				else {
					printf FH "%-1s", $averaged[$timepoint]{"time"};
				}

				for my $column (@columns_ref) {
					if ( $column eq "time" || $column eq $control_col ) {
						next;
					}
					else {
						printf FH "  %.8e", $averaged[$timepoint]{$column};
					}
				}
				printf FH "%s", "\n";
			}
			close FH;
		}
	}

	#print Dumper @averaged;

	# Now let's take care of our .failed files
	my @failed_runs = ();
	foreach (@directory_contents) {
		my $ext;
		($ext) = $_ =~ /(\.[^.]+)$/;

		if ( defined $ext && $ext eq ".failed" ) {
			$basename = $_;

			#$basename = ((split(/_/, $_))[0]) . "_" . ((split(/_/, $_))[1]);
			$basename =~ s/_\d+\.failed/\.failed/;
			push @failed_runs, $basename;
		}
	}

	my %dupes = ();
	map( $dupes{$_}++, @failed_runs );

	for (@failed_runs) {
		unlink glob "$current_directory/*.failed";
	}

	foreach my $key ( keys %dupes ) {
		if ( $dupes{$key} == ( $values_ref->{smoothing} ) ) {
			open( TMP, ">>$current_directory/$key" );
			close TMP;
		}
	}
}

sub zero_rand($) {
	my ($range) = @_;
	return 0 unless ($range);

	return rand($range);
}

sub rlnorm($$) {
	my ( $mean, $stddev ) = @_;

	my $r1 = sqrt( -2 * log( rand() ) ) * sin( 6.28318531 * rand() );

	return exp( log($mean) + $r1 * $stddev );
}

sub empty_sort(@) {
	return () unless (@_);

	return sort(@_);
}

sub numalpha_sort {
	my ( $adir, $anum, $aextra ) = $a =~ /k([mp])(\d+)(.*)/;
	my ( $bdir, $bnum, $bextra ) = $b =~ /k([mp])(\d+)(.*)/;

	if ( $adir && $bdir ) {

		# Titles are of the form kp17 or km14_blah.

		return ( $anum <=> $bnum ) unless ( $anum == $bnum );

		# kp* should come before km*
		return ( $bdir cmp $adir ) unless ( $adir eq $bdir );

		return ( $aextra cmp $bextra );
	}

	my ($numa) = $a =~ /(\d+)/;

	my ($numb) = $b =~ /(\d+)/;

	return $numa <=> $numb if ( $numa && $numb && ( $numa != $numb ) );
	return $a cmp $b;
}

sub combine_columns($$$@) {
	show_output "\n*** Function: combine_columns() ***\n"
	  if ( $verbosity >= 4 );

	my $scan_parameter    = shift @_;
	my $combined_filename = shift @_;

	# First file is experiment file.
	my $exp = shift @_;

	my @filenames = sort numalpha_sort @_;

	unshift @filenames, $exp
	  if ($exp);

	# Remove all extensions for use in column headers.
	my @fn_roots = @filenames;
	for ( my $i = 0 ; $i < @fn_roots ; $i++ ) {

		# Do not remove exp extension.
		$fn_roots[$i] =~ s/\.(gdat)$//;

		if ( $fn_roots[$i] =~ /_(perm\d+)/ ) {

			# Just use the perm designation.
			$fn_roots[$i] = $1;
		}
	}

	# Open all files.
	my @filehandles;
	my $header = "";
	my @used_columns;
	my @all_columns_lookup;
	my $filename;

	foreach $filename (@filenames) {
		my $fh;
		my $fhh;

		open $fh, $filename
		  or show_output "die", "Error: Unable to open $filename in combine_columns\n";

		# Get header.
		my $new_header = <$fh>;

		# Remove line endings.
		$new_header =~ tr/\n\r//d;

		# Remove leading space.
		$new_header =~ s/^\s*//;

		my @new_columns = split( /\s+/, $new_header );

		show_output "die", "Error: No header in $_\n"
		  unless @new_columns;

		# Remove hash from first column.
		my $hash = shift @new_columns;

		show_output "die", "Error: File $_ missing hash in first column of first line.\n"
		  unless ( $hash eq "#" );

		unless (@used_columns) {

			# Create initial list of columns to use.
			@used_columns = @new_columns;
		}

		# Only use columns that are present in all files.
		my @new_used_columns = ();
		my %new_columns      = ();

		for ( my $i = 0 ; $i < @new_columns ; $i++ ) {

			# Create column name lookup table for each file.
			$new_columns{ $new_columns[$i] } = $i;
			if ( grep /^$new_columns[$i]$/, @used_columns ) {

				# Only output columns that all files contain.
				push @new_used_columns, $new_columns[$i];
			}
		}

		@used_columns = @new_used_columns;

		push @all_columns_lookup, \%new_columns;
		push @filehandles,        $fh;
	}

	# Output combined header.
	open my $combined_fh, ">", $combined_filename
	  or show_output "die", "Error: Unable to open $combined_filename for writing in combine_columns\n";

	# Get rid of first column name if it is just a hashmark.
	shift @used_columns if ( $used_columns[0] eq "#" );

	# Real first column should be time.
	show_output "die", "Error: Time not found as first column in all files\n"
	  unless ( $used_columns[0] =~ /time/i || $used_columns[0] =~ /$scan_parameter/i );

	# Print time.
	printf $combined_fh "# %1s", shift @used_columns;

	# Print experiment file separately.
	my @non_exp_roots = @fn_roots;
	shift @non_exp_roots;

	foreach my $column (@used_columns) {
		printf $combined_fh " %1s", "${column}_exp";
		foreach my $fn (@non_exp_roots) {
			printf $combined_fh " %1s", "${column}_$fn";
		}
	}
	printf $combined_fh "\n";

	# Add time column back on.
	unshift @used_columns, "$scan_parameter";

	# Output each combined line in turn.
	my $exp_line_no = 0;
	for ( ; ; ) {
		my $file_end  = 0;
		my @all_lines = ();
		foreach (@filehandles) {
			my @line;
			for ( ; ; ) {
				my $line = <$_>;
				last unless ($line);
				$line =~ tr/\n\r//d;

				# Remove any leading spaces.
				$line =~ s/^\s+//;
				print $line;
				@line = split /\s+/, $line;

				# Don't try to read anymore if nothing left in file.
				last unless @line;

				# Read line from exp file no matter what so that we
				# have something to compare to.
				last unless ( exists( $all_lines[0] ) && exists( $all_lines[0][0] ) );

				# Use current line if timepoint exists in experiment file.
				# Read next line if timepoint does not exist in experiment file.
				last if ( $all_lines[0][0] == $line[0] );
			}

			if (@line) {
				print "\n";
				push @all_lines, \@line;

			}
			else {
				$file_end = 1;
			}
		}

		if ($file_end) {
			die "Error, .exp file and simulation output (.gdat files) are of different length. Please be sure both files start and end at the same timepoint.\n"
			  if (@all_lines);
			last;
		}

		# Make sure time agrees in all files.
		my $current_time = $all_lines[0][0];
		for ( my $i = 0 ; $i < @filehandles ; $i++ ) {
			#print "ct is $current_time and other is $all_lines[$i][0]\n";
			if ( $current_time != $all_lines[$i][0] ) {
				my $string = "Error: Time does not match in $filenames[$i]: $current_time != $all_lines[$i][0].\n";
				show_output "die", $string;
			}
		}

		printf $combined_fh " %.8e", $current_time;

		for ( my $i = 1 ; $i < @used_columns ; $i++ ) {
			for ( my $j = 0 ; $j < @filehandles ; $j++ ) {
				if ( $all_lines[$j][ $all_columns_lookup[$j]{ $used_columns[$i] } ] eq "NaN" ) {
					printf $combined_fh " %s", "NaN";
				}
				else {
					printf $combined_fh " %.12e", $all_lines[$j][ $all_columns_lookup[$j]{ $used_columns[$i] } ];
				}
			}
		}

		printf $combined_fh "\n";

		$exp_line_no++;
	}
}

sub finished($$$) {
	my ( $values_ref, $orig_config_file, $exit_message ) = @_;

	$exit_message .= "Results can be found in $values_ref->{output_dir}$values_ref->{job_name}/Results\n"
	  unless $values_ref->{bootstrap};

	my ( $best_run_filename, $best_chi, $best_perm ) = run_consolidate( $values_ref, $orig_config_file, $exit_message );

	unlink(".lock_$values_ref->{job_name}");

	my $string = "End timestamp: " . time() . "\n";
	show_output $string;

	exit 0
	  unless $values_ref->{bootstrap};

	return ( $best_run_filename, $best_chi );
}

sub graph_chi($) {

	my ($values_ref) = @_;
	my ( @x, @y );
	my @lines;
	my @graph_values;

	my $filename = "$values_ref->{output_dir}$values_ref->{job_name}/Results/sorted_params.txt";

	open SUMM, $filename
	  or show_output "die", "Error: Cannot open sorted_params.txt for graphing";

	my $line;
	while (<SUMM>) {
		next if $. == 1;    # skip header line

		$line = ( split /\s+/ )[0];
		$line .= " ";
		$line .= ( split /\s+/ )[1];
		$line =~ s/gen//;
		$line =~ s/perm/-/;

		push @lines, $line;
	}
	close SUMM;

	@lines = sort { ( split /[-\s]/, $a )[1] <=> ( split /[-,\s]/, $b )[1] } @lines;
	@lines = sort { ( split /[-\s]/, $a )[0] <=> ( split /[-,\s]/, $b )[0] } @lines;

	my $x;
	my $y;

	my $current_gen;

	foreach my $perm (@lines) {
		$x = ( split( /\s+/, $perm ) )[0];
		$y = ( split( /\s+/, $perm ) )[1];

		push @x, $x;
		push @y, $y;
	}

	@graph_values = ( \@x, \@y );
	my $max = max(@y);
	
	my $x_lbl_skip = scalar(@x)/10;

	# Create the graph and give it dimensions
	my $mygraph = GD::Graph::lines->new( 800, 600 );
	$mygraph->set_title_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_legend_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_x_label_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_y_label_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_x_axis_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_y_axis_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_values_font( 'fonts/OpenSans-Regular.ttf', 9 );
	$mygraph->set_text_clr("black");

	# Set the graph options
	$mygraph->set(
		x_label  => 'Permutation',
		y_label  => 'Goodness of Fit',
		title    => 'Goodness of Fit vs. Permutation',
		r_margin => 10,

		transparent => 0,
		bgclr       => 'white',
		fgclr       => 'black',

		#y_max_value => $height,
		labelclr     => 'black',
		legendclr    => 'black',
		axislabelclr => 'black',
		textclr      => 'black',

		#x_tick_number => 10,
		x_label_skip      => $x_lbl_skip,
		x_long_ticks      => 0,
		x_label_position  => 0.5,
		x_labels_vertical => 1,

		y_tick_number => 10,

		#y_label_skip => 100,
		#y_long_ticks => 1,
		y_label_position => 0.5,
		y_max_value      => $max * 1.01,
		y_number_format  => \&y_format,

		legend_placement => 'RC',

		#long_ticks => 1,

		line_width => 3,
	) or show_output $mygraph->error;

	# Make the image
	my $myimage = $mygraph->plot( \@graph_values ) or show_output "die", $mygraph->error;

	open OUT, '>', "$values_ref->{output_dir}$values_ref->{job_name}/Results/objfunc_$values_ref->{job_name}.png"
	  or show_output "die", "Error: Couldn't open for output: $!";

	binmode(OUT);
	print OUT $myimage->png();
	close OUT;
}

sub y_format {
	my $value = shift;
	my $ret;

	$ret = sprintf( "%.2e", $value );

	return $ret;
}

sub graph_outs($) {

	my ($values_ref) = @_;
	my @lines;
	my @graph_values;
	my $all_together = 0;

	my $all_summary_file = "$values_ref->{output_dir}$values_ref->{job_name}/Results/sorted_params.txt";
	my $best_file;
	my $best_gen;
	my $best_perm;

	my @outputs;

	my $model_root = $values_ref->{model};
	$model_root =~ s/(\.[^.]*)$//;
	my $stripped_model_root = fileparse($model_root);

	my @exp_files = @{ $values_ref->{exp_file} };

	my ( $expdata, $expcolumns ) = load_exp( \@exp_files, $values_ref );
	my @exp_roots = keys %$expdata;
	my %data      = ();

	foreach my $root_filename (@exp_roots) {
		my $exp_name = fileparse($root_filename);
		$exp_name =~ s/.exp$//;
		$exp_name =~ s/${stripped_model_root}_//;

		my $filetemp = fileparse( ${root_filename} );

		my $best_data_file = "$values_ref->{output_dir}$values_ref->{job_name}/Results/${exp_name}_bestfit.gdat";

		my ( $data_ref, @columns ) = load_data($best_data_file);
		$data{$root_filename} = $data_ref;
	}

	foreach my $exp_root (@exp_roots) {

		my $filename = fileparse($exp_root);
		$filename =~ s/${stripped_model_root}_//;

		# Control column will not be 'time' in dose-response data
		my $control_col = "time";
		if ( $values_ref->{scan_parameter} ) {
			foreach ( @{ $values_ref->{scan_parameter} } ) {
				if ( ( split( ' ', $_ ) )[0] eq $filename ) {
					$control_col = ( split( ' ', $_ ) )[1];
				}
			}
		}

		my @sim = @{ $data{$exp_root} };
		my @exp = @{ $expdata->{$exp_root} };

		my $num_exp_outputs  = scalar @exp;
		my $num_best_outputs = scalar @sim;
		my @x;

		foreach my $col ( @{ $expcolumns->{$exp_root} } ) {
			my @y_exp;
			my @y_best;

			next if ( $col =~ /_SD/ );

			for ( my $i = 0 ; $i < @sim ; $i++ ) {
				for ( my $j = 0 ; $j < @exp ; $j++ ) {
					if ( $col eq $control_col ) {
						next if ( defined $x[$i] && $x[$i] );
					}
					else {
						next if ( ( defined $y_exp[$i] && $y_exp[$i] ) );
					}
					if ( abs( $exp[$j]{$control_col} ) == abs( $sim[$i]{$control_col} ) ) {
						if ( $col eq $control_col ) {
							$x[$i] = sprintf( "%.2e", $sim[$i]{$col} );
						}
						else {
							if ( $exp[$j]{$col} eq "NaN" ) {
								$y_exp[$i] = undef;
							}
							else {
								$y_exp[$i] = abs( $exp[$j]{$col} );
							}
						}
					}
					else {
						if ( $col eq $control_col ) {
							$x[$i] = '';
						}
						else {
							$y_exp[$i] = undef;
						}
					}
				}
				push @y_best, abs( $sim[$i]{$col} );
			}
			next if ( $col eq $control_col );

			@graph_values = ( \@x, \@y_exp, \@y_best );

			my @legend_keys = ( "${col}_EXP", "${col}_BEST" );

			#my $height = max(max(@y_exp),max(@y_best));

			my $exp_max = 0;
			my $exp_min = 1000000000000;
			for ( my $i = 0 ; $i < @y_exp ; $i++ ) {
				next if ( !defined $y_exp[$i] );
				$exp_max = $y_exp[$i] if $y_exp[$i] > $exp_max;
				$exp_min = $y_exp[$i] if $y_exp[$i] < $exp_min;
			}
			my $best_max = 0;
			my $best_min = 1000000000000;
			for ( my $i = 0 ; $i < @y_best ; $i++ ) {
				next if ( !defined $y_best[$i] );
				$best_max = $y_best[$i] if $y_best[$i] > $best_max;
				$best_min = $y_best[$i] if $y_best[$i] < $best_min;
			}

			my $all_max = max( $exp_max, $best_max );
			my $all_min = min( $exp_min, $best_min );

			# Create the graph and give it dimensions
			my $mygraph = GD::Graph::mixed->new( 800, 600 );
			$mygraph->set_legend(@legend_keys);
			$mygraph->set_title_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_legend_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_x_label_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_y_label_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_x_axis_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_y_axis_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_values_font( 'fonts/OpenSans-Regular.ttf', 9 );
			$mygraph->set_text_clr("black");

			my $x_tick_num = scalar(@x);

			# Set the graph options
			$mygraph->set(
				x_label => $control_col,
				y_label => $col,
				title   => 'Best Fit vs Experimental Data',

				transparent => 0,
				bgclr       => 'white',

				y_max_value => $all_max * 1.01,
				y_min_value => $all_min,

				#x_tick_number => "auto",
				#x_label_skip => scalar,
				x_long_ticks      => 0,
				x_label_position  => 0.5,
				x_labels_vertical => 1,

				#y_tick_number => 10,
				#y_label_skip => 100,
				#y_long_ticks => 1,
				y_label_position => 0.5,
				y_number_format  => \&y_format,
				dclrs            => [qw(red black)],

				legend_placement => 'RC',

				#long_ticks => 1,

				line_width => 2,

				types => [qw( points lines )],

				marker_size => 3,
				markers     => [7],
			) or show_output $mygraph->error;

			# Make the image
			my $myimage = $mygraph->plot( \@graph_values )
			  or show_output "die", $mygraph->error;

			open OUT, '>', "$values_ref->{output_dir}$values_ref->{job_name}/Results/exp_vs_best_${col}_$values_ref->{job_name}.png"
			  or show_output "die", "Error: Couldn't open for output: $!";
			binmode(OUT);
			print OUT $myimage->png();
			close OUT;
		}
	}
}

sub get_arch_info($) {
	my ($values_ref) = @_;

	# Run uname in order to get OS and architecture information
	my $os = $Config{osname};

	# We'll assume we're using 32-bit by default
	my $using_64 = 0;

	# But if perl is able to do 8-bit integers we'll assume 64-bit
	if ( $Config{longsize} == 8 ) {
		$using_64 = 1;
	}

	# This runs once.  Determines proper command_success and nfsim_command based on OS and architecture

	# Using Cygwin
	if ( $os =~ /cygwin/i ) {

		# Different systems return different success codes when we run a command, so we need to store the
		# correct success code to tell if our commands are running correctly.
		if ( !exists $values_ref->{command_success} ) {
			$values_ref->{command_success} = 0;
		}

		# Similarly, different systems allow the use of different pipe methods for storing the output from a command..
		if ( !exists $values_ref->{pipe_char} ) {
			$values_ref->{pipe_char} = ">&";
		}
	}

	# Using Mac OS X
	elsif ( $os =~ /darwin/i ) {
		if ( !exists $values_ref->{command_success} ) {
			$values_ref->{command_success} = -1;
		}

		if ( !exists $values_ref->{pipe_char} ) {
			$values_ref->{pipe_char} = ">&";
		}
	}

	# Using Linux
	elsif ( $os =~ /linux/i ) {
		if ( !exists $values_ref->{command_success} ) {
			$values_ref->{command_success} = -1;
		}

		if ( !exists $values_ref->{pipe_char} ) {
			$values_ref->{pipe_char} = ">";
		}

	}

	# If all else fails, let's default to options that will work with x86 Linux
	else {
		if ( !exists $values_ref->{command_success} ) {
			$values_ref->{command_success} = -1;
		}

		if ( !exists $values_ref->{pipe_char} ) {
			$values_ref->{pipe_char} = ">";
		}
	}

	my $answer;

	if ( !-x "$values_ref->{bng_command}" ) {
		do {
			print "Your BioNetGen executable does not have execute permissions. Would you like me to fix this? ";
			$answer = <STDIN>;
			chomp $answer;
		} until ( ( $answer eq 'y' ) || ( $answer eq 'n' ) || ( $answer eq 'Y' ) || ( $answer eq 'N' ) );

		show_output "die", "Please set your BioNetGen as executable (chmod +x BNG2.pl).\n"
		  if ( ( $answer eq 'N' ) || ( $answer eq 'n' ) );

		chmod 0755, "$values_ref->{bng_command}"
		  or show_output "die", "Error: BioNetFit wasn't able to make your BNG2.pl executable. Please set your BioNetGen as executable (chmod +x BNG2.pl).\n";
	}

	if ( $values_ref->{use_cluster} && !-x "$0" ) {
		do {
			print "BioNetFit has permission issues that could prevent qsub from running it properly. Would you like me to fix this? ";
			$answer = <STDIN>;
			chomp $answer if $answer;
		} until ( ( $answer eq 'y' ) || ( $answer eq 'n' ) || ( $answer eq 'Y' ) || ( $answer eq 'N' ) );

		show_output "die", "Please set your BioNetFit binary as executable (chmod +x [BioNetFit.pl]).\n"
		  if ( ( $answer eq 'N' ) || ( $answer eq 'n' ) );

		chmod 0755, "$0"
		  or show_output "die", "Error: BioNetFit wasn't able to make itself executable. Please set your BioNetFit binary as executable (chmod +x [BioNetFit.pl]).\n";
	}
}

sub bootstrap_run($$$) {
	my ( $values_ref, $total_sum, $best_run_filename ) = @_;

	show_output "Preparing to bootstrap next run...\n"
	  if ( $verbosity >= 2 );

	my $param_filename = "$best_run_filename.bngl";

	unless ( -f $param_filename ) {
		$param_filename = "$best_run_filename.txt";
	}

	# Load parameter values from best sim
	my ( $names_ref, $vars_ref ) = get_model_params( $param_filename, $values_ref );

	# Set up directories and file names
	my $bootstrap_dir      = "$values_ref->{output_dir}$values_ref->{job_name}_bootstrap";
	my $bootstrap_filename = "$bootstrap_dir/params.txt";
	my $new_config_file    = "$bootstrap_dir/bootstrap.conf";
	my $bootstrap_num      = $values_ref->{bootstrap_num};

	show_output "Setting up new configuration file...\n"
	  if ( $verbosity >= 3 );

	# Fill values_ref with original config file so we can use it to output a new bootstrapping config file
	my $orig_config_file = $values_ref->{orig_config};
	my @config_files     = ($orig_config_file);
	my $values_ref_new   = parse_file(@config_files);
	update_values($values_ref_new);

	# Get out bootstrap_num in order
	$values_ref_new->{bootstrap_num} = $bootstrap_num;
	$values_ref_new->{bootstrap_num} = $values_ref_new->{bootstrap_num} + 1
	  unless ( $total_sum >= $values_ref_new->{bootstrap_chi} );

	if ( $total_sum >= $values_ref_new->{bootstrap_chi} ) {
		$values_ref_new->{bootstrap_retry} += 1;
	}
	else {
		$values_ref_new->{bootstrap_retry} = 0;
	}

	# Change a few values in new config file
	$values_ref_new->{show_welcome_message} = 0;
	$values_ref_new->{ask_create}           = 0;
	$values_ref_new->{ask_overwrite}        = 0;

	# Output new config file
	unless ( $values_ref_new->{bootstrap_num} > $values_ref_new->{bootstrap} ) {
		show_output "Outputting new config file...\n"
		  if ( $verbosity >= 3 );

		open NEW_CONFIG, '>', $new_config_file
		  or show_output "die", "Error: Unable to open $new_config_file for writing.\n";

		print NEW_CONFIG "# Automatically generated, job:$values_ref_new->{job_name}\n";

		foreach my $name ( sort keys %$values_ref_new ) {
			if ( grep( /^$name$/, @list_variables ) || grep( /^$name$/, @var_variables ) ) {
				for ( my $i = 0 ; $i < @{ $values_ref_new->{$name} } ; $i++ ) {
					print NEW_CONFIG "$name $values_ref_new->{$name}[$i]\n";
				}
			}
			else {
				print NEW_CONFIG "$name $values_ref_new->{$name}\n";
			}
		}
		close NEW_CONFIG;
	}

	# Output to params file

	unless ( $total_sum >= $values_ref->{bootstrap_chi} ) {
		show_output "Writing to params.txt (chi^2 of $total_sum is less than limit of $values_ref->{bootstrap_chi})\n"
		  if ( $verbosity >= 3 );

		my @var_names  = @$names_ref;
		my @var_values = @$vars_ref;
		my $exists     = 0;

		if ( -f $bootstrap_filename ) {
			$exists = 1;
		}

		open BOOTSTRAP, '>>', $bootstrap_filename
		  or show_output "die", "Error: Unable to open $bootstrap_filename for writing.\n";

		unless ($exists) {
			printf BOOTSTRAP "%-1s",  "Run";
			printf BOOTSTRAP " %-1s", "Chi-Sq";
			foreach (@var_names) {

				#printf BOOTSTRAP " %-14s", $_;
				printf BOOTSTRAP " %-1s", $_;
			}
			printf BOOTSTRAP "\n";
		}
		printf BOOTSTRAP "%-1s",   $values_ref_new->{bootstrap_num} - 1;
		printf BOOTSTRAP " %.12f", $total_sum;
		foreach (@var_values) {
			if ( looks_like_number($_) ) {
				printf BOOTSTRAP " %.8e", $_;
			}
			else {
				printf BOOTSTRAP " %.8e", 0;
			}
		}

		print BOOTSTRAP "\n";
		close BOOTSTRAP;
	}

	# Only save run results if our chi^2 is less than bootstrap_chi
	unless ( $total_sum >= $values_ref_new->{bootstrap_chi} ) {
		show_output "Moving results to bootstrap directory...\n"
		  if ( $verbosity >= 3 );

		my $current = $values_ref_new->{bootstrap_num} - 1;
		if ( -d "$bootstrap_dir/Results_$current" ) {
			rmtree("$bootstrap_dir/Results_$current");
		}
		verbose_mkdir("$bootstrap_dir/Results_$current");
		my @files_to_copy = glob "$values_ref_new->{output_dir}/$values_ref_new->{job_name}/Results/*";
		for (@files_to_copy) {
			copy( $_, "$bootstrap_dir/Results_$current" );
		}
	}
	else {
		my $current = $values_ref_new->{bootstrap_num};
		show_output "\nChi^2 of best run ($total_sum) was not below threshold of $values_ref->{bootstrap_chi}. Rerunning $current of $values_ref->{bootstrap}.\n";

		if ( $values_ref_new->{bootstrap_retry} > $values_ref->{bootstrap_retries} ) {
			show_output "die", "\nWe've retried this run the maximum number of times allowed by 'bootstrap_retries' in your .conf file. You can find your results in $values_ref->{output_dir}$values_ref->{job_name}.\n";
		}
	}

	# Done!
	if ( $values_ref_new->{bootstrap_num} > $values_ref_new->{bootstrap} ) {
		show_output "\nBootstrapped $values_ref_new->{bootstrap} runs. Bootstrapping complete. Results can be found in $bootstrap_dir.\n";
	}

	# More bootstrap runs necessary, do another one
	else {
		# Create the command to run the next generation
		my $command = "perl $0 $submit_command $new_config_file";

		show_output "\nStarting new bootstrap run ($values_ref_new->{bootstrap_num} of $values_ref_new->{bootstrap}).\n"
		  if ( $verbosity >= 1 );

		( $values_ref->{command_success} == system($command) )
		  or show_output "die", "Error: Unable to run command: $command\n";
		exit 0;
	}
}

sub show_output {
	my $action;
	my $message;

	if ( @_ == 1 ) {
		$message = $_[0];
	}
	else {
		$action  = $_[0];
		$message = $_[1];
	}

	# Save all output to file if we're doing a qsub run
	if ($qsub) {
		my $opened = 1;

		open( OUTFILE, '>>', $qsub )
		  or $opened = 0;
		print OUTFILE $message if $opened;
		close(OUTFILE) if $opened;
	}
	else {
		print "$message";
	}

	if ( $action && $action eq "die" ) {
		print "\n";
		unlink "$medir/.lock_$qsub";
		die;
	}
}

sub generate_bootstrap_file ($) {

	my ( $values_ref, $model_roots ) = @_;

	my @bootstrap = ( {} );

	my @exp_files = @{ $values_ref->{exp_file} };
	my ( $expdata, $expcolumns ) = load_exp( \@exp_files, $values_ref );
	my @exp_roots = keys %$expdata;

	my $model_name = fileparse( $values_ref->{model} );
	$model_name =~ s/.\w+$//;

	foreach my $root_filename (@exp_roots) {
		undef @bootstrap;
		my @exp = @{ $expdata->{$root_filename} };

		my $filename = fileparse($root_filename);
		$filename =~ s/${model_name}_//;

		# Control column will not be 'time' in dose-response data
		my $control_col = "time";
		if ( $values_ref->{scan_parameter} ) {
			foreach ( @{ $values_ref->{scan_parameter} } ) {
				if ( ( split( ' ', $_ ) )[0] eq $filename ) {
					$control_col = ( split( ' ', $_ ) )[1];
				}
			}
		}

		foreach my $col ( @{ $expcolumns->{$root_filename} } ) {
			next if ( $col eq $control_col );

			# Fill exp hash with 0s.
			for ( my $i = 0 ; $i < @exp ; $i++ ) {
				$bootstrap[$i]{$col} = 0;
			}

			# Add 1 to each timepoint selected randomly
			for ( my $i = 0 ; $i < @exp ; $i++ ) {
				my $timepoint = int( rand(@exp) );
				$bootstrap[$timepoint]{$col} = $bootstrap[$timepoint]{$col} + 1;
			}
		}

		my $bootstrap_dir = "$values_ref->{output_dir}$values_ref->{job_name}_bootstrap";

		open RAND_EXP, '>', "$bootstrap_dir/$filename.dat"
		  or die "Error: Can't open $bootstrap_dir/filename.dat for writing.\n";

		for ( my $i = 0 ; $i < @bootstrap ; $i++ ) {
			foreach my $col ( sort numalpha_sort @{ $expcolumns->{$root_filename} } ) {
				next if ( $col eq $control_col );
				print RAND_EXP "$bootstrap[$i]{$col} ";
			}
			print RAND_EXP "\n";
		}
		close RAND_EXP;
	}

	show_output "Created new random bootstrap map\n"
	  if ( $verbosity >= 3 );
}

sub run_monitor ($) {
	my ($values_ref) = @_;

	my $finished_flag = 0;
	my %finished_generations;
	my $current_gen   = 1;
	my $last_filesize = 0;
	my $last_line     = 0;
	my $filesize      = 0;
	my $bootstrap_num = $values_ref->{bootstrap_num};

	# Make sure $qsub exists. Otherwise we might run into an undef if the file hasn't been written to yet.
	unless ( -f $qsub ) {
		open( TMP, ">>$qsub" );
		close TMP;
	}

	if ( $values_ref->{bootstrap} ) {

		if ( -f "$values_ref->{output_dir}$values_ref->{job_name}_bootstrap/bootstrap.conf" ) {
			my @config_filenames;
			push @config_filenames, "$values_ref->{output_dir}/$values_ref->{job_name}_bootstrap/bootstrap.conf";

			# Create our $values_ref from the config file(s)
			my $values_ref = parse_file(@config_filenames);

			# Create defaults and modify values where appropriate.
			update_values($values_ref);

			$bootstrap_num = $values_ref->{bootstrap_num};
		}
	}

	while (1) {
		my $outdir = "$values_ref->{output_dir}$values_ref->{job_name}";
		$outdir .= "_bootstrap"
		  if ( $values_ref->{bootstrap} );

		opendir GEN_DIR, $outdir
		  or next;

		my @output_directories = readdir(GEN_DIR);
		closedir(GEN_DIR);

		for (@output_directories) {
			if ( $values_ref->{bootstrap} ) {
				if ( $_ eq "Results_$bootstrap_num" ) {
					$finished_flag = 1;
				}
			}
			else {
				if ( $_ eq "Results" ) {
					$finished_flag = 1;
				}
			}
		}

		sleep 2;

		$filesize = ( -s $qsub );
		if ( $filesize != $last_filesize ) {

			$last_filesize = -s $qsub;

			open JOB_OUTPUT, '<', $qsub;
			my @job_output = <JOB_OUTPUT>;

			for ( my $i = $last_line ; $i < @job_output ; $i++ ) {
				print $job_output[$i];
			}
			$last_line = @job_output;
		}

		last
		  if ( $finished_flag && !defined $values_ref->{bootstrap} );

		if ( $finished_flag && $values_ref->{bootstrap} ) {
			$finished_flag = 0;
			last
			  if ( $bootstrap_num > $values_ref->{bootstrap} );
		}
	}
	exit 0;
}

sub check_queue($$) {
	my ( $values_ref, $status ) = @_;

	my @flags = @$status;
	my $prefix = substr( $values_ref->{job_name}, 0, 4 );

	my $search_string = "qstat -u \$USER";
	if ( $values_ref->{cluster_software} eq "slurm" ) {
		$search_string = "squeue -u \$USER";
	}

	if (@flags) {
		$search_string .= " | grep";
	}
	for (@flags) {
		$search_string .= " -e ' $_ '";
	}
	$search_string .= " | wc -l";

	my $total_jobs = 0;
	$total_jobs = `$search_string`;

	chomp($total_jobs) if $total_jobs;

	return $total_jobs if $total_jobs;
	return 0;
}

sub submit_qsub_job($$) {
	my ( $values_ref, $command ) = @_;
	my @flags = qw(R Q qw r t PD);

	my $job         = 0;
	my $retry_limit = 5;

	show_output "Attempting to submit cluster job with command: $command\n"
	  if ( $verbosity >= 3 );

	for ( my $i = 0 ; $i < $retry_limit ; $i++ ) {
		if ( check_queue( $values_ref, \@flags ) <= $values_ref->{job_limit} ) {
			my $open_status = 1;
			open( QSUB, "$command |" )
			  or $open_status = 0;

			$job = <QSUB> if $open_status;
			close QSUB if $open_status;
		}
		else {
			show_output "die", "Error: You have too many cluster jobs running to start more. Please adjust your .conf file to lower the number of qsub jobs used, or change 'job_limit' to a higher value.\n";
		}

		last unless ( !defined $job || $job =~ /pbs_iff/ );
		sleep 5
		  unless $job;
	}
	chomp $job if $job;

	show_output "Successfully ran qsub job: $job\n"
	  if ( $verbosity >= 3 && $job );

	show_output "Couldn't run qsub job...\n"
	  if ( $verbosity >= 3 && !$job );

	return $job;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
