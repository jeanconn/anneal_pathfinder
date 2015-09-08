#! /usr/bin/env perl
use warnings;

##***************************************************************************
##  Inspect Aspect camera telemetry, possibly in real time
##***************************************************************************

use Getopt::Long;
use FileHandle;
use Time::HiRes qw ( time sleep );
use Bit::Vector;
use Tlm;
use Tk;
use Tk::Table;
use Chandra::Time;

get_options();

##***************************************************************************
## Initialization
##***************************************************************************
#
$aca_packet_vec = Bit::Vector->new(8*224);
$img_data_vec   = Bit::Vector->new((27*4+5)*8);

for $slot (0..7) {
    for $i (0..120) {
	$img_data[$slot][$i] = 0.0;
    }
}

@img_data_length = (32, 0, 5+2*27, 0, 0, 0, 0, 5+4*27);

$synfail = $chkfail = $calfail = $highbgd = $ramfail = $romfail = $pwrfail = 0;

$PERSIST = 1;
$PERSIST = 1;
$CURSOR_HOME  = ($opt_page ? `tput home` : '');
$CLEAR_SCREEN = ($opt_page ? `tput clear`: '');
$slot68 = 0;   # this is used to identify a slot having 6x6 or 8x8 readouts
$MIN_UPDATE_TIME = 0.5;
$last_update = 0;
$finished_skip = 0;      # only skip forward to a starting VCDU count once!

$time0 = time;

#  Get location of ACA pixel data within tlm stream
#
for $slot (0..7) {
    for $i (0..7) {
	for $j (0..7) {
	    $img_array[$slot][$i][$j] = 0.0;
	}
    }
}
@img_data = ('','','','','','','','');
@img_total = (0,0,0,0,0,0,0,0);
@img_min   = (0,0,0,0,0,0,0,0);
@img_max   = (10,10,10,10,10,10,10,10);
get_aca_tlm_pos();
get_aca_decom_info();
get_img_array_map();
$zoom = 10;

Tlm::open_telem (file_name => $opt_i,
		 read_size => 10290,
		 quiet     => 0,
		);

$aca_packet = '';		# 224 byte packet of ACA data

set_limits();
make_gui();			# Initialize graphics screen
open_bin_output()  if (defined $opt_out); # Open binary output file if needed
%WROTE_HEADER = {} if (defined $opt_out);

print $CLEAR_SCREEN;
$top->after (100, \&process_telem);

MainLoop();

##***************************************************************************
sub process_telem {
# 
# Main processing loop, which gets frames of telemetry, processes them,
# and updates the event log as necessary
#
##***************************************************************************

  # Get telemetry

  if ($b_pause->cget('-text') eq "Pause") {
#    &Tlm::skip_frames ($opt_vcdu - $Tlm::vcdu) if ($Tlm::vcdu < $opt_vcdu);
    if (($Tlm::vcdu < $opt_vcdu) && !$finished_skip) {
	&Tlm::skip_frames ($opt_vcdu - $Tlm::vcdu);
    };
    $finished_skip = 1;
    &Tlm::get_minor_frame() ?
	process_frame ($b_fast->cget('-text') eq "Normal") : end_splat();
  }

  end_splat() if ($b_pause->cget('-text') eq 'Done');


#  update_event_log();		# Update event log

  if (time - $last_update > $MIN_UPDATE_TIME) {
    $top->update();
    $last_update = time;
  }
  $top->after ($update_delay, \&process_telem);
}

##***************************************************************************
sub set_limits {
##***************************************************************************
  $limits{red}{'Syntax err'} = [-1, 1, 'Status flag'];
  $limits{red}{'Cal fail'} = [-1, 1, 'Status flag'];
  $limits{red}{'RAM fail'} = [-1, 1, 'Status flag'];
  $limits{red}{'ROM fail'} = [-1, 1, 'Status flag'];
  $limits{red}{'Power fail'} = [-1, 1, 'Status flag'];
  $limits{red}{'High bgd'} = [-1, 1, 'Status flag'];
  $limits{red}{'Checksum err'} = [-1, 1, 'Status flag'];
  $limits{red}{'Cmd count'} = [-1, 1, 'Command'];
  $limits{red}{'Cmd progress'} = [-1, 1, 'Command'];
#  $limits{red}{''} = [,];
  foreach (keys %{$limits{red}}) {
      $ {$limits{lastval}{$_}} = 0;
  }
}

##***************************************************************************
sub check_limits {
##***************************************************************************
    foreach (keys %{$limits{red}}) {
	if (($ {$limits{var}{$_}}  <= $limits{red}{$_}[0]
	     || $ {$limits{var}{$_}}  >= $limits{red}{$_}[1])) {
	    if ($ {$limits{lastval}{$_}} != $ {$limits{var}{$_}}) {
		$limits{widget}{$_}->configure(-background => 'red');
		update_event_log($limits{red}{$_}[2], $_, $ {$limits{var}{$_}});
	    }
	} else {
	    $limits{widget}{$_}->configure(-background => 'green');
	}
	$ {$limits{lastval}{$_}} = $ {$limits{var}{$_}};
    }
}

##***************************************************************************
sub update_event_log {
##***************************************************************************
#  $event_log->insert('end', "$vcdu: hello \n");
  my ($message, $var, $value) = @_;
  $event_log->insert('end', "$tlm_date $tlm_gmt :: $message");
  $event_log->insert('end', " : $var = $value") if (defined $var);
  $event_log->insert('end', "\n");

  $event_log->see('end');
  if ($event_log->index('end') > $MAX_EVENT_LOG_LINES) {
    $event_log->delete("1.0 linestart", "1.0 lineend + 1 c");
  }
}


##***************************************************************************
sub make_gui {
## Graphics initialization
##***************************************************************************
  $top = MainWindow->new();
  
  # Times and images

  $times_and_images = $top->Frame();
  $times_and_images->pack();

  @times = (['VCDU',   \$vcdu],
	 ['Format', \$format],
	 ['Date',   \$tlm_date],
	 ['GMT',   \$tlm_gmt],
	 );
  $times = $times_and_images->Table(-columns => 2,
				    -rows    => $#ts+1,
				    -scrollbars => 'o',
				   );
  $times->pack(-expand=> 1, -fill => 'both', -side => 'left');
  for $i (0 .. $#times) {
      my $l = $times->Label(-text => $times[$i][0], -relief => 'groove');
      $times->put($i, 0, $l);
      $l = $times->Label(-textvariable => $times[$i][1], -relief => 'groove');
      $times->put($i, 1, $l);
  }

  # Pixel images of 8 slots, at top

  $slot_img  = $times_and_images->Photo( 'slot_img' , -palette => 140);
  $all_img  = $times_and_images->Photo( 'all_img' , -palette => 140);
  $all_img->put (("#ffffff"), -to => (0, 0, 8*8*$zoom+9, 8*$zoom));
  $times_and_images->Label('-image'=> $all_img  )->pack;
  
  # Frame with all the other information and buttons

  $info_and_buttons = $top->Frame();
  $info_and_buttons->pack();

  @ts = (
	 ['Integ time', \$integ_time],
	 ['Cmd count', \$cmd_cnt],
	 ['Cmd progress', \$cmd_prog],
	 ['Syntax err', \$synfail],
	 ['Checksum err', \$chkfail],
	 ['Cal fail', \$calfail],
	 ['RAM fail', \$ramfail],
	 ['ROM fail', \$romfail],
	 ['Power fail', \$pwrfail],
	 ['High bgd', \$highbgd],
	);
	 
  $times_and_status = $info_and_buttons->Table(-columns => 2,
					       -rows    => $#ts+1,
					       -scrollbars => 'o',
					       );
  $times_and_status->pack(-expand=> 1, -fill => 'both', -side => 'left');
  for $i (0 .. $#ts) {
    my $l = $times_and_status->Label(-text => $ts[$i][0], -relief => 'groove');
    $times_and_status->put($i, 0, $l);
    $l = $times_and_status->Label(-textvariable => $ts[$i][1], -relief => 'groove');
    $times_and_status->put($i, 1, $l);
    if (exists $limits{red}{$ts[$i][0]}) {
      $limits{var}{$ts[$i][0]} = $ts[$i][1]; 
      $limits{widget}{$ts[$i][0]} = $l;
    }
  }

  # ACA info table for each of eight slots

  @info_msid   = qw(IMGNUM1 IMGFID1 IMGFUNC1 IMGSTAT IMGTYPE IMGROW0
		    IMGCOL0 IMGSCALE BGDAVG BGDSTAT);
  @info_format = qw(%i %i %i %i %i %5i %5i %6.1f  %6.1f %i);
  
  $info_table  = $info_and_buttons->Table(-columns => 9,
					  -rows => $#info_msid+1,
					  -scrollbars => 'o',
					 );
  $info_table->pack(-expand=> 1, -fill => 'both', -side => 'left');
  
  for $j (0..$#info_msid) {
    for $i (0..8) {
      $info_table_txt[$j][$i] = ($i == 0) ? $info_msid[$j] : '      ';
      my $l = $info_table->Label(-textvariable => \$info_table_txt[$j][$i],
				 -relief => 'groove');
      $info_table->put($j, $i, $l);
    }
  }
  
  $b_pause = $info_and_buttons->Button(-text => $opt_pause ? 'Resume' : 'Pause',
				       -command => sub {
					 $b_pause->cget('-text') eq "Pause" ?
					   $b_pause->configure(-text => 'Resume') :
					     $b_pause->configure(-text => 'Pause');
#					 $info_and_buttons->update();
					 $top->update();
					 $last_update=time;
				       }
				      );
  $b_pause->pack;
  
  $b_fast = $info_and_buttons->Button(-text => 'Fast',
				       -command => sub {
					 $b_fast->cget('-text') eq "Fast" ?
					   $b_fast->configure(-text => 'Normal') :
					     $b_fast->configure(-text => 'Fast');
#					 $info_and_buttons->update();
					 $top->update();
					 $last_update=time;
				       }
				      );
  $b_fast->pack;
  
  $b_quit = $info_and_buttons->Button(-text => 'Quit',
				      -command => sub {
					$b_pause->configure(-text => 'Done');
				      });
  $b_quit->pack;
  
  #$show = $info_and_buttons->Entry (width   =>  20)->pack;
  #$show->bind('<KeyPress-Return>', sub {
  #	      $update_delay = $show->get();
  #	    });
  
  $update_delay = 0;
  $delay = $info_and_buttons->Scale(-orient => 'vertical',
				    -from   => 0,
				    -to     => 2000,
				    -tickinterval => 500,
				    -label  => 'Delay',
				    -variable => \$update_delay
		    );
  $delay->pack;
  
  # Event log

  $top->Label(-text => 'Event log')->pack();
  $event_log = $top->Scrolled('Text', -height => 12, -width => 100);
  $event_log->pack(-side => 'top', -anchor => 'w');

}


##***************************************************************************
sub end_splat {
#  Close telemetry file and output binary file if necessary
##***************************************************************************
    &Tlm::close_telem();
    
    if (defined $opt_out) {
	foreach (0..7) {
	    if (defined $DAT_fh[$_]) {
		close ($DAT_fh[$_]) or die "Couldn't close DAT file $_\n";
	    }
	}
    }

    print STDERR "SPLAT% Finished. \n";
    $top->withdraw();
    exit(0);
}

##***************************************************************************
sub open_bin_output {
##  Open output binary file if needed
##***************************************************************************
  eval "use FileHandle";
  
  @out_num = split(/,/ , $opt_out);
  foreach (@out_num) {
    $DAT_fh[$_] = new FileHandle "> acaslot$_.dat";
    die "Couldn't open acaslot$_.dat\n" unless defined $DAT_fh[$_];
    $DAT_fh[$_]->autoflush(1);
  }
}

##***************************************************************************
sub get_options {
##***************************************************************************
  $Usage = <<'USAGE';
  Usage: splat.pl [-h] [-i tlm_file] [-image] [-limits <file>] [-logsc] 
                  [-out <file>] [-page] [-pea (a|b)] [-raw] [-time <msec>] [-tctm]
     -h     : Print help
     -i     : Optional input telemetry file
     -image : Show images (default, use -noimage otherwise)
     -limits: Limits file (not implemented yet)
     -logsc : Display log-scale images instead of linear-scale images
     -out   : Write ACA data from slots to ASCII files (e.g. -out=1,2,5,6)
     -page  : Run in display page mode
     -pea   : PEA unit ('a' or 'b')
     -raw   : Show raw data 
     -time  : Minimum time (msec) between screen updates
     -vcdu  : VCDU count value for start of processing
     -pause : Start out paused (default=false)
     -scale : Apply pixel signal scaling: DN = (scale_factor/32 * pixsig - 50) 
  
Examples:

Note the use of splat instead of splat.pl.  This is required to grab the
necessary ACA EGSE perl libraries.

From solaris in ASCDS environment (to get correct GMT):

  /proj/telmon/aca/telem_tools/dump2occ.pl telem.sto.gz > telem.occ

From linux:

  splat -i telem.occ

  gunzip --stdout telem.sto.gz \
  | $SKA/data/aca_egse/telem_convert/ehs5tdsn.pl \
  | splat

USAGE
  ; 

  $ACA_DECOM_HOME = "$ENV{SKA}/data/aca_egse/decom/aca";
  $MAX_EVENT_LOG_LINES = 10000;
  $opt_page = 1;
  $opt_time = 0;
  $opt_vcdu = 0;
  $opt_raw  = 1;
  $opt_i    = undef;
  $opt_image= 1;
  $opt_h    = 0;
  $opt_tctm = 0;
  $opt_logsc= 1;
  $opt_pea  = 'a';
  $opt_pause= 0;
  $opt_scale  = 1;
  
  &GetOptions ('h',		# Help
	       'page!',		# Run in display page mode
	       'image!',	# Display images
	       'raw!',		# Display raw data
	       'i=s',		# Optional input file
	       'time=i',	# Min time between updates
	       'vcdu=i',	# Minimum VCDU counter value
	       'tctm!',		# TCTM format?
	       'logsc!',	# Display log-scale grey-scale images
	       'out=s',		# Write ACA data from slots to an ASCII file
	       'pea=s',		# PEA ("a" or "b")
               'pause!',        # Start out paused
               'scale!',        # Apply pixel signal scaling: DN = (scale_factor/32 * pixsig - 50)
	      );
  
  die $Usage if ($opt_h);
  
  $opt_time = $opt_time / 1000.0; # convert to seconds
}

##***************************************************************************
sub set_info_table {
##***************************************************************************
    my $slt;
    my $j;
    my $txt;
    for $slt (1 .. 8) {
	for $j (0 .. $#info_msid) {
	    if (defined $cal{$info_msid[$j]}[$slt-1]) {
		$txt = sprintf($info_format[$j], $cal{$info_msid[$j]}[$slt-1]);
		if ($txt ne $info_table_txt[$j][$slt]) {
		  $info_table_txt[$j][$slt] = $txt;
		}
	    }
	}
    }
}

##***************************************************************************
sub process_frame {
##***************************************************************************
#   Get vcdu at *beginning* of 4 ACA minor frames 
  my $fast = shift;

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime($tlm_time);
  $tlm_gmt = sprintf "%02d:%02d:%02d",$hour,$min,$sec;
  $tlm_date = sprintf "%4d:%03d", $year+1900,$yday+1;
 
  return if $fast;
   
    # check if min VCDU reached

#    if ($vcdu < $opt_vcdu) {
#	print $CURSOR_HOME, "VCDU = $vcdu\n";
#	return;
#    }

    if ($vcdu % 4 == 0) {	# should have a complete aca packet (224 bytes)
	if ($opt_time > 0) {
	    $dt = time - $time0;
	    sleep ($opt_time - $dt) if ($dt < $opt_time);
	    $time0 = time;
	}
	print $CURSOR_HOME;
        printf "Format:%2d  VCDU:%9d (%6d,%3.3d)  ",$format,$vcdu,$vcdu/128,$vcdu%128; 
        if ($opt_tctm) {
	    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday) = gmtime($utc1970);
            printf "UTC:%14.3f = %4d %3.3d %2.2d:%2.2d:%2.2d",
	    $utc1970,$year+1900,$yday+1,$hour,$min,$sec; 
	}
        print "\n\n";

	if (length $aca_packet == 224) {
	    $aca_packet_vec->Block_Store(scalar reverse $aca_packet);
	    decom_aca_packet();
	    if ($pea_on) {
	      check_limits();	# Check limits for various values
	      process_mem_dump();	# Look for mem dump data, log to file if there
	      $out = '';
	      print_raw() if ($opt_raw);
	      print_img_array() if ($opt_image);
	      print $out;
	    }
	}

	$aca_packet = '';
    }

    for $i (0 .. $#aca_tlm_pos) {
	$aca_packet .= substr $mf, $aca_tlm_pos[$i], $aca_tlm_len[$i];
    }
}

##***************************************************************************
sub print_img_array {
##***************************************************************************
    for $slot (0, 2, 4, 6) {
	if (defined $cal{IMGCOL0}[$slot]) { 
	    $col0 = $cal{IMGCOL0}[$slot];
	    $row0 = $cal{IMGROW0}[$slot];
	} else {$col0 = 0; $row0 = 0};
	$out .= join '',"-" x 11, "SLOT ",$slot;
	$out .= sprintf (" (%4.0f,%4.0f)", $col0,$row0);
	$out .= sprintf (" (%7.0f %3s)", $img_total[$slot], $opt_scale ? 'DN' : 'raw');
	$out .= join '', "-"x18, "SLOT ",$slot+1;
	if (defined $cal{IMGCOL0}[$slot+1]) { 
	    $col0 = $cal{IMGCOL0}[$slot+1];
	    $row0 = $cal{IMGROW0}[$slot+1];
	} else {$col0 = 0; $row0 = 0};
	$out .= sprintf (" (%4.0f,%4.0f)", $col0,$row0);
	$out .= sprintf (" (%7.0f %3s)", $img_total[$slot+1], $opt_scale ? 'DN' : 'raw');
	$out .= join '', "-"x8,"\n";
	for $i (7, 6, 5, 4, 3, 2, 1, 0) {
	    $out .= join '', "| ";
	    for $j (0 .. 7) {
		$out .= sprintf ("%5.0f ", $img_array[$slot][$i][$j]);
	    }
	    $out .= join '', "| ";
	    for $j (0 .. 7) {
		$out .= sprintf ("%5.0f ", $img_array[$slot+1][$i][$j]);
	    }
	    $out .= join '', "|\n";
	}
    }
    $out .= join '', "-" x 101, "\n";
}

##***************************************************************************
sub decom_aca_packet {
##***************************************************************************
    my $slot;
    my $img_type;
    my $offset;

    # Extract global status stuff
    $startbit = 0;
    $integ_time = extract ($aca_packet_vec, 16) * 0.016;

    $startbit = 16;
    $highbgd = extract ($aca_packet_vec, 1);
    $ramfail = extract ($aca_packet_vec, 1);
    $romfail = extract ($aca_packet_vec, 1);
    $pwrfail = extract ($aca_packet_vec, 1);
    $calfail = extract ($aca_packet_vec, 1);
    $chkfail = extract ($aca_packet_vec, 1);
    $startbit = 23;
    $synfail = extract ($aca_packet_vec, 1);
    $cmd_cnt = extract($aca_packet_vec, 6);
    extract ($aca_packet_vec, 2);
    $cmd_prog= extract ($aca_packet_vec, 6);

    # If the PEA is off, then clear each the image slot, redrawing if needed

    if (substr ($aca_packet, 3, 40) eq "\377"x40) {  # PEA off
        $pea_on = 0;
	for $slot (0..7) {
	    my $update_slot = 0;
	    for $i (0..7) {
		for $j (0..7) {
		    $update_slot = 1 if ($img_array[$slot][$i][$j] != 0.0);
		    $img_array[$slot][$i][$j] = 0.0;
		}
	    }
	    draw_slot($slot);
	}
	return;
    }

    $pea_on = 1;
    $startbit = 40;		# start of "image types" in aca_packet_bits
    for $slot (0 .. 7) {
#       next unless (extr_val ($aca_packet_vec, 68+$slot*216, 2));    # idle slot?
	$img_type = extract ($aca_packet_vec, 3);
	$cal{IMGTYPE}[$slot] = $img_type;

	# Put header data into $img_data[$slot] if needed
	if ($img_type == 0 || $img_type == 1 || $img_type == 4) {
	    $img_data[$slot] = substr $aca_packet, 0, 5;
	}

	# identify a 6x6 or 8x8 readout slot, for printing out temperatures
	if ($img_type == 2 || $img_type == 5) {$slot68 = $slot};

	# Copy data from the aca_packet (224 bytes) to img_data[slot]
	$img_data[$slot] .= substr $aca_packet, $slot*27+8, 27;

	# Decom data in $img_data[$slot] if needed
	if ($img_type == 0 || $img_type == 2 || $img_type == 7) {
	    $len = $img_data_length[$img_type];
	    $lenbuf = 27*4+5;	# length of $img_data_vec in bytes
	    if (length $img_data[$slot] < $lenbuf) {
		$img_data[$slot] .= "\000" x ($lenbuf - length $img_data[$slot]);
	    }

	    $img_data_vec->Block_Store(scalar reverse $img_data[$slot]);

	    # Extract the data

	    @{$fields[$slot]} = ();
	    for $i (0 .. $#d_field) {
		$img_data_vec_startbit = $d_byte[$i] * 8 + $d_bit[$i];
		last if ($img_data_vec_startbit + $d_len[$i] > $len * 8);
		push @{$fields[$slot]}, $d_field[$i];
		$raw{$d_field[$i]}[$slot] =
		    extr_val ($img_data_vec, $img_data_vec_startbit, $d_len[$i]);
	    }

	    apply_cal ($slot);
	    set_image_array ($slot, $img_type);
	    write_slot_header($img_type, $slot)
		if (defined $opt_out
                        && defined $DAT_fh[$slot]
                            && not defined $WROTE_HEADER[$slot]
                    );
	    write_slot_data($img_type, $slot)
		if (defined $opt_out && defined $DAT_fh[$slot]);
	}
    }
    $top->update();
}    


##***************************************************************************
sub write_slot_header {
##***************************************************************************
# write slot data to a file, for the specified slot
# if temperature etc. data is not available, put out dummy data = -99.
    my $img_type = shift;
    my $slot = shift;

    print {$DAT_fh[$slot]} "VCDU time DUTC SLOT IMGTYPE ";
    for $i (0 .. 11) {
	print {$DAT_fh[$slot]} $fields[$slot][$i] , " ";
    }
    print {$DAT_fh[$slot]} "BGDRMS TEMPCD TEMPHOUS TEMPPRIM TEMPSEC BGDSTAT ";
    for $i (0 .. 7) {
	for $j (0 .. 7) {
	    print {$DAT_fh[$slot]} "r${i}_c${j} ";
	}
    }
    print {$DAT_fh[$slot]} "\n";
    $WROTE_HEADER[$slot] = 1;
}


##***************************************************************************
sub write_slot_data {
##***************************************************************************
# write slot data to a file, for the specified slot
# if temperature etc. data is not available, put out dummy data = -99.
    my $img_type = shift;
    my $slot = shift;

    $dutc = $opt_tctm ? $utc1970 : -99;
    my $ctime = Chandra::Time->new(time, {format => 'unix'})->secs;
    print {$DAT_fh[$slot]} $vcdu," ",$ctime," ",$dutc," ",$slot," ",$img_type," ";
    for $i (0 .. 11) {
	print {$DAT_fh[$slot]} $cal{ $fields[$slot][$i] }[$slot], " ";
    }
    if (defined $fields[$slot][36]) {
	for $i (31 .. 36) {
	    print {$DAT_fh[$slot]} $cal{ $fields[$slot][$i] }[$slot], " ";
	}
    }
    else {
	print {$DAT_fh[$slot]} "-99 -99 -99 -99 -99 -99 ";
    } 
    for $i (0 .. 7) {
	for $j (0 .. 7) {
	    print {$DAT_fh[$slot]} $img_array[$slot][$i][$j]," ";
	}
    }
    print {$DAT_fh[$slot]} "\n";
}

##***************************************************************************
sub set_image_array {
##***************************************************************************
    # Set fixed size (8x8) image arrays with pixel data
    my ($slot, $img_type) = @_;
    $img_total[$slot] = 0.0;
    my $n_pix = 0;
    $img_max[$slot] = 10;
    $img_min[$slot] = 1000;
    for $i (0 .. 7) {
	for $j (0 .. 7) {
	    if (defined $img_array_map[$img_type][$i][$j]) {
		$img_array[$slot][$i][$j] = $cal{$img_array_map[$img_type][$i][$j]}[$slot];
		$img_array[$slot][$i][$j] = 99999.0 if ($img_array[$slot][$i][$j] > 99999.);
#		$img_total[$slot] += $img_array[$slot][$i][$j] - $cal{BGDAVG}[$slot];
		$n_pix += 1;
		$img_total[$slot] += $img_array[$slot][$i][$j];
	    } else {
		$img_array[$slot][$i][$j] = 0.0;
	    }
	    $img_max[$slot] = $img_array[$slot][$i][$j] 
		if ($img_array[$slot][$i][$j] > $img_max[$slot]);
	    $img_min[$slot] = $img_array[$slot][$i][$j] 
		if ($img_array[$slot][$i][$j] < $img_min[$slot])
	}
    }
    $cal{MYBGD}[$slot] = ($img_array[$slot][0][0] + $img_array[$slot][0][1]
			  + $img_array[$slot][0][6] + $img_array[$slot][0][7]
			  + $img_array[$slot][7][0] + $img_array[$slot][7][1]
			  + $img_array[$slot][7][6] + $img_array[$slot][7][7]) / 8.0;
    $img_total[$slot] -= $n_pix * $cal{MYBGD}[$slot];

    draw_slot($slot);
    set_info_table();
}
    

##***************************************************************************
sub get_img_array_map {
##***************************************************************************

#
#  Get image array mapping
#
    my $map_info = "$ACA_DECOM_HOME/img_array_map.dat";
    open (MAP, $map_info) or die "Couldn't open $map_info\n";
    while (<MAP>) {
	next if (/^;/ || /^\s+$/);
	my ($img_type, $i, $j, $field) = split;
	$img_array_map[$img_type][$i][$j] = "SIGPIX" . $field;
    }
    close MAP;
}

##***************************************************************************
sub apply_cal {
##***************************************************************************
    # Apply calibrations
    my $slot = shift;
    foreach (@{$fields[$slot]}) {
	if (/SIGPIX/) {  # BGDAVG and BGDRMS are not scaled
            $cal{$_}[$slot] =  $opt_scale ?
              $raw{$_}[$slot] * $raw{IMGSCALE}[$slot] / 32.0 - 50.0 : $raw{$_}[$slot];
	} elsif (/INTEG/) {
	    $cal{$_}[$slot] = $raw{$_}[$slot] * 0.016;
	} elsif (/TEMPCD/ || /TEMPHOUS/ || /TEMPPRIM/ || /TEMPSEC/) {
	    $cal{$_}[$slot] = $raw{$_}[$slot] * 0.4;
	    if ($raw{$_}[$slot] > 127) { $cal{$_}[$slot] = ($raw{$_}[$slot] - 256) * 0.4 };
	} elsif (/IMGCOL/ || /IMGROW/) {
	    $cal{$_}[$slot] = $raw{$_}[$slot];
	    if ($raw{$_}[$slot] > 511) { $cal{$_}[$slot] = $raw{$_}[$slot] - 1024 };
	} else {
	    $cal{$_}[$slot] = $raw{$_}[$slot];
	}
#	print "$_ = $cal{$_}[$slot]\n";
    }
}


##***************************************************************************
sub print_raw {
##***************************************************************************
    $out .= join '', "Header: ";
    for $i (0..7) {
      $out .= sprintf ("%02x ", unpack("C", substr $aca_packet, $i, 1));
    }
    
    $startbit = 0;
    $inttime = extract ($aca_packet_vec, 16, ' int_tim:%04x');
    
    $startbit = 20;
    extract ($aca_packet_vec, 1, ' cal:%01x');
    extract ($aca_packet_vec, 1, 'chk:%01x');
    $startbit = 23;
    extract ($aca_packet_vec, 1, 'syn:%01x');
    extract($aca_packet_vec, 6, 'cnt:%02x  ');
    extract ($aca_packet_vec, 2);
    $cmd_prog= extract ($aca_packet_vec, 6);
    
    $out .= join '', "\n";
    $out .= sprintf ( "INT_TIME=%8.3f sec  ", $inttime*0.016);
    if (defined $cal{TEMPCD}[$slot68]) {
	$out .= join '', "CCD_Temp= ",$cal{TEMPCD}[$slot68],"  ";
    }
    if (defined $cal{TEMPHOUS}[$slot68]) {
	$out .= join '', "AChous_Temp= ",$cal{TEMPHOUS}[$slot68],"  ";
    }

    if ($cmd_cnt > 0) {
	$out .= sprintf ("CMD_CNT=%2d  ", $cmd_cnt);
    } else {
	$out .= sprintf ("            ");
    }

    if ($cmd_prog > 0) {
	$out .= sprintf ("CMD_PROG=%2d  ", $cmd_prog);
    } else {
	$out .= sprintf ("             ");
    }
    
    $out .= "\n";
    
    # Raw bytes
    
    foreach $slot (0 .. 7) {
      $out .= join '', "Raw $slot: ";
      @aca_slot = unpack("C*", substr $aca_packet, $slot*27+8, 27);
      foreach (@aca_slot) {
	$out .= sprintf ("%02x ", $_);
      }
      $out .= join '', "\n";
    }

    $out .= "\n";
}


##***************************************************************************
sub process_mem_dump {
##***************************************************************************
  return if ($dump_file_open_failed);
  $dumping = 0;
  foreach $slot (0..7) {
    # Memory dump has IMGTYPE == 3
    next unless (defined $cal{IMGTYPE}[$slot] && $cal{IMGTYPE}[$slot] == 3); 
    $dumping = 1;

    # Is this the beginning of memory dump?
    if (! $dump_FH) {
      my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday) = gmtime($tlm_time);
      $sec = sprintf "%2.2d", $sec;
      $min = sprintf "%2.2d", $min;
      $hour = sprintf "%2.2d", $hour;
      $year = sprintf "%4.4d", $year+1900;
      $yday = sprintf "%3.3d", $yday+1;
      my $root_dir = "$ENV{SKA}/data/aca_egse";
      my $dump_file = "memdump.$year$yday.$hour$min$sec";
      unless ($dump_FH = new FileHandle "> $root_dir/$dump_file") {
	$dump_file_open_failed = 1;
	update_event_log("Could not open $dump_file in $root_dir");
	return;
      } 
      update_event_log("Memory dump file is $root_dir/$dump_file");
    }
    
    my $p_addr = unpack("S", substr $aca_packet, $slot*27+8+0,  2);
    my $prom = $p_addr & 0x8000;
    my $addr = ($p_addr & 0x7fff) << 1;
    my @data = unpack("S*", substr $aca_packet, $slot*27+8+2, 24);
    my $chksum = unpack("C*", substr $aca_packet, $slot*27+8+26, 1);
    
    printf $dump_FH "%08x ", $addr;
    foreach (@data) { printf $dump_FH "%04x ", $_ }
    print $dump_FH $prom ? "PROM " : "RAM  ";
    printf $dump_FH "%02x\n", $chksum;
  }

  if ($dump_FH && !$dumping) { # Dump has stopped, so close file
    close $dump_FH ? update_event_log("Memory dump file closed")
      : update_event_log("Memory dump file could not be closed");
    undef $dump_FH;
  }
}

##***************************************************************************
sub extract {
##***************************************************************************
    my ($vector, $len, $format) = @_;
    $val = $vector->Chunk_Read($len, $vector->Size()-$startbit - $len);
    $startbit += $len;

    if (defined $format) {
	$format .= " ";
	$out .= sprintf ($format, $val);
    }
    return $val;
}

##***************************************************************************
sub extr_val {
##***************************************************************************
    my ($vector, $startbit, $len) = @_;
    $val = $vector->Chunk_Read($len, $vector->Size() - $startbit - $len);
}


##***************************************************************************
sub get_aca_tlm_pos {
##***************************************************************************
#
##  Get location of ACA pixel data within tlm stream
#
    $stripform = "$ACA_DECOM_HOME/aids_pea_${opt_pea}.strip";
    @aca_tlm_pos = ();
    open (STRIP, $stripform) or die "Couldn't open $stripform\n";
    while (<STRIP>) {
	next if /^;/;
	next if /^\s+$/;
	@vals = split;
	next unless ($vals[0] == 2);
	push @aca_tlm_pos, $vals[1];
	push @aca_tlm_len, $vals[2];
    }
    close STRIP;    
}

##***************************************************************************
sub get_aca_decom_info {
##***************************************************************************
#
#  Get location of telemetry items in ACA image data
#
    $msid = 0;
    $decom_info = "$ACA_DECOM_HOME/aise_decom_info.dat";
    open (DECOM, $decom_info) or die "Couldn't open $decom_info\n";
    while (<DECOM>) {
	next if /^;/;
	next if /^\s+$/;
	($msid, $field, $byte, $bit, $len) = split;
	push @d_field, $field;
	push @d_byte, $byte;
	push @d_bit, $bit;
	push @d_len, $len;
    }
    close DECOM;
}

##***************************************************************************
sub draw_slot {
##***************************************************************************
    my ($slot) = @_;
    @carr = ([]);
    $pixmin = $img_min[$slot];
    if ($pixmin < 9) {$pixmin = 9};
    $pixmax = $img_max[$slot] - $pixmin;
    if ($opt_logsc) { $pixmin = log($pixmin); $pixmax = log($img_max[$slot]) -$pixmin };
    for $i (0 .. 7) {
	for $j (0 .. 7) {
	    $img_scale = $img_array[$slot][7-$i][$j];
	    if ($img_scale < 9) {$img_scale = 9};
	    if ($opt_logsc) { $img_scale = log($img_scale) };
	    $img_scale = int (($img_scale  - $pixmin) /$pixmax * 255);
	    $img_scale = 0 if ($img_scale < 0);
	    $img_scale = 255 if ($img_scale > 255);
	    $color = sprintf ("%02x", $img_scale);
	    push @{$carr[$i]}, "#$color$color$color";
	}
    }
    $slot_img->put(\@carr);
    $all_img->copy ($slot_img, -zoom => $zoom, -to => ($slot*8*$zoom+$slot+1, 0));
#    $top->update();
#    $last_update = time;
}
    
