#!/bin/sh 
set fp [open "~/Documents/CS4701/color.txt" r]
while {[gets $fp line] >= 0} {
	set slen [string length $line]
	if {$slen > 0} {
		set ldata [split $line " "]
		set fr [lindex $ldata 0]
		set molI [lindex $ldata 1]
		set cInter [split [lindex $ldata 2] ";"]
		set c [lindex $cInter 0]
		animate goto $fr
		set sel [atomselect 0 [concat "index " $molI]]
		$sel frame $fr
		$sel set user $c
		$sel delete
	}
}
mol modcolor 0 0 User 
mol colupdate 0 0 1 
mol scaleminmax 0 0 0.0 100
animate forward
