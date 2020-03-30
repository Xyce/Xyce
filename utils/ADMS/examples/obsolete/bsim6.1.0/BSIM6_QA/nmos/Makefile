
#
#	Example Makefile to run tests and check results.
#
#	This is an example makefile for running QA tests on a
#	model and then checking the simulated results against
#	reference results. A separate target is defined for each
#	variant of the model. The program runQaTests.pl runs the
#	tests, and that program expects a perl module SIMULATOR.pm
#	to be provided for each simulator that is tested.
#	Examples of these are provided.
#

qaSpecFile		=	qaSpec
qaResultsDirectory	=	results
testProgramName		=	../../../../lib/runQaTests.pl
#testProgramName		=	runQaTests.pl
testProgramFlags        =       -nw
#testProgramFlags        =       -d

help:
	@echo "" ; \
	echo "Valid targets are:" ; \
	echo "" ; \
	echo "all                run tests and compare results for all simulators" ; \
	echo "" ; \
	echo "spice              run tests and compare results spice" ; \
	echo "hspice             run tests and compare results hspice" ; \
	echo "spectre            run tests and compare results spectre" ; \
	echo "ads                run tests and compare results ads" ; \
	echo "eldo               run tests and compare results eldo" ; \
	echo "nexxim             run tests and compare results nexxim" ; \
	echo "smartspice         run tests and compare results smartspice" ; \
	echo "" ; \
	echo "clean              remove all previously generated simulation results"; \
	echo "" ; \
	echo "NOTE: if test results exist they are not resimulated" ; \
	echo "NOTE: to force resimulation run \"make clean\" first" ; \
	echo ""

all:	spice hspice spectre ads eldo nexxim smartspice

#####
##### spice tests
#####

spice:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s spice $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for spice"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s spice $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (spice): $$test" ; \
	    for variant in `$(testProgramName) -lv -s spice $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s spice -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/spice` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/spice/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for spice"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s spice $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (spice): $$test" ; \
		    for variant in `$(testProgramName) -lv -s spice $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s spice -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### hspice tests
#####

hspice:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s hspice $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for hspice"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s hspice $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (hspice): $$test" ; \
	    for variant in `$(testProgramName) -lv -s hspice $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s hspice -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/hspice` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/hspice/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for hspice"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s hspice $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (hspice): $$test" ; \
		    for variant in `$(testProgramName) -lv -s hspice $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s hspice -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### spectre tests
#####

spectre:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s spectre $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for spectre"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s spectre $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (spectre): $$test" ; \
	    for variant in `$(testProgramName) -lv -s spectre $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s spectre -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/spectre` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/spectre/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for spectre"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s spectre $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (spectre): $$test" ; \
		    for variant in `$(testProgramName) -lv -s spectre $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s spectre -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### ads tests
#####

ads:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s ads $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for ads"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s ads $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (ads): $$test" ; \
	    for variant in `$(testProgramName) -lv -s ads $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s ads -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/ads` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/ads/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for ads"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s ads $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (ads): $$test" ; \
		    for variant in `$(testProgramName) -lv -s ads $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s ads -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### eldo tests
#####

eldo:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s eldo $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for eldo"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s eldo $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (eldo): $$test" ; \
	    for variant in `$(testProgramName) -lv -s eldo $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s eldo -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/eldo` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/eldo/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for eldo"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s eldo $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (eldo): $$test" ; \
		    for variant in `$(testProgramName) -lv -s eldo $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s eldo -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### nexxim tests
#####

nexxim:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s nexxim $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for nexxim"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s nexxim $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (nexxim): $$test" ; \
	    for variant in `$(testProgramName) -lv -s nexxim $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s nexxim -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/nexxim` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/nexxim/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for nexxim"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s nexxim $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (nexxim): $$test" ; \
		    for variant in `$(testProgramName) -lv -s nexxim $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s nexxim -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done

#####
##### smartspice tests
#####

smartspice:
	@-echo ""; \
	localPlatform=`$(testProgramName) -platform` ; \
	localVersion=`$(testProgramName) -sv -s smartspice $(qaSpecFile)` ; \
	localVersionAndPlatform=$$localVersion._.$$localPlatform ; \
	echo "******"; \
	echo "****** $(qaSpecFile) tests for smartspice"; \
	echo "****** (for version $$localVersion on platform $$localPlatform)"; \
	echo "******"; \
	for test in `$(testProgramName) -lt -s smartspice $(qaSpecFile)` ; \
	do \
	    echo ""; \
	    echo "****** Checking test (smartspice): $$test" ; \
	    for variant in `$(testProgramName) -lv -s smartspice $(qaSpecFile)` ; \
	    do \
		$(testProgramName) $(testProgramFlags) -s smartspice -r -t $$test -var $$variant $(qaSpecFile) ; \
	    done ; \
	done ; \
	for version in `ls -C1 $(qaResultsDirectory)/smartspice` ; \
	do \
	    for platform in `ls -C1 $(qaResultsDirectory)/smartspice/$$version` ; \
	    do \
		versionAndPlatform=$$version._.$$platform ; \
		if [ $$versionAndPlatform = $$localVersionAndPlatform ] ; \
		then \
		    break ; \
		fi ; \
		echo "" ; \
		echo "******"; \
		echo "****** Comparing previously run $(qaSpecFile) tests for smartspice"; \
		echo "****** (for version $$version on platform $$platform)"; \
		echo "******"; \
		for test in `$(testProgramName) -lt -s smartspice $(qaSpecFile)` ; \
		do \
		    echo ""; \
		    echo "****** Checking test (smartspice): $$test" ; \
		    for variant in `$(testProgramName) -lv -s smartspice $(qaSpecFile)` ; \
		    do \
			$(testProgramName) -c $$version $$platform -s smartspice -t $$test -var $$variant $(qaSpecFile) ; \
		    done ; \
		done ; \
	    done ; \
	done


clean:
	@/bin/rm -rf $(qaResultsDirectory)/spice spiceCkt* b3v3check.log
	@/bin/rm -rf $(qaResultsDirectory)/hspice hspiceCkt*
	@/bin/rm -rf $(qaResultsDirectory)/spectre spectreCkt* *ahdlcmi
	@/bin/rm -rf $(qaResultsDirectory)/ads adsCkt* spectra.raw
	@/bin/rm -rf $(qaResultsDirectory)/eldo eldoCkt*
	@/bin/rm -rf $(qaResultsDirectory)/nexxim *.sdf *.csv *.csv.*
	@/bin/rm -rf $(qaResultsDirectory)/smartspice smartspiceCkt*
