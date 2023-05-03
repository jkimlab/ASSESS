include inc/common.mk


all: topLibs 

alpha: clean utils-alpha cgi-alpha blatSuite

beta: check-beta clean cgi-beta

# do a git update and clean
update:
	git pull
	${MAKE} clean

topLibs:
#@./machTest.sh
	cd lib && ${MAKE}
#	cd jkOwnLib && ${MAKE}

hgLib:
	@./hg/sqlEnvTest.sh
	cd hg/lib && ${MAKE}

libs: topLibs hgLib
	cd hg/protein/lib && ${MAKE}

cgi: libs
	cd hg && ${MAKE} cgi

cgi-alpha: libs
	cd hg && ${MAKE} alpha

cgi-beta: check-beta libs
	cd hg && ${MAKE} beta

check-beta:
	# this will fail if we are not in a beta checkout:
	git branch | egrep 'v[0-9]+_branch' > /dev/null	

blatSuite: topLibs 
	cd blat && ${MAKE}
	cd gfClient && ${MAKE}
	cd gfServer && ${MAKE}
	cd webBlat && ${MAKE}
	cd hg/pslPretty && ${MAKE}
	cd hg/pslReps && ${MAKE}
	cd hg/pslSort && ${MAKE}
	cd utils/nibFrag && ${MAKE}
	cd utils/faToNib && ${MAKE}
	cd utils/faToTwoBit && ${MAKE}
	cd utils/twoBitToFa && ${MAKE}
	cd utils/twoBitInfo && ${MAKE}
	cd isPcr && ${MAKE}
	cd blatz && ${MAKE}

# all of these application makefiles have been updated to include use
#	of DESTDIR and BINDIR

UTILS_APPLIST = bedClip bedGraphToBigWig bedToBigBed \
        bigBedInfo bigBedSummary bigBedToBed bigWigAverageOverBed  \
	bigWigInfo bigWigSummary bigWigToBedGraph bigWigToWig\
	faCount faFrag faOneRecord faPolyASizes faRandomize faSize \
	faSomeRecords faToNib \
	faToTwoBit htmlCheck nibFrag paraFetch paraSync sizeof stringify \
	textHistogram twoBitInfo twoBitToFa userApps wigCorrelate wigToBigWig

HG_APPLIST = bedSort liftUp liftOver bedItemOverlapCount encode/validateFiles \
	 ratStuff/mafsInRegion ratStuff/mafSpeciesSubset \
	 pslCDnaFilter genePredToGtf pslPretty pslReps pslSort 


HG_UTILS_APPLIST = bedExtendRanges gff3ToGenePred gtfToGenePred \
	hubCheck overlapSelect makeTableList

HG_HGWIGGLE_APP = hgWiggle

BLAT_APPLIST = blat gfClient gfServer 

userApps: topLibs hgLib
	@for P in ${UTILS_APPLIST}; do \
	    ( cd utils/$${P} && echo utils/$${P} && ${MAKE} ) ; \
	done
	@for P in ${HG_APPLIST}; do \
	    ( cd hg/$${P} && echo hg/$${P} && ${MAKE} ) ; \
	done
	@for P in ${HG_UTILS_APPLIST}; do \
	    ( cd hg/utils/$${P} && echo hg/utils/$${P} && ${MAKE} ) ; \
	done
	@for P in ${HG_HGWIGGLE_APP}; do \
	    ( cd hg/makeDb/hgLoadWiggle && echo hg/makeDb/hgLoadWiggle && ${MAKE} $${P} ) ; \
	done
	@for P in ${BLAT_APPLIST}; do \
	    ( cd $${P} && echo $${P} && ${MAKE} ) ; \
	done

clean:
	cd lib && ${MAKE} clean 

testDirs = lib blat gfServer hg
test:: $(testDirs:%=%.test)
%.test:
	cd $* && ${MAKE} test

LIB_TAGS_IN = lib/*.[hc] hg/lib/*.[hc] jkOwnLib/*.c inc/*.h hg/inc/*.h

# build tags for libraries
.PHONY: tags
tags:
	ctags ${LIB_TAGS_IN}

# build emacs tags for libraries
.PHONY: etags
etags:
	etags ${LIB_TAGS_IN}

# build tags for all files
.PHONY: tags-all
tags-all:
	find . -name '*.[ch]' | ctags -L -

# build emacs tags for all files
.PHONY: etags-all
etags-all:
	find . -name '*.[ch]' | etags -

