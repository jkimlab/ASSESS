all:
	cd lib/kent/src/lib && ${MAKE}
	cd util/samtools-1.5 && ./configure --without-curses --disable-bz2 --disable-lzma && ${MAKE}
	cd util/bedtools2 && ${MAKE}
	cd code/makeBlocks && ${MAKE}
	cd code && ${MAKE}

clean:
	cd lib/kent/src/lib && ${MAKE} clean
	cd util/samtools-1.5 && ${MAKE} clean
	cd util/bedtools2 && ${MAKE} clean
	cd code/makeBlocks && ${MAKE} clean
	cd code && ${MAKE} clean
