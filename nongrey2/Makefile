#* Makefile for CEPAM -Code d'Evolution Planetaire Adaptatif et Modulaire
#* (c) T. Guillot 2004
#* make : compiles all programs and creates bin/cepam
#* make test : run the code and compare output with the Test.csv file
#* The compilation uses either g95, gfortan, xlf, or g77, or f77, in that order
#* In order to force the use of a given compiler, use (e.g.):
#*	make -e force "FOR=/sw/bin/g77"

#------------------------------

# PARAMETERS (Can be redefined)
#REP = Cepam
#SAVED = ${REP}/Makefile ${REP}/README ${REP}/MODIFS ${REP}/DATE ${REP}/src ${REP}/expls ${REP}/graph ${REP}/input ${REP}/bin ${REP}/data
FEXEC = NonGrey
FFLAGS = 
#FFLAGS = -O2 -I. -I./src -I./src/Communs -module ./src -parallel
SHELL = /bin/bash
#------------------------------
OBJS := $(shell echo `ls -m *.f90`| sed 's/\.f/\.o/g' | sed 's/\,//g')

#INCLUDES = src/cepam.inc
FOR := $(shell if which ifort; then :;\
	elif which gfortran; then :;\
	elif which g95; then :;\
	elif which xlf; then :;\
	elif which f77; then :;\
	elif which g77; then :;\
	else echo 'I couldn''t find gfortran, ifort, g95, xlf, g77, f77'; exit 11; fi)
THISDIR := $(shell pwd)
#DATAS = ./test_expls/evolutions/CoRoT-Exo-3_std/data ./test_expls/evolutions/HD149026_std/data ./test_expls/evolutions/HD209458_diss/data 
TESTFILE:=./PTprofile-test.csv

FEXEC : *.f90 
	${FOR} *.f90 -o ${FEXEC}

test  : ${FEXEC}
	@./${FEXEC}
	@ echo "Compare the begining and end of the downloaded file Test.csv and the newly created file PTprofile.csv. This comparison make sense only if you kept the same values of the coefficients as in the downloaded code."
	@ echo ""
	@ echo "For the head"
	@ head -2 Test.csv PTprofile.csv
	@ echo ""
	@ echo "Now for the end of the file"
	@ echo ""
	@ tail -2  Test.csv PTprofile.csv

force : 
	/bin/rm -f *.o ${FEXEC}
	 ${MAKE}
	
clean :
	/bin/rm -f *.o
