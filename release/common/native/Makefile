#
# command-line: make STATE_COUNT=<#> mac/linux
#

STATE_COUNT = 4

INNAME 		:= NewNativeLikelihoodCore.c
OUTNAME		:= NativeLikelihoodCore-$(STATE_COUNT)

############

MAC_JAVA_HOME	:= /Library/Java/Home
LINUX_JAVA_HOME	:= /usr/java/latest

MAC_INCLUDES	+= -I$(MAC_JAVA_HOME)/include
LINUX_INCLUDES	+= -I$(LINUX_JAVA_HOME)/include -I$(LINUX_JAVA_HOME)/include/linux

ARCH		:= native

MAC_LINK	:= bundle  # Can also be 'dynamiclib'
LINUX_LINK	:= shared

OPTIONS		:= -funroll-loops -ffast-math -fstrict-aliasing

############################################################


mac :
	cc -c -arch $(ARCH) -O3 -fast $(OPTIONS) $(MAC_INCLUDES) -std=c99 \
	   -DSTATE_COUNT=$(STATE_COUNT) \
	   -o $(OUTNAME).$(ARCH).o $(INNAME)
	cc -o lib$(OUTNAME).$(ARCH).jnilib -framework JavaVM -arch $(ARCH) \
	   -$(MAC_LINK) $(OUTNAME).$(ARCH).o
	lipo -create lib$(OUTNAME).$(ARCH).jnilib \
	     -output lib$(OUTNAME).jnilib

linux :
	gcc -c -O4 $(OPTIONS) $(LINUX_INCLUDES) -c $(INNAME) -std=c99 -DSTATE_COUNT=$(STATE_COUNT)  -o lib$(OUTNAME).o
	ld -$(LINUX_LINK) -o lib$(OUTNAME).so lib$(OUTNAME).o

