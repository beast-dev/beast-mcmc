#cc -c -O3 -funroll-loops -I/Library/Java/Home/include NucleotideLikelihoodCore.c
#cc -bundle -o libNucleotideLikelihoodCore.jnilib -framework JavaVM NucleotideLikelihoodCore.o


# This now creates a 'fat' dylib which runs on intel & PPC macs:
cc -c -arch ppc -O3 -funroll-loops -mtune=G5 -I/Library/Java/Home/include -o NucleotideLikelihoodCore.PPC.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.PPC.jnilib -framework JavaVM -arch ppc -dynamiclib NucleotideLikelihoodCore.PPC.o

cc -c -arch i386 -O3 -funroll-loops -I/Library/Java/Home/include -o NucleotideLikelihoodCore.i386.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib NucleotideLikelihoodCore.PPC.o

lipo -create libNucleotideLikelihoodCore.PPC.jnilib libNucleotideLikelihoodCore.i386.jnilib -output libNucleotideLikelihoodCore.jnilib
