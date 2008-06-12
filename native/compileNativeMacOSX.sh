#cc -c -O3 -funroll-loops -I/Library/Java/Home/include NucleotideLikelihoodCore.c
#cc -bundle -o libNucleotideLikelihoodCore.jnilib -framework JavaVM NucleotideLikelihoodCore.o


# This now creates a 'fat' dylib which runs on intel & PPC macs:
cc -c -arch ppc -O2 -funroll-loops -I/Library/Java/Home/include -o NucleotideLikelihoodCore.PPC.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.PPC.jnilib -framework JavaVM -arch ppc -dynamiclib NucleotideLikelihoodCore.PPC.o

cc -c -arch i386 -O2 -funroll-loops -I/Library/Java/Home/include -o NucleotideLikelihoodCore.i386.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib NucleotideLikelihoodCore.i386.o

lipo -create libNucleotideLikelihoodCore.PPC.jnilib libNucleotideLikelihoodCore.i386.jnilib -output libNucleotideLikelihoodCore.jnilib

##AminoAcidLikelihoodCore
cc -c -arch ppc -O2 -funroll-loops -I/Library/Java/Home/include -o AminoAcidLikelihoodCore.PPC.o AminoAcidLikelihoodCore.c
cc -o libAminoAcidLikelihoodCore.PPC.jnilib -framework JavaVM -arch ppc -dynamiclib AminoAcidLikelihoodCore.PPC.o

cc -c -arch i386 -O2 -funroll-loops -I/Library/Java/Home/include -o AminoAcidLikelihoodCore.i386.o AminoAcidLikelihoodCore.c
cc -o libAminoAcidLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib AminoAcidLikelihoodCore.i386.o

lipo -create libAminoAcidLikelihoodCore.PPC.jnilib libAminoAcidLikelihoodCore.i386.jnilib -output libAminoAcidLikelihoodCore.jnilib
