#cc -c -O3 -funroll-loops -I/Library/Java/Home/include NucleotideLikelihoodCore.c
#cc -bundle -o libNucleotideLikelihoodCore.jnilib -framework JavaVM NucleotideLikelihoodCore.o


# This now creates a 'fat' dylib which runs on intel & PPC macs:
cc -c -arch ppc -O3 -fast -funroll-loops -I/Library/Java/Home/include -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers -o NucleotideLikelihoodCore.PPC.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.PPC.jnilib -framework JavaVM -arch ppc -dynamiclib NucleotideLikelihoodCore.PPC.o

cc -c -arch i386 -O3 -fast -funroll-loops -I/Library/Java/Home/include -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers -o NucleotideLikelihoodCore.i386.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib NucleotideLikelihoodCore.i386.o

cc -c -arch x86_64 -O3 -fast -funroll-loops -I/Library/Java/Home/include -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers -o NucleotideLikelihoodCore.x86_64.o NucleotideLikelihoodCore.c
cc -o libNucleotideLikelihoodCore.x86_64.jnilib -framework JavaVM -arch x86_64 -dynamiclib NucleotideLikelihoodCore.x86_64.o

lipo -create libNucleotideLikelihoodCore.PPC.jnilib libNucleotideLikelihoodCore.i386.jnilib libNucleotideLikelihoodCore.x86_64.jnilib -output libNucleotideLikelihoodCore.jnilib

##AminoAcidLikelihoodCore
#cc -c -arch ppc -O2 -funroll-loops -I/Library/Java/Home/include -o AminoAcidLikelihoodCore.PPC.o AminoAcidLikelihoodCore.c
#cc -o libAminoAcidLikelihoodCore.PPC.jnilib -framework JavaVM -arch ppc -dynamiclib AminoAcidLikelihoodCore.PPC.o

#cc -c -arch i386 -O2 -funroll-loops -I/Library/Java/Home/include -o AminoAcidLikelihoodCore.i386.o AminoAcidLikelihoodCore.c
#cc -o libAminoAcidLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib AminoAcidLikelihoodCore.i386.o

#lipo -create libAminoAcidLikelihoodCore.PPC.jnilib libAminoAcidLikelihoodCore.i386.jnilib -output libAminoAcidLikelihoodCore.jnilib

##GeneralLikelihoodCore
#cc -c -arch i386  -O2 -funroll-loops -I/Library/Java/Home/include \
#   -o GeneralLikelihoodCore.i386.o GeneralLikelihoodCore.c

#cc -o libGeneralLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 -dynamiclib GeneralLikelihoodCore.i386.o

#lipo -create libGeneralLikelihoodCore.i386.jnilib -output libGeneralLikelihoodCore.jnilib

#NativeMemoryLikelihoodCore -- Generates pointer cast warnings
#cc -c -arch i386  -O3 -fast -funroll-loops -I/Library/Java/Home/include \
#   -o NativeMemoryLikelihoodCore.i386.o NativeMemoryLikelihoodCore.c 

#cc -c -arch x86_64  -O3 -fast -funroll-loops -I/Library/Java/Home/include \
#   -o NativeMemoryLikelihoodCore.x86_64.o NativeMemoryLikelihoodCore.c

#cc -c -arch i386  -O2 -funroll-loops -I/Library/Java/Home/include \
#   -o NativeSubstitutionModel.i386.o NativeSubstitutionModel.c

#cc -o libNativeMemoryLikelihoodCore.i386.jnilib -framework JavaVM -arch i386 \
#    -dynamiclib NativeMemoryLikelihoodCore.i386.o NativeSubstitutionModel.i386.o


#cc -o libNativeSubstitutionModel.i386.jnilib -framework JavaVM -arch i386 -dynamiclib NativeSubstitutionModel.i386.o

#lipo -create libNativeMemoryLikelihoodCore.i386.jnilib -output libNativeMemoryLikelihoodCore.jnilib 
    




