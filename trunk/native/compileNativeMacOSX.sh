cc -c -O3 -funroll-loops -I/Library/Java/Home/include NucleotideLikelihoodCore.c
cc -bundle -o libNucleotideLikelihoodCore.jnilib -framework JavaVM NucleotideLikelihoodCore.o
