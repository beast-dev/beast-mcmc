gcc -O4 -march=core2 -funroll-loops -ffast-math -fPIC -fstrict-aliasing -c -I/usr/java/default/include/ -I/usr/java/default/include/linux NucleotideLikelihoodCore.c -o libNucleotideLikelihoodCore.o
ld -shared -o libNucleotideLikelihoodCore.so libNucleotideLikelihoodCore.o

gcc -O4 -march=pentiumpro -mcpu=pentiumpro -funroll-loops -ffast-math -fstrict-aliasing -c -I/usr/java/j2sdk1.4.1_01/include/ -I/usr/java/j2sdk1.4.1_01/include/linux NucleotideLikelihoodCore.c -o libNucleotideLikelihoodCore.o
ld -shared -o libNucleotideLikelihoodCore.so libNucleotideLikelihoodCore.o 

gcc -O4 -march=pentiumpro -mcpu=pentiumpro -funroll-loops -ffast-math -fstrict-aliasing -c -I/usr/java/j2sdk1.4.1_01/include/ -I/usr/java/j2sdk1.4.1_01/include/linux AminoAcidLikelihoodCore.c -o libAminoAcidLikelihoodCore.o
ld -shared -o libAminoAcidLikelihoodCore.so libAminoAcidLikelihoodCore.o 
