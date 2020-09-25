All Fortran code is in the folder source. There is a Makefile, but it is specific to my desktop. I couldn't figure out how to compile this code on Mio or AuN.
I set up the makefile to put the executables in the folders 200*ex*.  Those folders all correspond to examples in the papers. Each contains a folder in and a 
folder out.  The folder in contains files with parameters that the codes read.  The usage of the codes is ./200*ex* file.txt, where file.txt contains parameters.
See the included files or the source code for usage.  Also in each example folder is a shell script that runs the code for all files in the in directory, and puts
the output in a file of the same name in the out directory.  There are also matlab scripts check.m that read output and make the tables I reported.

The third inversion method in 2D doesn't seem to work.  The other 2 2D methods seem to work but take a very long time.