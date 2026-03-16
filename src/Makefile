#
CC=gcc
COPTS= -O2 -Wall


mesh:	mesh.o main_mesh.o eigen.o libmesh6.o  
	$(CC) $(COPTS) -o mesh mesh.o  main_mesh.o eigen.o libmesh6.o -lm

libmesh6.o :	libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS)  -I. libmesh6.c  

mesh.o :	mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h 
	$(CC) -c $(COPTS)  -I. eigen.c

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c  

clean :
	-rm mesh mesh.o main_mesh.o eigen.o libmesh6.o 
