all: ring broadcast gather alltoall

ring: ring.c
	mpicc -Wall -o ring ring.c

broadcast: broadcast.c
	mpicc -Wall -o broadcast broadcast.c

gather: gather.c
	mpicc -Wall -o gather gather.c

alltoall: alltoall.c
	mpicc -Wall -o alltoall alltoall.c

clean:
	rm -rf ring broadcast gather alltoall slurm*