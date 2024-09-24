all: ring broadcast gather

ring: ring.c
	mpicc -Wall -o ring ring.c

broadcast: broadcast.c
	mpicc -Wall -o broadcast broadcast.c

gather: gather.c
	mpicc -Wall -o gather gather.c

clean:
	rm -rf ring broadcast gather slurm*