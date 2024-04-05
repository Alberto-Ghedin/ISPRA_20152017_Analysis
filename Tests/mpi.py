from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    c = 0
    t=time.time()
    for i in range(size):
        for _ in range(10**7):
            c+=1
    t=time.time()-t
    print(t)
    comm.barrier()
else:
    comm.barrier()

c=0
lista_indici=list(range(8))

t=time.time()
lista=[]
for indice in lista_indici[rank::size]:
    lista.append(indice)
t=time.time()-t
global_list=comm.gather(lista, root=0)
print(f'io sono il rank {rank} e ho impiegato {t} secondi')
if rank==0:
    print(f'Io sono il rank {rank} e ho raccolto questi dati: {global_list}')
    risultati=np.array(global_list)
    print(risultati)
    print(risultati.transpose().flatten())

   

