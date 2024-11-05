import matplotlib.pyplot as plt
import numpy as np

input_file = "28000.txt"
output_file = "28000.jpg"

processes = []
times = []

with open(input_file, 'r') as file:
    for line in file:
        if line == "\n":
            break

        time, proc = map(float, line.split())
        times.append(time)
        processes.append(int(proc))

initial_time = times[0]
speedup = [initial_time / t for t in times]

plt.figure(figsize=(10, 6))
plt.plot(processes, speedup, marker='o', linestyle='-', color='b', label="Ускорение")
plt.plot(range(2,30),range(2,30),'-',c="blue",linewidth=0.5,label="Linear speedup")

plt.xlabel('Количество процессов')
plt.ylabel('Ускорение')
plt.title('График строгой масштабируемости')

plt.grid(True)
plt.legend()

plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
