import matplotlib.pyplot as plt

input_file = "data.txt"
output_file = "graph.jpg"

processes = []
times1 = []
times2 = []

with open(input_file, 'r') as file:
    for i, line in enumerate(file):
        time, proc = map(float, line.split())
        if i < 4: 
            times1.append(time)
            processes.append(int(proc)) 
        else: 
            times2.append(time)

initial_time1 = times1[0]
initial_time2 = times2[0]

speedup1 = [initial_time1 / t * 4 for t in times1]
speedup2 = [initial_time2 / t * 4 for t in times2]

plt.figure(figsize=(10, 6))

plt.plot(processes, speedup1, marker='o', linestyle='-', color='b', label="28000")

plt.plot(processes, speedup2, marker='s', linestyle='--', color='r', label="45000")

plt.plot(range(1,32),range(1,32),'-',c="blue",linewidth=0.5,label="Linear speedup")


plt.xlabel('Number of process')
plt.ylabel('speedup')

plt.grid(True)
plt.legend()

plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
