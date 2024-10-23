import matplotlib.pyplot as plt
import numpy as np

# Считываем имя входного и выходного файлов от пользователя
input_file = input("Введите имя файла с данными: ")
output_file = input("Введите имя выходного файла для сохранения графика: ")

# Чтение данных из файла
processes = []
times = []

with open(input_file, 'r') as file:
    for line in file:
        time, proc = map(float, line.split())
        times.append(time)
        processes.append(int(proc))

# Вычисление ускорения
initial_time = times[0]
speedup = [initial_time / t * 4 for t in times]

# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(processes, speedup, marker='o', linestyle='-', color='b', label="Ускорение")
plt.plot(range(2,30),range(2,30),'-',c="blue",linewidth=0.5,label="Linear speedup")


# Подписи осей
plt.xlabel('Количество процессов')
plt.ylabel('Ускорение')
plt.title('График строгой масштабируемости')

# Добавление сетки и легенды
plt.grid(True)
plt.legend()

# Сохранение графика в файл
plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
