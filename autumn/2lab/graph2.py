import matplotlib.pyplot as plt

# Считываем имя входного и выходного файлов от пользователя
# input_file = input("Введите имя файла с данными: ")
# output_file = input("Введите имя выходного файла для сохранения графика: ")

input_file = "monte_result.txt"
output_file = "graph2.jpg"


# Чтение данных из файла
processes = []
times1 = []
times2 = []

with open(input_file, 'r') as file:
    for i, line in enumerate(file):
        time, proc = map(float, line.split())
        if i < 4:  # Первые 4 строки для первого графика
            times1.append(time)
            processes.append(int(proc))  # Процессы достаточно сохранить только один раз
        else:  # Последующие строки для второго графика
            times2.append(time)

# Вычисление ускорения для обоих наборов данных
initial_time1 = times1[0]
initial_time2 = times2[0]

speedup1 = [initial_time1 / t * 4 for t in times1]
speedup2 = [initial_time2 / t * 4 for t in times2]

# Построение графиков
plt.figure(figsize=(10, 6))

# Первый график
plt.plot(processes, speedup1, marker='o', linestyle='-', color='b', label="10^7")

# Второй график
plt.plot(processes, speedup2, marker='s', linestyle='--', color='r', label="10^8")

plt.plot(range(2,30),range(2,30),'-',c="blue",linewidth=0.5,label="Linear speedup")


# Подписи осей
plt.xlabel('Количество процессов')
plt.ylabel('Ускорение')
plt.title('Графики строгой масштабируемости')

# Добавление сетки и легенды
plt.grid(True)
plt.legend()

# Сохранение графика в файл
plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
