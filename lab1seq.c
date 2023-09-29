#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Функция для выделения памяти под трехмерный массив
double*** create3DArray(int xSize, int ySize, int zSize) {
    double*** arr = (double***)malloc(xSize * sizeof(double**));
    for (int i = 0; i < xSize; i++) {
        arr[i] = (double**)malloc(ySize * sizeof(double*));
        for (int j = 0; j < ySize; j++) {
            arr[i][j] = (double*)malloc(zSize * sizeof(double));
        }
    }
    return arr;
}

// Функция для освобождения памяти, выделенной под трехмерный массив
void free3DArray(double*** arr, int xSize, int ySize) {
    for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

// Функция для определения температуры в заданной точке (x, y, z)
double calculateTemperature(int x, int y, int z, int k, int m, int n, double*** temperature, int SphereCenterX, int SphereCenterY, int SphereCenterZ, double SphereRadius, double alphaInside, double alphaOutside) {
    double epsilon = 0.001; // Порог сходимости
    double maxDiff;
    double startTime = omp_get_wtime(); // Засекаем начальное время
    do {
        maxDiff = 0.0;
        for (int i = 1; i < k - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                for (int l = 1; l < n - 1; l++) {
                    double alpha = alphaOutside; // Коэффициент теплопроводности по умолчанию (вне сферы)
                    double distanceToSphereCenter = sqrt((i - SphereCenterX) * (i - SphereCenterX) +
                                                        (j - SphereCenterY) * (j - SphereCenterY) +
                                                        (l - SphereCenterZ) * (l - SphereCenterZ));
                    if (distanceToSphereCenter <= SphereRadius) {
                        alpha = alphaInside; // Коэффициент теплопроводности внутри сферы
                    }
                    double prevTemperature = temperature[i][j][l];
                    temperature[i][j][l] = temperature[i][j][l] + alpha * (temperature[i + 1][j][l] + temperature[i - 1][j][l]
                                                   + temperature[i][j + 1][l] + temperature[i][j - 1][l]
                                                   + temperature[i][j][l + 1] + temperature[i][j][l - 1]
                                                   - 6 * temperature[i][j][l]);
                    double diff = fabs(temperature[i][j][l] - prevTemperature);
                    if (diff > maxDiff) {
                        maxDiff = diff;
                    }
                }
            }
        }
    } while (maxDiff > epsilon);
    double endTime = omp_get_wtime(); // Засекаем конечное время
    printf("Время выполнения: %.2f секунд\n", endTime - startTime);
    return temperature[x][y][z];
}

int main(int argc, char *argv[]) {
    int x = atoi(argv[1]); // Координата точки для поиска по оси X
    int y = atoi(argv[2]); // Координата точки для поиска по оси Y
    int z = atoi(argv[3]); // Координата точки для поиска по оси Z
    int k = atoi(argv[4]); // Размер по оси X
    int m = atoi(argv[5]); // Размер по оси Y
    int n = atoi(argv[6]); // Размер по оси Z
    double T1 = atof(argv[7]); // Температура верхней грани
    double T2 = atof(argv[8]); // Температура сторон верхней грани
    double T3 = atof(argv[9]); // Температура начального состояния
    double alphaInside = atof(argv[10]); // Коэффициент теплопроводности внутри сферы
    double alphaOutside = atof(argv[11]); // Коэффициент теплопроводности вне сферы
    double SphereRadius = atof(argv[12]); // Радиус сферы
    int SphereCenterX = atoi(argv[13]); // Координата центра сферы по оси X
    int SphereCenterY = atoi(argv[14]); // Координата центра сферы по оси Y
    int SphereCenterZ = atoi(argv[15]); // Координата центра сферы по оси Z
    if (x < 0 || x >= k || y < 0 || y >= m || z < 0 || z >= n) {
        printf("Неверные координаты точки.\n");
    } else {
        double*** temperature = create3DArray(k, m, n);
        // Инициализация начального состояния
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < m; j++) {
                for (int l = 0; l < n; l++) {
                    temperature[i][j][l] = T3;
                }
            }
        }
        // Установка граничных условий
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < m; j++) {
                temperature[i][j][0] = 0; // Нижняя грань
                temperature[i][j][n - 1] = T1; // Верхняя грань
            }
        }
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < m; j++) {
                temperature[i][j][n - 2] = T2; // Стороны верхней грани
            }
        }
        double result = calculateTemperature(x, y, z, k, m, n, temperature, SphereCenterX, SphereCenterY, SphereCenterZ, SphereRadius, alphaInside, alphaOutside);
        printf("Температура в точке (%d, %d, %d) = %.2f\n", x, y, z, result);
        free3DArray(temperature, k, m);
    }
    return 0;
}