#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#define EPS 10E-9

class S21Matrix {
 public:
  // CONSTRUCTORS:
  S21Matrix();  // Базовый конструктор, инициализирующий матрицу некоторой
                // заранее заданной размерностью
  S21Matrix(int rows, int cols);  // Параметризированный конструктор с
                                  // количеством строк и столбцов
  S21Matrix(const S21Matrix& other);  // Конструктор копирования
  S21Matrix(S21Matrix&& other);  // Конструктор переноса
  ~S21Matrix();                  // Деструктор

  // OPERATORS:
  double& operator()(
      int i, int j);  // Индексация по элементам матрицы (строка, колонка)
  S21Matrix& operator=(
      const S21Matrix& other);  // Присвоение матрице значений другой матрицы
  S21Matrix operator+(const S21Matrix& other);  // Сложение двух матриц
  S21Matrix operator-(
      const S21Matrix& other);  // Вычитание одной матрицы из другой
  S21Matrix operator*(const S21Matrix& other);  // Умножение матриц
  S21Matrix operator*(const double& num);  // Умножение матрицы на число
  bool operator==(
      const S21Matrix& other) const;  // Проверка на равенство матриц (EqMatrix)
  void operator+=(const S21Matrix& other);  // Присвоение сложения (SumMatrix)
  void operator-=(const S21Matrix& other);  // Присвоение разности (SubMatrix)
  void operator*=(const S21Matrix& other);  // Присвоение умножения (MulMatrix)
  void operator*=(const double& num);  // Присвоение умножения (MulNumber)

  // OPERATIONS:
  bool EqMatrix(const S21Matrix& other)
      const;  // Проверяет матрицы на равенство между собой
  void SumMatrix(
      const S21Matrix& other);  // Прибавляет вторую матрицы к текущей
  void SubMatrix(const S21Matrix& other);  // Вычитает из текущей матрицы другую
  void MulNumber(const double num);  // Умножает текущую матрицу на число
  void MulMatrix(const S21Matrix& other);  // Умножает текущую матрицу на вторую
  S21Matrix Transpose();  // Создает новую транспонированную матрицу из текущей
                          // и возвращает ее
  S21Matrix Minor(int iRow,
                  int jCol);  // Второстепенный метод позволяет применять
                              // матрицу всех размеров, а не только квадратную
  S21Matrix CalcComplements();  // Вычисляет матрицу алгебраических дополнений
                                // текущей матрицы и возвращает ее
  S21Matrix InverseMatrix();  // Вычисляет и возвращает обратную матрицу
  double Determinant();  // Вычисляет и возвращает определитель текущей матрицы

  // ADDITIONAL METHODS:
  void AlocMatrix(
      double*** matrix, int rows,
      int cols);  // Выделяет память под двумерный массив (матрицу) типа double
  void DelMatrix(double** matrix);  // Освобождает память
  void CheckIndexes(int i, int j);  // Проверяет индексы
  void SwapMaxRows(double** matrix, int maxRow,
                   int j);  // Меняет местами максимальное количество строк
  double GaussDet();  // Вычисляет определитель квадратной матрицы
  int CheckZero();  // Проверяет содержит ли матрица нули

  // ACCESSOR AND MUTATOR:
  int GetRows() const;  // Возвращает количество строк в матрице
  int GetCols() const;  // Возвращает количество столбцов в матрице
  void SetRows(int rows);  // Изменяет количество строк в матрице
  void SetCols(int cols);  // изменяет количество столбцов в матрице

 private:
  int rows_, cols_;
  double** matrix_;
};

#endif  // S21_MATRIX_OOP_H