%Rafik Badilah_L0325009_Praktikum 5_Metode Numerik

% Matriks koefisien A dan vektor b dari Soal
A = [ 2  1 -1;
      4  3  1;
     -2  1  2];

b = [3; 9; 4];

disp('=================== 1. ELIMINASI GAUSS ===================');
x_gauss = EliminasiGauss(A, b);
disp('Solusi x (Gauss):');
disp(x_gauss);

disp('================ 2. ELIMINASI GAUSS-JORDAN ================');
x_jordan = EliminasiGaussJordan(A, b);
disp('Solusi x (Gauss-Jordan):');
disp(x_jordan);

disp('===================== 3. DEKOMPOSISI LU ===================');
x_lu = LU_Solusi(A, b);
disp('Solusi x (LU Decomposition):');
disp(x_lu);

