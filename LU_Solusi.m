%Rafik Badilah_L0325009_Praktikum 5_Metode Numerik

function x = LU_Solusi(A, b)
    tic;
    fprintf('Hasil dekomposisi LU\n');
    [L, U, b_mod] = LU_fwdelim(A, b);

    disp('Matriks L ='); disp(L);
    disp('Matriks U ='); disp(U);

    [n, ~] = size(L);
    y = zeros(n, 1);
    x = zeros(n, 1);

    % Forward substitution: L * y = b_mod
    y(1) = b_mod(1) / L(1, 1);
    for i = 2:n
        y(i) = (b_mod(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end

    % Backward substitution: U * x = y
    x(n) = y(n) / U(n, n);
    for i = n-1:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
    toc;
end

function [L, U, b] = LU_fwdelim(A, b)
    [n, ~] = size(A);
    L = eye(n); % Inisialisasi matriks L sebagai matriks identitas

    for i = 1:n-1
        [pivot, k] = max(abs(A(i:n, i)));
        if k ~= 1
            % Penukaran baris pada A dan b
            temp1 = A(i, :);
            temp2 = b(i, :);
            A(i, :) = A(i+k-1, :);
            b(i, :) = b(i+k-1, :);
            A(i+k-1, :) = temp1;
            b(i+k-1, :) = temp2;

            % Penukaran elemen L yang sudah tercatat sebelumnya
            if i > 1
                tempL = L(i, 1:i-1);
                L(i, 1:i-1) = L(i+k-1, 1:i-1);
                L(i+k-1, 1:i-1) = tempL;
            end
        end

        for h = i+1:n
            m = A(h, i) / A(i, i);
            L(h, i) = m; % Simpan pengali m ke matriks L
            A(h, :) = A(h, :) - m * A(i, :); % Eliminasi matriks A saja
        end
    end

    U = A; % Matriks A hasil eliminasi menjadi matriks U
end
