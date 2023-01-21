clear; clc

% input
n = input("격자 수 : ");
a = input("열 확산도 : ");
delx = input("격자 간격 : ");
delt = input("시간 간격 : ");
t = input("경과 시간 : ");

% lamda
lamda = a * delt / (delx ^ 2);

% preallocation
plate = zeros(n + 2, n + 2);
L = zeros(n + 2, n + 2);
U = zeros(n + 2, n + 2);
temp = zeros(n + 2, 1);
z = zeros(n + 2, 1);

% boundary condition
plate(1, :) = 100;
plate(end, :) = 0;
plate(:, 1) = 75;
plate(:, end) = 50;

% LU decomposition of Tridiagonal matrix
U(2, 2) = 2 * (1 + lamda); U(2, 3) = - lamda; L(2, 2) = 1; L(3, 2) = - lamda / (2 * (1 + lamda));
for i = 3 : length(plate) - 1
    L(i, i) = 1;
    L(i, i - 1) = - lamda / U(i - 1, i - 1);
    U(i, i) = 2 * (1 + lamda) - L(i, i - 1) * - lamda;
    U(i, i + 1) = - lamda;
end
U(:, end) = 0;

% iterate until t(sec)
for tt = 1 : t / delt
    % step 1(l -> l + 1/2)
    for i = 2 : length(plate) - 1
        for j = length(plate) - 1 : -1 : 2
            temp(j) = lamda * plate(length(plate) - j + 1, i + 1) + 2 * (1 - lamda) * plate(length(plate) - j + 1, i) + lamda * plate(length(plate) - j + 1, i - 1);
            if j == length(plate) - 1
                temp(j) = temp(j) + lamda * plate(1, j);
            end
            if j == 2
                temp(j) = temp(j) + lamda * plate(end, j);
            end
        end
        % update plate by using LU decomposition
        for j = 2 : length(plate) - 1
            z(j) = temp(j) - L(j, j - 1) * z(j - 1);
        end
        for j = length(plate) - 1 : -1 : 2
            plate(length(plate) - j + 1, i) = (z(j) - U(j, j + 1) * plate(length(plate) - j, i)) / U(j, j);
        end
    end
    
    % step 2(l + 1/2 -> l + 1)
    for j = length(plate) - 1 : -1 : 2
        for i = 2 : length(plate) - 1
            temp(i) = lamda * plate(j + 1, i) + 2 * (1 - lamda) * plate(j, i) + lamda * plate(j - 1, i);
            if i == 2
                temp(i) = temp(i) + lamda * plate(j, 1);
            end
            if i == length(plate) - 1
                temp(i) = temp(i) + lamda * plate(j, end);
            end
        end
        % update plate by using LU decomposition
        for i = 2 : length(plate) - 1
            z(i) = temp(i) - L(i, i - 1) * z(i - 1);
        end
        for i = length(plate) - 1 : -1 : 2
            plate(j, i) = (z(i) - U(i, i + 1) * plate(j, i + 1)) / U(i, i);
        end
    end
end