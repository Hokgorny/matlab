function [f0_c, f1_c, f2_c] = Blasius_RK4(f2_init, h, lenii) % f'' 초기값, h, length(ii)
    Fd3 = @(x, y) -1/2 * x * y;     % Blasius 방정식

    f2 = zeros(1, lenii);           % 사전할당
    f1 = zeros(1, lenii);
    f0 = zeros(1, lenii);
    
    f0(1) = 0; f1(1) = 0; f2(1) = f2_init;  % 초기값 지정

    for i = 1 : (lenii - 1)         % 4차 Runge-Kutta
        k1_f2 = Fd3(f0(i), f2(i));  
        k1_f1 = f2(i);
        k1_f0 = f1(i);
    
        k2_f2 = Fd3(f0(i) + 1/2 * k1_f0 * h, f2(i) + 1/2 * k1_f2 * h);
        k2_f1 = f2(i) + 1/2 * k1_f1 * h;
        k2_f0 = f1(i) + 1/2 * k1_f0 * h;
    
        k3_f2 = Fd3(f0(i) + 1/2 * k2_f0 * h, f2(i) + 1/2 * k2_f2 * h);
        k3_f1 = f2(i) + 1/2 * k2_f1 * h;
        k3_f0 = f1(i) + 1/2 * k2_f0 * h;
    
        k4_f2 = Fd3(f0(i) + k3_f0 * h, f2(i) + k3_f2 * h);
        k4_f1 = f2(i) + k3_f1 * h;
        k4_f0 = f1(i) + k3_f0 * h;
    
        f2(i + 1) = f2(i) + 1/6 * (k1_f2 + 2 * k2_f2 + 2 * k3_f2 + k4_f2) * h;
        f1(i + 1) = f1(i) + 1/6 * (k1_f1 + 2 * k2_f1 + 2 * k3_f1 + k4_f1) * h;
        f0(i + 1) = f0(i) + 1/6 * (k1_f0 + 2 * k2_f0 + 2 * k3_f0 + k4_f0) * h;
    end

    f0_c = f0;
    f1_c = f1;
    f2_c = f2;
end