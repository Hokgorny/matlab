clear; clc                                  % 초기화

h = 0.00001;                                % η 간격
ii = 0 : h : 8;                             % η 범위
lenii = length(ii);                         % η 배열 길이
ll = 0; rl = 8; f2_init = (rl + ll) / 2;    % 이분법 초기값

% 사격법
while true
    [~, f1, ~] = Blasius_RK4(f2_init, h, lenii);
    % 이분법
    if f1(end) < 1
        ll = f2_init;
        f2_init = (f2_init + rl) / 2;
    elseif f1(end) > 1
        rl = f2_init;
        f2_init = (f2_init + ll) / 2;
    elseif f1(end) - 1 < 0.001
        break;
    end
end

[f0, f1, f2] = Blasius_RK4(f2_init, h, lenii);

plot(ii, f0, 'r');
hold on
plot(ii, f1, 'b');
plot(ii, f2, 'g');
grid on
title(['f''''(0) : ', num2str(f2_init)]);
axis([0 8 0 2]);
xlabel('η');
ylabel('f, f'' and f''''');
legend('f(η)', 'f''(η)', 'f''''(η)');