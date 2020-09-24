% Verify nonnegativity

tmp = Hq(:)' * X(:) + 2 * gq' * x + fq;
nost = solvesdp(con, tmp, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end
fprintf('verify nonneg = %f\n', double(tmp))

tmp = 2 * gl' * x + fl;
nost = solvesdp(con, tmp, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end
fprintf('verify nonneg = %f\n', double(tmp))

tmp = Hq(:)' * X(:) + 2 * (gq + gl)' * x + (fq + fl) - qlmin;
nost = solvesdp(con, tmp, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end
fprintf('verify nonneg = %f\n', double(tmp))

con = [con; ...
    (r + R) * R * (Hq(:)' * X(:) + 2 * gq' * x + fq) + ...
    (r + R) * (2 * gl' * X * b + (fl * b - 2 * a * gl)' * x - a * fl) >= ...
    qlmin * trace(X) + r * R * (Hq(:)' * X(:) + 2 * (gq + gl)' * x + (fq + fl)) ...
    - (2 * gl' * X * c + fl * c' * x) - cmax * R * (2 * gl' * x + fl)
];
