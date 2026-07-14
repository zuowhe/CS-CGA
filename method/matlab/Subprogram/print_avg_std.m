function [] = print_avg_std(eoer_avg,eoer_std)
fprintf('-----------------------------------------------------------------------------------\n');
fprintf('Avg   %9.5f    %9.5f    %9.5f   %4.2f %11.3f  %11.3f     %7.3f\n',...
        eoer_avg(1),eoer_avg(2),eoer_avg(3),eoer_avg(7),eoer_avg(4),eoer_avg(5),eoer_avg(6));
fprintf('Std   %9.5f    %9.5f    %9.5f   %4.2f  %10.3f  %11.3f     %7.3f\n',...
        eoer_std(1),eoer_std(2),eoer_std(3),eoer_std(7),eoer_std(4),eoer_std(5),eoer_std(6));
fprintf('===================================================================================\n');
end
