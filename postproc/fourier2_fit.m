function [x,y,fitresult, gof] = fourier2_fit(x0, y0)

[xData, yData] = prepareCurveData( x0, y0 );

ft = fittype( 'fourier2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.StartPoint = [0 0 0 0 0 0.0556034097980494];
[fitresult, gof] = fit( xData, yData, ft, opts );

x = x0;
y = fitresult.a0 + fitresult.a1*cos(x*fitresult.w) + fitresult.b1*sin(x*fitresult.w) + fitresult.a2*cos(2*x*fitresult.w) + fitresult.b2*sin(2*x*fitresult.w);

% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y5 vs. x5', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% xlabel( 'x5', 'Interpreter', 'none' );
% ylabel( 'y5', 'Interpreter', 'none' );
% grid on

end

