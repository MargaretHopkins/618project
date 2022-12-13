function nz = addpipe( oz, f, l )
%   NZ = ADDPIPE( OZ, F, L ) takes frequency-domain input impedance
%   values provided in vector OZ at frequencies F and returns a new
%   theoretical input impedance in vector NZ with an additional
%   cylindrical pipe of length L added to the input.
%
%   This is used to simulate the effect of a mouthpiece at the input
%   of an instrument when the input impedance is measured without a
%   mouthpiece present.  The impedance data is assumed to be
%   normalized with respect to the characteristic impedance at the
%   entrance.
%
%   By Gary P. Scavone, McGill University, 21 May 2008.

if nargin~=3
  error('addpipe: Number of arguments is incorrect.');
end

% Physical constants and evaluation frequencies

omega = 2 * pi * f;     % radian frequencies
c = 347.23;             % speed of sound (m/sec) %
%rho = 1.1769;           % density of air (kg/m^3) %
k = omega/c;

% Pipe characteristic impedance

%Zo = rho * c ./ ( pi * ra .^ 2 );
Zo = 1;

% Lossless Cylindrical Transmission Matrix Elements

sinL = sin(l*k);
cosL = cos(l*k);
A = cosL;
B = 1i * (Zo * ones(size(k))) .* sinL ;
C = 1i * sinL ./ (Zo * ones(size(k)));
D = cosL;

nz = ( B + A.*oz ) ./ ( D + C .* oz );
