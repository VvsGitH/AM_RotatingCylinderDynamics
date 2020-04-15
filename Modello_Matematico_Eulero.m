%  PROGETTO DI MECCANICA APPLICATA II - parte2
%  Modello matematico di un cilindro rotante nel vuoto, colpito all'istante
%  0 da una forza impulsiva. 

%  Il sdr solidale pone l'asse z allineato con l'asse principale del
%  cilindro. Il sdr inerziale è allineato in t=0 al sdr solidale.

%  Questo script, a partiredalla matrice di rotazione R, calcolata nella
%  parte1, calcolata l'andamento nel tempo degli angoli di Eulero 
%  (precessione, nutazione e rollio).

%% Generazione delle equazioni
% A partire dagli angoli di Eulero è possibile definire tre matrici
%   R1 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]; %Da XYZ(SDR FISSO) a x'y'z'
%   R2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]; %Da x'y'z' a x''y''z''
%   R3 = [cos(xi) sin(xi) 0; -sin(xi) cos(xi) 0; 0 0 1]; %Da x''y''z'' a xyz(SDR SOLIDALE)
% La matrice di rotazione R, calcolata nella parte 1 è uguale a:
%   R = R3*R2*R1
% Ho dunque a disposizione 9 equazioni (una per ogni elemento di R) in tre
% incognite (phi, theta e xi).
% Ho scelto le tre equazioni più semplici:
%   sin(phi)*sin(theta) = R31
%   sin(xi)*sin(theta) = R13
%   cos(theta) = R33
% Ho poi isolato phi, theta e xi, prima di utilizzarle di seguito.

%% Carico i dati
if exist('Ws','var')== 0
    load('cylinder.mat');
end
t1 = t;
R13 = R(:,1,3);
R33 = R(:,3,3);
R31 = R(:,3,1);

%% Calcolo gli angoli di Eulero
phi = asin(R31./sqrt(1.-R33.^2));
xi = asin(R13./sqrt(1.-R33.^2));
theta = acos(R33);

%% Plotto gli angoli di Eulero
disp('Plotto degli angoli di Eulero..')

figure(3);
subplot(3,1,1)
plot(t1,phi), grid on, title('Precessione - Asse Z')
subplot(3,1,2)
plot(t1,theta), grid on, title('Nutazione - Asse X1')
subplot(3,1,3)
plot(t1,xi), grid on, title('Rollio - Asse Z2')

%% Calcolo le velocità di rotazione approssimate
disp('Calcolo le velocità di rotazione approssimate..')

Dt = diff(t1); % Vettore [t2-t1, .., ti+1 - t1, .., tn-tn-1]
Dphi = diff(phi);
Dtheta = diff(theta);
Dxi = diff(xi);

Wphi = Dphi./Dt; % Calcolo la velocità angolare come delta_phi/delta_t
Wtheta = Dtheta./Dt;
Wxi = Dxi./Dt;

%% Plotto le velocità di rotazione di Eulero
disp('Plotto le velocità di rotazione..')

figure(4);
subplot(3,1,1)
plot(t1(1:end-1),Wphi), grid on, title('Velocità di Precessione - Asse Z')
ylim([-20, 20])
subplot(3,1,2)
plot(t1(1:end-1),Wtheta), grid on, title('Velocità di Nutazione - Asse X1')
subplot(3,1,3)
plot(t1(1:end-1),Wxi), grid on, title('Velocità di Rollio - Asse Z2')

% Fine dello script
disp('Script end')