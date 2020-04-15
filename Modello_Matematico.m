%  PROGETTO DI MECCANICA APPLICATA II
%  Modello matematico di un cilindro rotante nel vuoto, colpito all'istante
%  0 da una forza impulsiva.

%  Il sdr solidale pone l'asse z allineato con l'asse principale del
%  cilindro. Il sdr inerziale è allineato in t=0 al sdr solidale.

%  Questo script risolve un sistema di equazioni differenziali e calcola la
%  velocità angolare del cilindro (SDR solidale) e la sua matrice di rotazione.  

clearvars
global I

%% Definizione parametri
disp('Definizione Parametri..')

r= 0.5; %raggio(m)
l= 1; %lunghezza(m)
dens = 2700; %alluminio(kg/m^3)
%m= pi*r^2*l*dens; %massa(kg)
m = 2000;
Vin= 5; %velocità iniziale(m/s)
Fy= 2000; %forza impulsiva(N)

%% Definizione Matrici
disp('Definizione del Tensore d''Inerzia e calcolo delle condizioni iniziali..')

I = m*[r^2/4+l^2/12 0 0;0 r^2/4+l^2/12 0; 0 0 r^2/2]; %Tensore d'inerzia
W0 = [0; 0; Vin/r] + (I^-1)*[Fy*l; 0; 0]; %Velocità angolare in 0+

% %% Creazione equazioni
% syms t Ixx Iyy Izz x11 x12 x13 x21 x22 x23 x31 x32 x33 wx wy wz
% R = [x11 x12 x13; x21 x22 x23; x31 x32 x33];
% W = [wx wy wz];
% I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
% eq1 = diff(W,t) = -(I^-1)*cross(W, I*W);
% eq2 = diff(R,t) == R*W;

%% Risolvo le equazioni differenziali
disp('Risolzione delle equazioni differenziali..')

CI = [W0;1;0;0;0;1;0;0;0;1];
Tf = 10;
% { W_dot = -(I^-1)*cross(W,I*W);
% { R_dot = R*W_skew;
[t,Results] = ode45(@Rotation3D, [0 Tf], CI);

%% Creo il vettore delle velocità angolari
disp('Creazione del vettore delle velocità angolari nel sistema solidale..')

Ws = zeros(3,length(t));
for i = 1:length(t)
    Ws(:,i) = [Results(i,1), Results(i,2), Results(i,3)];
end

%% Creo la matrice di rotazione
disp('Creazione della matrice di rotazione..')

R = zeros(length(t),3,3);
for i = 1:length(t)
    R(i,:,:) = [Results(i,4), Results(i,5), Results(i,6); 
                Results(i,7), Results(i,8), Results(i,9);
                Results(i,10), Results(i,11), Results(i,12)];
end

%% Creo un vettore dei tempi equispaziato 

FPS = 40; % Frame per Second desiderati
Frame_Duration = 1/FPS;
Total_Frames = Tf*FPS;
index = zeros(Total_Frames,1); % Vettore che conterrà gli indici degli elementi del vettore t equispaziati
for i = 1:Total_Frames
    val = Frame_Duration*i; % Calcolo l'i-esimo instante di tempo per ottenere gli FPS desiderati  
    [~,index(i)] = min( abs( t-val ) ); % Estraggo dal vettore t l'indice dell'elemento che meglio approssima l'istante di tempo desiderato
end

%% PLOTTING
disp('Plotting dei risultati: ')
disp('  * Per i grafici delle velocità angolari, premere 1.')
disp('  * Per la rotazione del sistema di riferimento, premere 2.')
disp('  * Per uscire dallo script premere 3.')

while(true)
   
   x = input('Inserire un numero: ');
    
   switch x
       case 1
    
        % Plotto le velocità angolari
        disp('Grafici delle velocità angolari')
        
        figure(1);
        subplot(3,1,1)
        plot(t,Ws(1,:)), grid on, title('Wx')
        subplot(3,1,2)
        plot(t,Ws(2,:)), grid on, title('Wy')
        subplot(3,1,3)
        plot(t,Ws(3,:)), grid on, title('Wz')
        
       case 2
        
        % Plotto la variazione della matrice di rotazione
        disp('Plotto la rotazione del sistema di riferimento')
        
        figure(2);
        axis([-1,1,-1,1,-1,1]) % Configuro gli assi del grafico
        grid on
        hold on
        % Creo delle frecce per rappresentare i tre assi
        a1 = quiver3(0,0,0,R(1,1,1),R(1,1,2),R(1,1,3),'r'); % AsseX solidale
        a2 = quiver3(0,0,0,R(1,2,1),R(1,2,2),R(1,2,3),'b'); % AsseY solidale
        a3 = quiver3(0,0,0,R(1,3,1),R(1,3,2),R(1,3,3),'g'); % AsseZ solidale
        
        for i=1:Total_Frames
            tic % Inizio a misurare il tempo di esecuzione
            a1.UData = R(index(i),1,1); % Aggiorno l'asse X in base alle tre componenti della prima colonna di R
            a1.VData = R(index(i),1,2);
            a1.WData = R(index(i),1,3);
            a2.UData = R(index(i),2,1); % Aggiorno l'asse Y in base alle tre componenti della seconda colonna di R
            a2.VData = R(index(i),2,2);
            a2.WData = R(index(i),2,3);
            a3.UData = R(index(i),3,1); % Aggiorno l'asse Z in base alle tre componenti della terza colonna di R
            a3.VData = R(index(i),3,2);
            a3.WData = R(index(i),3,3);
            pause(Frame_Duration - toc); % Imposto una pausa pari alla durata desiderata del frame, meno il tempo di esecuzione delle istruzioni
        end
        
       otherwise
        
        break;
        
   end
    
end

%% Salvo i risultati su file
save('cylinder.mat','t','R','Ws');

% Fine dello script
disp('Script end')

%% Funzione: Rotation3D_Velocity
function [Xdot] = Rotation3D(~,x)
global I

% { W_dot = -(I^-1)*cross(W,I*W);
% { R_dot = R*W_skew;

Xdot = [ (I(2,2)*x(2)*x(3) - I(3,3)*x(2)*x(3))/I(1,1);
        -(I(1,1)*x(1)*x(3) - I(3,3)*x(1)*x(3))/I(2,2);
         (I(1,1)*x(1)*x(2) - I(2,2)*x(1)*x(2))/I(3,3);
          x(3)*x(5) - x(2)*x(6);
          x(1)*x(6) - x(3)*x(4);
          x(2)*x(4) - x(1)*x(5);
          x(3)*x(8) - x(2)*x(9);
          x(1)*x(9) - x(3)*x(7);
          x(2)*x(7) - x(1)*x(8);
          x(3)*x(11) - x(2)*x(12);
          x(1)*x(12) - x(3)*x(10);
          x(2)*x(10) - x(1)*x(11) ];  
end
