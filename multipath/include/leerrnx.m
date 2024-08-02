function [L1,L2,P1,P2,C1,temps,TOW,sat_PRNs]=leerrnx(fobservacion)
% LEERRNX  lee un fichero RINEX de observaciones
%          y formatea los datos en matrices donde
%          cada fila contiene las observaciones de una �poca
%          y las columnas corresponden a los sat�lites
%
% [L1,L2,P1,P2,C1,temps]=leerrnx('ebre1050.00o');
%
% L1 y L2: unidades en ciclos completos de portadora
% P1 y P2: unidades en metros
% temps es una estructura con los campos:
% anyo mes dia hora minutos segundos
% Primero comprueba que no exista nigun fichero nombrearchivo.mat
% Si existe lo lee directamente, si no existe lo crea.
% Carles

disp('Iniciando el proceso de lectura de observaciones.');

   fid = fopen(fobservacion);
   if fid==-1, disp('Error abriendo fichero de datos'), 
      disp ('El proceso de lectura de datos no se ha podido completar correctamente')
      L1=0;L2=0;P1=0;P2=0;C1=0;temps=0;
      return; %salimos del programa
   end;
   
   % Buscamos informaci�n de los tipos de observables en la cabecera
   lineas_cabecera=0;
   while 1
      lineas_cabecera=lineas_cabecera+1;
      linea = fgetl(fid);
      %
      observables = findstr(linea,'# / TYPES OF OBSERV');
      if ~isempty(observables),
         [numtip,resto]=strtok(linea); 
         numtipnom=str2num(numtip); %numtipnom: n�mero de tipos de observables nominal
         if numtipnom==6
            disp('Este fichero contiene demasiados tipos de datos')
            ok=0;
            break; %salimos del while
         end %numtipnom
         
         orden=zeros(1,numtipnom);
         for ntipo=1:numtipnom
            [tip,resto]=strtok(resto);
            switch tip,
            case 'L1',
               orden(ntipo)=1;
            case 'L2',
               orden(ntipo)=2;
            case 'C1',
               orden(ntipo)=3;
            case 'P1',
               orden(ntipo)=4;
            case 'P2',
               orden(ntipo)=5;
            case 'D1',
               disp('El fichero contiene datos Doppler en L1')
               orden(ntipo)=6;
            case 'D2',
               disp('El fichero contiene datos Doppler en L2')
               orden(ntipo)=6;
            case 'T1',
               disp('El fichero contiene datos Transit (Doppler en 150 MHz)')
               orden(ntipo)=6;
            case 'T2',
               disp('El fichero contiene datos Transit (Doppler en 400 MHz)')
               orden(ntipo)=6;
   
            otherwise break;
            end %switch tip
         end % bucle ntipo
      end %linea de tipos de observables
      %
      final = findstr(linea,'END OF HEADER');
      if ~isempty(final), break; end;
   end; %while
   % Si salimos del bucle es porque ya hemos acabado la cabecera
   
   % Contamos el numero de observaciones
   numepocas=0;
   while 1
      % Esto es muy lento
      linea = fgetl(fid);
      if ~isstr(linea), 
         break;
      end;
      final = findstr(linea,'END OF FILE');
      if ~isempty(final),
         break; 
      end;
      if isempty(str2num(linea(10))), numepocas = numepocas+1; end;
   end;
   
   disp(strcat('N�mero de �pocas observadas: ',num2str(numepocas)))
   % num=input('�Cuantas �pocas quiere procesar?');
   num=numepocas;
   frewind(fid); %Nos situamos al principio del fichero
   for i=1:lineas_cabecera, linea=fgetl(fid); end;  %Saltamos la cabecera
   
   
   for epoca=1:num
      % Leer linea de cabecera de �poca
      linea = fgetl(fid);
      if ~isstr(linea) break; end % a veces los ficheros se acaban "de golpe"
      [kk,longlinea]=size(linea);
      if longlinea<32 break; end
      numobs=str2num(linea(31:32)); %N�mero de sat�lites observados
      if longlinea<(32+3*numobs) break; end
      for n=1:numobs
         satelite(n)=str2num(linea(31+3*n:32+3*n));
      end;
      
      % Estructura temps
      temps(epoca,1)= struct('anyo',str2num(linea(2:3)),'mes',str2num(linea(5:6)),'dia',str2num(linea(8:9)),...
         'hora',str2num(linea(11:12)),'minutos',str2num(linea(14:15)),'segundos',str2num(linea(17:26)));
      
        %computation of the corresponding julian day
        jd = julday(temps(epoca,1).anyo+2000, temps(epoca,1).mes, temps(epoca,1).dia, 0);
        %computation of the GPS time in weeks and seconds of week
        [week, sec_of_week] = gps_time(jd); %#ok<ASGLU>
        TOW(epoca) = sec_of_week + temps(epoca,1).hora*3600+temps(epoca,1).minutos*60+temps(epoca,1).segundos;
        
      % Leer lineas de datos
      for sat=1:numobs,
         datos=fgetl(fid);
         % L1 es una matriz donde cada fila es el conjunto
         % de observaciones L1 durante una �poca, y la columna
         % representa el n�mero del satelite
         
         [kk,numcar]=size(datos);
         clear d;
         if exist('C1')==0
             L1=zeros(numepocas,numobs);
             L2=zeros(numepocas,numobs);
             C1=zeros(numepocas,numobs);
             P1=zeros(numepocas,numobs);
             P2=zeros(numepocas,numobs);
         end
         if numcar==80
            d=[str2num(datos(2:14)) str2num(datos(18:30))...
                     str2num(datos(34:46)) str2num(datos(50:62)) str2num(datos(66:78))];

            if length(d)>numobs
                disp('Error reading observations');
                break;
            end
         end
     
         for index=1:numtipnom, %para cada uno de los tipos de datos
            switch orden(index),
            case 1, %medida de fase en L1
               L1(epoca,sat)=d(index);
            case 2, %medida de fase en L2
               L2(epoca,sat)=d(index);
            case 3, %pseudodistancia utilizando c�digo C/A en L1
               C1(epoca,sat)=d(index);
            case 4, %pseudodistancia utilizando c�digo P en L1
               P1(epoca,sat)=d(index);
            case 5, %pseudodistancia utilizando c�digo P en L2
               P2(epoca,sat)=d(index);
            case 6,
               % se trata de T1, T2, D1, D2
               % no hacemos nada
               break;
            otherwise,
               disp('Error leyendo los datos');
               break;
            end %switch
         end %bucle index
      end; %bucle sat
   end; %bucle epoca
   
   sat_PRNs=satelite;
   fclose all; % Cerramos el fichero original de datos

   
