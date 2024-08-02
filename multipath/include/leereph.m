function [eph]=leereph(archivoeph)
%LEEREPH Lee un fichero de mensaje de navegación RINEX
%   	   formatea los datos en una matriz de 35 filas
%        y tantas columnas como número de observaciones
%        contenga el fichero de efemérides.
%        Primero comprueba que no exista nigun fichero archivoeph.mat
%        Si existe lo lee directamente, si no existe lo crea.
% Ej: eph=leereph('95oct18casa____r0.eph');
% Filas de eph:
% 1:#sat 2:año 3:mes 4:dia 5:hora 6:minutos 7:segundos 8:a0 9:a1 10:a2 11:aode 12:crs
% 13:deltan 14:Mo 15:cuc 16:ecc 17:cus 18:raiza 19:toe 20:cic 21:Omega 22:cis 23:io
% 24:crc 25:omega 26:omegadot 27:idot 28:codes 29:weekno 30:L2flag 31:svaccuracy
% 32:svhealth 33:tgd 34:aodc 35:txtime
%
% Carles

% Unidades en segundos, metros o radianes
disp('Iniciando proceso de lectura de efemerides...')

% Quitamos la extensión al nombre de archivo
nombrearchivo=strtok(archivoeph,'.');
mat='.mat';
nombremat=strcat(nombrearchivo,mat);
% En nombremat tenemos el mismo nombre de archivo con la extension .mat

fidmat = fopen(nombremat,'r');
if fidmat == -1,
   % El correspondiente archivo .mat no existe
   % Realizamos la lectura del fichero y creación de variable eph
   
   fide = fopen(archivoeph);
   if fide ==-1, disp('Error abriendo fichero de efemérides'); eph=0; break; end;
   
   % Saltamos la cabecera
   head_lines = 0;
   while 1			    
      head_lines = head_lines+1;
      line = fgetl(fide);
      answer = findstr(line,'END OF HEADER');
      if  ~isempty(answer), break;  end;
   end; %end while
   
   noeph = -1;
   while 1
      noeph = noeph+1;
      line = fgetl(fide);
      if line == -1, break;  end
   end; %end while
   noeph = noeph/8;
   fprintf('Hay un total de %3.0f efemerides disponibles\n',noeph)
   frewind(fide);
   for i = 1:head_lines, line = fgetl(fide); end; % Saltamos cabecera

   % Guardar memoria para los datos
   svprn = zeros(1,noeph);
   anyo = zeros(1,noeph);
   mes = zeros(1,noeph);
   dia = zeros(1,noeph);
   hora = zeros(1,noeph);
   minutos = zeros(1,noeph);
   segundos = zeros(1,noeph);
   a0 = zeros(1,noeph);
   a1 = zeros(1,noeph);
   a2 = zeros(1,noeph);
   weekno = zeros(1,noeph);
   t0c   = zeros(1,noeph);
   tgd   = zeros(1,noeph);
   aodc  = zeros(1,noeph);
   toe   = zeros(1,noeph);
   af2   = zeros(1,noeph);
   af1   = zeros(1,noeph);
   af0   = zeros(1,noeph);
   aode  = zeros(1,noeph);
   deltan = zeros(1,noeph);
   Mo    = zeros(1,noeph);
   ecc   = zeros(1,noeph);
   raiza = zeros(1,noeph);
   toe   = zeros(1,noeph);
   cic   = zeros(1,noeph);
   crc   = zeros(1,noeph);
   cis   = zeros(1,noeph);
   crs   = zeros(1,noeph);
   cuc   = zeros(1,noeph);
   cus   = zeros(1,noeph);
   Omega = zeros(1,noeph);
   omega = zeros(1,noeph);
   io = zeros(1,noeph);
   omegadot = zeros(1,noeph);
   idot  = zeros(1,noeph);
   svaccuracy = zeros(1,noeph);
   svhealth = zeros(1,noeph);
   codes = zeros(1,noeph);
   L2flag = zeros(1,noeph);
   tgd = zeros(1,noeph);
   iodc = zeros(1,noeph);
   txtime = zeros(1,noeph);
   svaccuracy = zeros(1,noeph);
   
   for i = 1:noeph
      line = fgetl(fide);	  %
      svprn(i) = str2num(line(1:2));
      anyo(i) = str2num(line(3:6));
      mes(i) = str2num(line(7:9));
      dia(i) = str2num(line(10:12));
      hora(i) = str2num(line(13:15));
      minutos(i) = str2num(line(16:18));
      segundos(i) = str2num(line(19:22));
      a0(i) = str2num(line(23:41));
      a1(i) = str2num(line(42:60));
      a2(i) = str2num(line(61:79));
      line = fgetl(fide);	  %
      aode(i) = str2num(line(5:22));
      crs(i) = str2num(line(23:41));
      deltan(i) = str2num(line(42:60));
      Mo(i) = str2num(line(61:79));
      line = fgetl(fide);	  %
      cuc(i) = str2num(line(5:22));
      ecc(i) = str2num(line(23:41));
      cus(i) = str2num(line(42:60));
      raiza(i) = str2num(line(61:79));
      line=fgetl(fide);       %
      toe(i) = str2num(line(5:22));
      cic(i) = str2num(line(23:41));
      Omega(i) = str2num(line(42:60));
      cis(i) = str2num(line(61:79));
      line = fgetl(fide);	    %
      io(i) =  str2num(line(5:22));
      crc(i) = str2num(line(23:41));
      omega(i) = str2num(line(42:60));
      omegadot(i) = str2num(line(61:79));
      line = fgetl(fide);	    %
      idot(i) = str2num(line(5:22));
      codes(i) = str2num(line(23:41));
      weekno(i) = str2num(line(42:60));
      L2flag(i) = str2num(line(61:79));
      line = fgetl(fide);	    %
      svaccuracy(i) = str2num(line(5:22));
      svhealth(i) = str2num(line(23:41));
      tgd(i) = str2num(line(42:60));
      iodc(i) = str2num(line(61:79));
      line = fgetl(fide);	    %
      txtime(i) = str2num(line(5:22));
      %spare = line(23:41);
      %spare = line(42:60);
      %spare = line(61:79);
   end %bucle i
   
   estado = fclose(fide);
   if estado==-1, disp('Error cerrando el fichero de efemérides'), end;
   
   %  Descripción de la variable eph.
   eph(1,:) = svprn;    % numero de satelite
   eph(2,:) = anyo;
   eph(3,:) = mes;
   eph(4,:) = dia;
   eph(5,:) = hora;
   eph(6,:) = minutos;
   eph(7,:) = segundos;
   eph(8,:) = a0;       % coeficiente polinomico de correccion del reloj: sesgo(segundos)
   eph(9,:) = a1;       % coeficiente polinomico de correccion del reloj: drift (seg/seg)
   eph(10,:) = a2;      % coeficiente polinomico de correccion del reloj: drift-rate (seg/seg^2)
   eph(11,:) = aode;    % Age Of Data Ephemeris, en segundos
   eph(12,:) = crs;     % Amplitud de la correccion armonica del termino del seno al radio de la orbita (metros)
   eph(13,:) = deltan;  % diferencia media de movimiento del valor calculado (rad/seg)
   eph(14,:) = Mo;      % Anomalia media en el tiempo de referencia (rad)    
   eph(15,:) = cuc;     % Amplitud de la correccion armonica del termino del coseno al argumento de latitud (rad)
   eph(16,:) = ecc;     % Eccentricidad
   eph(17,:) = cus;     % Amplitud de la correccion armonica del termino del seno al argumento de latitud (rad)
   eph(18,:) = raiza;   % raiz cuadrada del semieje mayor (m^1/2)
   eph(19,:) = toe;     % tiempo de referencia,segundos en la semana GPS (Time of Ephemeris)
   eph(20,:) = cic;     % Amplitud de la correccion armonica del termino del coseno al angulo de inclinacion (rad)
   eph(21,:) = Omega;   % Longitud del nodo ascendente en el tiempo de referencia
   eph(22,:) = cis;     % Amplitud de la correccion armonica del termino del seno al angulo de inclinacion (rad)
   eph(23,:) = io;      % angulo de inclinacion en el tiempo de referencia (rad)
   eph(24,:) = crc;     % Amplitud de la correccion armonica del termino del coseno al radio de la orbita (metros)
   eph(25,:) = omega;   % argumento del perigeo (rad)
   eph(26,:) = omegadot;% velocidad de cambio de ascencion derecha (rad/s)
   eph(27,:) = idot;    % velocidad de cambio de inclinacion (rad/s)
   eph(28,:) = codes;
   eph(29,:) = weekno;
   eph(30,:) = L2flag;
   eph(31,:) = svaccuracy;
   eph(32,:) = svhealth;
   eph(33,:) = tgd;
   eph(34,:) = aodc;
   eph(35,:) = txtime;
   
   disp('Guardando archivo .mat ...')
   % Guardamos variable eph en un fichero .mat con el mismo nombre
   save(nombremat,'eph');
      
else
   % El correspondiente archivo .mat ya se ha creado en una lectura anterior
   fclose('all'); % Lo cerramos; lo leeremos con load
   disp('Leyendo archivo .mat ...')
   load(nombrearchivo); % ahora tenemos en memoria la variable eph
   [treintaycinco,noeph]=size(eph);
   fprintf('Hay un total de %3.0f efemerides disponibles\n',noeph)
end

disp('El proceso de lectura de efemerides se ha completado correctamente')