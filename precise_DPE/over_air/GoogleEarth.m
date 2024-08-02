function GoogleEarth(navSolutions,settings)
% GoogleEarth(navSolutions,settings)
% Creates KML file for viewing in Google Earth

long = sum(navSolutions.longitude)/length( navSolutions.longitude );
lat = sum(navSolutions.latitude)/length( navSolutions.latitude );
height = sum(navSolutions.height)/length( navSolutions.height );

% Create a KML file (Keyhole Markup Language)
% http://code.google.com/apis/kml/
fname = settings.fileName;
% Extract .bin
lf=length(fname);
if ~isunix,
    ind = regexp(fname,'\');
    if isempty(ind),
        ind = regexp(fname,'/');
    end

else
    ind = regexp(fname,'/');
end
if isempty(ind),
    ind=0;
end
fname_no_extension = fname((ind(end)+1):lf-4);
fileNameStr = strcat(fname_no_extension,'.kml');
[fid, message] = fopen(fileNameStr, 'wt');
if (fid > 0),
    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid,'<kml xmlns="http://earth.google.com/kml/2.2">');
    fprintf(fid,'<Document>');
    fprintf(fid,'	<name>%s</name>',fileNameStr);
    fprintf(fid,'       	<StyleMap id="msn_ylw-pushpin">\n');
    fprintf(fid,'		<Pair>\n');
    fprintf(fid,'			<key>normal</key>\n');
    fprintf(fid,'			<styleUrl>#sn_ylw-pushpin</styleUrl>\n');
    fprintf(fid,'		</Pair>\n');
    fprintf(fid,'		<Pair>\n');
    fprintf(fid,'			<key>highlight</key>\n');
    fprintf(fid,'			<styleUrl>#sh_ylw-pushpin</styleUrl>\n');
    fprintf(fid,'		</Pair>\n');
    fprintf(fid,'	</StyleMap>\n');
    fprintf(fid,'	<Style id="sh_ylw-pushpin">\n');
    fprintf(fid,'		<IconStyle>\n');
    fprintf(fid,'			<scale>1.3</scale>\n');
    fprintf(fid,'			<Icon>\n');
    fprintf(fid,'				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n');
    fprintf(fid,'			</Icon>\n');
    fprintf(fid,'			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>\n');
    fprintf(fid,'		</IconStyle>\n');
    fprintf(fid,'	</Style>\n');
    fprintf(fid,'	<Style id="sn_ylw-pushpin">\n');
    fprintf(fid,'		<IconStyle>\n');
    fprintf(fid,'			<scale>1.1</scale>\n');
    fprintf(fid,'			<Icon>\n');
    fprintf(fid,'				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n');
    fprintf(fid,'			</Icon>\n');
    fprintf(fid,'			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>\n');
    fprintf(fid,'		</IconStyle>\n');
    fprintf(fid,'	</Style>\n');
    fprintf(fid,'	<Placemark>\n');
    fprintf(fid,'		<name>%s</name>\n',fname_no_extension);
    fprintf(fid,'		<description>Signal captured by the GN3S Sampler. Processed by Carles Fernandez at CTTC http://www.cttc.cat</description>\n');
    fprintf(fid,'		<LookAt>\n');
    fprintf(fid,'			<longitude>%17.15f</longitude>\n',long);
    fprintf(fid,'			<latitude>%16.14f</latitude>\n',lat);
    fprintf(fid,'			<altitude>%4.2f</altitude>\n',height);
    fprintf(fid,'			<tilt>0</tilt>\n');
    fprintf(fid,'			<heading>0</heading>\n');
    fprintf(fid,'		</LookAt>\n');
    fprintf(fid,'		<styleUrl>#msn_ylw-pushpin</styleUrl>\n');
    fprintf(fid,'		<Point>\n');
    fprintf(fid,'			<coordinates>%17.15f,%16.14f,%4.2f</coordinates>\n',long,lat,height);
    fprintf(fid,'		</Point>\n');
    fprintf(fid,'    </Placemark>\n');
    
    fprintf(fid,'	 <Placemark>\n');
    fprintf(fid,'	  <name>GNSS-SDR PVT</name>\n');
    fprintf(fid,'	  <description>GNSS-SDR position log</description>\n');
    fprintf(fid,'	  <styleUrl>#yellowLineGreenPoly</styleUrl>\n');
    fprintf(fid,'	  <LineString>\n');
    fprintf(fid,'	  <extrude>0</extrude>\n');
    fprintf(fid,'	  <tessellate>1</tessellate>\n');
    fprintf(fid,'	  <altitudeMode>absolute</altitudeMode>\n');
    fprintf(fid,'			<coordinates>\n');
    for ind=1:1:length(navSolutions.longitude)
       fprintf(fid,'			%17.15f,%16.14f,%4.2f\n',navSolutions.longitude(ind),navSolutions.latitude(ind),navSolutions.height(ind));
    end
    fprintf(fid,'			</coordinates>\n');
    fprintf(fid,'	  </LineString>\n');
    fprintf(fid,'	  </Placemark>\n');
    fprintf(fid,'</Document>\n');
    fprintf(fid,'</kml>\n');
    fclose(fid);

    if not(isunix), % Windows, execute Google Earth
        webtosee=['http://maps.google.com/maps?q=' num2str(lat) ',' num2str(long)];
        web(webtosee)
        current_dir=pwd;
        file_kml=[current_dir '\' fileNameStr];
        cd('C:\Program Files\Google\Google Earth Pro\client');
        toexecute=['googleearth.exe ' '"' file_kml '"'  ];
        system(toexecute);
        cd(current_dir);
    else %we are in Unix / Linux environment
%         webtosee=['http://maps.google.com/maps?q=' num2str(lat) ',' num2str(long)];
%         web(webtosee)
        %disp('There is no Google Earth for Linux yet...')

    end
else
    %=== Error while opening the data file ================================
    error('Unable to read file %s: %s.', fileNameStr, message);
end % if fid
